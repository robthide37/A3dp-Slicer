///|/ Copyright (c) Prusa Research 2022 Enrico Turri @enricoturri1966, Vojtěch Bubník @bubnikv
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#include "Thumbnails.hpp"
#include "../miniz_extension.hpp"
#include "../format.hpp"

#include <qoi/qoi.h>
#include <jpeglib.h>
#include <jerror.h>

#include <boost/algorithm/string.hpp>
#include <string>

namespace Slic3r::GCodeThumbnails {

using namespace std::literals;

struct CompressedPNG : CompressedImageBuffer 
{
    ~CompressedPNG() override { if (data) mz_free(data); }
    std::string_view tag() const override { return "thumbnail_PNG"sv; }
};

struct CompressedJPG : CompressedImageBuffer
{
    ~CompressedJPG() override { free(data); }
    std::string_view tag() const override { return "thumbnail_JPG"sv; }
};

struct CompressedQOI : CompressedImageBuffer
{
    ~CompressedQOI() override { free(data); }
    std::string_view tag() const override { return "thumbnail_QOI"sv; }
};

struct CompressedBIQU : CompressedImageBuffer
{
    ~CompressedBIQU() override { free(data); }
    std::string_view tag() const override { return "thumbnail_BIQU"sv; }
};

std::unique_ptr<CompressedImageBuffer> compress_thumbnail_png(const ThumbnailData& data)
{
    auto out = std::make_unique<CompressedPNG>();
    out->data = tdefl_write_image_to_png_file_in_memory_ex((const void*)data.pixels.data(), data.width, data.height, 4, &out->size, MZ_DEFAULT_LEVEL, 1);
    return out;
}

std::unique_ptr<CompressedImageBuffer> compress_thumbnail_biqu(const ThumbnailData& data)
{
    // Take vector of RGBA pixels and flip the image vertically
    std::vector<uint8_t> rgba_pixels(data.pixels.size());
    const size_t row_size = data.width * 4;
    for (size_t y = 0; y < data.height; ++y)
        ::memcpy(rgba_pixels.data() + (data.height - y - 1) * row_size, data.pixels.data() + y * row_size, row_size);

    auto out = std::make_unique<CompressedBIQU>();
    //size: height is number of lines. Add 2 byte to each line for the ';' and '\n'. Each pixel is 4 byte, +1 for the 0 of the c_str
    out->size = data.height * (2 + data.width * 4) + 1;
    out->data = malloc(out->size);

    int idx = 0;
    std::stringstream tohex;
    tohex << std::setfill('0') << std::hex;
    for (size_t y = 0; y < data.height; ++y) {
        tohex << ";";
        for (size_t x = 0; x < data.width; ++x) {
            uint16_t pixel = 0;
            //r
            pixel |= uint16_t((rgba_pixels[y * row_size + x * 4 + 0 ] & 0x000000F8) >> 3);
            //g
            pixel |= uint16_t((rgba_pixels[y * row_size + x * 4 + 1 ] & 0x000000FC) << 3);
            //b
            pixel |= uint16_t((rgba_pixels[y * row_size + x * 4 + 2 ] & 0x000000F8) << 8);
            tohex << std::setw(4) << pixel;
        }
        tohex << "\n";
    }
    std::string str = tohex.str();
    assert(str.size() + 1 == out->size);
    ::memcpy(out->data, (const void*)str.c_str(), out->size);
    return out;
}

std::unique_ptr<CompressedImageBuffer> compress_thumbnail_jpg(const ThumbnailData& data)
{
    // Take vector of RGBA pixels and flip the image vertically
    std::vector<unsigned char> rgba_pixels(data.pixels.size());
    const unsigned int row_size = data.width * 4;
    for (unsigned int y = 0; y < data.height; ++y) {
        ::memcpy(rgba_pixels.data() + (data.height - y - 1) * row_size, data.pixels.data() + y * row_size, row_size);
    }

    // Store pointers to scanlines start for later use
    std::vector<unsigned char*> rows_ptrs;
    rows_ptrs.reserve(data.height);
    for (unsigned int y = 0; y < data.height; ++y) {
        rows_ptrs.emplace_back(&rgba_pixels[y * row_size]);
    }

    std::vector<unsigned char> compressed_data(data.pixels.size());
    unsigned char* compressed_data_ptr = compressed_data.data();
    unsigned long compressed_data_size = data.pixels.size();

    jpeg_error_mgr err;
    jpeg_compress_struct info;
    info.err = jpeg_std_error(&err);
    jpeg_create_compress(&info);
    jpeg_mem_dest(&info, &compressed_data_ptr, &compressed_data_size);

    info.image_width = data.width;
    info.image_height = data.height;
    info.input_components = 4;
    info.in_color_space = JCS_EXT_RGBA;

    jpeg_set_defaults(&info);
    jpeg_set_quality(&info, 85, TRUE);
    jpeg_start_compress(&info, TRUE);

    jpeg_write_scanlines(&info, rows_ptrs.data(), data.height);
    jpeg_finish_compress(&info);
    jpeg_destroy_compress(&info);

    // FIXME -> Add error checking

    auto out = std::make_unique<CompressedJPG>();
    out->data = malloc(compressed_data_size);
    out->size = size_t(compressed_data_size);
    ::memcpy(out->data, (const void*)compressed_data.data(), out->size);
    return out;
}

std::unique_ptr<CompressedImageBuffer> compress_thumbnail_qoi(const ThumbnailData &data)
{
    qoi_desc desc;
    desc.width      = data.width;
    desc.height     = data.height;
    desc.channels   = 4;
    desc.colorspace = QOI_SRGB;

    // Take vector of RGBA pixels and flip the image vertically
    std::vector<uint8_t> rgba_pixels(data.pixels.size() * 4);
    size_t row_size = data.width * 4;
    for (size_t y = 0; y < data.height; ++ y)
        memcpy(rgba_pixels.data() + (data.height - y - 1) * row_size, data.pixels.data() + y * row_size, row_size);
    
    auto out = std::make_unique<CompressedQOI>();
    int  size;
    out->data = qoi_encode((const void*)rgba_pixels.data(), &desc, &size);
    out->size = size;
    return out;
}

std::unique_ptr<CompressedImageBuffer> compress_thumbnail(const ThumbnailData &data, GCodeThumbnailsFormat format)
{
    switch (format) {
        case GCodeThumbnailsFormat::PNG:
        default:
            return compress_thumbnail_png(data);
        case GCodeThumbnailsFormat::JPG:
            return compress_thumbnail_jpg(data);
        case GCodeThumbnailsFormat::QOI:
            return compress_thumbnail_qoi(data);
        case GCodeThumbnailsFormat::BIQU:
            return compress_thumbnail_biqu(data);
    }
}

std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list_from_prusa(const std::string& thumbnails_string, const std::string_view def_ext /*= "PNG"sv*/)
{
    if (thumbnails_string.empty())
        return {};

    std::istringstream is(thumbnails_string);
    std::string point_str;

    ThumbnailErrors errors;

    // generate thumbnails data to process it

    GCodeThumbnailDefinitionsList thumbnails_list;
    while (std::getline(is, point_str, ',')) {
        Vec2d point(Vec2d::Zero());
        GCodeThumbnailsFormat format;
        std::istringstream iss(point_str);
        std::string coord_str;
        if (std::getline(iss, coord_str, 'x') && !coord_str.empty()) {
            std::istringstream(coord_str) >> point(0);
            if (std::getline(iss, coord_str, '/') && !coord_str.empty()) {
                std::istringstream(coord_str) >> point(1);

                if (0 < point(0) && point(0) < 1000 && 0 < point(1) && point(1) < 1000) {
                    std::string ext_str;
                    std::getline(iss, ext_str, '/');

                    if (ext_str.empty())
                        ext_str = def_ext.empty() ? "PNG"sv : def_ext;

                    // check validity of extention
                    boost::to_upper(ext_str);
                    if (!ConfigOptionEnum<GCodeThumbnailsFormat>::from_string(ext_str, format)) {
                        format = GCodeThumbnailsFormat::PNG;
                        errors = enum_bitmask(errors | ThumbnailError::InvalidExt);
                    }

                    thumbnails_list.emplace_back(std::make_pair(format, point));
                }
                else
                    errors = enum_bitmask(errors | ThumbnailError::OutOfRange);
                continue;
            }
        }
        errors = enum_bitmask(errors | ThumbnailError::InvalidVal);
    }

    return std::make_pair(std::move(thumbnails_list), errors);
}

std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list(const std::vector<Vec2d>& thumbnails, GCodeThumbnailsFormat format)
{
    if (thumbnails.empty())
        return {};

    ThumbnailErrors errors;
    GCodeThumbnailDefinitionsList thumbnails_list;

    // generate thumbnails data to process it
    for (const Vec2d &point : thumbnails) {
        thumbnails_list.emplace_back(std::make_pair(format, point));
    }

    return std::make_pair(std::move(thumbnails_list), errors);
}

std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list(const ConfigBase& config)
{
    // ??? Unit tests or command line slicing may not define "thumbnails" or "thumbnails_format".
    // ??? If "thumbnails_format" is not defined, export to PNG.

    // generate thumbnails data to process it
    assert(config.option<ConfigOptionPoints>("thumbnails"));
    assert(config.option<ConfigOptionEnum<GCodeThumbnailsFormat>>("thumbnails_format"));
    if (const auto thumbnails_value = config.option<ConfigOptionPoints>("thumbnails"))
        return make_and_check_thumbnail_list(thumbnails_value->get_values(), config.option<ConfigOptionEnum<GCodeThumbnailsFormat>>("thumbnails_format")->value);

    return {};
}

std::string get_error_string(const ThumbnailErrors& errors)
{
    std::string error_str;

    if (errors.has(ThumbnailError::InvalidVal))
        error_str += "\n - " + format("Invalid input format. Expected vector of dimensions in the following format: \"%1%\"", "XxY/EXT, XxY/EXT, ...");
    if (errors.has(ThumbnailError::OutOfRange))
        error_str += "\n - Input value is out of range";
    if (errors.has(ThumbnailError::InvalidExt))
        error_str += "\n - Some extension in the input is invalid";

    return error_str;
}

} // namespace Slic3r::GCodeThumbnails
