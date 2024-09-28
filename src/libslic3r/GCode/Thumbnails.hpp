///|/ Copyright (c) Prusa Research 2022 Lukáš Hejl @hejllukas, Vojtěch Bubník @bubnikv
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
#ifndef slic3r_GCodeThumbnails_hpp_
#define slic3r_GCodeThumbnails_hpp_

#include "../Point.hpp"
#include "../PrintConfig.hpp"
#include "ThumbnailData.hpp"

#include <vector>
#include <memory>
#include <string_view>

#include <LibBGCode/binarize/binarize.hpp>

#include <boost/beast/core/detail/base64.hpp>

#include "../libslic3r/enum_bitmask.hpp"

namespace Slic3r {
    enum class ThumbnailError : int { InvalidVal, OutOfRange, InvalidExt };
    using ThumbnailErrors = enum_bitmask<ThumbnailError>;
    ENABLE_ENUM_BITMASK_OPERATORS(ThumbnailError);
}

namespace Slic3r::GCodeThumbnails {

constexpr std::string_view EMPTY_TAG = "thumbnail";

struct CompressedImageBuffer
{
    void       *data { nullptr };
    size_t      size { 0 };
    virtual ~CompressedImageBuffer() {}
    virtual std::string_view tag() const = 0;
};

std::unique_ptr<CompressedImageBuffer> compress_thumbnail(const ThumbnailData &data, GCodeThumbnailsFormat format);

typedef std::vector<std::pair<GCodeThumbnailsFormat, Vec2d>> GCodeThumbnailDefinitionsList;

using namespace std::literals;
std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list_from_prusa(const std::string& thumbnails_string, const std::string_view def_ext = "PNG"sv);
std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list(const std::vector<Vec2d>& thumbnails, GCodeThumbnailsFormat format);
std::pair<GCodeThumbnailDefinitionsList, ThumbnailErrors> make_and_check_thumbnail_list(const ConfigBase &config);

std::string get_error_string(const ThumbnailErrors& errors);

template<typename WriteToOutput, typename ThrowIfCanceledCallback>
inline void export_thumbnails_to_file(ThumbnailsGeneratorCallback &thumbnail_cb,
                                      const std::vector<Vec2d> &   sizes,
                                      bool                         with_bed,
                                      GCodeThumbnailsFormat        format,
                                      bool                         with_tag_format,
                                      WriteToOutput                output,
                                      ThrowIfCanceledCallback      throw_if_canceled)
{
    // Write thumbnails using base64 encoding
    if (thumbnail_cb != nullptr) {
        // 0-size is the same as no sizes.
        std::vector<Vec2d> good_sizes;
        for (const Vec2d& size : sizes)
            if (size.x() > 0 && size.y() > 0)
                good_sizes.push_back(size);
        if (good_sizes.empty()) return;

        //Create the thumbnails
        static constexpr const size_t max_row_length = 78;
        ThumbnailsList thumbnails = thumbnail_cb(ThumbnailsParams{ good_sizes, true, true, with_bed, true });
        for (const ThumbnailData& data : thumbnails)
            if (data.is_valid()) {
                auto compressed = compress_thumbnail(data, format);
                if (compressed->data && compressed->size) {
                    if (format == GCodeThumbnailsFormat::BIQU) {
                        // BIQU firmware need to have nothing before the thumbnail
                        //output((boost::format("\n;\n; %s begin %dx%d %d\n") 
                        //    % (with_tag_format ? compressed->tag() : EMPTY_TAG)
                        //    % data.width % data.height % (compressed->size - 1)).str().c_str());
                        //print size in hex
                        std::stringstream ss;
                        ss << std::setfill('0') << std::hex;
                        //biqu header
                        ss << ";" << std::setw(4) << data.width << std::setw(4) << data.height << "\n";
                        output(ss.str().c_str());
                        if (((char*)compressed->data)[compressed->size -1] == '\0')
                            output((char*)(compressed->data));
                        else
                            assert(false);
                        output("; bigtree thumbnail end\n");
                    } else {
                        std::string encoded;
                        encoded.resize(boost::beast::detail::base64::encoded_size(compressed->size));
                        encoded.resize(boost::beast::detail::base64::encode((void*)encoded.data(), (const void*)compressed->data, compressed->size));

                        output((boost::format("\n;\n; %s begin %dx%d %d\n") 
                            % (with_tag_format ? compressed->tag() : EMPTY_TAG)
                            % data.width % data.height % encoded.size()).str().c_str());
                        while (encoded.size() > max_row_length) {
                            output((boost::format("; %s\n") % encoded.substr(0, max_row_length)).str().c_str());
                            encoded = encoded.substr(max_row_length);
                        }

                        if (encoded.size() > 0)
                            output((boost::format("; %s\n") % encoded).str().c_str());
                    }
                    output((boost::format("; %s end\n;\n") % (with_tag_format ? compressed->tag() : EMPTY_TAG)).str().c_str());
                }
                throw_if_canceled();
            }
    }
}

template<typename ThrowIfCanceledCallback>
inline void generate_binary_thumbnails(ThumbnailsGeneratorCallback& thumbnail_cb, std::vector<bgcode::binarize::ThumbnailBlock>& out_thumbnails,
    const std::vector<std::pair<GCodeThumbnailsFormat, Vec2d>> &thumbnails_list, ThrowIfCanceledCallback throw_if_canceled)
{
    using namespace bgcode::core;
    using namespace bgcode::binarize;
    out_thumbnails.clear();
    assert(thumbnail_cb != nullptr);
    if (thumbnail_cb != nullptr) {
        for (const auto& [format, size] : thumbnails_list) {
            ThumbnailsList thumbnails = thumbnail_cb(ThumbnailsParams{ {size}, true, true, true, true });
            for (const ThumbnailData &data : thumbnails)
                if (data.is_valid()) {
                    auto compressed = compress_thumbnail(data, format);
                    if (compressed->data != nullptr && compressed->size > 0) {
                        ThumbnailBlock& block = out_thumbnails.emplace_back(ThumbnailBlock());
                        block.params.width = (uint16_t)data.width;
                        block.params.height = (uint16_t)data.height;
                        switch (format) {
                        case GCodeThumbnailsFormat::PNG: { block.params.format = (uint16_t)EThumbnailFormat::PNG; break; }
                        case GCodeThumbnailsFormat::JPG: { block.params.format = (uint16_t)EThumbnailFormat::JPG; break; }
                        case GCodeThumbnailsFormat::QOI: { block.params.format = (uint16_t)EThumbnailFormat::QOI; break; }
                        }
                        block.data.resize(compressed->size);
                        memcpy(block.data.data(), compressed->data, compressed->size);
                    }
                }
        }
    }
}

} // namespace Slic3r::GCodeThumbnails

#endif // slic3r_GCodeThumbnails_hpp_
