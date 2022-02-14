#ifndef SLARASTER_CPP
#define SLARASTER_CPP

#include <functional>

#include <libslic3r/SLA/RasterBase.hpp>
#include <libslic3r/SLA/AGGRaster.hpp>

// minz image write:
#include <miniz.h>

namespace Slic3r { namespace sla {

const RasterBase::TMirroring RasterBase::NoMirror = {false, false};
const RasterBase::TMirroring RasterBase::MirrorX  = {true, false};
const RasterBase::TMirroring RasterBase::MirrorY  = {false, true};
const RasterBase::TMirroring RasterBase::MirrorXY = {true, true};



static void pwx_get_pixel_span(const std::uint8_t* ptr, const std::uint8_t* end,
                         std::uint8_t& pixel, size_t& span_len)
{
    size_t max_len;

    span_len = 0;
    pixel = (*ptr) & 0xF0;
    // the maximum length of the span depends on the pixel color
    max_len = (pixel == 0 || pixel == 0xF0) ? 0xFFF : 0xF;
    while (((*ptr) & 0xF0) == pixel && ptr < end && span_len < max_len) {
        span_len++;
        ptr++;
    }
}

EncodedRaster PWXRasterEncoder::operator()(const void *ptr, size_t w, size_t h,
                                           size_t      num_components)
{
    std::vector<uint8_t> dst;
    size_t span_len;
    std::uint8_t pixel;    
    auto size = w * h * num_components;   
    dst.reserve(size);

    const std::uint8_t* src = reinterpret_cast<const std::uint8_t*>(ptr);
    const std::uint8_t* src_end = src + size;
    while (src < src_end) {
        pwx_get_pixel_span(src, src_end, pixel, span_len);
        src += span_len;
        // fully transparent of fully opaque pixel
        if (pixel == 0 || pixel == 0xF0) {
            pixel = pixel  | (span_len >> 8);
            std::copy(&pixel, (&pixel) + 1, std::back_inserter(dst));
            pixel = span_len & 0xFF;
            std::copy(&pixel, (&pixel) + 1, std::back_inserter(dst));
        } 
        // antialiased pixel
        else {
            pixel = pixel | span_len;
            std::copy(&pixel, (&pixel) + 1, std::back_inserter(dst));
        }
    }
  
    return EncodedRaster(std::move(dst), "pwx");
}

EncodedRaster PNGRasterEncoder::operator()(const void *ptr, size_t w, size_t h,
                                           size_t      num_components)
{
    std::vector<uint8_t> buf;
    size_t s = 0;
    
    void *rawdata = tdefl_write_image_to_png_file_in_memory(
        ptr, int(w), int(h), int(num_components), &s);
    
    // On error, data() will return an empty vector. No other info can be
    // retrieved from miniz anyway...
    if (rawdata == nullptr) return EncodedRaster({}, "png");
    
    auto pptr = static_cast<std::uint8_t*>(rawdata);
    
    buf.reserve(s);
    std::copy(pptr, pptr + s, std::back_inserter(buf));
    
    MZ_FREE(rawdata);
    return EncodedRaster(std::move(buf), "png");
}

std::ostream &operator<<(std::ostream &stream, const EncodedRaster &bytes)
{
    stream.write(reinterpret_cast<const char *>(bytes.data()),
                 std::streamsize(bytes.size()));
    
    return stream;
}

EncodedRaster PPMRasterEncoder::operator()(const void *ptr, size_t w, size_t h,
                                           size_t      num_components)
{
    std::vector<uint8_t> buf;
    
    auto header = std::string("P5 ") +
            std::to_string(w) + " " +
            std::to_string(h) + " " + "255 ";
    
    auto sz = w * h * num_components;
    size_t s = sz + header.size();
    
    buf.reserve(s);

    auto buff = reinterpret_cast<const std::uint8_t*>(ptr);
    std::copy(header.begin(), header.end(), std::back_inserter(buf));
    std::copy(buff, buff+sz, std::back_inserter(buf));
    
    return EncodedRaster(std::move(buf), "ppm");
}

std::unique_ptr<RasterBase> create_raster_grayscale_aa(
    const Resolution        &res,
    const PixelDim          &pxdim,
    double                   gamma,
    const RasterBase::Trafo &tr)
{
    std::unique_ptr<RasterBase> rst;
    
    if (gamma > 0)
        rst = std::make_unique<RasterGrayscaleAAGammaPower>(res, pxdim, tr, gamma);
    else if (std::abs(gamma - 1.) < 1e-6)
        rst = std::make_unique<RasterGrayscaleAA>(res, pxdim, tr, agg::gamma_none());
    else
        rst = std::make_unique<RasterGrayscaleAA>(res, pxdim, tr, agg::gamma_threshold(.5));
    
    return rst;
}

} // namespace sla
} // namespace Slic3r

#endif // SLARASTER_CPP
