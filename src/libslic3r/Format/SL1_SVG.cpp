#include "SL1_SVG.hpp"
#include "SLA/RasterBase.hpp"
#include "libslic3r/LocalesUtils.hpp"

namespace Slic3r {

namespace {

void transform(ExPolygon &ep, const sla::RasterBase::Trafo &tr, const BoundingBox &bb)
{
    if (tr.flipXY) {
        for (auto &p : ep.contour.points) std::swap(p.x(), p.y());
        for (auto &h : ep.holes)
            for (auto &p : h.points) std::swap(p.x(), p.y());
    }

    if (tr.mirror_x){
        for (auto &p : ep.contour.points) p.x() = bb.max.x() - p.x() + bb.min.x();
        for (auto &h : ep.holes)
            for (auto &p : h.points) p.x() = bb.max.x() - p.x() + bb.min.x();
    }

    if (tr.mirror_y){
        for (auto &p : ep.contour.points) p.y() = bb.max.y() - p.y() + bb.min.y();
        for (auto &h : ep.holes)
            for (auto &p : h.points) p.y() = bb.max.y() - p.y() + bb.min.y();
    }
}

void append_svg(std::string &buf, const Polygon &poly)
{
    buf += "<path d=\"M";
    for (auto &p : poly) {
        buf += " ";
        buf += float_to_string_decimal_point(unscaled<float>(p.x()));
        buf += " ";
        buf += float_to_string_decimal_point(unscaled<float>(p.y()));
    }
    buf += " z\""; // mark path as closed
    buf += " />\n";
}

} // namespace

// A fake raster from SVG
class SVGRaster : public sla::RasterBase {
    // Resolution here will be used for svg boundaries
    BoundingBox m_bb;
    Trafo       m_trafo;

    std::string m_svg;

public:
    SVGRaster(Resolution res = {}, Trafo tr = {})
        : m_bb{BoundingBox{{0, 0}, Vec2crd{res.width_px, res.height_px}}}
        , m_trafo{tr}
    {
        std::string w = float_to_string_decimal_point(unscaled<float>(res.width_px));
        std::string h = float_to_string_decimal_point(unscaled<float>(res.height_px));
        // Add svg header.
        m_svg =
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n"
            "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
            "<svg height=\"" + h + "mm" + "\" width=\"" + w + "mm" + "\" viewBox=\"0 0 " + w + " " + h +
            "\" style=\"fill: white; stroke: none; fill-rule: nonzero\" "
            "xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

        // Add black background;
        m_svg += "<rect fill='black' stroke='none' x='0' y='0' width='" + w + "' height='" + h + "'/>\n";
    }

    void draw(const ExPolygon& poly) override
    {
        auto cpoly = poly;

        transform(cpoly, m_trafo, m_bb);
        append_svg(m_svg, cpoly.contour);
        for (auto &h : cpoly.holes)
            append_svg(m_svg, h);
    }

    Resolution resolution() const override
    {
        return {size_t(m_bb.size().x()), size_t(m_bb.size().y())};
    }

    // Pixel dimension is undefined in this case.
    PixelDim pixel_dimensions() const override { return {0, 0}; }

    Trafo      trafo() const override { return m_trafo; }

    sla::EncodedRaster encode(sla::RasterEncoder /*encoder*/) const override
    {
        std::vector<uint8_t> data;
        constexpr const char finish[] = "</svg>\n";

        data.reserve(m_svg.size() + std::size(finish));

        std::copy(m_svg.begin(), m_svg.end(), std::back_inserter(data));
        std::copy(finish, finish + std::size(finish) - 1, std::back_inserter(data));

        return sla::EncodedRaster{std::move(data), "svg"};
    }
};

std::unique_ptr<sla::RasterBase> SL1_SVGArchive::create_raster() const
{
    auto w = scaled<size_t>(cfg().display_width.getFloat());
    auto h = scaled<size_t>(cfg().display_height.getFloat());

    std::array<bool, 2>         mirror;

    mirror[X] = cfg().display_mirror_x.getBool();
    mirror[Y] = cfg().display_mirror_y.getBool();

    auto ro = cfg().display_orientation.getInt();
    sla::RasterBase::Orientation orientation =
        ro == sla::RasterBase::roPortrait ? sla::RasterBase::roPortrait :
                                            sla::RasterBase::roLandscape;

    if (orientation == sla::RasterBase::roPortrait) {
        std::swap(w, h);
    }

    sla::RasterBase::Resolution res{w, h};
    sla::RasterBase::Trafo tr{orientation, mirror};

    // Gamma does not really make sense in an svg, right?
    // double gamma = cfg().gamma_correction.getFloat();

    return std::make_unique<SVGRaster>(SVGRaster::Resolution{w, h}, tr);
}

sla::RasterEncoder SL1_SVGArchive::get_encoder() const
{
    return nullptr;
}

} // namespace Slic3r
