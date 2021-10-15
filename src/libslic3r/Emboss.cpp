#include "Emboss.hpp"
#include <stdio.h>
#include <cstdlib>
#include <boost/nowide/convert.hpp>
#include <ClipperUtils.hpp> // union_ex

#define STB_TRUETYPE_IMPLEMENTATION // force following include to generate implementation
#include "imgui/imstb_truetype.h" // stbtt_fontinfo

#include <Triangulation.hpp> // CGAL project

using namespace Slic3r;

// do not expose out of this file stbtt_ data types
class Private
{
public: 
    Private() = delete;

    static std::optional<stbtt_fontinfo> load_font_info(const Emboss::Font &font);
    static std::optional<Emboss::Glyph> get_glyph(stbtt_fontinfo &font_info, int unicode_letter, float flatness = 2.f);
    static std::optional<Emboss::Glyph> get_glyph(int unicode, const Emboss::Font &font, const FontProp &font_prop, 
        Emboss::Glyphs &cache, std::optional<stbtt_fontinfo> &font_info_opt);

    static FontItem create_font_item(std::wstring name, std::wstring path);

    /// <summary>
    /// TODO: move to ExPolygon utils
    /// Remove multi points. When exist multi point dilate it by rect 3x3 and union result.
    /// </summary>
    /// <param name="expolygons">Shape which can contain same point, will be extended by dilatation rects</param>
    /// <returns>ExPolygons with only unique points</returns>
    static ExPolygons dilate_to_unique_points(ExPolygons &expolygons);
};

std::optional<stbtt_fontinfo> Private::load_font_info(const Emboss::Font &font)
{
    int font_offset = stbtt_GetFontOffsetForIndex(font.buffer.data(), font.index);
    if (font_offset < 0) {
        std::cerr << "Font index("<<font.index<<") doesn't exist.";
        return {};        
    }
    stbtt_fontinfo font_info;
    if (stbtt_InitFont(&font_info, font.buffer.data(), font_offset) == 0) {
        std::cerr << "Can't initialize font.";
        return {};
    }
    return font_info;
}

std::optional<Emboss::Glyph> Private::get_glyph(stbtt_fontinfo &font_info, int unicode_letter, float flatness)
{
    int glyph_index = stbtt_FindGlyphIndex(&font_info, unicode_letter);
    if (glyph_index == 0) { 
        std::cerr << "Character codepoint(" << unicode_letter 
            << " = '" << (char) unicode_letter << "') is not defined in the font.";
        return {};
    }

    Emboss::Glyph glyph;
    stbtt_GetGlyphHMetrics(&font_info, glyph_index, &glyph.advance_width, &glyph.left_side_bearing);

    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&font_info, glyph_index, &vertices);
    if (num_verts <= 0) return glyph; // no shape

    int *contour_lengths = NULL;
    int  num_countour_int = 0;

    stbtt__point *points = stbtt_FlattenCurves(vertices, num_verts,
        flatness, &contour_lengths, &num_countour_int, font_info.userdata);

    size_t   num_contour = static_cast<size_t>(num_countour_int);
    Polygons glyph_polygons;
    glyph_polygons.reserve(num_contour);
    size_t pi = 0; // point index
    for (size_t ci = 0; ci < num_contour; ++ci) {
        int length = contour_lengths[ci];
        // check minimal length for triangle
        if (length < 4) {
            // weird font
            pi+=length;
            continue;
        }
        // last point is first point
        --length;
        Points pts;
        pts.reserve(length);
        for (int i = 0; i < length; ++i) {
            const stbtt__point &point = points[pi++];
            pts.emplace_back(point.x, point.y);
        }
        // last point is first point
        assert(pts.front() == Point(points[pi].x, points[pi].y));
        ++pi;

        // change outer cw to ccw and inner ccw to cw order
        std::reverse(pts.begin(), pts.end());
        glyph_polygons.emplace_back(pts);
    }
    // fix for bad defined fonts
    glyph.shape = Slic3r::union_ex(glyph_polygons);
    // inner cw - hole
    // outer ccw - contour
    return glyph;
}

std::optional<Emboss::Glyph> Private::get_glyph(
    int                            unicode,
    const Emboss::Font &           font,
    const FontProp &               font_prop,
    Emboss::Glyphs &               cache,
    std::optional<stbtt_fontinfo> &font_info_opt)
{
    auto glyph_item = cache.find(unicode);
    if (glyph_item != cache.end())
        return glyph_item->second;

    
    if (!font_info_opt.has_value()) {
        font_info_opt = Private::load_font_info(font);
        // can load font info?
        if (!font_info_opt.has_value()) return {};
    }
    std::optional<Emboss::Glyph> glyph_opt =
        Private::get_glyph(*font_info_opt, unicode, font_prop.flatness);

    // IMPROVE: multiple loadig glyph without data
    // has definition inside of font?
    if (!glyph_opt.has_value()) return {};
    cache[unicode] = *glyph_opt;
    return glyph_opt;
}

FontItem Private::create_font_item(std::wstring name, std::wstring path) {
    return FontItem(boost::nowide::narrow(name.c_str()),
                    boost::nowide::narrow(path.c_str()));
}

ExPolygons Private::dilate_to_unique_points(ExPolygons &expolygons)
{   
    std::set<Point> points;
    std::set<Point> multi_points;
    auto find_multipoint = [&points, &multi_points](const Points &pts) {
        for (const Point &p : pts) {
            auto it = points.find(p);
            if (it != points.end())
                multi_points.insert(p);
            else
                points.insert(p);
        }
    };
    for (const ExPolygon &expolygon : expolygons) {
        find_multipoint(expolygon.contour.points);
        for (const Slic3r::Polygon &hole : expolygon.holes)
            find_multipoint(hole.points);
    }
    // speed up, no multipoints
    if (multi_points.empty()) return expolygons;

    // CCW rectangle around zero with size 3*3 px for dilatation
    const Points rect_3_3{Point(1, 1), Point(-1, 1), Point(-1, -1), Point(1, -1)};
    const Points rect_side{Point(1, 0), Point(0, 1), Point(-1, 0), Point(0, -1)};

    // all new added points for reduction
    std::set<Point> rects_points;

    // extends expolygons with dilatation rectangle
    expolygons.reserve(expolygons.size() + multi_points.size());
    for (const Point &multi_point : multi_points) {
        Slic3r::Polygon rect(rect_3_3); // copy points
        rect.translate(multi_point);
        for (const Point p : rect.points) rects_points.insert(p);
        // add side point to be sure with result
        for (const Point p : rect_side) rects_points.insert(p + multi_point);
        expolygons.emplace_back(rect);
    }
    ExPolygons result = union_ex(expolygons);

    // reduce new created close points
    auto reduce_close_points = [&rects_points](Points &pts) {
        bool is_first = false;
        size_t offset = 0;
        bool is_prev_rect = false;
        for (size_t i = 0; i < pts.size(); i++) { 
            const Point &p = pts[i];
            bool is_rect = (rects_points.find(p) != rects_points.end());
            if (is_prev_rect && is_rect) ++offset;
            if (offset != 0) pts[i - offset] = p;
            if (i == 0 && is_rect) is_first = true;
            is_prev_rect = is_rect;
        }
        // remove last
        if (is_first && is_prev_rect) ++offset;
        if (offset != 0)
            pts.erase(pts.begin() + (pts.size() - offset), pts.end());        
    };
    for (ExPolygon &expolygon : result) {
        reduce_close_points(expolygon.contour.points);
        for (Slic3r::Polygon &hole : expolygon.holes)
            reduce_close_points(hole.points);
    }
    return result;
}

#ifdef _WIN32
#include <windows.h>
#include <wingdi.h>
#include <windef.h>
#include <WinUser.h>

// Get system font file path
std::optional<std::wstring> Emboss::get_font_path(const std::wstring &font_face_name)
{
    static const LPWSTR fontRegistryPath = L"Software\\Microsoft\\Windows NT\\CurrentVersion\\Fonts";
    HKEY hKey;
    LONG result;

    // Open Windows font registry key
    result = RegOpenKeyEx(HKEY_LOCAL_MACHINE, fontRegistryPath, 0, KEY_READ, &hKey);
    if (result != ERROR_SUCCESS) return {};    

    DWORD maxValueNameSize, maxValueDataSize;
    result = RegQueryInfoKey(hKey, 0, 0, 0, 0, 0, 0, 0, &maxValueNameSize, &maxValueDataSize, 0, 0);
    if (result != ERROR_SUCCESS) return {};

    DWORD valueIndex = 0;
    LPWSTR valueName = new WCHAR[maxValueNameSize];
    LPBYTE valueData = new BYTE[maxValueDataSize];
    DWORD valueNameSize, valueDataSize, valueType;
    std::wstring wsFontFile;

    // Look for a matching font name
    do {
        wsFontFile.clear();
        valueDataSize = maxValueDataSize;
        valueNameSize = maxValueNameSize;

        result = RegEnumValue(hKey, valueIndex, valueName, &valueNameSize, 0, &valueType, valueData, &valueDataSize);

        valueIndex++;
        if (result != ERROR_SUCCESS || valueType != REG_SZ) {
            continue;
        }

        std::wstring wsValueName(valueName, valueNameSize);

        // Found a match
        if (_wcsnicmp(font_face_name.c_str(), wsValueName.c_str(), font_face_name.length()) == 0) {

            wsFontFile.assign((LPWSTR)valueData, valueDataSize);
            break;
        }
    }while (result != ERROR_NO_MORE_ITEMS);

    delete[] valueName;
    delete[] valueData;

    RegCloseKey(hKey);

    if (wsFontFile.empty()) return {};
    
    // Build full font file path
    WCHAR winDir[MAX_PATH];
    GetWindowsDirectory(winDir, MAX_PATH);

    std::wstringstream ss;
    ss << winDir << "\\Fonts\\" << wsFontFile;
    wsFontFile = ss.str();

    return wsFontFile;
}

FontList Emboss::get_font_list()
{
    //FontList list1 = get_font_list_by_enumeration();
    //FontList list2 = get_font_list_by_register();
    //FontList list3 = get_font_list_by_folder();
    return get_font_list_by_register();
}

FontList Emboss::get_font_list_by_register() {
    static const LPWSTR fontRegistryPath = L"Software\\Microsoft\\Windows NT\\CurrentVersion\\Fonts";
    HKEY hKey;
    LONG result;

    // Open Windows font registry key
    result = RegOpenKeyEx(HKEY_LOCAL_MACHINE, fontRegistryPath, 0, KEY_READ, &hKey);
    if (result != ERROR_SUCCESS) { 
        std::wcerr << L"Can not Open register key (" << fontRegistryPath << ")" 
            << L", function 'RegOpenKeyEx' return code: " << result <<  std::endl;
        return {}; 
    }

    DWORD maxValueNameSize, maxValueDataSize;
    result = RegQueryInfoKey(hKey, 0, 0, 0, 0, 0, 0, 0, &maxValueNameSize,
                             &maxValueDataSize, 0, 0);
    if (result != ERROR_SUCCESS) { 
        std::cerr << "Can not earn query key, function 'RegQueryInfoKey' return code: " 
            << result <<  std::endl;
        return {}; 
    }

    // Build full font file path
    WCHAR winDir[MAX_PATH];
    GetWindowsDirectory(winDir, MAX_PATH);
    std::wstring font_path = std::wstring(winDir) + L"\\Fonts\\";

    FontList font_list;
    DWORD    valueIndex = 0;
    // Look for a matching font name
    LPWSTR font_name = new WCHAR[maxValueNameSize];
    LPBYTE fileTTF_name = new BYTE[maxValueDataSize];
    DWORD  font_name_size, fileTTF_name_size, valueType;
    do {
        fileTTF_name_size = maxValueDataSize;
        font_name_size = maxValueNameSize;

        result = RegEnumValue(hKey, valueIndex, font_name, &font_name_size, 0,
                              &valueType, fileTTF_name, &fileTTF_name_size);
        valueIndex++;
        if (result != ERROR_SUCCESS || valueType != REG_SZ) continue;
        std::wstring font_name_w(font_name, font_name_size);
        std::wstring file_name_w((LPWSTR) fileTTF_name, fileTTF_name_size);
        std::wstring path_w = font_path + file_name_w;

        // filtrate .fon from lists
        size_t pos = font_name_w.rfind(L" (TrueType)");
        if (pos >= font_name_w.size()) continue;
        // remove TrueType text from name
        font_name_w = std::wstring(font_name_w, 0, pos);
        font_list.emplace_back(Private::create_font_item(font_name_w, path_w));
    } while (result != ERROR_NO_MORE_ITEMS);
    delete[] font_name;
    delete[] fileTTF_name;

    RegCloseKey(hKey);
    return font_list;
}

// TODO: Fix global function
bool CALLBACK EnumFamCallBack(LPLOGFONT       lplf,
                              LPNEWTEXTMETRIC lpntm,
                              DWORD           FontType,
                              LPVOID          aFontList)
{
    std::vector<std::wstring> *fontList =
        (std::vector<std::wstring> *) (aFontList);
    if (FontType & TRUETYPE_FONTTYPE) {
        std::wstring name = lplf->lfFaceName;
        fontList->push_back(name);
    }
    return true;
    // UNREFERENCED_PARAMETER(lplf);
    UNREFERENCED_PARAMETER(lpntm);
}

FontList Emboss::get_font_list_by_enumeration() {   

    HDC                       hDC = GetDC(NULL);
    std::vector<std::wstring> font_names;
    EnumFontFamilies(hDC, (LPCTSTR) NULL, (FONTENUMPROC) EnumFamCallBack,
                     (LPARAM) &font_names);

    FontList font_list;
    for (const std::wstring &font_name : font_names) {
        font_list.emplace_back(Private::create_font_item(font_name, L""));
    }    
    return font_list;
}

FontList Emboss::get_font_list_by_folder() {
    FontList result;
    WCHAR winDir[MAX_PATH];
    UINT winDir_size = GetWindowsDirectory(winDir, MAX_PATH);
    std::wstring search_dir = std::wstring(winDir, winDir_size) + L"\\Fonts\\";
    WIN32_FIND_DATA fd;
    HANDLE          hFind;
    // By https://en.wikipedia.org/wiki/TrueType has also suffix .tte
    std::vector<std::wstring> suffixes = {L"*.ttf", L"*.ttc", L"*.tte"};
    for (const std::wstring &suffix : suffixes) {
        hFind = ::FindFirstFile((search_dir + suffix).c_str(), &fd);
        if (hFind == INVALID_HANDLE_VALUE) continue;
        do {
            // skip folder . and ..
            if (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) continue;
            std::wstring file_name(fd.cFileName);
            // TODO: find font name instead of filename
            result.emplace_back(Private::create_font_item(file_name, search_dir + file_name));
        } while (::FindNextFile(hFind, &fd));
        ::FindClose(hFind);
    }
    return result;
}

#else
FontList Emboss::get_font_list() { 
    // not implemented
    return {}; 
}

std::optional<std::wstring> Emboss::get_font_path(const std::wstring &font_face_name){
    // not implemented
    return {};
}
#endif

std::optional<Emboss::Font> Emboss::load_font(std::vector<unsigned char> data)
{
    Font res;
    res.buffer                = std::move(data);

    unsigned int index       = 0;
    int          font_offset = 0;
    while (font_offset >= 0) {
        font_offset = stbtt_GetFontOffsetForIndex(res.buffer.data(), index++);
    }
    --index; // last one is bad
    // at least one font must be inside collection
    if (index < 1) {
        std::cerr << "There is no font collection inside data.";
        return {};
    }
    // select default font on index 0
    res.index = 0;
    res.count = index;

    auto font_info = Private::load_font_info(res);
    if (!font_info.has_value()) return {};
    const stbtt_fontinfo *info = &(*font_info);
    // load information about line gap
    stbtt_GetFontVMetrics(info, &res.ascent, &res.descent, &res.linegap);

    return res;
}

std::optional<Emboss::Font> Emboss::load_font(const char *file_path)
{
    FILE *file = fopen(file_path, "rb");
    if (file == nullptr) {
        std::cerr << "Couldn't open " << file_path << " for reading.";
        return {};
    }

    // find size of file
    if (fseek(file, 0L, SEEK_END) != 0) {
        std::cerr << "Couldn't fseek file " << file_path << " for size measure.";
        return {};
    }
    size_t size = ftell(file);
    if (size == 0) {
        std::cerr << "Size of font file is zero. Can't read.";
        return {};    
    }
    rewind(file);

    std::vector<unsigned char> buffer(size);
    size_t count_loaded_bytes = fread((void *) &buffer.front(), 1, size, file);
    if (count_loaded_bytes != size) {
        std::cerr << "Different loaded(from file) data size.";
        return {};
    }
    return load_font(std::move(buffer));
}


#ifdef _WIN32
std::optional<Emboss::Font> Emboss::load_font(HFONT hfont) 
{
    HDC hdc = ::CreateCompatibleDC(NULL);
    if (hdc == NULL) {
        std::cerr << "Can't create HDC by CreateCompatibleDC(NULL).";
        return {};
    }

    // To retrieve the data from the beginning of the file for TrueType
    // Collection files specify 'ttcf' (0x66637474).
    DWORD dwTable = 0x66637474;
    DWORD dwOffset = 0;

    ::SelectObject(hdc, hfont);
    size_t size = ::GetFontData(hdc, dwTable, dwOffset, NULL, 0);
    if (size == GDI_ERROR) { 
        // HFONT is NOT TTC(collection)
        dwTable = 0;
        size = ::GetFontData(hdc, dwTable, dwOffset, NULL, 0);
    }

    if (size == 0 || size == GDI_ERROR) {
        std::cerr << "HFONT doesn't have size.";
        ::DeleteDC(hdc);
        return {};    
    }

    std::vector<unsigned char> buffer(size);
    size_t loaded_size = ::GetFontData(hdc, dwTable, dwOffset, buffer.data(), size);
    ::DeleteDC(hdc);
    if (size != loaded_size) {
        std::cerr << "Different loaded(from HFONT) data size.";
        return {};    
    }

    return load_font(std::move(buffer));
}
#endif // _WIN32

std::optional<Emboss::Glyph> Emboss::letter2glyph(const Font &font,
                                                  int         letter,
                                                  float       flatness)
{
    auto font_info_opt = Private::load_font_info(font);
    if (!font_info_opt.has_value()) return {};
    return Private::get_glyph(*font_info_opt, (int) letter, flatness);
}

ExPolygons Emboss::text2shapes(Font &          font,
                               const char *    text,
                               const FontProp &font_prop)
{
    std::optional<stbtt_fontinfo> font_info_opt;
    
    Point    cursor(0, 0);
    ExPolygons result;

    std::wstring ws = boost::nowide::widen(text);
    for (wchar_t wc: ws){
        if (wc == '\n') { 
            cursor.x() = 0;
            cursor.y() -= font.ascent - font.descent + font.linegap + font_prop.line_gap;
            continue;
        } 
        if (wc == '\t') {
            // '\t' = 4*space => same as imgui
            const int count_spaces = 4;
            std::optional<Glyph> space_opt = Private::get_glyph(int(' '), font, font_prop, font.cache, font_info_opt);
            if (!space_opt.has_value()) continue;
            cursor.x() += count_spaces *(space_opt->advance_width + font_prop.char_gap);
            continue;
        }

        int unicode = static_cast<int>(wc);
        std::optional<Glyph> glyph_opt = Private::get_glyph(unicode, font, font_prop, font.cache, font_info_opt);
        if (!glyph_opt.has_value()) continue;
        
        // move glyph to cursor position
        ExPolygons expolygons = glyph_opt->shape; // copy
        for (ExPolygon &expolygon : expolygons) 
            expolygon.translate(cursor);
        cursor.x() += glyph_opt->advance_width + font_prop.char_gap;
        expolygons_append(result, expolygons);
    }
    result = Slic3r::union_ex(result);
    return Private::dilate_to_unique_points(result);
}

indexed_triangle_set Emboss::polygons2model(const ExPolygons &shape2d,
                                            const IProject &projection)
{
    indexed_triangle_set result;
    size_t count_point = count_points(shape2d);
    result.vertices.reserve(2 * count_point);

    std::vector<Vec3f> &front_points = result.vertices;
    std::vector<Vec3f>  back_points;
    back_points.reserve(count_point);

    auto insert_point = [&projection, &front_points,
                         &back_points](const Polygon& polygon) {
        for (const Point& p : polygon.points) {
            auto p2 = projection.project(p);
            front_points.emplace_back(p2.first);
            back_points.emplace_back(p2.second);
        }
    };
    for (const ExPolygon &expolygon : shape2d) {
        insert_point(expolygon.contour);
        for (const Polygon &hole : expolygon.holes) insert_point(hole);
    }
    // insert back points, front are already in
    result.vertices.insert(result.vertices.end(),
                           std::make_move_iterator(back_points.begin()),
                           std::make_move_iterator(back_points.end()));

    // CW order of triangle indices
    std::vector<Vec3i> shape_triangles = Triangulation::triangulate(shape2d);
    result.indices.reserve(shape_triangles.size() * 2 + count_point * 2);
    // top triangles - change to CCW
    for (const Vec3i &t : shape_triangles)
        result.indices.emplace_back(t.x(), t.z(), t.y());
    // bottom triangles - use CW
    for (const Vec3i &t : shape_triangles)
        result.indices.emplace_back(t.x() + count_point, t.y() + count_point,
                                    t.z() + count_point);

    // quads around - zig zag by triangles
    size_t polygon_offset = 0;
    auto add_quads = [&result,&polygon_offset, count_point](const Polygon& polygon) {
        uint32_t polygon_points = polygon.points.size();
        for (uint32_t p = 0; p < polygon_points; p++) { 
            uint32_t i = polygon_offset + p;
            // previous index
            uint32_t ip = (p == 0) ? (polygon_offset + polygon_points - 1) : (i - 1);
            // bottom indices
            uint32_t i2  = i + count_point;
            uint32_t ip2 = ip + count_point;

            result.indices.emplace_back(i, i2, ip);
            result.indices.emplace_back(ip2, ip, i2);
        }
        polygon_offset += polygon_points;
    };
    
    for (const ExPolygon &expolygon : shape2d) {
        add_quads(expolygon.contour);
        for (const Polygon &hole : expolygon.holes) add_quads(hole);
    }
    return result;
}

std::pair<Vec3f, Vec3f> Emboss::ProjectZ::project(const Point &p) const
{
    Vec3f front(p.x(),p.y(),0.f);
    Vec3f back = front; // copy
    back.z() = m_depth;
    return std::make_pair(front, back);
}
