#include "Emboss.hpp"
#include <stdio.h>
#include <cstdlib>
#include <boost/nowide/convert.hpp>
#include <ClipperUtils.hpp> // union_ex

#define STB_TRUETYPE_IMPLEMENTATION // force following include to generate implementation
#include "imgui/imstb_truetype.h" // stbtt_fontinfo
#include "Utils.hpp" // ScopeGuard

#include <Triangulation.hpp> // CGAL project
#include "libslic3r.h"

#include "ClipperUtils.hpp" // for boldness - polygon extend(offset)

using namespace Slic3r;

double Emboss::SHAPE_SCALE = 0.001;//SCALING_FACTOR;

// do not expose out of this file stbtt_ data types
class Private
{
public: 
    Private() = delete;
    static bool is_valid(const Emboss::FontFile &font, unsigned int index);
    static std::optional<stbtt_fontinfo> load_font_info(const unsigned char *data, unsigned int index = 0);
    static std::optional<Emboss::Glyph> get_glyph(const stbtt_fontinfo &font_info, int unicode_letter, float flatness);

    // take glyph from cache
    static const Emboss::Glyph* get_glyph(int unicode, const Emboss::FontFile &font, const FontProp &font_prop, 
        Emboss::Glyphs &cache, std::optional<stbtt_fontinfo> &font_info_opt);

    static FontItem create_font_item(std::wstring name, std::wstring path);

    /// <summary>
    /// TODO: move to ExPolygon utils
    /// Remove multi points. When exist multi point dilate it by rect 3x3 and union result.
    /// </summary>
    /// <param name="expolygons">Shape which can contain same point, will be extended by dilatation rects</param>
    /// <returns>ExPolygons with only unique points</returns>
    static ExPolygons dilate_to_unique_points(ExPolygons &expolygons);

    // scale and convert float to int coordinate
    static Point to_point(const stbtt__point &point);
};

bool Private::is_valid(const Emboss::FontFile &font, unsigned int index) {
    if (font.data == nullptr) return false;
    if (font.data->empty()) return false;
    if (index >= font.infos.size()) return false;
    return true;
}

std::optional<stbtt_fontinfo> Private::load_font_info(
    const unsigned char *data, unsigned int index)
{
    int font_offset = stbtt_GetFontOffsetForIndex(data, index);
    if (font_offset < 0) {
        assert(false);
        // "Font index(" << index << ") doesn't exist.";
        return {};        
    }
    stbtt_fontinfo font_info;
    if (stbtt_InitFont(&font_info, data, font_offset) == 0) {
        // Can't initialize font.
        assert(false);
        return {};
    }
    return font_info;
}

std::optional<Emboss::Glyph> Private::get_glyph(const stbtt_fontinfo &font_info, int unicode_letter, float flatness)
{
    int glyph_index = stbtt_FindGlyphIndex(&font_info, unicode_letter);
    if (glyph_index == 0) {
        //wchar_t wchar = static_cast<wchar_t>(unicode_letter); 
        //<< "Character unicode letter ("
        //<< "decimal value = " << std::dec << unicode_letter << ", "
        //<< "hexadecimal value = U+" << std::hex << unicode_letter << std::dec << ", "
        //<< "wchar value = " << wchar
        //<< ") is NOT defined inside of the font. \n";
        return {};
    }

    Emboss::Glyph glyph;
    stbtt_GetGlyphHMetrics(&font_info, glyph_index, &glyph.advance_width, &glyph.left_side_bearing);

    stbtt_vertex *vertices;
    int num_verts = stbtt_GetGlyphShape(&font_info, glyph_index, &vertices);
    if (num_verts <= 0) return glyph; // no shape
    ScopeGuard sg1([&vertices]() { free(vertices); });

    int *contour_lengths = NULL;
    int  num_countour_int = 0;
    stbtt__point *points = stbtt_FlattenCurves(vertices, num_verts,
        flatness, &contour_lengths, &num_countour_int, font_info.userdata);
    if (!points) return glyph; // no valid flattening
    ScopeGuard sg2([&contour_lengths, &points]() {
        free(contour_lengths); 
        free(points); 
    });

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
        for (int i = 0; i < length; ++i) 
            pts.emplace_back(to_point(points[pi++]));
        
        // last point is first point --> closed contour
        assert(pts.front() == to_point(points[pi]));
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

const Emboss::Glyph* Private::get_glyph(
    int                            unicode,
    const Emboss::FontFile &       font,
    const FontProp &               font_prop,
    Emboss::Glyphs &               cache,
    std::optional<stbtt_fontinfo> &font_info_opt)
{
    const double RESOLUTION = 0.0125; // TODO: read from printer configuration
    auto glyph_item = cache.find(unicode);
    if (glyph_item != cache.end()) return &glyph_item->second;

    unsigned int font_index = font_prop.collection_number.has_value()?
            *font_prop.collection_number : 0;
    if (!is_valid(font, font_index)) return nullptr;

    if (!font_info_opt.has_value()) {
        
        font_info_opt  = Private::load_font_info(font.data->data(), font_index);
        // can load font info?
        if (!font_info_opt.has_value()) return nullptr;
    }

    float flatness = static_cast<float>(font.infos[font_index].ascent * RESOLUTION / font_prop.size_in_mm);
    std::optional<Emboss::Glyph> glyph_opt =
        Private::get_glyph(*font_info_opt, unicode, flatness);

    // IMPROVE: multiple loadig glyph without data
    // has definition inside of font?
    if (!glyph_opt.has_value()) return nullptr;

    if (font_prop.char_gap.has_value()) 
        glyph_opt->advance_width += *font_prop.char_gap;

    // scale glyph size
    glyph_opt->advance_width = 
        static_cast<int>(glyph_opt->advance_width / Emboss::SHAPE_SCALE);
    glyph_opt->left_side_bearing = 
        static_cast<int>(glyph_opt->left_side_bearing / Emboss::SHAPE_SCALE);

    if (font_prop.boldness.has_value()) {
        float delta = *font_prop.boldness / Emboss::SHAPE_SCALE / font_prop.size_in_mm;
        glyph_opt->shape = offset_ex(glyph_opt->shape, delta);
    }

    if (font_prop.skew.has_value()) {
        const float &ratio = *font_prop.skew;
        auto skew = [&ratio](Slic3r::Polygon &polygon) {
            for (Slic3r::Point &p : polygon.points) { p.x() += p.y() * ratio; }
        };
        for (ExPolygon &expolygon : glyph_opt->shape) {
            skew(expolygon.contour);
            for (Slic3r::Polygon &hole : expolygon.holes) skew(hole);
        }
    }

    // union of shape
    // (for sure) I do not believe in font corectness
    // modification like bold or skew could create artefacts
    glyph_opt->shape = Slic3r::union_ex(glyph_opt->shape);
    // unify multipoints with similar position. Could appear after union
    dilate_to_unique_points(glyph_opt->shape); 
    auto it = cache.insert({unicode, std::move(*glyph_opt)});
    assert(it.second);
    return &it.first->second;
}

FontItem Private::create_font_item(std::wstring name, std::wstring path) {
    return { boost::nowide::narrow(name.c_str()),
             boost::nowide::narrow(path.c_str()),
             FontItem::Type::file_path, FontProp() };
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
        for (const Point& p : rect.points) rects_points.insert(p);
        // add side point to be sure with result
        for (const Point& p : rect_side) rects_points.insert(p + multi_point);
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

Point Private::to_point(const stbtt__point &point) {
    return Point(static_cast<int>(std::round(point.x / Emboss::SHAPE_SCALE)),
                 static_cast<int>(std::round(point.y / Emboss::SHAPE_SCALE)));
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
        assert(false);
        //std::wcerr << L"Can not Open register key (" << fontRegistryPath << ")" 
        //    << L", function 'RegOpenKeyEx' return code: " << result <<  std::endl;
        return {}; 
    }

    DWORD maxValueNameSize, maxValueDataSize;
    result = RegQueryInfoKey(hKey, 0, 0, 0, 0, 0, 0, 0, &maxValueNameSize,
                             &maxValueDataSize, 0, 0);
    if (result != ERROR_SUCCESS) {
        assert(false);
        // Can not earn query key, function 'RegQueryInfoKey' return code: result
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

std::unique_ptr<Emboss::FontFile> Emboss::create_font_file(
    std::unique_ptr<std::vector<unsigned char>> data)
{
    int collection_size = stbtt_GetNumberOfFonts(data->data());
    // at least one font must be inside collection
    if (collection_size < 1) {
        assert(false);
        // There is no font collection inside font data
        return nullptr;
    }

    unsigned int c_size = static_cast<unsigned int>(collection_size);
    std::vector<FontFile::Info> infos;
    infos.reserve(c_size);
    for (unsigned int i = 0; i < c_size; ++i) {
        auto font_info = Private::load_font_info(data->data(), i);
        if (!font_info.has_value()) return nullptr;

        const stbtt_fontinfo *info = &(*font_info);
        // load information about line gap
        int ascent, descent, linegap;
        stbtt_GetFontVMetrics(info, &ascent, &descent, &linegap);

        float pixels       = 1000.; // value is irelevant
        float em_pixels    = stbtt_ScaleForMappingEmToPixels(info, pixels);
        int   units_per_em = static_cast<int>(std::round(pixels / em_pixels));

        infos.emplace_back(FontFile::Info{ascent, descent, linegap, units_per_em});
    }
    return std::make_unique<Emboss::FontFile>(std::move(data), std::move(infos));
}

std::unique_ptr<Emboss::FontFile> Emboss::create_font_file(const char *file_path)
{
    FILE *file = fopen(file_path, "rb");
    if (file == nullptr) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Couldn't open " << file_path << " for reading." << std::endl;
        return nullptr;
    }

    // find size of file
    if (fseek(file, 0L, SEEK_END) != 0) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Couldn't fseek file " << file_path << " for size measure." << std::endl;
        return nullptr;
    }
    size_t size = ftell(file);
    if (size == 0) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Size of font file is zero. Can't read." << std::endl;
        return nullptr;    
    }
    rewind(file);
    auto buffer = std::make_unique<std::vector<unsigned char>>(size);
    size_t count_loaded_bytes = fread((void *) &buffer->front(), 1, size, file);
    if (count_loaded_bytes != size) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Different loaded(from file) data size." << std::endl;
        return nullptr;
    }
    return create_font_file(std::move(buffer));
}


#ifdef _WIN32
static bool load_hfont(void* hfont, DWORD &dwTable, DWORD &dwOffset, size_t& size, HDC hdc = nullptr){
    bool del_hdc = false;
    if (hdc == nullptr) { 
        del_hdc = true;
        hdc = ::CreateCompatibleDC(NULL);
        if (hdc == NULL) return false;
    }
    
    // To retrieve the data from the beginning of the file for TrueType
    // Collection files specify 'ttcf' (0x66637474).
    dwTable  = 0x66637474;
    dwOffset = 0;

    ::SelectObject(hdc, hfont);
    size = ::GetFontData(hdc, dwTable, dwOffset, NULL, 0);
    if (size == GDI_ERROR) {
        // HFONT is NOT TTC(collection)
        dwTable = 0;
        size    = ::GetFontData(hdc, dwTable, dwOffset, NULL, 0);
    }

    if (size == 0 || size == GDI_ERROR) {
        if (del_hdc) ::DeleteDC(hdc);
        return false;
    }
    return true;
}

void * Emboss::can_load(HFONT hfont)
{
    DWORD dwTable=0, dwOffset=0;
    size_t size = 0;
    if (!load_hfont(hfont, dwTable, dwOffset, size)) return nullptr;
    return hfont;
}

std::unique_ptr<Emboss::FontFile> Emboss::create_font_file(HFONT hfont)
{
    HDC hdc = ::CreateCompatibleDC(NULL);
    if (hdc == NULL) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Can't create HDC by CreateCompatibleDC(NULL)." << std::endl;
        return nullptr;
    }

    DWORD dwTable=0,dwOffset = 0;
    size_t size;
    if (!load_hfont(hfont, dwTable, dwOffset, size, hdc)) {
        ::DeleteDC(hdc);
        return nullptr;
    }
    auto buffer = std::make_unique<std::vector<unsigned char>>(size);
    size_t loaded_size = ::GetFontData(hdc, dwTable, dwOffset, buffer->data(), size);
    ::DeleteDC(hdc);
    if (size != loaded_size) {
        assert(false);
        // BOOST_LOG_TRIVIAL(error) << "Different loaded(from HFONT) data size." << std::endl;
        return nullptr;    
    }
    return create_font_file(std::move(buffer));
}
#endif // _WIN32

std::optional<Emboss::Glyph> Emboss::letter2glyph(const FontFile &font,
                                                  unsigned int    font_index,
                                                  int             letter,
                                                  float           flatness)
{
    if (!Private::is_valid(font, font_index)) return {};
    auto font_info_opt = Private::load_font_info(font.data->data(), font_index);
    if (!font_info_opt.has_value()) return {};
    return Private::get_glyph(*font_info_opt, letter, flatness);
}

ExPolygons Emboss::text2shapes(FontFileWithCache &font_with_cache,
                               const char *    text,
                               const FontProp &font_prop)
{
    assert(font_with_cache.has_value());

    std::optional<stbtt_fontinfo> font_info_opt;    
    Point    cursor(0, 0);
    ExPolygons result;
    const FontFile& font = *font_with_cache.font_file;
    unsigned int font_index = font_prop.collection_number.has_value()?
        *font_prop.collection_number : 0;
    if (!Private::is_valid(font, font_index)) return {};
    const FontFile::Info& info = font.infos[font_index];
    Emboss::Glyphs& cache = *font_with_cache.cache;
    std::wstring ws = boost::nowide::widen(text);
    for (wchar_t wc: ws){
        if (wc == '\n') { 
            int line_height = info.ascent - info.descent + info.linegap;
            if (font_prop.line_gap.has_value())
                line_height += *font_prop.line_gap;
            line_height = static_cast<int>(line_height / SHAPE_SCALE);

            cursor.x() = 0;
            cursor.y() -= line_height;
            continue;
        } 
        if (wc == '\t') {
            // '\t' = 4*space => same as imgui
            const int count_spaces = 4;
            const Glyph* space = Private::get_glyph(int(' '), font, font_prop, cache, font_info_opt);
            if (space == nullptr) continue;
            cursor.x() += count_spaces * space->advance_width;
            continue;
        }
        if (wc == '\r') continue;

        int unicode = static_cast<int>(wc);
        const Glyph* glyph_ptr = Private::get_glyph(unicode, font, font_prop, cache, font_info_opt);
        if (glyph_ptr == nullptr) continue;
        
        // move glyph to cursor position
        ExPolygons expolygons = glyph_ptr->shape; // copy
        for (ExPolygon &expolygon : expolygons) 
            expolygon.translate(cursor);

        cursor.x() += glyph_ptr->advance_width;
        expolygons_append(result, std::move(expolygons));
    }
    result = Slic3r::union_ex(result);
    return Private::dilate_to_unique_points(result);
}

void Emboss::apply_transformation(const FontProp &font_prop,
                                  Transform3d    &transformation)
{
    if (font_prop.angle.has_value()) {
        double angle_z = *font_prop.angle;
        transformation *= Eigen::AngleAxisd(angle_z, Vec3d::UnitZ());
    }
    if (font_prop.distance.has_value()) {
        Vec3d translate = Vec3d::UnitZ() * (*font_prop.distance);
        transformation.translate(translate);
    }
}

bool Emboss::is_italic(const FontFile &font, unsigned int font_index)
{
    if (font_index >= font.infos.size()) return false;
    std::optional<stbtt_fontinfo> font_info_opt = Private::load_font_info(font.data->data(), font_index);

    if (!font_info_opt.has_value()) return false;
    stbtt_fontinfo *info = &(*font_info_opt);

    // https://docs.microsoft.com/cs-cz/typography/opentype/spec/name
    // https://developer.apple.com/fonts/TrueType-Reference-Manual/RM06/Chap6name.html
    // 2 ==> Style / Subfamily name
    int name_id = 2;
    int length;
    const char* value = stbtt_GetFontNameString(info, &length,
                                               STBTT_PLATFORM_ID_MICROSOFT,
                                               STBTT_MS_EID_UNICODE_BMP,
                                               STBTT_MS_LANG_ENGLISH,                            
                                               name_id);

    // value is big endian utf-16 i need extract only normal chars
    std::string value_str;
    value_str.reserve(length / 2);
    for (int i = 1; i < length; i += 2)
        value_str.push_back(value[i]);

    // lower case
    std::transform(value_str.begin(), value_str.end(), value_str.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    const std::vector<std::string> italics({"italic", "oblique"});
    for (const std::string &it : italics) { 
        if (value_str.find(it) != std::string::npos) { 
            return true; 
        }
    }
    return false; 
}

std::string Emboss::create_range_text(const std::string &text,
                                      const FontFile    &font,
                                      unsigned int       font_index,
                                      bool              *exist_unknown)
{
    if (!Private::is_valid(font, font_index)) return {};
            
    std::wstring ws = boost::nowide::widen(text);

    // need remove symbols not contained in font
    std::sort(ws.begin(), ws.end());

    auto font_info_opt = Private::load_font_info(font.data->data(), 0);
    if (!font_info_opt.has_value()) return {};
    const stbtt_fontinfo *font_info = &(*font_info_opt);

    if (exist_unknown != nullptr) *exist_unknown = false;
    int prev_unicode = -1;
    ws.erase(std::remove_if(ws.begin(), ws.end(),
        [&prev_unicode, font_info, exist_unknown](wchar_t wc) -> bool {
            int unicode = static_cast<int>(wc);

            // skip white spaces
            if (unicode == '\n' || 
                unicode == '\r' || 
                unicode == '\t') return true;

            // is duplicit?
            if (prev_unicode == unicode) return true;
            prev_unicode = unicode;

            // can find in font?
            bool is_unknown = !stbtt_FindGlyphIndex(font_info, unicode);
            if (is_unknown && exist_unknown != nullptr)
                *exist_unknown = true;
            return is_unknown;
        }), ws.end());

    return boost::nowide::narrow(ws);
}

double Emboss::get_shape_scale(const FontProp &fp, const FontFile &ff)
{
    const auto  &cn          = fp.collection_number;
    unsigned int font_index  = (cn.has_value()) ? *cn : 0;
    int          unit_per_em = ff.infos[font_index].unit_per_em;
    double       scale       = fp.size_in_mm / unit_per_em;
    // Shape is scaled for store point coordinate as integer
    return scale * Emboss::SHAPE_SCALE;
}

indexed_triangle_set Emboss::polygons2model(const ExPolygons &shape2d,
                                            const IProjection &projection)
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
            auto p2 = projection.create_front_back(p);
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

std::pair<Vec3f, Vec3f> Emboss::ProjectZ::create_front_back(const Point &p) const
{
    Vec3f front(
        static_cast<float>(p.x() * SHAPE_SCALE), 
        static_cast<float>(p.y() * SHAPE_SCALE),
        0.f);
    return std::make_pair(front, project(front));
}

Vec3f Emboss::ProjectZ::project(const Vec3f &point) const 
{
    Vec3f res = point; // copy
    res.z() = m_depth;
    return res;
}

Transform3d Emboss::create_transformation_onto_surface(const Vec3f &position,
                                                       const Vec3f &normal,
                                                       float        up_limit)
{
    // up and emboss direction for generated model
    Vec3d text_up_dir     = Vec3d::UnitY();
    Vec3d text_emboss_dir = Vec3d::UnitZ();

    // wanted up direction of result
    Vec3d wanted_up_side = Vec3d::UnitZ();
    if (std::fabs(normal.z()) > up_limit) wanted_up_side = Vec3d::UnitY();

    Vec3d wanted_emboss_dir = normal.cast<double>();
    // after cast from float it needs to be normalized again
    wanted_emboss_dir.normalize(); 

    // create perpendicular unit vector to surface triangle normal vector
    // lay on surface of triangle and define up vector for text
    Vec3d wanted_up_dir = wanted_emboss_dir
        .cross(wanted_up_side)
        .cross(wanted_emboss_dir);
    // normal3d is NOT perpendicular to normal_up_dir
    wanted_up_dir.normalize(); 

    // perpendicular to emboss vector of text and normal
    Vec3d  axis_view  = text_emboss_dir.cross(wanted_emboss_dir);
    double angle_view = std::acos(text_emboss_dir.dot(wanted_emboss_dir)); // in rad
    axis_view.normalize();

    Eigen::AngleAxis view_rot(angle_view, axis_view);
    Vec3d wanterd_up_rotated = view_rot.matrix().inverse() * wanted_up_dir;
    wanterd_up_rotated.normalize();
    double angle_up = std::acos(text_up_dir.dot(wanterd_up_rotated));

    // text_view and text_view2 should have same direction
    Vec3d text_view2 = text_up_dir.cross(wanterd_up_rotated);
    Vec3d diff_view  = text_emboss_dir - text_view2;
    if (std::fabs(diff_view.x()) > 1. ||
        std::fabs(diff_view.y()) > 1. ||
        std::fabs(diff_view.z()) > 1.) // oposit direction
        angle_up *= -1.;

    Eigen::AngleAxis up_rot(angle_up, text_emboss_dir);

    Transform3d transform = Transform3d::Identity();
    transform.translate(position.cast<double>());
    transform.rotate(view_rot);
    transform.rotate(up_rot);
    return transform;
}


// OrthoProject

std::pair<Vec3f, Vec3f> Emboss::OrthoProject::create_front_back(const Point &p) const {
    Vec3d front(p.x(), p.y(), 0.);
    Vec3f front_tr = (m_matrix * front).cast<float>();
    return std::make_pair(front_tr, project(front_tr));
}

Vec3f Emboss::OrthoProject::project(const Vec3f &point) const
{
    return point + m_direction;
}
