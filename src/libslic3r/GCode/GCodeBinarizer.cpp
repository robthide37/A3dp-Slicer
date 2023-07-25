#include "GCodeBinarizer.hpp"

#if ENABLE_BINARIZED_GCODE_DEBUG
#define NOMINMAX
#include <windows.h>
#include <debugapi.h>
#endif // ENABLE_BINARIZED_GCODE_DEBUG

extern "C" {
#include <heatshrink/heatshrink_encoder.h>
#include <heatshrink/heatshrink_decoder.h>
}

#include <algorithm>
#include <cassert>

namespace BinaryGCode {

static size_t g_checksum_max_cache_size = 65536;
static constexpr const size_t MAX_GCODE_CACHE_SIZE = 65536;

namespace MeatPack {
static constexpr const uint8_t Command_None{ 0 };
//#Command_TogglePacking = 253 -- Currently unused, byte 253 can be reused later.
static constexpr const uint8_t Command_EnablePacking{ 251 };
static constexpr const uint8_t Command_DisablePacking{ 250 };
static constexpr const uint8_t Command_ResetAll{ 249 };
static constexpr const uint8_t Command_QueryConfig{ 248 };
static constexpr const uint8_t Command_EnableNoSpaces{ 247 };
static constexpr const uint8_t Command_DisableNoSpaces{ 246 };
static constexpr const uint8_t Command_SignalByte{ 0xFF };

static constexpr const uint8_t Flag_OmitWhitespaces{ 0x01 };
static constexpr const uint8_t Flag_RemoveComments{ 0x02 };

static constexpr const uint8_t BothUnpackable{ 0b11111111 };
static constexpr const char SpaceReplacedCharacter{ 'E' };

static const std::unordered_map<char, uint8_t> ReverseLookupTbl = {
    { '0',  0b00000000 },
    { '1',  0b00000001 },
    { '2',  0b00000010 },
    { '3',  0b00000011 },
    { '4',  0b00000100 },
    { '5',  0b00000101 },
    { '6',  0b00000110 },
    { '7',  0b00000111 },
    { '8',  0b00001000 },
    { '9',  0b00001001 },
    { '.',  0b00001010 },
    { ' ',  0b00001011 },
    { '\n', 0b00001100 },
    { 'G',  0b00001101 },
    { 'X',  0b00001110 },
    { '\0', 0b00001111 } // never used, 0b1111 is used to indicate the next 8-bits is a full character
};

class MPBinarizer
{
public:
    explicit MPBinarizer(uint8_t flags = 0) : m_flags(flags) {}

    void initialize(std::vector<uint8_t>& dst) {
        initialize_lookup_tables();
        append_command(Command_EnablePacking, dst);
        m_binarizing = true;
    }

    void finalize(std::vector<uint8_t>& dst) {
        if ((m_flags & Flag_RemoveComments) != 0) {
            assert(m_binarizing);
            append_command(Command_ResetAll, dst);
            m_binarizing = false;
        }
    }

    void binarize_line(const std::string& line, std::vector<uint8_t>& dst) {
        auto unified_method = [this](const std::string& line) {
            const std::string::size_type g_idx = line.find('G');
            if (g_idx != std::string::npos) {
                if (g_idx + 1 < line.size() && line[g_idx + 1] >= '0' && line[g_idx + 1] <= '9') {
                    if ((m_flags & Flag_OmitWhitespaces) != 0) {
                        std::string result = line;
                        std::replace(result.begin(), result.end(), 'e', 'E');
                        std::replace(result.begin(), result.end(), 'x', 'X');
                        std::replace(result.begin(), result.end(), 'g', 'G');
                        result.erase(std::remove(result.begin(), result.end(), ' '), result.end());
                        if (result.find('*') != std::string::npos) {
                            size_t checksum = 0;
                            result.erase(std::remove(result.begin(), result.end(), '*'), result.end());
                            for (size_t i = 0; i < result.size(); ++i) {
                                checksum ^= static_cast<uint8_t>(result[i]);
                            }
                            result += "*" + std::to_string(checksum);
                        }
                        result += '\n';
                        return result;
                    }
                    else {
                        std::string result = line;
                        std::replace(result.begin(), result.end(), 'x', 'X');
                        std::replace(result.begin(), result.end(), 'g', 'G');
                        result.erase(std::remove(result.begin(), result.end(), ' '), result.end());
                        if (result.find('*') != std::string::npos) {
                            size_t checksum = 0;
                            result.erase(std::remove(result.begin(), result.end(), '*'), result.end());
                            for (size_t i = 0; i < result.size(); ++i) {
                                checksum ^= static_cast<uint8_t>(result[i]);
                            }
                            result += "*" + std::to_string(checksum);
                        }
                        result += '\n';
                        return result;
                    }
                }
            }
            return line;
        };
        auto is_packable = [](char c) {
            return (s_lookup_tables.packable[static_cast<uint8_t>(c)] != 0);
        };
        auto pack_chars = [](char low, char high) {
            return (((s_lookup_tables.value[static_cast<uint8_t>(high)] & 0xF) << 4) |
                (s_lookup_tables.value[static_cast<uint8_t>(low)] & 0xF));
        };

        if (!line.empty()) {
            if ((m_flags & Flag_RemoveComments) == 0) {
                if (line[0] == ';') {
                    if (m_binarizing) {
                        append_command(Command_DisablePacking, dst);
                        m_binarizing = false;
                    }

                    dst.insert(dst.end(), line.begin(), line.end());
                    return;
                }
            }

            if (line[0] == ';' ||
                line[0] == '\n' ||
                line[0] == '\r' ||
                line.size() < 2)
                return;

            std::string modifiedLine = line.substr(0, line.find(';'));
            if (modifiedLine.empty())
                return;
            auto trim_right = [](const std::string& str) {
                if (str.back() != ' ')
                    return str;
                auto bit = str.rbegin();
                while (bit != str.rend() && *bit == ' ') {
                    ++bit;
                }
                return str.substr(0, std::distance(str.begin(), bit.base()));
            };
            modifiedLine = trim_right(modifiedLine);
            modifiedLine = unified_method(modifiedLine);
            if (modifiedLine.back() != '\n')
                modifiedLine.push_back('\n');
            const size_t line_len = modifiedLine.size();
            std::vector<uint8_t> temp_buffer;
            temp_buffer.reserve(line_len);

            for (size_t line_idx = 0; line_idx < line_len; line_idx += 2) {
                const bool skip_last = line_idx == (line_len - 1);
                const char char_1 = modifiedLine[line_idx];
                const char char_2 = (skip_last ? '\n' : modifiedLine[line_idx + 1]);
                const bool c1_p = is_packable(char_1);
                const bool c2_p = is_packable(char_2);

                if (c1_p) {
                    if (c2_p)
                        temp_buffer.emplace_back(static_cast<uint8_t>(pack_chars(char_1, char_2)));
                    else {
                        temp_buffer.emplace_back(static_cast<uint8_t>(pack_chars(char_1, '\0')));
                        temp_buffer.emplace_back(static_cast<uint8_t>(char_2));
                    }
                }
                else {
                    if (c2_p) {
                        temp_buffer.emplace_back(static_cast<uint8_t>(pack_chars('\0', char_2)));
                        temp_buffer.emplace_back(static_cast<uint8_t>(char_1));
                    }
                    else {
                        temp_buffer.emplace_back(static_cast<uint8_t>(BothUnpackable));
                        temp_buffer.emplace_back(static_cast<uint8_t>(char_1));
                        temp_buffer.emplace_back(static_cast<uint8_t>(char_2));
                    }
                }
            }

            if (!m_binarizing && !temp_buffer.empty()) {
                append_command(Command_EnablePacking, dst);
                m_binarizing = true;
            }

            dst.insert(dst.end(), temp_buffer.begin(), temp_buffer.end());
        }
    }

private:
    unsigned char m_flags{ 0 };
    bool m_binarizing{ false };

    struct LookupTables
    {
        std::array<uint8_t, 256> packable;
        std::array<uint8_t, 256> value;
        bool initialized;
        unsigned char flags;
    };

    static LookupTables s_lookup_tables;

    void append_command(unsigned char cmd, std::vector<uint8_t>& dst) {
        dst.emplace_back(Command_SignalByte);
        dst.emplace_back(Command_SignalByte);
        dst.emplace_back(cmd);
    }

    void initialize_lookup_tables() {
        if (s_lookup_tables.initialized && m_flags == s_lookup_tables.flags)
            return;

        for (const auto& [c, value] : ReverseLookupTbl) {
            const int index = static_cast<int>(c);
            s_lookup_tables.packable[index] = 1;
            s_lookup_tables.value[index] = value;
        }

        if ((m_flags & Flag_OmitWhitespaces) != 0) {
            s_lookup_tables.value[static_cast<uint8_t>(SpaceReplacedCharacter)] = ReverseLookupTbl.at(' ');
            s_lookup_tables.packable[static_cast<uint8_t>(SpaceReplacedCharacter)] = 1;
            s_lookup_tables.packable[static_cast<uint8_t>(' ')] = 0;
        }
        else {
            s_lookup_tables.packable[static_cast<uint8_t>(SpaceReplacedCharacter)] = 0;
            s_lookup_tables.packable[static_cast<uint8_t>(' ')] = 1;
        }

        s_lookup_tables.initialized = true;
        s_lookup_tables.flags = m_flags;
    }
};

MPBinarizer::LookupTables MPBinarizer::s_lookup_tables = { { 0 }, { 0 }, false, 0 };

static constexpr const unsigned char SecondNotPacked{ 0b11110000 };
static constexpr const unsigned char FirstNotPacked{ 0b00001111 };
static constexpr const unsigned char NextPackedFirst{ 0b00000001 };
static constexpr const unsigned char NextPackedSecond{ 0b00000010 };

// See for reference: https://github.com/scottmudge/Prusa-Firmware-MeatPack/blob/MK3_sm_MeatPack/Firmware/meatpack.cpp
static void unbinarize(const std::vector<uint8_t>& src, std::string& dst) {
    bool unbinarizing = false;
    bool cmd_active = false;             // Is a command pending
    uint8_t char_buf = 0;                // Buffers a character if dealing with out-of-sequence pairs
    size_t cmd_count = 0;                // Counts how many command bytes are received (need 2)
    size_t full_char_queue = 0;          // Counts how many full-width characters are to be received
    std::array<uint8_t, 2> char_out_buf; // Output buffer for caching up to 2 characters
    size_t char_out_count = 0;           // Stores number of characters to be read out

    auto handle_command = [&](uint8_t c) {
        switch (c)
        {
        case Command_EnablePacking:  { unbinarizing = true; break; }
        case Command_DisablePacking: { unbinarizing = false; break; }
        case Command_ResetAll:       { unbinarizing = false; break; }
        default:
        case Command_QueryConfig:    { break; }
        }
    };

    auto handle_output_char = [&](uint8_t c) {
        char_out_buf[char_out_count++] = c;
    };

    auto get_char = [](uint8_t c) {
        switch (c)
        {
        case 0b0000: { return '0'; }
        case 0b0001: { return '1'; }
        case 0b0010: { return '2'; }
        case 0b0011: { return '3'; }
        case 0b0100: { return '4'; }
        case 0b0101: { return '5'; }
        case 0b0110: { return '6'; }
        case 0b0111: { return '7'; }
        case 0b1000: { return '8'; }
        case 0b1001: { return '9'; }
        case 0b1010: { return '.'; }
        case 0b1011: { return 'E'; }
        case 0b1100: { return '\n'; }
        case 0b1101: { return 'G'; }
        case 0b1110: { return 'X'; }
        }
        return '\0';
    };

    auto unpack_chars = [&](uint8_t pk, std::array<uint8_t, 2>& chars_out) {
        uint8_t out = 0;

        // If lower 4 bytes is 0b1111, the higher 4 are unused, and next char is full.
        if ((pk & FirstNotPacked) == FirstNotPacked)
            out |= NextPackedFirst;
        else
          chars_out[0] = get_char(pk & 0xF); // Assign lower char

        // Check if upper 4 bytes is 0b1111... if so, we don't need the second char.
        if ((pk & SecondNotPacked) == SecondNotPacked)
            out |= NextPackedSecond;
        else
            chars_out[1] = get_char((pk >> 4) & 0xF); // Assign upper char

        return out;
    };

    auto handle_rx_char = [&](uint8_t c) {
        if (unbinarizing) {
            if (full_char_queue > 0) {
                handle_output_char(c);
                if (char_buf > 0) {
                    handle_output_char(char_buf);
                    char_buf = 0;
                }
                --full_char_queue;
            }
            else {
                std::array<uint8_t, 2> buf = { 0, 0 };
                const uint8_t res = unpack_chars(c, buf);

                if ((res & NextPackedFirst) != 0) {
                    ++full_char_queue;
                    if ((res & NextPackedSecond) != 0)
                        ++full_char_queue;
                    else
                        char_buf = buf[1];
                }
                else {
                    handle_output_char(buf[0]);
                    if (buf[0] != '\n') {
                        if ((res & NextPackedSecond) != 0)
                            ++full_char_queue;
                        else
                            handle_output_char(buf[1]);
                    }
                }
            }
        }
        else // Packing not enabled, just copy character to output
            handle_output_char(c);
    };

    auto get_result_char = [&](std::array<char, 2>& chars_out) {
        if (char_out_count > 0) {
            const size_t res = char_out_count;
            for (uint8_t i = 0; i < char_out_count; ++i) {
                chars_out[i] = (char)char_out_buf[i];
            }
            char_out_count = 0;
            return res;
        }
        return (size_t)0;
    };

    std::vector<uint8_t> unbin_buffer(2 * src.size(), 0);
    auto it_unbin_end = unbin_buffer.begin();

#if ENABLE_BINARIZED_GCODE_DEBUG
    size_t line_start = 0;
#endif // ENABLE_BINARIZED_GCODE_DEBUG
    bool add_space = false;

    auto begin = src.begin();
    auto end = src.end();

    auto it_bin = begin;
    while (it_bin != end) {
        uint8_t c_bin = *it_bin;
        if (c_bin == Command_SignalByte) {
            if (cmd_count > 0) {
                cmd_active = true;
                cmd_count = 0;
            }
            else
              ++cmd_count;
        }
        else {
            if (cmd_active) {
                handle_command(c_bin);
                cmd_active = false;
            }
            else {
                if (cmd_count > 0) {
                    handle_rx_char(Command_SignalByte);
                    cmd_count = 0;
                }

                handle_rx_char(c_bin);
            }
        }

        std::array<char, 2> c_unbin{ 0, 0 };
        const size_t char_count = get_result_char(c_unbin);
        for (size_t i = 0; i < char_count; ++i) {
            // GCodeReader::parse_line_internal() is unable to parse a G line where the data are not separated by spaces
            // so we add them where needed
            if (c_unbin[i] == 'G' && std::distance(unbin_buffer.begin(), it_unbin_end) > 0 && *std::prev(it_unbin_end, 1) == '\n')
                add_space = true;
            else if (c_unbin[i] == '\n')
                add_space = false;

            if (add_space && *std::prev(it_unbin_end, 1) != ' ' &&
                (c_unbin[i] == 'X' || c_unbin[i] == 'Y' || c_unbin[i] == 'Z' || c_unbin[i] == 'E' || c_unbin[i] == 'F')) {
                *it_unbin_end = ' ';
                ++it_unbin_end;
            }

            if (c_unbin[i] != '\n' || std::distance(unbin_buffer.begin(), it_unbin_end) == 0 || *std::prev(it_unbin_end, 1) != '\n') {
                *it_unbin_end = c_unbin[i];
                ++it_unbin_end;
            }

#if ENABLE_BINARIZED_GCODE_DEBUG
            if (c_unbin[i] == '\n') {
                const std::string out(unbin_buffer.begin() + line_start, it_unbin_end);
                if (!out.empty()) {
                    OutputDebugStringA(out.c_str());
                    line_start = std::distance(unbin_buffer.begin(), it_unbin_end);
                }
            }
#endif // ENABLE_BINARIZED_GCODE_DEBUG
        }

        ++it_bin;
    }

    dst.insert(dst.end(), unbin_buffer.begin(), it_unbin_end);
}
} // namespace MeatPack

std::string translate_result(BinaryGCode::EResult result)
{
    switch (result)
    {
    case BinaryGCode::EResult::Success:                     { return "Success"; }
    case BinaryGCode::EResult::ReadError:                   { return "Read error"; }
    case BinaryGCode::EResult::WriteError:                  { return "Write error"; }
    case BinaryGCode::EResult::InvalidMagicNumber:          { return "Invalid magic number"; }
    case BinaryGCode::EResult::InvalidVersionNumber:        { return "Invalid version number"; }
    case BinaryGCode::EResult::InvalidChecksumType:         { return "Invalid checksum type"; }
    case BinaryGCode::EResult::InvalidBlockType:            { return "Invalid block type"; }
    case BinaryGCode::EResult::InvalidCompressionType:      { return "Invalid compression type"; }
    case BinaryGCode::EResult::InvalidMetadataEncodingType: { return "Invalid metadata encoding type"; }
    case BinaryGCode::EResult::InvalidGCodeEncodingType:    { return "Invalid gcode encoding type"; }
    case BinaryGCode::EResult::DataCompressionError:        { return "Data compression error"; }
    case BinaryGCode::EResult::DataUncompressionError:      { return "Data uncompression error"; }
    case BinaryGCode::EResult::MetadataEncodingError:       { return "Data encoding error"; }
    case BinaryGCode::EResult::MetadataDecodingError:       { return "Data decoding error"; }
    case BinaryGCode::EResult::GCodeEncodingError:          { return "GCode encoding error"; }
    case BinaryGCode::EResult::GCodeDecodingError:          { return "GCode decoding error"; }
    case BinaryGCode::EResult::BlockNotFound:               { return "Block not found"; }
    case BinaryGCode::EResult::InvalidChecksum:             { return "Invalid checksum"; }
    case BinaryGCode::EResult::InvalidThumbnailFormat:      { return "Invalid thumbnail format"; }
    case BinaryGCode::EResult::InvalidThumbnailWidth:       { return "Invalid thumbnail width"; }
    case BinaryGCode::EResult::InvalidThumbnailHeight:      { return "Invalid thumbnail height"; }
    case BinaryGCode::EResult::InvalidThumbnailDataSize:    { return "Invalid thumbnail data size"; }
    }
    return std::string();
}

size_t get_checksum_max_cache_size() { return g_checksum_max_cache_size; }
void set_checksum_max_cache_size(size_t size) { g_checksum_max_cache_size = size; }

static uint16_t checksum_types_count()          { return 1 + (uint16_t)EChecksumType::CRC32; }
static uint16_t block_types_count()             { return 1 + (uint16_t)EBlockType::Thumbnail; }
static uint16_t compression_types_count()       { return 1 + (uint16_t)ECompressionType::Heatshrink_12_4; }
static uint16_t thumbnail_formats_count()       { return 1 + (uint16_t)EThumbnailFormat::QOI; }
static uint16_t metadata_encoding_types_count() { return 1 + (uint16_t)EMetadataEncodingType::INI; }
static uint16_t gcode_encoding_types_count()    { return 1 + (uint16_t)EGCodeEncodingType::MeatPackComments; }

static bool write_to_file(FILE& file, const void* data, size_t data_size)
{
    fwrite(data, 1, data_size, &file);
    return !ferror(&file);
}

static bool read_from_file(FILE& file, void* data, size_t data_size)
{
    fread(data, 1, data_size, &file);
    return !ferror(&file);
}

static bool encode_metadata(const std::vector<std::pair<std::string, std::string>>& src, std::vector<uint8_t>& dst,
    EMetadataEncodingType encoding_type)
{
    for (const auto& [key, value] : src) {
        switch (encoding_type)
        {
        case EMetadataEncodingType::INI:
        {
            dst.insert(dst.end(), key.begin(), key.end());
            dst.emplace_back('=');
            dst.insert(dst.end(), value.begin(), value.end());
            dst.emplace_back('\n');
            break;
        }
        }
    }
    return true;
}

static bool decode_metadata(const std::vector<uint8_t>& src, std::vector<std::pair<std::string, std::string>>& dst,
    EMetadataEncodingType encoding_type)
{
    switch (encoding_type)
    {
    case EMetadataEncodingType::INI:
    {
        auto begin_it = src.begin();
        auto end_it = src.begin();
        while (end_it != src.end()) {
            while (end_it != src.end() && *end_it != '\n') {
                ++end_it;
            }
            const std::string item(begin_it, end_it);
            const size_t pos = item.find_first_of('=');
            if (pos != std::string::npos) {
                dst.emplace_back(std::make_pair(item.substr(0, pos), item.substr(pos + 1)));
                begin_it = ++end_it;
            }
        }
        break;
    }
    }

    return true;
}

static bool encode_gcode(const std::string& src, std::vector<uint8_t>& dst, EGCodeEncodingType encoding_type)
{
    switch (encoding_type)
    {
    case EGCodeEncodingType::None:
    {
        dst.insert(dst.end(), src.begin(), src.end());
        break;
    }
    case EGCodeEncodingType::MeatPack:
    case EGCodeEncodingType::MeatPackComments:
    {
        uint8_t binarizer_flags = (encoding_type == EGCodeEncodingType::MeatPack) ? MeatPack::Flag_RemoveComments : 0;
        binarizer_flags |= MeatPack::Flag_OmitWhitespaces;
        MeatPack::MPBinarizer binarizer(binarizer_flags);
        binarizer.initialize(dst);
        auto begin_it = src.begin();
        auto end_it = src.begin();
        while (end_it != src.end()) {
            while (end_it != src.end() && *end_it != '\n') {
                ++end_it;
            }
            const std::string line(begin_it, ++end_it);
            binarizer.binarize_line(line, dst);
            begin_it = end_it;
        }
        binarizer.finalize(dst);
        break;
    }
    }
    return true;
}

static bool decode_gcode(const std::vector<uint8_t>& src, std::string& dst, EGCodeEncodingType encoding_type)
{
    switch (encoding_type)
    {
    case EGCodeEncodingType::None:
    {
        dst.insert(dst.end(), src.begin(), src.end());
        break;
    }
    case EGCodeEncodingType::MeatPack:
    case EGCodeEncodingType::MeatPackComments:
    {
        MeatPack::unbinarize(src, dst);
        break;
    }
    }
    return true;
}

static bool compress(const std::vector<uint8_t>& src, std::vector<uint8_t>& dst, ECompressionType compression_type)
{
    switch (compression_type)
    {
    case ECompressionType::Heatshrink_11_4:
    case ECompressionType::Heatshrink_12_4:
    {
        const uint8_t window_sz = (compression_type == ECompressionType::Heatshrink_11_4) ? 11 : 12;
        const uint8_t lookahead_sz = 4;
        heatshrink_encoder* encoder = heatshrink_encoder_alloc(window_sz, lookahead_sz);
        if (encoder == nullptr)
            return false;

        // calculate the maximum compressed size (assuming a conservative estimate)
        const size_t src_size = src.size();
        const size_t max_compressed_size = src_size + (src_size >> 2);
        dst.resize(max_compressed_size);

        uint8_t* buf = const_cast<uint8_t*>(src.data());
        uint8_t* outbuf = dst.data();

        // compress data
        size_t tosink = src_size;
        size_t output_size = 0;
        while (tosink > 0) {
            size_t sunk = 0;
            const HSE_sink_res sink_res = heatshrink_encoder_sink(encoder, buf, tosink, &sunk);
            if (sink_res != HSER_SINK_OK) {
                heatshrink_encoder_free(encoder);
                return false;
            }
            if (sunk == 0)
                // all input data processed
                break;

            tosink -= sunk;
            buf += sunk;

            size_t polled = 0;
            const HSE_poll_res poll_res = heatshrink_encoder_poll(encoder, outbuf + output_size, max_compressed_size - output_size, &polled);
            if (poll_res < 0) {
                heatshrink_encoder_free(encoder);
                return false;
            }
            output_size += polled;
        }

        // input data finished
        const HSE_finish_res finish_res = heatshrink_encoder_finish(encoder);
        if (finish_res < 0) {
            heatshrink_encoder_free(encoder);
            return false;
        }

        // poll for final output
        size_t polled = 0;
        const HSE_poll_res poll_res = heatshrink_encoder_poll(encoder, outbuf + output_size, max_compressed_size - output_size, &polled);
        if (poll_res < 0) {
            heatshrink_encoder_free(encoder);
            return false;
        }
        dst.resize(output_size + polled);
        heatshrink_encoder_free(encoder);
        break;
    }
    case ECompressionType::None:
    default:
    {
        break;
    }
    }

    return true;
}

static bool uncompress(const std::vector<uint8_t>& src, std::vector<uint8_t>& dst, ECompressionType compression_type, size_t uncompressed_size)
{
    switch (compression_type)
    {
    case ECompressionType::Heatshrink_11_4:
    case ECompressionType::Heatshrink_12_4:
    {
        const uint8_t window_sz = (compression_type == ECompressionType::Heatshrink_11_4) ? 11 : 12;
        const uint8_t lookahead_sz = 4;
        const uint16_t input_buffer_size = 2048;
        heatshrink_decoder* decoder = heatshrink_decoder_alloc(input_buffer_size, window_sz, lookahead_sz);
        if (decoder == nullptr)
            return false;

        dst.resize(uncompressed_size);

        uint8_t* buf = const_cast<uint8_t*>(src.data());
        uint8_t* outbuf = dst.data();

        uint32_t sunk = 0;
        uint32_t polled = 0;

        const size_t compressed_size = src.size();
        while (sunk < compressed_size) {
            size_t count = 0;
            const HSD_sink_res sink_res = heatshrink_decoder_sink(decoder, &buf[sunk], compressed_size - sunk, &count);
            if (sink_res < 0) {
                heatshrink_decoder_free(decoder);
                return false;
            }

            sunk += (uint32_t)count;

            HSD_poll_res poll_res;
            do {
                poll_res  = heatshrink_decoder_poll(decoder, &outbuf[polled], uncompressed_size - polled, &count);
                if (poll_res < 0) {
                    heatshrink_decoder_free(decoder);
                    return false;
                }
                polled += (uint32_t)count;
            } while (polled < uncompressed_size && poll_res == HSDR_POLL_MORE);
        }

        const HSD_finish_res finish_res = heatshrink_decoder_finish(decoder);
        if (finish_res < 0) {
            heatshrink_decoder_free(decoder);
            return false;
        }

        heatshrink_decoder_free(decoder);
        break;
    }
    case ECompressionType::None:
    default:
    {
        break;
    }
    }

    return true;
}

static uint32_t crc32_sw(const uint8_t* buffer, uint32_t length, uint32_t crc)
{
    uint32_t value = crc ^ 0xFFFFFFFF;
    while (length--) {
        value ^= (uint32_t)*buffer++;
        for (int bit = 0; bit < 8; bit++) {
            if (value & 1)
                value = (value >> 1) ^ 0xEDB88320;
            else
                value >>= 1;
        }
    }
    value ^= 0xFFFFFFFF;
    return value;
}

std::vector<uint8_t> encode(const void* data, size_t data_size)
{
    std::vector<uint8_t> ret(data_size);
    memcpy(ret.data(), data, data_size);
    return ret;
}

Checksum::Checksum(EChecksumType type)
: m_type(type)
{
    if (m_type != EChecksumType::None)
        m_checksum = std::vector<uint8_t>(checksum_size(m_type), '\0');
}

EChecksumType Checksum::get_type() const
{
    return m_type;
}

void Checksum::append(const std::vector<uint8_t>& data)
{
    size_t remaining_data_size = std::distance(data.begin(), data.end());
    auto it_begin = data.begin();
    while (remaining_data_size + m_cache.size() > g_checksum_max_cache_size) {
        update();
        if (remaining_data_size > g_checksum_max_cache_size) {
            m_cache.insert(m_cache.end(), it_begin, it_begin + g_checksum_max_cache_size);
            it_begin += g_checksum_max_cache_size;
            remaining_data_size -= g_checksum_max_cache_size;
        }
    }

    m_cache.insert(m_cache.end(), it_begin, data.end());
}

bool Checksum::matches(Checksum& other)
{
    update();
    other.update();
    return m_checksum == other.m_checksum;
}

EResult Checksum::write(FILE& file)
{
    if (m_type != EChecksumType::None) {
        update();    
        if (!write_to_file(file, (const void*)m_checksum.data(), m_checksum.size()))
            return EResult::WriteError;
    }
    return EResult::Success;
}

EResult Checksum::read(FILE& file)
{
    if (m_type != EChecksumType::None) {
        if (!read_from_file(file, (void*)m_checksum.data(), m_checksum.size()))
            return EResult::ReadError;
    }
    return EResult::Success;
}

void Checksum::update()
{
    if (m_cache.empty())
        return;

    switch (m_type)
    {
    case EChecksumType::None:
    {
        break;
    }
    case EChecksumType::CRC32:
    {
        const uint32_t old_crc = *(uint32_t*)m_checksum.data();
        const uint32_t new_crc = crc32_sw(m_cache.data(), (uint32_t)m_cache.size(), old_crc);
        *(uint32_t*)m_checksum.data() = new_crc;
        break;
    }
    }

    m_cache.clear();
}

EResult FileHeader::write(FILE& file) const
{
    if (magic != *(uint32_t*)(MAGIC.data()))
        return EResult::InvalidMagicNumber;
    if (checksum_type >= checksum_types_count())
        return EResult::InvalidChecksumType;

    if (!write_to_file(file, (const void*)&magic, sizeof(magic)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&version, sizeof(version)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&checksum_type, sizeof(checksum_type)))
        return EResult::WriteError;

    return EResult::Success;
}

EResult FileHeader::read(FILE& file, const uint32_t* const max_version)
{
    if (!read_from_file(file, (void*)&magic, sizeof(magic)))
        return EResult::ReadError;
    if (magic != *(uint32_t*)(MAGIC.data()))
        return EResult::InvalidMagicNumber;

    if (!read_from_file(file, (void*)&version, sizeof(version)))
        return EResult::ReadError;
    if (max_version != nullptr && version > *max_version)
        return EResult::InvalidVersionNumber;

    if (!read_from_file(file, (void*)&checksum_type, sizeof(checksum_type)))
        return EResult::ReadError;
    if (checksum_type >= checksum_types_count())
        return EResult::InvalidChecksumType;

    return EResult::Success;
}

void BlockHeader::update_checksum(Checksum& checksum) const
{
    checksum.append(encode((const void*)&type, sizeof(type)));
    checksum.append(encode((const void*)&compression, sizeof(compression)));
    checksum.append(encode((const void*)&uncompressed_size, sizeof(uncompressed_size)));
    if (compression != (uint16_t)ECompressionType::None)
        checksum.append(encode((const void*)&compressed_size, sizeof(compressed_size)));
}

EResult BlockHeader::write(FILE& file) const
{
    if (!write_to_file(file, (const void*)&type, sizeof(type)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&compression, sizeof(compression)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&uncompressed_size, sizeof(uncompressed_size)))
        return EResult::WriteError;
    if (compression != (uint16_t)ECompressionType::None) {
        if (!write_to_file(file, (const void*)&compressed_size, sizeof(compressed_size)))
            return EResult::WriteError;
    }
    return EResult::Success;
}

EResult BlockHeader::read(FILE& file)
{
    if (!read_from_file(file, (void*)&type, sizeof(type)))
        return EResult::ReadError;
    if (type >= block_types_count())
        return EResult::InvalidBlockType;

    if (!read_from_file(file, (void*)&compression, sizeof(compression)))
        return EResult::ReadError;
    if (compression >= compression_types_count())
        return EResult::InvalidCompressionType;

    if (!read_from_file(file, (void*)&uncompressed_size, sizeof(uncompressed_size)))
        return EResult::ReadError;
    if (compression != (uint16_t)ECompressionType::None) {
        if (!read_from_file(file, (void*)&compressed_size, sizeof(compressed_size)))
            return EResult::ReadError;
    }

    return EResult::Success;
}

EResult BaseMetadataBlock::write(FILE& file, EBlockType block_type, ECompressionType compression_type, Checksum& checksum) const
{
    if (encoding_type > metadata_encoding_types_count())
        return EResult::InvalidMetadataEncodingType;

    BlockHeader block_header = { (uint16_t)block_type, (uint16_t)compression_type, (uint32_t)0 };
    std::vector<uint8_t> out_data;
    if (!raw_data.empty()) {
        // process payload encoding
        std::vector<uint8_t> uncompressed_data;
        if (!encode_metadata(raw_data, uncompressed_data, (EMetadataEncodingType)encoding_type))
            return EResult::MetadataEncodingError;
        // process payload compression
        block_header.uncompressed_size = (uint32_t)uncompressed_data.size();
        std::vector<uint8_t> compressed_data;
        if (compression_type != ECompressionType::None) {
            if (!compress(uncompressed_data, compressed_data, compression_type))
                return EResult::DataCompressionError;
            block_header.compressed_size = (uint32_t)compressed_data.size();
        }
        out_data.swap((compression_type == ECompressionType::None) ? uncompressed_data : compressed_data);
    }

    // write block header
    EResult res = block_header.write(file);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block payload
    if (!write_to_file(file, (const void*)&encoding_type, sizeof(encoding_type)))
        return EResult::WriteError;
    if (!out_data.empty()) {
        if (!write_to_file(file, (const void*)out_data.data(), out_data.size()))
            return EResult::WriteError;
    }

    if (checksum.get_type() != EChecksumType::None) {
        // update checksum with block header
        block_header.update_checksum(checksum);
        // update checksum with block payload
        checksum.append(encode((const void*)&encoding_type, sizeof(encoding_type)));
        if (!out_data.empty())
            checksum.append(out_data);
    }

    return EResult::Success;
}

EResult BaseMetadataBlock::read_data(FILE& file, const BlockHeader& block_header)
{
    const ECompressionType compression_type = (ECompressionType)block_header.compression;

    if (!read_from_file(file, (void*)&encoding_type, sizeof(encoding_type)))
        return EResult::ReadError;
    if (encoding_type > metadata_encoding_types_count())
        return EResult::InvalidMetadataEncodingType;

    std::vector<uint8_t> data;
    const size_t data_size = (compression_type == ECompressionType::None) ? block_header.uncompressed_size : block_header.compressed_size;
    if (data_size > 0) {
        data.resize(data_size);
        if (!read_from_file(file, (void*)data.data(), data_size))
            return EResult::ReadError;
    }

    std::vector<uint8_t> uncompressed_data;
    if (compression_type != ECompressionType::None) {
        if (!uncompress(data, uncompressed_data, compression_type, block_header.uncompressed_size))
            return EResult::DataUncompressionError;
    }

    if (!decode_metadata((compression_type == ECompressionType::None) ? data : uncompressed_data, raw_data, (EMetadataEncodingType)encoding_type))
        return EResult::MetadataDecodingError;

    return EResult::Success;
}

EResult FileMetadataBlock::write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const
{
    Checksum cs(checksum_type);

    // write block header, payload
    EResult res = BaseMetadataBlock::write(file, EBlockType::FileMetadata, compression_type, cs);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block checksum
    if (checksum_type != EChecksumType::None)
        return cs.write(file);

    return EResult::Success;
}

EResult FileMetadataBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // read block payload
    EResult res = BaseMetadataBlock::read_data(file, block_header);
    if (res != EResult::Success)
        // propagate error
        return res;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult ThumbnailBlock::write(FILE& file, EChecksumType checksum_type) const
{
    if (format >= thumbnail_formats_count())
        return EResult::InvalidThumbnailFormat;
    if (width == 0)
        return EResult::InvalidThumbnailWidth;
    if (height == 0)
        return EResult::InvalidThumbnailHeight;
    if (data.size() == 0)
        return EResult::InvalidThumbnailDataSize;

    // write block header
    const BlockHeader block_header = { (uint16_t)EBlockType::Thumbnail, (uint16_t)ECompressionType::None, (uint32_t)data.size() };
    EResult res = block_header.write(file);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block payload
    if (!write_to_file(file, (const void*)&format, sizeof(format)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&width, sizeof(width)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)&height, sizeof(height)))
        return EResult::WriteError;
    if (!write_to_file(file, (const void*)data.data(), data.size()))
        return EResult::WriteError;

    if (checksum_type != EChecksumType::None) {
        Checksum cs(checksum_type);
        // update checksum with block header
        block_header.update_checksum(cs);
        // update checksum with block payload
        update_checksum(cs);
        // write block checksum
        res = cs.write(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult ThumbnailBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // read block payload
    if (!read_from_file(file, (void*)&format, sizeof(format)))
        return EResult::ReadError;
    if (format >= thumbnail_formats_count())
        return EResult::InvalidThumbnailFormat;
    if (!read_from_file(file, (void*)&width, sizeof(width)))
        return EResult::ReadError;
    if (width == 0)
        return EResult::InvalidThumbnailWidth;
    if (!read_from_file(file, (void*)&height, sizeof(height)))
        return EResult::ReadError;
    if (height == 0)
        return EResult::InvalidThumbnailHeight;
    if (block_header.uncompressed_size == 0)
        return EResult::InvalidThumbnailDataSize;

    data.resize(block_header.uncompressed_size);
    if (!read_from_file(file, (void*)data.data(), block_header.uncompressed_size))
        return EResult::ReadError;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        const EResult res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

void ThumbnailBlock::update_checksum(Checksum& checksum) const
{
    checksum.append(encode((const void*)&format, sizeof(format)));
    checksum.append(encode((const void*)&width, sizeof(width)));
    checksum.append(encode((const void*)&height, sizeof(height)));
    checksum.append(data);
}

EResult PrinterMetadataBlock::write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const
{
    Checksum cs(checksum_type);

    // write block header, payload
    EResult res = BaseMetadataBlock::write(file, EBlockType::PrinterMetadata, compression_type, cs);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block checksum
    if (checksum_type != EChecksumType::None)
        return cs.write(file);

    return EResult::Success;
}

EResult PrinterMetadataBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // read block payload
    EResult res = BaseMetadataBlock::read_data(file, block_header);
    if (res != EResult::Success)
        // propagate error
        return res;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult PrintMetadataBlock::write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const
{
    Checksum cs(checksum_type);

    // write block header, payload
    EResult res = BaseMetadataBlock::write(file, EBlockType::PrintMetadata, compression_type, cs);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block checksum
    if (checksum_type != EChecksumType::None)
        return cs.write(file);

    return EResult::Success;
}

EResult PrintMetadataBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // read block payload
    EResult res = BaseMetadataBlock::read_data(file, block_header);
    if (res != EResult::Success)
        // propagate error
        return res;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult SlicerMetadataBlock::write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const
{
    Checksum cs(checksum_type);

    // write block header, payload
    EResult res = BaseMetadataBlock::write(file, EBlockType::SlicerMetadata, compression_type, cs);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block checksum
    if (checksum_type != EChecksumType::None)
        return cs.write(file);

    return EResult::Success;
}

EResult SlicerMetadataBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // read block payload
    EResult res = BaseMetadataBlock::read_data(file, block_header);
    if (res != EResult::Success)
        // propagate error
        return res;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult GCodeBlock::write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const
{
    if (encoding_type > gcode_encoding_types_count())
        return EResult::InvalidGCodeEncodingType;

    BlockHeader block_header = { (uint16_t)EBlockType::GCode, (uint16_t)compression_type, (uint32_t)0 };
    std::vector<uint8_t> out_data;
    if (!raw_data.empty()) {
        // process payload encoding
        std::vector<uint8_t> uncompressed_data;
        if (!encode_gcode(raw_data, uncompressed_data, (EGCodeEncodingType)encoding_type))
            return EResult::GCodeEncodingError;
        // process payload compression
        block_header.uncompressed_size = (uint32_t)uncompressed_data.size();
        std::vector<uint8_t> compressed_data;
        if (compression_type != ECompressionType::None) {
            if (!compress(uncompressed_data, compressed_data, compression_type))
                return EResult::DataCompressionError;
            block_header.compressed_size = (uint32_t)compressed_data.size();
        }
        out_data.swap((compression_type == ECompressionType::None) ? uncompressed_data : compressed_data);
    }

    // write block header
    EResult res = block_header.write(file);
    if (res != EResult::Success)
        // propagate error
        return res;

    // write block payload
    if (!write_to_file(file, (const void*)&encoding_type, sizeof(encoding_type)))
        return EResult::WriteError;
    if (!out_data.empty()) {
#if ENABLE_BINARIZED_GCODE_DEBUG
        const std::string out = "GCodeBlock data size:" + std::to_string(out_data.size()) + "\n";
        OutputDebugStringA(out.c_str());
#endif // ENABLE_BINARIZED_GCODE_DEBUG
        if (!write_to_file(file, (const void*)out_data.data(), out_data.size()))
            return EResult::WriteError;
    }

    // write checksum
    if (checksum_type != EChecksumType::None) {
        Checksum cs(checksum_type);
        // update checksum with block header
        block_header.update_checksum(cs);
        // update checksum with block payload
        cs.append(encode((const void*)&encoding_type, sizeof(encoding_type)));
        if (!out_data.empty())
            cs.append(out_data);
        res = cs.write(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

EResult GCodeBlock::read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
  const ECompressionType compression_type = (ECompressionType)block_header.compression;

    if (!read_from_file(file, (void*)&encoding_type, sizeof(encoding_type)))
        return EResult::ReadError;
    if (encoding_type > gcode_encoding_types_count())
        return EResult::InvalidGCodeEncodingType;

    std::vector<uint8_t> data;
    const size_t data_size = (compression_type == ECompressionType::None) ? block_header.uncompressed_size : block_header.compressed_size;
    if (data_size > 0) {
        data.resize(data_size);
        if (!read_from_file(file, (void*)data.data(), data_size))
            return EResult::ReadError;
    }

    std::vector<uint8_t> uncompressed_data;
    if (compression_type != ECompressionType::None) {
        if (!uncompress(data, uncompressed_data, compression_type, block_header.uncompressed_size))
            return EResult::DataUncompressionError;
    }

    if (!decode_gcode((compression_type == ECompressionType::None) ? data : uncompressed_data, raw_data, (EGCodeEncodingType)encoding_type))
        return EResult::GCodeDecodingError;

    const EChecksumType checksum_type = (EChecksumType)file_header.checksum_type;
    if (checksum_type != EChecksumType::None) {
        // read block checksum
        Checksum cs(checksum_type);
        const EResult res = cs.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

    return EResult::Success;
}

#if ENABLE_CHECKSUM_BLOCK
EResult ChecksumBlock::write(FILE& file) const
{
    if (!data.empty()) {
        const BlockHeader block_header = { (uint16_t)EBlockType::Checksum, (uint16_t)ECompressionType::None, (uint32_t)data.size() };
        // write block header
        const EResult res = block_header.write(file);
        if (res != EResult::Success)
            // propagate error
            return res;
        // write block payload
        if (!write_to_file(file, (const void*)data.data(), data.size()))
            return EResult::WriteError;
    }

    return EResult::Success;
}

EResult ChecksumBlock::read_data(FILE& file, const BlockHeader& block_header)
{
    if (block_header.uncompressed_size > 0) {
        data.resize(block_header.uncompressed_size);
        if (!read_from_file(file, (void*)data.data(), block_header.uncompressed_size))
            return EResult::ReadError;
    }
    else
        data.clear();

    return EResult::Success;
}
#endif // ENABLE_CHECKSUM_BLOCK

EResult Binarizer::initialize(FILE& file, const BinarizerConfig& config)
{
    if (!m_enabled)
        return EResult::Success;

    m_file = &file;

    m_config = config;
#if ENABLE_CHECKSUM_BLOCK
    // initialize checksum
    m_checksum = ChecksumBlock();
#endif // ENABLE_CHECKSUM_BLOCK

    // save header
    FileHeader file_header;
    file_header.checksum_type = (uint16_t)m_config.checksum;
    EResult res = file_header.write(*m_file);
    if (res != EResult::Success)
        return res;

    // save file metadata block
    res = m_binary_data.file_metadata.write(*m_file, m_config.compression.file_metadata, m_config.checksum);
    if (res != EResult::Success)
        return res;

    // save printer metadata block
    res = m_binary_data.printer_metadata.write(*m_file, m_config.compression.printer_metadata, m_config.checksum);
    if (res != EResult::Success)
        return res;

    // save thumbnail blocks
    for (const ThumbnailBlock& block : m_binary_data.thumbnails) {
        res = block.write(*m_file, m_config.checksum);
        if (res != EResult::Success)
            return res;
    }

    // save print metadata block
    res = m_binary_data.print_metadata.write(*m_file, m_config.compression.print_metadata, m_config.checksum);
    if (res != EResult::Success)
        return res;

    // save slicer metadata block
    res = m_binary_data.slicer_metadata.write(*m_file, m_config.compression.slicer_metadata, m_config.checksum);
    if (res != EResult::Success)
        return res;

    return EResult::Success;
}

static EResult write_gcode_block(FILE& file, const std::string& raw_data, const BinarizerConfig& config)
{
    GCodeBlock block;
    block.encoding_type = (uint16_t)config.gcode_encoding;
    block.raw_data = raw_data;
    return block.write(file, config.compression.gcode, config.checksum);
}

EResult Binarizer::append_gcode(const std::string& gcode)
{
    if (gcode.empty())
        return EResult::Success;

    assert(m_file != nullptr);
    if (m_file == nullptr)
        return EResult::WriteError;

    auto it_begin = gcode.begin();
    do {
        const size_t begin_pos = std::distance(gcode.begin(), it_begin);
        const size_t end_line_pos = gcode.find_first_of('\n', begin_pos);
        if (end_line_pos == std::string::npos)
            return EResult::WriteError;

        const size_t line_size = 1 + end_line_pos - begin_pos;
        if (line_size + m_gcode_cache.length() > MAX_GCODE_CACHE_SIZE) {
            if (!m_gcode_cache.empty()) {
                const EResult res = write_gcode_block(*m_file, m_gcode_cache, m_config);
                if (res != EResult::Success)
                    // propagate error
                    return res;
                m_gcode_cache.clear();
            }
        }

        if (line_size > MAX_GCODE_CACHE_SIZE)
            return EResult::WriteError;

        m_gcode_cache.insert(m_gcode_cache.end(), it_begin, it_begin + line_size);
        it_begin += line_size;
    }
    while (it_begin != gcode.end());

    return EResult::Success;
}

EResult Binarizer::finalize()
{
    if (!m_enabled)
        return EResult::Success;

    // save gcode cache, if not empty
    if (!m_gcode_cache.empty()) {
        const EResult res = write_gcode_block(*m_file, m_gcode_cache, m_config);
        if (res != EResult::Success)
            // propagate error
            return res;
    }

#if ENABLE_CHECKSUM_BLOCK
    if (m_checksum_type != EChecksumType::None) {
        // save checksum
        // dummy checksum until it is not properly implemented
        switch (m_checksum_type)
        {
        case EChecksumType::CRC32:
        case EChecksumType::MD5:
        {
            m_checksum.data.clear();
            break;
        }
        }

        res = m_checksum.write(file);
        if (res != EResult::Success)
            return res;
    }
#endif // ENABLE_CHECKSUM_BLOCK

    return EResult::Success;
}

bool is_valid_binary_gcode(FILE& file)
{
    // cache file position
    const long curr_pos = ftell(&file);
    rewind(&file);

    std::array<uint8_t, 4> magic;
    fread((void*)magic.data(), 1, magic.size(), &file);
    if (ferror(&file))
        return false;
    else {
        // restore file position
        fseek(&file, curr_pos, SEEK_SET);
        return magic == MAGIC;
    }
}

EResult read_header(FILE& file, FileHeader& header, const uint32_t* const max_version)
{
    rewind(&file);
    return header.read(file, max_version);
}

static EResult checksums_match(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    // cache file position
    const long curr_pos = ftell(&file);

    Checksum curr_cs((EChecksumType)file_header.checksum_type);
    // update block checksum block header
    block_header.update_checksum(curr_cs);

    // read block payload
    size_t remaining_payload_size = block_payload_size(block_header);
    while (remaining_payload_size > 0) {
        const size_t size_to_read = std::min(remaining_payload_size, g_checksum_max_cache_size);
        std::vector<uint8_t> payload(size_to_read);
        if (!read_from_file(file, payload.data(), payload.size()))
            return EResult::ReadError;
        curr_cs.append(payload);
        remaining_payload_size -= size_to_read;
    }

    // read checksum
    Checksum read_cs((EChecksumType)file_header.checksum_type);
    EResult res = read_cs.read(file);
    if (res != EResult::Success)
        // propagate error
        return res;

    // Verify checksum 
    if (!curr_cs.matches(read_cs))
        return EResult::InvalidChecksum;

    // restore file position
    fseek(&file, curr_pos, SEEK_SET);

    return EResult::Success;
}

EResult read_next_block_header(FILE& file, const FileHeader& file_header, BlockHeader& block_header, bool verify_checksum)
{
    if (verify_checksum && (EChecksumType)file_header.checksum_type != EChecksumType::None) {
        const EResult res = block_header.read(file);
        if (res != EResult::Success)
            // propagate error
            return res;

        return checksums_match(file, file_header, block_header);
    }
    else
        return block_header.read(file);
}

EResult read_next_block_header(FILE& file, const FileHeader& file_header, BlockHeader& block_header, EBlockType type, bool verify_checksum)
{
    // cache file position
    const long curr_pos = ftell(&file);

    do {
        EResult res = read_next_block_header(file, file_header, block_header, false);
        if (res != EResult::Success)
            // propagate error
            return res;
        else if (feof(&file)) {
            // block not found
            // restore file position
            fseek(&file, curr_pos, SEEK_SET);
            return EResult::BlockNotFound;
        }
        else if ((EBlockType)block_header.type == type) {
            // block found
            if (verify_checksum) {
                res = checksums_match(file, file_header, block_header);
                if (res != EResult::Success)
                    // propagate error
                    return res;
                else
                    break;
            }
        }

        if (!feof(&file)) {
            res = skip_block_content(file, file_header, block_header);
            if (res != EResult::Success)
                // propagate error
                return res;
        }
    } while (true);

    return EResult::Success;
}

EResult skip_block_payload(FILE& file, const BlockHeader& block_header)
{
    fseek(&file, (long)block_payload_size(block_header), SEEK_CUR);
    return ferror(&file) ? EResult::ReadError : EResult::Success;
}

EResult skip_block_content(FILE& file, const FileHeader& file_header, const BlockHeader& block_header)
{
    fseek(&file, (long)block_content_size(file_header, block_header), SEEK_CUR);
    return ferror(&file) ? EResult::ReadError : EResult::Success;
}

size_t block_parameters_size(EBlockType type)
{
    switch (type)
    {
    case EBlockType::FileMetadata:    { return FileMetadataBlock::get_parameters_size(); }
    case EBlockType::GCode:           { return GCodeBlock::get_parameters_size(); }
    case EBlockType::SlicerMetadata:  { return SlicerMetadataBlock::get_parameters_size(); }
    case EBlockType::PrinterMetadata: { return PrinterMetadataBlock::get_parameters_size(); }
    case EBlockType::PrintMetadata:   { return PrintMetadataBlock::get_parameters_size(); }
    case EBlockType::Thumbnail:       { return ThumbnailBlock::get_parameters_size(); }
    }
    return 0;
}

size_t block_payload_size(const BlockHeader& block_header)
{
    size_t ret = block_parameters_size((EBlockType)block_header.type);
    ret += ((ECompressionType)block_header.compression == ECompressionType::None) ?
        block_header.uncompressed_size : block_header.compressed_size;
    return ret;
}

size_t checksum_size(EChecksumType type)
{
    switch (type)
    {
    case EChecksumType::None:  { return 0; }
    case EChecksumType::CRC32: { return 4; }
    }
    return 0;
}

extern size_t block_content_size(const FileHeader& file_header, const BlockHeader& block_header)
{
#if ENABLE_CHECKSUM_BLOCK
    return ((EBlockType)block_header.type == EBlockType::Checksum) ?
        block_payload_size(block_header) : block_payload_size(block_header) + checksum_size((EChecksumType)file_header.checksum_type);
#else
    return block_payload_size(block_header) + checksum_size((EChecksumType)file_header.checksum_type);
#endif // ENABLE_CHECKSUM_BLOCK
}

} // namespace BinaryGCode

