#include "GCodeBinarizer.hpp"

#if ENABLE_BINARIZED_GCODE_DEBUG
#define NOMINMAX
#include <windows.h>
#include <debugapi.h>
#endif // ENABLE_BINARIZED_GCODE_DEBUG

#include <algorithm>
#include <cassert>

namespace BinaryGCode {

static size_t g_checksum_max_cache_size = 65536;
static const size_t MAX_GCODE_CACHE_SIZE = 65536;

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
static uint16_t compression_types_count()       { return 1 + (uint16_t)ECompressionType::None; }
static uint16_t thumbnail_formats_count()       { return 1 + (uint16_t)EThumbnailFormat::QOI; }
static uint16_t metadata_encoding_types_count() { return 1 + (uint16_t)EMetadataEncodingType::INI; }
static uint16_t gcode_encoding_types_count()    { return 1 + (uint16_t)EGCodeEncodingType::MeatPack; }

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
    {
        // TODO
        break;
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
        auto start_it = src.begin();
        auto end_it = src.begin();
        while (end_it != src.end()) {
            while (end_it != src.end() && *end_it != '\n') {
                ++end_it;
            }
            const std::string item(start_it, end_it);
            const size_t pos = item.find_first_of('=');
            if (pos != std::string::npos) {
                dst.emplace_back(std::make_pair(item.substr(0, pos), item.substr(pos + 1)));
                start_it = ++end_it;
            }
        }
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
    {
        // TODO
        break;
    }
    }
    return true;
}

static bool compress(const std::vector<uint8_t>& src, std::vector<uint8_t>& data_out, ECompressionType compression_type)
{
    return true;
}

static bool uncompress(const std::vector<uint8_t>& src, std::vector<uint8_t>& data_out, ECompressionType compression_type)
{
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
        if (!uncompress(data, uncompressed_data, compression_type))
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
        if (!uncompress(data, uncompressed_data, compression_type))
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

EResult Binarizer::initialize(FILE& file, EGCodeEncodingType gcode_encoding_type, EChecksumType checksum_type)
{
    if (!m_enabled)
        return EResult::Success;

    m_file = &file;

    m_gcode_encoding_type = gcode_encoding_type;
    m_checksum_type = checksum_type;
#if ENABLE_CHECKSUM_BLOCK
    // initialize checksum
    m_checksum = ChecksumBlock();
#endif // ENABLE_CHECKSUM_BLOCK

    // save header
    FileHeader file_header;
    file_header.checksum_type = (uint16_t)m_checksum_type;
    EResult res = file_header.write(*m_file);
    if (res != EResult::Success)
        return res;

    // save file metadata block
    res = m_binary_data.file_metadata.write(*m_file, m_compression_type, m_checksum_type);
    if (res != EResult::Success)
        return res;

    // save printer metadata block
    res = m_binary_data.printer_metadata.write(*m_file, m_compression_type, m_checksum_type);
    if (res != EResult::Success)
        return res;

    // save thumbnail blocks
    for (const ThumbnailBlock& block : m_binary_data.thumbnails) {
        res = block.write(*m_file, m_checksum_type);
        if (res != EResult::Success)
            return res;
    }

    // save print metadata block
    res = m_binary_data.print_metadata.write(*m_file, m_compression_type, m_checksum_type);
    if (res != EResult::Success)
        return res;

    // save slicer metadata block
    res = m_binary_data.slicer_metadata.write(*m_file, m_compression_type, m_checksum_type);
    if (res != EResult::Success)
        return res;

    return EResult::Success;
}

static EResult write_gcode_block(FILE& file, const std::string& raw_data, EGCodeEncodingType encoding_type, ECompressionType compression_type,
    EChecksumType checksum_type)
{
    GCodeBlock block;
    block.encoding_type = (uint16_t)encoding_type;
    block.raw_data = raw_data;
    return block.write(file, compression_type, checksum_type);
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
                const EResult res = write_gcode_block(*m_file, m_gcode_cache, m_gcode_encoding_type, m_compression_type, m_checksum_type);
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
        const EResult res = write_gcode_block(*m_file, m_gcode_cache, m_gcode_encoding_type, m_compression_type, m_checksum_type);
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

