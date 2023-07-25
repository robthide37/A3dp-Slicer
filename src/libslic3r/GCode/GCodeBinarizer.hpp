#ifndef slic3r_GCode_GCodeBinarizer_hpp_
#define slic3r_GCode_GCodeBinarizer_hpp_

#ifdef _WIN32
#define ENABLE_BINARIZED_GCODE_DEBUG 1
#endif // _WIN32

#define ENABLE_CHECKSUM_BLOCK 0

#include <array>
#include <vector>
#include <string>
#include <cstdio>
#include <functional>

namespace BinaryGCode {

static const std::array<uint8_t, 4> MAGIC{ 'G', 'C', 'D', 'E' };
static const uint32_t VERSION = 1;

enum class EResult : uint16_t
{
    Success,
    ReadError,
    WriteError,
    InvalidMagicNumber,
    InvalidVersionNumber,
    InvalidChecksumType,
    InvalidBlockType,
    InvalidCompressionType,
    InvalidMetadataEncodingType,
    InvalidGCodeEncodingType,
    DataCompressionError,
    DataUncompressionError,
    MetadataEncodingError,
    MetadataDecodingError,
    GCodeEncodingError,
    GCodeDecodingError,
    BlockNotFound,
    InvalidChecksum,
    InvalidThumbnailFormat,
    InvalidThumbnailWidth,
    InvalidThumbnailHeight,
    InvalidThumbnailDataSize
};

// Returns a string description of the given result
extern std::string translate_result(BinaryGCode::EResult result);

enum class EChecksumType : uint16_t
{
    None,
    CRC32
};

class Checksum
{
public:
    // Constructs a checksum of the given type.
    // The checksum data are sized accordingly.
    explicit Checksum(EChecksumType type);

    EChecksumType get_type() const;

    // Appends the given data to the cache and performs a checksum update if 
    // the size of the cache exceeds the max checksum cache size.
    void append(const std::vector<uint8_t>& data);
    // Returns true if the given checksum is equal to this one
    bool matches(Checksum& other);

    EResult write(FILE& file);
    EResult read(FILE& file);

private:
    EChecksumType m_type;
    std::vector<uint8_t> m_cache;
    std::vector<uint8_t> m_checksum;

    void update();
};

struct FileHeader
{    
    uint32_t magic{ *(uint32_t*)(MAGIC.data()) };
    uint32_t version{ VERSION };
    uint16_t checksum_type{ (uint16_t)EChecksumType::None };

    EResult write(FILE& file) const;
    EResult read(FILE& file, const uint32_t* const max_version);
};

enum class EBlockType : uint16_t
{
#if ENABLE_CHECKSUM_BLOCK
    Checksum,
#endif // ENABLE_CHECKSUM_BLOCK
    FileMetadata,
    GCode,
    SlicerMetadata,
    PrinterMetadata,
    PrintMetadata,
    Thumbnail
};

enum class ECompressionType : uint16_t
{
    None,
    Heatshrink_11_4,
    Heatshrink_12_4,
};

struct BlockHeader
{
    uint16_t type{ 0 };
    uint16_t compression{ 0 };
    uint32_t uncompressed_size{ 0 };
    uint32_t compressed_size{ 0 };

    // Updates the given checksum with the data of this BlockHeader
    void update_checksum(Checksum& checksum) const;

    EResult write(FILE& file) const;
    EResult read(FILE& file);
};

enum class EMetadataEncodingType : uint16_t
{
    INI,
};

struct BaseMetadataBlock
{
    // type of data encoding
    uint16_t encoding_type{ 0 };
    // data in key/value form
    std::vector<std::pair<std::string, std::string>> raw_data;

    // write block header and data in encoded format
    EResult write(FILE& file, EBlockType block_type, ECompressionType compression_type, Checksum& checksum) const;
    // read block data in encoded format
    EResult read_data(FILE& file, const BlockHeader& block_header);

    static size_t get_parameters_size() { return sizeof(encoding_type); }
};

struct FileMetadataBlock : public BaseMetadataBlock
{
    // write block header and data
    EResult write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);
};

enum class EThumbnailFormat : uint16_t
{
    PNG,
    JPG,
    QOI
};

struct ThumbnailBlock
{
    uint16_t format{ 0 };
    uint16_t width{ 0 };
    uint16_t height{ 0 };
    std::vector<uint8_t> data;

    // write block header and data
    EResult write(FILE& file, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);

    static size_t get_parameters_size() { return sizeof(format) + sizeof(width) + sizeof(height); }

private:
    void update_checksum(Checksum& checksum) const;
};

struct PrinterMetadataBlock : public BaseMetadataBlock
{
    // write block header and data
    EResult write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);
};

struct PrintMetadataBlock : public BaseMetadataBlock
{
    // write block header and data
    EResult write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);
};

struct SlicerMetadataBlock : public BaseMetadataBlock
{
    // write block header and data
    EResult write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);
};

enum class EGCodeEncodingType : uint16_t
{
    None,
    MeatPack,
    MeatPackComments
};

struct GCodeBlock
{
    uint16_t encoding_type{ 0 };
    std::string raw_data;

    // write block header and data
    EResult write(FILE& file, ECompressionType compression_type, EChecksumType checksum_type) const;
    // read block data
    EResult read_data(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);

    static size_t get_parameters_size() { return sizeof(encoding_type); }
};

#if ENABLE_CHECKSUM_BLOCK
struct ChecksumBlock
{
    std::vector<uint8_t> data;

    // write block header and data
    EResult write(FILE& file) const;
    // read block data
    EResult read_data(FILE& file, const BlockHeader& block_header);
};
#endif // ENABLE_CHECKSUM_BLOCK

//=====================================================================================================================================
//
// PRUSASLICER INTERFACE
//  
//=====================================================================================================================================

struct BinaryData
{
    FileMetadataBlock file_metadata;
    PrinterMetadataBlock printer_metadata;
    std::vector<ThumbnailBlock> thumbnails;
    SlicerMetadataBlock slicer_metadata;
    PrintMetadataBlock print_metadata;

    void reset() {
        file_metadata.raw_data.clear();
        printer_metadata.raw_data.clear();
        thumbnails.clear();
        slicer_metadata.raw_data.clear();
        print_metadata.raw_data.clear();
    }
};

struct BinarizerConfig
{
    struct Compression
    {
        ECompressionType file_metadata{ ECompressionType::None };
        ECompressionType printer_metadata{ ECompressionType::None };
        ECompressionType thumbnail{ ECompressionType::None };
        ECompressionType print_metadata{ ECompressionType::None };
        ECompressionType slicer_metadata{ ECompressionType::None };
        ECompressionType gcode{ ECompressionType::None };
    };
    Compression compression;
    EGCodeEncodingType gcode_encoding{ EGCodeEncodingType::None };
    EMetadataEncodingType metadata_encoding{ EMetadataEncodingType::INI };
    EChecksumType checksum{ EChecksumType::CRC32 };
};

class Binarizer
{
public:
    bool is_enabled() const { return m_enabled; }
    void set_enabled(bool enable) { m_enabled = enable; }

    BinaryData& get_binary_data() { return m_binary_data; }
    const BinaryData& get_binary_data() const { return m_binary_data; }

    EResult initialize(FILE& file, const BinarizerConfig& config);
    EResult append_gcode(const std::string& gcode);
    EResult finalize();

private:
    bool m_enabled{ false };

    BinarizerConfig m_config;
    FILE* m_file{ nullptr };
    BinaryData m_binary_data;
    std::string m_gcode_cache;
#if ENABLE_CHECKSUM_BLOCK
    ChecksumBlock m_checksum;
#endif // ENABLE_CHECKSUM_BLOCK
};

//=====================================================================================================================================
//
// FIRMWARE INTERFACE
//  
//=====================================================================================================================================

// Get the max size of the cache used to calculate checksums, in bytes
size_t get_checksum_max_cache_size();
// Set the max size of the cache used to calculate checksums, in bytes
void set_checksum_max_cache_size(size_t size);

// Returns true if the given file is a valid binary gcode
// Does not modify the file position
extern bool is_valid_binary_gcode(FILE& file);

// Reads the file header.
// If max_version is not null, version is checked against the passed value
// If return == EResult::Success:
// - header will contain the file header
// - file position will be set at the start of the 1st block header
extern EResult read_header(FILE& file, FileHeader& header, const uint32_t* const max_version);

// Reads next block header from the current file position.
// File position must be at the start of a block header.
// If return == EResult::Success:
// - block_header will contain the header of the block
// - file position will be set at the start of the block parameters data
extern EResult read_next_block_header(FILE& file, const FileHeader& file_header, BlockHeader& block_header, bool verify_checksum);

// Searches and reads next block header with the given type from the current file position.
// File position must be at the start of a block header.
// If return == EResult::Success:
// - block_header will contain the header of the block with the required type
// - file position will be set at the start of the block parameters data
// otherwise:
// - file position will keep the current value
extern EResult read_next_block_header(FILE& file, const FileHeader& file_header, BlockHeader& block_header, EBlockType type, bool verify_checksum);

// Skips the payload (parameters + data) of the block with the given block header.
// File position must be at the start of the block parameters.
// If return == EResult::Success:
// - file position will be set at the start of the block checksum, if present, or of next block header
extern EResult skip_block_payload(FILE& file, const BlockHeader& block_header);

// Skips the content (parameters + data + checksum) of the block with the given block header.
// File position must be at the start of the block parameters.
// If return == EResult::Success:
// - file position will be set at the start of the next block header
extern EResult skip_block_content(FILE& file, const FileHeader& file_header, const BlockHeader& block_header);

// Returns the size of the parameters of the given block type, in bytes.
extern size_t block_parameters_size(EBlockType type);

// Returns the size of the payload (parameters + data) of the block with the given header, in bytes.
extern size_t block_payload_size(const BlockHeader& block_header);

// Returns the size of the checksum of the given type, in bytes.
extern size_t checksum_size(EChecksumType type);

// Returns the size of the content (parameters + data + checksum) of the block with the given header, in bytes.
extern size_t block_content_size(const FileHeader& file_header, const BlockHeader& block_header);

} // namespace BinaryGCode

#endif // slic3r_GCode_GCodeBinarizer_hpp_
