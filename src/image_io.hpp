#pragma once

#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

struct Gray {
    int w = 0, h = 0;
    std::vector<float> px;  // row-major, 0..255 scale whatever the source bit depth
};

// First page of a TIFF: 8/16-bit gray or RGB(A), uncompressed or LZW strips. Throws on anything else.
Gray readTiff(const std::filesystem::path& path);
// 8-bit RGB PNG bytes, uncompressed deflate blocks (no zlib needed).
std::vector<uint8_t> encodePng(int w, int h, const std::vector<uint8_t>& rgb);
// An Excel workbook (.xlsx) with one sheet per Sheet; a cell that reads as a number is stored as one, any
// other as text; a cell with a link is a hyperlink, styled blue and underlined. Sheet names are made valid
// for Excel (at most 31 characters, unique). Written as an uncompressed ZIP archive (no zip library needed).
struct Sheet {
    std::string name;
    std::vector<std::vector<std::string>> rows;
    std::map<std::pair<size_t, size_t>, std::string> links;  // (row, column) -> URL of a hyperlink in that cell
};
void writeXlsx(const std::filesystem::path& path, const std::vector<Sheet>& sheets);
// The running program's own file (absolute), or empty if the system cannot tell.
std::filesystem::path executablePath();
// Native multi-select "open" dialog starting in `startDir`; empty if cancelled or unavailable.
// macOS lets you pick TIFF files and folders together; Windows and Linux dialogs pick files only.
std::vector<std::filesystem::path> pickFiles(const std::filesystem::path& startDir);
