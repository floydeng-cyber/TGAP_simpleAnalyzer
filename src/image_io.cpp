#include "image_io.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <windows.h>
#include <commdlg.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

namespace fs = std::filesystem;

namespace {

std::runtime_error tiffError(const std::string& msg) { return std::runtime_error("TIFF: " + msg); }

// TIFF LZW: MSB-first codes of 9..12 bits; the code width grows one entry early.
std::vector<uint8_t> lzwDecode(const uint8_t* src, size_t n, size_t expected) {
    std::vector<std::vector<uint8_t>> table;
    std::vector<uint8_t> out;
    size_t bit = 0;
    int width = 9, prev = -1;
    while (out.size() < expected && bit + width <= n * 8) {
        int code = 0;
        for (int i = 0; i < width; ++i, ++bit) code = code << 1 | ((src[bit >> 3] >> (7 - (bit & 7))) & 1);
        if (code == 257) break;
        if (code == 256 || table.empty()) {  // clear code (the stream always starts with one)
            table.assign(258, {});
            for (int i = 0; i < 256; ++i) table[i] = {static_cast<uint8_t>(i)};
            width = 9;
            prev = -1;
            if (code == 256) continue;
        }
        std::vector<uint8_t> entry;
        if (code < static_cast<int>(table.size())) entry = table[code];
        else if (prev >= 0 && code == static_cast<int>(table.size())) entry = table[prev], entry.push_back(table[prev][0]);
        else throw tiffError("corrupt LZW data");
        if (prev >= 0 && table.size() < 4096) {
            table.push_back(table[prev]);
            table.back().push_back(entry[0]);
        }
        out.insert(out.end(), entry.begin(), entry.end());
        prev = code;
        if (table.size() + 1 >= (1u << width) && width < 12) ++width;
    }
    return out;
}

void putU32(std::vector<uint8_t>& out, uint32_t v) {
    for (int s = 24; s >= 0; s -= 8) out.push_back(static_cast<uint8_t>(v >> s));
}

uint32_t crc32(const uint8_t* p, size_t n) {
    uint32_t crc = ~0u;
    for (size_t i = 0; i < n; ++i) {
        crc ^= p[i];
        for (int k = 0; k < 8; ++k) crc = crc & 1 ? 0xEDB88320u ^ (crc >> 1) : crc >> 1;
    }
    return ~crc;
}

std::string xmlEscape(const std::string& s) {
    std::string out;
    for (char c : s) out += c == '&' ? "&amp;" : c == '<' ? "&lt;" : c == '>' ? "&gt;" : c == '"' ? "&quot;" : std::string(1, c);
    return out;
}

}  // namespace

Gray readTiff(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open " + path.string());
    const std::vector<uint8_t> b((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    if (b.size() < 8 || !((b[0] == 'I' && b[1] == 'I') || (b[0] == 'M' && b[1] == 'M'))) throw tiffError("not a TIFF file");
    const bool le = b[0] == 'I';
    auto u = [&](size_t off, int bytes) -> uint32_t {
        if (off + bytes > b.size()) throw tiffError("unexpected end of file");
        uint32_t v = 0;
        for (int i = 0; i < bytes; ++i) v |= static_cast<uint32_t>(b[off + i]) << 8 * (le ? i : bytes - 1 - i);
        return v;
    };
    if (u(2, 2) != 42) throw tiffError("not a classic TIFF (BigTIFF is not supported)");

    std::map<uint16_t, std::vector<uint32_t>> tags;  // SHORT and LONG tags only
    const size_t ifd = u(4, 4);
    for (uint32_t i = 0, n = u(ifd, 2); i < n; ++i) {
        const size_t e = ifd + 2 + 12 * i;
        const uint32_t type = u(e + 2, 2), count = u(e + 4, 4), size = type == 3 ? 2 : type == 4 ? 4 : 0;
        if (!size) continue;
        const size_t at = size * count <= 4 ? e + 8 : u(e + 8, 4);
        for (uint32_t k = 0; k < count; ++k) tags[u(e, 2)].push_back(u(at + size * k, size));
    }
    auto tag = [&](uint16_t id, uint32_t fallback) { return tags.count(id) ? tags[id][0] : fallback; };
    const uint32_t w = tag(256, 0), h = tag(257, 0), spp = tag(277, 1), bps = tag(258, 1), compression = tag(259, 1),
                   photometric = tag(262, 1), predictor = tag(317, 1);
    if (!w || !h || !tags.count(273) || tags.count(322)) throw tiffError("missing strips (tiled images are not supported)");
    if ((bps != 8 && bps != 16) || tag(339, 1) != 1) throw tiffError("only 8/16-bit unsigned samples are supported");
    if (compression != 1 && compression != 5) throw tiffError("only uncompressed or LZW images are supported");
    if (photometric > 2 || (photometric == 2 && spp < 3) || (spp > 1 && tag(284, 1) != 1))
        throw tiffError("only gray or interleaved RGB images are supported");

    const size_t bytesPer = bps / 8, rowBytes = w * spp * bytesPer, total = rowBytes * h;
    const size_t stripBytes = std::min(tag(278, h), h) * rowBytes;
    std::vector<uint8_t> raw;
    for (size_t s = 0; s < tags[273].size() && raw.size() < total; ++s) {
        const size_t off = tags[273][s], cnt = tags[279].at(s), want = std::min(stripBytes, total - raw.size());
        if (off + cnt > b.size()) throw tiffError("strip outside the file");
        std::vector<uint8_t> strip = compression == 5 ? lzwDecode(&b[off], cnt, want)
                                                       : std::vector<uint8_t>(b.begin() + off, b.begin() + off + std::min(cnt, want));
        if (strip.size() < want) throw tiffError("truncated strip");
        raw.insert(raw.end(), strip.begin(), strip.begin() + want);
    }
    if (raw.size() < total) throw tiffError("truncated image");

    auto sample = [&](size_t i) -> uint32_t {
        return bytesPer == 1 ? raw[i] : le ? raw[2 * i] | raw[2 * i + 1] << 8 : raw[2 * i] << 8 | raw[2 * i + 1];
    };
    const uint32_t maxValue = (1u << bps) - 1;
    if (predictor == 2) {  // horizontal differencing
        for (size_t i = 0; i < size_t(w) * h * spp; ++i) {
            if (i % (w * spp) < spp) continue;
            const uint32_t v = (sample(i) + sample(i - spp)) & maxValue;
            if (bytesPer == 1) raw[i] = static_cast<uint8_t>(v);
            else raw[2 * i + !le] = static_cast<uint8_t>(v), raw[2 * i + le] = static_cast<uint8_t>(v >> 8);
        }
    }

    Gray img{static_cast<int>(w), static_cast<int>(h), std::vector<float>(size_t(w) * h)};
    for (size_t p = 0; p < img.px.size(); ++p) {
        const size_t i = p * spp;
        const float v = photometric == 2 ? 0.299f * sample(i) + 0.587f * sample(i + 1) + 0.114f * sample(i + 2) : sample(i);
        img.px[p] = (photometric == 0 ? maxValue - v : v) * 255.0f / maxValue;
    }
    return img;
}

std::vector<uint8_t> encodePng(int w, int h, const std::vector<uint8_t>& rgb) {
    std::vector<uint8_t> raw;
    for (int y = 0; y < h; ++y) {
        raw.push_back(0);  // filter type: none
        raw.insert(raw.end(), rgb.begin() + size_t(y) * w * 3, rgb.begin() + size_t(y + 1) * w * 3);
    }
    std::vector<uint8_t> z = {0x78, 0x01};
    uint32_t a = 1, s = 0;
    for (size_t pos = 0; pos < raw.size(); pos += 65535) {
        const uint16_t len = static_cast<uint16_t>(std::min<size_t>(65535, raw.size() - pos));
        z.insert(z.end(), {static_cast<uint8_t>(pos + len == raw.size()), static_cast<uint8_t>(len), static_cast<uint8_t>(len >> 8),
                           static_cast<uint8_t>(~len), static_cast<uint8_t>(~len >> 8)});
        z.insert(z.end(), raw.begin() + pos, raw.begin() + pos + len);
    }
    for (uint8_t v : raw) a = (a + v) % 65521, s = (s + a) % 65521;
    putU32(z, s << 16 | a);

    std::vector<uint8_t> png = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1a, '\n'};
    auto chunk = [&](const char* type, std::vector<uint8_t> data) {
        std::vector<uint8_t> c;
        putU32(c, static_cast<uint32_t>(data.size()));
        c.insert(c.end(), type, type + 4);
        c.insert(c.end(), data.begin(), data.end());
        putU32(c, crc32(c.data() + 4, c.size() - 4));
        png.insert(png.end(), c.begin(), c.end());
    };
    std::vector<uint8_t> header;
    putU32(header, w);
    putU32(header, h);
    header.insert(header.end(), {8, 2, 0, 0, 0});  // 8-bit RGB
    chunk("IHDR", header);
    chunk("IDAT", z);
    chunk("IEND", {});
    return png;
}

void writeXlsx(const fs::path& path, const std::vector<Sheet>& sheets) {
    // The workbook's parts, as XML.
    const std::string head = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
    const std::string main = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
    const std::string rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";
    std::vector<std::pair<std::string, std::string>> parts;
    std::string types = head + "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
                        "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
                        "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
                        "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>";
    std::string book = head + "<workbook xmlns=\"" + main + "\" xmlns:r=\"" + rel + "\"><sheets>";
    std::string bookRels = head + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">";
    std::set<std::string> used;  // sheet names so far, upper case
    for (size_t i = 0; i < sheets.size(); ++i) {
        const std::string n = std::to_string(i + 1);
        // Excel sheet names: at most 31 characters, none of :\/?*[], unique ignoring case. A longer name
        // keeps its start and end (which tell folders apart) around a '~'; a repeat gets "~2", "~3", ...
        std::string base = sheets[i].name, name;
        for (char& c : base) c = std::string(":\\/?*[]").find(c) == std::string::npos ? c : '_';
        if (base.empty()) base = "Sheet" + n;
        for (int repeat = 1;; ++repeat) {
            const std::string suffix = repeat > 1 ? "~" + std::to_string(repeat) : "";
            const size_t max = 31 - suffix.size();
            name = (base.size() <= max ? base : base.substr(0, max / 2) + "~" + base.substr(base.size() - (max - 1 - max / 2))) + suffix;
            std::string key = name;
            for (char& c : key) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (used.insert(key).second) break;
        }
        types += "<Override PartName=\"/xl/worksheets/sheet" + n + ".xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>";
        book += "<sheet name=\"" + xmlEscape(name) + "\" sheetId=\"" + n + "\" r:id=\"rId" + n + "\"/>";
        bookRels += "<Relationship Id=\"rId" + n + "\" Type=\"" + rel + "/worksheet\" Target=\"worksheets/sheet" + n + ".xml\"/>";
        std::string xml = head + "<worksheet xmlns=\"" + main + "\" xmlns:r=\"" + rel + "\"><cols>";
        std::vector<size_t> widths;  // each column as wide as its longest entry
        for (const auto& row : sheets[i].rows)
            for (size_t c = 0; c < row.size(); ++c) widths.resize(std::max(widths.size(), c + 1)), widths[c] = std::max(widths[c], row[c].size());
        for (size_t c = 0; c < widths.size(); ++c)
            xml += "<col min=\"" + std::to_string(c + 1) + "\" max=\"" + std::to_string(c + 1) + "\" width=\"" + std::to_string(widths[c] + 2) + "\" customWidth=\"1\"/>";
        xml += "</cols><sheetData>";
        std::string hyperlinks, links;  // the sheet's <hyperlinks> and its relationships (their targets)
        int linkCount = 0;
        for (size_t r = 0; r < sheets[i].rows.size(); ++r) {
            xml += "<row r=\"" + std::to_string(r + 1) + "\">";
            for (size_t c = 0; c < sheets[i].rows[r].size(); ++c) {
                const std::string& v = sheets[i].rows[r][c];
                const std::string ref = std::string(c >= 26 ? 1 : 0, char('A' + c / 26 - 1)) + char('A' + c % 26) + std::to_string(r + 1);
                char* end = nullptr;
                std::strtod(v.c_str(), &end);
                const bool number = !v.empty() && !*end && v.find_first_not_of("0123456789.-") == std::string::npos;  // plain decimal
                const auto link = sheets[i].links.find({r, c});
                const std::string style = link != sheets[i].links.end() ? " s=\"1\"" : "";  // the hyperlink style
                xml += !number ? "<c r=\"" + ref + "\"" + style + " t=\"inlineStr\"><is><t>" + xmlEscape(v) + "</t></is></c>"  // text
                            : "<c r=\"" + ref + "\"" + style + "><v>" + v + "</v></c>";                                    // a number
                if (link != sheets[i].links.end()) {
                    const std::string id = "rId" + std::to_string(++linkCount);
                    hyperlinks += "<hyperlink ref=\"" + ref + "\" r:id=\"" + id + "\"/>";
                    links += "<Relationship Id=\"" + id + "\" Type=\"" + rel + "/hyperlink\" Target=\"" + xmlEscape(link->second) + "\" TargetMode=\"External\"/>";
                }
            }
            xml += "</row>";
        }
        xml += "</sheetData>";
        if (linkCount) {
            xml += "<hyperlinks>" + hyperlinks + "</hyperlinks>";
            parts.push_back({"xl/worksheets/_rels/sheet" + n + ".xml.rels",
                             head + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">" + links + "</Relationships>"});
        }
        parts.push_back({"xl/worksheets/sheet" + n + ".xml", xml + "</worksheet>"});
    }
    // Two cell styles: 0 plain, 1 a hyperlink (blue, underlined).
    types += "<Override PartName=\"/xl/styles.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>";
    bookRels += "<Relationship Id=\"rIdStyles\" Type=\"" + rel + "/styles\" Target=\"styles.xml\"/>";
    parts.push_back({"xl/styles.xml", head + "<styleSheet xmlns=\"" + main + "\"><fonts count=\"2\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font>"
                                       "<font><u/><sz val=\"11\"/><color rgb=\"FF0563C1\"/><name val=\"Calibri\"/></font></fonts>"
                                       "<fills count=\"2\"><fill><patternFill patternType=\"none\"/></fill><fill><patternFill patternType=\"gray125\"/></fill></fills>"
                                       "<borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders>"
                                       "<cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>"
                                       "<cellXfs count=\"2\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>"
                                       "<xf numFmtId=\"0\" fontId=\"1\" fillId=\"0\" borderId=\"0\" xfId=\"0\" applyFont=\"1\"/></cellXfs></styleSheet>"});
    parts.push_back({"[Content_Types].xml", types + "</Types>"});
    parts.push_back({"_rels/.rels", head + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
                                           "<Relationship Id=\"rId1\" Type=\"" + rel + "/officeDocument\" Target=\"xl/workbook.xml\"/></Relationships>"});
    parts.push_back({"xl/workbook.xml", book + "</sheets></workbook>"});
    parts.push_back({"xl/_rels/workbook.xml.rels", bookRels + "</Relationships>"});

    // A ZIP archive of the parts, stored (uncompressed): local headers with the data, then the central directory.
    std::vector<uint8_t> zip, dir;
    auto put = [](std::vector<uint8_t>& out, uint32_t v, int bytes) {  // little-endian
        for (int b = 0; b < bytes; ++b) out.push_back(static_cast<uint8_t>(v >> (8 * b)));
    };
    for (const auto& [name, data] : parts) {
        const uint32_t crc = crc32(reinterpret_cast<const uint8_t*>(data.data()), data.size()), size = static_cast<uint32_t>(data.size());
        const uint32_t offset = static_cast<uint32_t>(zip.size());
        for (std::vector<uint8_t>* out : {&zip, &dir}) {
            const bool central = out == &dir;
            put(*out, central ? 0x02014b50 : 0x04034b50, 4);
            if (central) put(*out, 20, 2);                      // made by
            put(*out, 20, 2), put(*out, 0, 2), put(*out, 0, 2);  // needed version, flags, stored
            put(*out, 0, 2), put(*out, 0x21, 2);                // time, date (1980-01-01)
            put(*out, crc, 4), put(*out, size, 4), put(*out, size, 4);
            put(*out, static_cast<uint32_t>(name.size()), 2), put(*out, 0, 2);  // name length, extra length
            if (central) put(*out, 0, 2), put(*out, 0, 2), put(*out, 0, 2), put(*out, 0, 4), put(*out, offset, 4);
            out->insert(out->end(), name.begin(), name.end());
        }
        zip.insert(zip.end(), data.begin(), data.end());
    }
    const uint32_t dirOffset = static_cast<uint32_t>(zip.size()), count = static_cast<uint32_t>(parts.size());
    zip.insert(zip.end(), dir.begin(), dir.end());
    put(zip, 0x06054b50, 4), put(zip, 0, 2), put(zip, 0, 2), put(zip, count, 2), put(zip, count, 2);
    put(zip, static_cast<uint32_t>(dir.size()), 4), put(zip, dirOffset, 4), put(zip, 0, 2);
    std::ofstream f(path, std::ios::binary);
    if (!f.write(reinterpret_cast<const char*>(zip.data()), static_cast<std::streamsize>(zip.size()))) throw std::runtime_error("cannot write " + path.string());
}

fs::path executablePath() {
#ifdef _WIN32
    std::wstring buf(32768, L'\0');
    buf.resize(GetModuleFileNameW(nullptr, buf.data(), static_cast<DWORD>(buf.size())));
    return fs::path(buf);
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::string buf(size, '\0');
    if (_NSGetExecutablePath(buf.data(), &size) != 0) return {};
    std::error_code ec;
    const fs::path p = fs::canonical(fs::path(buf.c_str()), ec);  // resolves symlinks and "./"
    return ec ? fs::path(buf.c_str()) : p;
#else
    std::error_code ec;
    const fs::path p = fs::read_symlink("/proc/self/exe", ec);
    return ec ? fs::path() : p;
#endif
}

#ifdef _WIN32
std::vector<fs::path> pickFiles(const fs::path& startDir) {
    std::vector<wchar_t> buf(1 << 16, L'\0');
    const std::wstring dir = startDir.wstring();
    OPENFILENAMEW ofn{};
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = L"TIFF images\0*.tif;*.tiff\0All files\0*.*\0";
    ofn.lpstrInitialDir = dir.c_str();
    ofn.lpstrFile = buf.data();
    ofn.nMaxFile = static_cast<DWORD>(buf.size());
    ofn.Flags = OFN_ALLOWMULTISELECT | OFN_EXPLORER | OFN_FILEMUSTEXIST;
    if (!GetOpenFileNameW(&ofn)) return {};
    // One file: a full path. Several: the folder, then each name; NUL-separated.
    std::vector<std::wstring> parts;
    for (const wchar_t* p = buf.data(); *p; p += parts.back().size() + 1) parts.emplace_back(p);
    if (parts.size() == 1) return {parts[0]};
    std::vector<fs::path> files;
    for (size_t i = 1; i < parts.size(); ++i) files.push_back(fs::path(parts[0]) / parts[i]);
    return files;
}
#else
std::vector<fs::path> pickFiles(const fs::path& startDir) {
    auto quote = [](const std::string& s) {  // single-quote for the shell
        std::string q = "'";
        for (char c : s) q += c == '\'' ? std::string("'\\''") : std::string(1, c);
        return q + "'";
    };
#ifdef __APPLE__
    // NSOpenPanel (via JavaScript for Automation) is the macOS dialog that accepts files and folders together.
    const std::string command = "osascript -l JavaScript -e " + quote(
        "function run(argv) { ObjC.import('Cocoa');"
        // osascript must act as an (accessory) app, or the panel cancels itself as soon as it is clicked.
        " var app = $.NSApplication.sharedApplication;"
        " app.setActivationPolicy($.NSApplicationActivationPolicyAccessory); app.activateIgnoringOtherApps(true);"
        " var p = $.NSOpenPanel.openPanel;"
        " p.canChooseFiles = true; p.canChooseDirectories = true; p.allowsMultipleSelection = true;"
        " p.allowedFileTypes = $(['tif', 'tiff']);"
        " p.directoryURL = $.NSURL.fileURLWithPath(argv[0]);"
        " p.message = 'Select TIFF images and/or folders';"
        " p.level = $.NSModalPanelWindowLevel;"  // float above Terminal, which otherwise hides it
        " if (p.runModal != 1) return '';"
        " var out = ''; for (var i = 0; i < p.URLs.count; i++) out += p.URLs.objectAtIndex(i).path.js + '\\n';"
        " return out; }") + " " + quote(startDir.string()) + " 2>/dev/null";
#else
    const std::string command = "zenity --file-selection --multiple --separator='\n' --file-filter='*.tif *.tiff' --filename=" +
                                quote(startDir.string() + "/") + " 2>/dev/null";
#endif
    std::vector<fs::path> files;
    if (FILE* pipe = popen(command.c_str(), "r")) {
        char line[4096];
        while (fgets(line, sizeof(line), pipe)) {
            std::string s(line);
            while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
            if (!s.empty()) files.emplace_back(s);
        }
        pclose(pipe);
    }
    return files;
}
#endif
