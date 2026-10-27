// trace_boxes: finds the horizontal and vertical trace in TIFF images, boxes each (axis-aligned), flags a
// neck connecting them, and plots each trace's cross-section. It writes <name>_traces.svg (the annotated
// image: boxes, long axes, ids, a neck's box and narrowest cut in red; next to one cross-section subplot
// per trace and one for the neck) for every image with a neck. The measurements of every image without a
// neck go to analysisResult.xlsx, one sheet per folder; images with a neck or without both traces are listed
// on its last sheet, "Special cases". similarityResult_<type>.svg plots each image's HT and VT similarity
// factors (width x length x depth).

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "traces.hpp"

namespace fs = std::filesystem;

namespace {

const char* kColors[] = {"#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6"};
const char* color(int id) { return kColors[(id - 1) % 7]; }
const char* kNeckColor = "#e34948";
constexpr int kMaxImageHeight = 800;  // the embedded image is downsampled to at most this height

// The image with boxes, dashed long axes, ids and the neck's narrowest cut (red), as PNG bytes.
std::vector<uint8_t> annotate(const Gray& img, const std::vector<Trace>& traces, int& outW, int& outH) {
    std::vector<uint8_t> rgb(size_t(img.w) * img.h * 3);
    for (size_t i = 0; i < img.px.size(); ++i) rgb[3 * i] = rgb[3 * i + 1] = rgb[3 * i + 2] = uint8_t(std::clamp(img.px[i], 0.f, 255.f));
    const int t = std::max(1, std::min(img.w, img.h) / 350);  // line thickness
    auto dot = [&](double x, double y, int size, const char* hex) {
        const unsigned c = std::stoul(hex + 1, nullptr, 16);
        for (int dy = 0; dy < size; ++dy)
            for (int dx = 0; dx < size; ++dx) {
                const int px = int(x) + dx - size / 2, py = int(y) + dy - size / 2;
                if (px < 0 || py < 0 || px >= img.w || py >= img.h) continue;
                uint8_t* p = &rgb[(size_t(py) * img.w + px) * 3];
                p[0] = uint8_t(c >> 16), p[1] = uint8_t(c >> 8), p[2] = uint8_t(c);
            }
    };
    auto line = [&](double x0, double y0, double x1, double y1, int size, const char* hex, bool dashed = false) {
        const int steps = std::max(1, int(std::hypot(x1 - x0, y1 - y0)));
        for (int i = 0; i <= steps; ++i)
            if (!dashed || (i / (4 * t)) % 2 == 0) dot(x0 + (x1 - x0) * i / steps, y0 + (y1 - y0) * i / steps, size, hex);
    };
    static const unsigned short digits[] = {075557, 026227, 071747, 071717, 055711, 074717, 074757, 071111, 075757, 075717};
    for (const Trace& tr : traces) {
        const double xs[] = {tr.x0(), tr.x1(), tr.x1(), tr.x0()}, ys[] = {tr.y0(), tr.y0(), tr.y1(), tr.y1()};
        const double cx = (tr.x0() + tr.x1()) / 2, cy = (tr.y0() + tr.y1()) / 2;
        for (int i = 0; i < 4; ++i) line(xs[i], ys[i], xs[(i + 1) % 4], ys[(i + 1) % 4], t, color(tr.id));
        if (tr.vertical) line(cx, tr.y0(), cx, tr.y1(), t, color(tr.id), true);
        else line(tr.x0(), cy, tr.x1(), cy, t, color(tr.id), true);
        // Id in a 3x5 pixel font on a black tile at the box centre.
        const std::string id = std::to_string(tr.id);
        const int s = 3 * t, x0 = int(cx) - int(id.size()) * 2 * s, y0 = int(cy) - 3 * s;
        for (int y = -1; y < 6; ++y)
            for (int x = -1; x < int(id.size()) * 4; ++x) dot(x0 + x * s, y0 + y * s, s, "#000000");
        for (size_t c = 0; c < id.size(); ++c)
            for (int row = 0; row < 5; ++row)
                for (int col = 0; col < 3; ++col)
                    if (digits[id[c] - '0'] >> ((4 - row) * 3 + 2 - col) & 1) dot(x0 + (4 * int(c) + col) * s + s / 2, y0 + row * s + s / 2, s, color(tr.id));
        if (tr.connected && tr.neckWidth > 0) line(tr.neckAt, tr.neckFrom, tr.neckAt, tr.neckTo, 3 * t, kNeckColor);  // the neck's narrowest cut
        if (tr.neckX1 >= tr.neckX0 && tr.neckY1 >= tr.neckY0) {  // the neck's box, dashed
            const double nx[] = {tr.neckX0, tr.neckX1, tr.neckX1, tr.neckX0}, ny[] = {tr.neckY0, tr.neckY0, tr.neckY1, tr.neckY1};
            for (int i = 0; i < 4; ++i) line(nx[i], ny[i], nx[(i + 1) % 4], ny[(i + 1) % 4], 2 * t, kNeckColor, true);
        }
    }

    // Box-average down to at most kMaxImageHeight rows to keep the SVG small.
    const int f = (img.h + kMaxImageHeight - 1) / kMaxImageHeight;
    outW = img.w / f, outH = img.h / f;
    std::vector<uint8_t> small(size_t(outW) * outH * 3);
    for (int y = 0; y < outH * f; ++y)
        for (int x = 0; x < outW * f; ++x)
            for (int c = 0; c < 3; ++c) small[(size_t(y / f) * outW + x / f) * 3 + c] += uint8_t(rgb[(size_t(y) * img.w + x) * 3 + c] / (f * f));
    return encodePng(outW, outH, small);
}

std::string base64(const std::vector<uint8_t>& d) {
    static const char* t = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    for (size_t i = 0; i < d.size(); i += 3) {
        const uint32_t v = d[i] << 16 | (i + 1 < d.size() ? d[i + 1] << 8 : 0) | (i + 2 < d.size() ? d[i + 2] : 0);
        out += {t[v >> 18 & 63], t[v >> 12 & 63], i + 1 < d.size() ? t[v >> 6 & 63] : '=', i + 2 < d.size() ? t[v & 63] : '='};
    }
    return out;
}

double niceStep(double range) {
    const double raw = range / 5, mag = std::pow(10, std::floor(std::log10(raw)));
    return raw <= mag ? mag : raw <= 2 * mag ? 2 * mag : raw <= 5 * mag ? 5 * mag : 10 * mag;
}

// One SVG: annotated image on the left, a grid of cross-section subplots on the right: one per trace,
// then one for the neck if there is one. Dotted lines and dots mark each profile's two edges (ridges),
// which are its boundaries (x = 0 and x = width).
void writeSvg(const fs::path& path, const std::string& name, const Gray& img, const std::vector<Trace>& traces) {
    struct Plot {
        std::string title;
        const char* color;
        const std::vector<double>* profile;
        double from, edge0, edge1;  // profile[i] is at position from + i; its edges
    };
    std::vector<Plot> plots;
    char buf[512];
    for (const Trace& t : traces) {
        std::snprintf(buf, sizeof buf, "Trace %d (%s): width %.0f px, length %.0f px", t.id, t.vertical ? "vertical" : "horizontal", t.width(), t.length());
        plots.push_back({buf, color(t.id), &t.crossSection, t.crossFrom, t.across0, t.across1});
    }
    for (const Trace& t : traces)
        if (t.connected && !t.neckProfile.empty()) {
            std::snprintf(buf, sizeof buf, "Neck (trace 1 into 2): narrowest width %.0f px (red line)", t.neckWidth);
            plots.push_back({buf, kNeckColor, &t.neckProfile, t.neckFrom0, t.neckEdge0, t.neckEdge1});
        }

    int iw, ih;
    const std::string png = base64(annotate(img, traces, iw, ih));
    const double imgH = 520, imgW = imgH * iw / ih, cw = 330, ch = 260, L = 52, R = 16, T = 54, B = 40, top = 34;
    const int n = int(plots.size()), cols = n > 1 ? 2 : 1, rows = std::max(1, (n + cols - 1) / cols);
    double lo = 255, hi = 0;  // shared intensity axis so profiles can be compared
    for (const Plot& p : plots)
        for (double v : *p.profile)
            if (!std::isnan(v)) lo = std::min(lo, v), hi = std::max(hi, v);
    if (hi <= lo) lo = 0, hi = 255;
    const double yStep = niceStep(hi - lo);
    lo = std::floor(lo / yStep) * yStep, hi = std::ceil(hi / yStep) * yStep;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path.string());
    std::snprintf(buf, sizeof buf, "<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' width='%g' height='%g' "
                  "font-family='Helvetica,Arial,sans-serif' font-size='11'>\n<rect width='100%%' height='100%%' fill='white'/>\n"
                  "<text x='12' y='22' font-size='14'>%s: %d trace(s); cross-sections are averaged along each trace (0-255)</text>\n",
                  12 + imgW + 12 + cols * cw, top + std::max(imgH + 12, rows * ch), name.c_str(), int(traces.size()));
    f << buf << "<image x='12' y='" << top + 10 << "' width='" << imgW << "' height='" << imgH << "' xlink:href='data:image/png;base64," << png << "'/>\n";
    for (int i = 0; i < n; ++i) {
        const Plot& p = plots[i];
        const std::vector<double>& prof = *p.profile;
        const double width = p.edge1 - p.edge0, x0 = p.from - p.edge0, x1 = p.from + prof.size() - 1 - p.edge0;
        const double ox = 24 + imgW + (i % cols) * cw + L, oy = top + (i / cols) * ch + T, pw = cw - L - R, ph = ch - T - B;
        auto X = [&](double s) { return ox + (s - x0) / std::max(1.0, x1 - x0) * pw; };
        auto Y = [&](double v) { return oy + (hi - v) / (hi - lo) * ph; };
        std::snprintf(buf, sizeof buf, "<text x='%g' y='%g' font-size='13' fill='%s'>%s</text>\n", ox, oy - 18, p.color, p.title.c_str());
        f << buf;
        for (double v = lo; v <= hi + 1e-9; v += yStep) {
            std::snprintf(buf, sizeof buf, "<line x1='%g' x2='%g' y1='%g' y2='%g' stroke='#e1e0d9'/><text x='%g' y='%g' text-anchor='end' fill='#898781'>%g</text>\n",
                          ox, ox + pw, Y(v), Y(v), ox - 6, Y(v) + 4, v);
            f << buf;
        }
        const double xStep = niceStep(std::max(1.0, x1 - x0));
        for (double s = std::ceil(x0 / xStep) * xStep; s <= x1 + 1e-9; s += xStep) {
            std::snprintf(buf, sizeof buf, "<text x='%g' y='%g' text-anchor='middle' fill='#898781'>%g</text>\n", X(s), oy + ph + 16, s + 0.0);
            f << buf;
        }
        std::snprintf(buf, sizeof buf, "<text x='%g' y='%g' text-anchor='middle' fill='#52514e'>across (px from left/top boundary)</text>\n"
                      "<polyline fill='none' stroke='%s' stroke-width='2' points='", ox + pw / 2, oy + ph + 32, p.color);
        f << buf;
        for (size_t s = 0; s < prof.size(); ++s)
            if (!std::isnan(prof[s])) f << X(p.from + s - p.edge0) << ',' << Y(prof[s]) << ' ';
        f << "'/>\n";
        for (double edge : {0.0, width}) {
            const int s = std::clamp(int(std::lround(p.edge0 + edge - p.from)), 0, int(prof.size()) - 1);
            std::snprintf(buf, sizeof buf, "<line x1='%g' x2='%g' y1='%g' y2='%g' stroke='#52514e' stroke-dasharray='2 3'/>"
                          "<circle cx='%g' cy='%g' r='4' fill='%s' stroke='white' stroke-width='1.5'/>\n",
                          X(edge), X(edge), oy, oy + ph, X(edge), Y(prof[s]), p.color);
            f << buf;
        }
    }
    f << "</svg>\n";
}

// A trace's depth: the mean height of its two ridges in its cross-section over the lowest point between them.
double depth(const Trace& t) {
    const std::vector<double>& c = t.crossSection;
    const int last = static_cast<int>(c.size()) - 1;
    const int a = std::clamp(static_cast<int>(std::lround(t.across0 - t.crossFrom)), 0, last);
    const int b = std::clamp(static_cast<int>(std::lround(t.across1 - t.crossFrom)), a, last);
    return (c[a] + c[b]) / 2 - *std::min_element(c.begin() + a, c.begin() + b + 1);
}

std::string number(double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return buf;
}

// An image name in three parts: the prefix up to its first '_' ("E06_"), the name from there to the next
// '-' ("M0002") and the suffix, the rest ("-02MS.tif").
std::vector<std::string> splitName(const std::string& file) {
    const size_t a = file.find('_') == std::string::npos ? 0 : file.find('_') + 1;
    size_t b = file.find('-', a);
    if (b == std::string::npos) b = std::min(file.size(), file.find('.', a));
    return {file.substr(0, a), file.substr(a, b - a), file.substr(b)};
}

// Appends a row starting with an image name's three parts; the prefix and suffix are shown only where they
// differ from the row above.
void addImageRow(std::vector<std::vector<std::string>>& rows, std::vector<std::string>& last, const std::string& file,
                 std::vector<std::string> rest, std::vector<std::string> before = {}) {
    std::vector<std::string> parts = splitName(file), row = std::move(before);
    row.push_back(last.empty() || parts[0] != last[0] ? parts[0] : "");
    row.push_back(parts[1]);
    row.push_back(last.empty() || parts[2] != last[2] ? parts[2] : "");
    row.insert(row.end(), rest.begin(), rest.end());
    rows.push_back(std::move(row));
    last = std::move(parts);
}

// One image without a neck, for the type's analysis plot: its type (the selected folder's name up to its first
// '_', e.g. "TGAP"), its folder, its name and its traces' similarity factors (width x length x depth).
struct Sample {
    std::string type, folder, image;
    double ht, vt;
    std::vector<std::string> measurements;  // its row in analysisResult.xlsx, after the name
    fs::path file;                          // the raw image
};

// A file:// URL for a local path (absolute; anything but unreserved characters and '/' percent-encoded).
std::string fileUrl(const fs::path& p) {
    std::string path = fs::absolute(p).lexically_normal().generic_string(), url = "file://";
    if (path.empty() || path[0] != '/') url += '/';  // "C:/..." on Windows
    for (unsigned char c : path) {
        if (std::isalnum(c) || std::string("-._~/:").find(static_cast<char>(c)) != std::string::npos) url += static_cast<char>(c);
        else {
            char hex[4];
            std::snprintf(hex, sizeof hex, "%%%02X", c);
            url += hex;
        }
    }
    return url;
}

const std::vector<std::string> kResultHeader = {"Prefix", "Image", "Suffix", "HT width (px)", "HT length (px)", "HT depth (0-255)",
                                                 "VT width (px)", "VT length (px)", "VT depth (0-255)"};

// Mean and sample standard deviation.
std::pair<double, double> meanSd(const std::vector<double>& v) {
    double mean = 0, var = 0;
    for (double x : v) mean += x / std::max<size_t>(1, v.size());
    for (double x : v) var += (x - mean) * (x - mean) / std::max<size_t>(1, v.size() - 1);
    return {mean, std::sqrt(var)};
}

// similarityResult_<type>.svg: "<type> Analysis Result", every image of the type as a dot at its HT and VT
// similarity factors, coloured by folder. Both axes are normalized over all of the type's images: (value -
// mean) / standard deviation, so the mean is 0 and mean +- SD is +-1 (the shaded band). Each folder's mean
// +- SD is listed in a table (also the legend) with all images' raw mean and SD.
void writeSimilarity(const fs::path& path, const std::string& type, const std::vector<Sample>& samples) {
    std::vector<std::string> folders;  // in order of appearance; a folder's label is the last word of its name
    for (const Sample& s : samples)
        if (std::find(folders.begin(), folders.end(), s.folder) == folders.end()) folders.push_back(s.folder);
    auto lastWord = [](const std::string& f) {
        const std::string name = fs::path(f).filename().string();
        const size_t cut = name.find_last_of("_- ");
        return cut == std::string::npos ? name : name.substr(cut + 1);
    };
    std::set<std::string> labels;
    for (const std::string& f : folders) labels.insert(lastWord(f));
    auto label = [&](const std::string& f) { return labels.size() == folders.size() ? lastWord(f) : fs::path(f).filename().string(); };

    // Mean and SD of (x, y) over the samples of one folder (or all of them).
    struct Stats {
        double mx = 0, my = 0, sx = 0, sy = 0;
        size_t n = 0;
    };
    auto statsOf = [&](const std::vector<double>& xs, const std::vector<double>& ys, const std::string* folder) {
        Stats st;
        for (size_t i = 0; i < samples.size(); ++i)
            if (!folder || samples[i].folder == *folder) ++st.n, st.mx += xs[i], st.my += ys[i];
        st.mx /= std::max<size_t>(1, st.n), st.my /= std::max<size_t>(1, st.n);
        for (size_t i = 0; i < samples.size(); ++i)
            if (!folder || samples[i].folder == *folder) st.sx += (xs[i] - st.mx) * (xs[i] - st.mx), st.sy += (ys[i] - st.my) * (ys[i] - st.my);
        st.sx = std::sqrt(st.sx / std::max<size_t>(1, st.n - 1)), st.sy = std::sqrt(st.sy / std::max<size_t>(1, st.n - 1));
        return st;
    };
    std::vector<double> ht, vt;
    for (const Sample& s : samples) ht.push_back(s.ht), vt.push_back(s.vt);
    const Stats raw = statsOf(ht, vt, nullptr);
    std::vector<double> zx, zy;  // normalized
    for (const Sample& s : samples) zx.push_back((s.ht - raw.mx) / (raw.sx > 0 ? raw.sx : 1)), zy.push_back((s.vt - raw.my) / (raw.sy > 0 ? raw.sy : 1));

    double x0 = -1, x1 = 1, y0 = -1, y1 = 1;  // at least +-1 shown
    for (size_t i = 0; i < samples.size(); ++i) x0 = std::min(x0, zx[i]), x1 = std::max(x1, zx[i]), y0 = std::min(y0, zy[i]), y1 = std::max(y1, zy[i]);
    x0 = std::floor(x0 * 2) / 2, x1 = std::ceil(x1 * 2) / 2, y0 = std::floor(y0 * 2) / 2, y1 = std::ceil(y1 * 2) / 2;  // to 0.5
    const double L = 90, T = 60, pw = 640, ph = 480, W = L + pw + 330, H = T + ph + 60;
    auto X = [&](double v) { return L + (v - x0) / (x1 - x0) * pw; };
    auto Y = [&](double v) { return T + (y1 - v) / (y1 - y0) * ph; };
    auto big = [](double v) {  // a whole number with thousands separators, e.g. 2,000,000
        std::string digits = std::to_string(static_cast<long long>(std::llround(v))), out;
        for (size_t i = 0; i < digits.size(); ++i) out += (i && (digits.size() - i) % 3 == 0 && digits[i - 1] != '-' ? "," : "") + std::string(1, digits[i]);
        return out;
    };

    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path.string());
    char buf[1024];
    std::snprintf(buf, sizeof buf, "<svg xmlns='http://www.w3.org/2000/svg' width='%g' height='%g' font-family='Helvetica,Arial,sans-serif' font-size='12'>\n"
                  "<rect width='100%%' height='100%%' fill='white'/>\n<text x='%g' y='34' font-size='18' text-anchor='middle'>%s Analysis Result</text>\n",
                  W, H, L + pw / 2, type.c_str());
    f << buf;
    // Grid every 0.5; the mean (0) solid, mean +- SD (+-1) dashed, the band between shaded.
    std::snprintf(buf, sizeof buf, "<rect x='%.1f' y='%.1f' width='%.1f' height='%.1f' fill='#f2f1ec'/>\n", X(-1), Y(1), X(1) - X(-1), Y(-1) - Y(1));
    f << buf;
    for (double v = x0; v <= x1 + 1e-9; v += 0.5) {
        const bool unit = std::abs(std::abs(v) - 1) < 1e-9, zero = std::abs(v) < 1e-9;
        std::snprintf(buf, sizeof buf, "<line x1='%.1f' x2='%.1f' y1='%g' y2='%g' stroke='%s'%s/><text x='%.1f' y='%g' text-anchor='middle' fill='#898781'>%g</text>\n",
                      X(v), X(v), T, T + ph, zero || unit ? "#898781" : "#eeede7", unit ? " stroke-dasharray='5 4'" : "", X(v), T + ph + 16, zero ? 0.0 : v);
        f << buf;
    }
    for (double v = y0; v <= y1 + 1e-9; v += 0.5) {
        const bool unit = std::abs(std::abs(v) - 1) < 1e-9, zero = std::abs(v) < 1e-9;
        std::snprintf(buf, sizeof buf, "<line x1='%g' x2='%g' y1='%.1f' y2='%.1f' stroke='%s'%s/><text x='%g' y='%.1f' text-anchor='end' fill='#898781'>%g</text>\n",
                      L, L + pw, Y(v), Y(v), zero || unit ? "#898781" : "#eeede7", unit ? " stroke-dasharray='5 4'" : "", L - 6, Y(v) + 4, zero ? 0.0 : v);
        f << buf;
    }
    std::snprintf(buf, sizeof buf, "<rect x='%g' y='%g' width='%g' height='%g' fill='none' stroke='#c3c2b7'/>\n"
                  "<text x='%g' y='%g' text-anchor='middle' fill='#52514e' font-size='13'>HT similarity factor</text>\n"
                  "<text x='%g' y='%g' text-anchor='middle' fill='#52514e' font-size='13' transform='rotate(-90 %g %g)'>VT similarity factor</text>\n",
                  L, T, pw, ph, L + pw / 2, T + ph + 42, L - 50, T + ph / 2, L - 50, T + ph / 2);
    f << buf;
    for (size_t i = 0; i < samples.size(); ++i) {
        const Sample& s = samples[i];
        const size_t k = std::find(folders.begin(), folders.end(), s.folder) - folders.begin();
        std::snprintf(buf, sizeof buf, "<circle cx='%.1f' cy='%.1f' r='3.5' fill='%s' fill-opacity='0.6'><title>%s %s: HT %s (%+.2f), VT %s (%+.2f)</title></circle>\n",
                      X(zx[i]), Y(zy[i]), kColors[k % 7], label(s.folder).c_str(), s.image.c_str(), big(s.ht).c_str(), zx[i], big(s.vt).c_str(), zy[i]);
        f << buf;
    }

    // Each folder's normalized mean and SD, for the table.
    std::vector<Stats> per;
    for (const std::string& folder : folders) per.push_back(statsOf(zx, zy, &folder));

    // Table (also the legend): each folder's normalized mean +- SD; below, all images' raw mean and SD.
    const double tx = L + pw + 24;
    std::snprintf(buf, sizeof buf, "<text x='%g' y='%g' fill='#52514e'>Mean &#177; SD (normalized)</text>"
                  "<text x='%g' y='%g' text-anchor='end' fill='#52514e'>HT</text><text x='%g' y='%g' text-anchor='end' fill='#52514e'>VT</text>\n",
                  tx, T + 10, tx + 200, T + 10, tx + 290, T + 10);
    f << buf;
    for (size_t k = 0; k < folders.size(); ++k) {
        const double y = T + 34 + 22.0 * k;
        std::snprintf(buf, sizeof buf, "<circle cx='%g' cy='%g' r='5' fill='%s'/><text x='%g' y='%g'>%s (%zu)</text>"
                      "<text x='%g' y='%g' text-anchor='end'>%+.2f &#177; %.2f</text><text x='%g' y='%g' text-anchor='end'>%+.2f &#177; %.2f</text>\n",
                      tx + 5, y - 4, kColors[k % 7], tx + 14, y, label(folders[k]).c_str(), per[k].n, tx + 200, y, per[k].mx, per[k].sx, tx + 290, y, per[k].my, per[k].sy);
        f << buf;
    }
    const double ry = T + 34 + 22.0 * folders.size() + 20;
    std::snprintf(buf, sizeof buf, "<text x='%g' y='%g' font-weight='bold'>All %zu images: 0 &#177; 1</text>\n"
                  "<text x='%g' y='%g' fill='#52514e'>0 is the mean, &#177;1 the mean &#177; SD:</text>\n"
                  "<text x='%g' y='%g' fill='#52514e'>HT: %s &#177; %s</text>\n<text x='%g' y='%g' fill='#52514e'>VT: %s &#177; %s</text>\n",
                  tx, ry, raw.n, tx, ry + 24, tx, ry + 44, big(raw.mx).c_str(), big(raw.sx).c_str(), tx, ry + 64, big(raw.my).c_str(), big(raw.sy).c_str());
    f << buf << "</svg>\n";
}

// The folder holding the floAnalyzer folder, wherever it is: from the program's own location, up to the
// first folder holding this project (CMakeLists.txt and src/main.cpp), then its parent. If the program runs
// outside the project, the folder it runs in.
fs::path workspaceFolder() {
    const fs::path exe = executablePath();
    for (fs::path dir = exe.parent_path(); !exe.empty() && dir.has_relative_path(); dir = dir.parent_path())
        if (fs::exists(dir / "CMakeLists.txt") && fs::exists(dir / "src" / "main.cpp")) return dir.parent_path();
    return fs::current_path();
}

bool isTiff(const fs::path& p) {
    std::string ext = p.extension().string();
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext == ".tif" || ext == ".tiff";
}

}  // namespace

int main(int argc, char** argv) {
    fs::path outDir;
    std::vector<fs::path> inputs;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if ((arg == "-o" || arg == "--out") && i + 1 < argc) outDir = argv[++i];
        else if (arg[0] != '-') inputs.emplace_back(arg);
        else {
            std::cout << "Usage: trace_boxes [-o OUTPUT_DIR] [image.tif | folder ...]\n"
                         "With no inputs a dialog opens. Folders are searched recursively.\n"
                         "Output: <folder holding floAnalyzer>/analysisOutput (or OUTPUT_DIR)/<image's or selected folder's name>/<same sub-folders>.\n";
            return arg == "-h" || arg == "--help" ? 0 : 2;
        }
    }
    const fs::path workspace = workspaceFolder();  // the file dialog starts here; analysisOutput goes here
    if (inputs.empty() && (inputs = pickFiles(workspace)).empty()) {
        std::cerr << "Nothing selected.\n";
        return 0;  // cancelling the dialog is not an error
    }

    // Each image, its output folder and its folder as listed in the workbook (from the selected folder down).
    // All output goes under analysisOutput in the folder holding floAnalyzer, in a folder named after the selected
    // folder (or a selected image's folder), keeping the sub-folder structure.
    struct Job {
        fs::path image, outDir;
        std::string folder;
    };
    std::vector<Job> jobs;
    if (outDir.empty()) outDir = workspace / "analysisOutput";
    for (fs::path in : inputs) {
        in = fs::absolute(in).lexically_normal();
        if (!in.has_filename()) in = in.parent_path();  // a folder given with a trailing separator
        if (!fs::is_directory(in)) {
            const fs::path dir = outDir / in.parent_path().filename();
            jobs.push_back({in, dir, in.parent_path().filename().generic_string()});
            continue;
        }
        const fs::path base = outDir / in.filename();
        std::vector<fs::path> found;
        for (const auto& e : fs::recursive_directory_iterator(in, fs::directory_options::skip_permission_denied))
            if (e.is_regular_file() && isTiff(e.path())) found.push_back(e.path());
        std::sort(found.begin(), found.end());
        for (const fs::path& p : found)
            jobs.push_back({p, (base / p.parent_path().lexically_relative(in)).lexically_normal(),
                            p.parent_path().lexically_relative(in.parent_path()).generic_string()});
        if (found.empty()) std::cerr << in.string() << ": no TIFF images found\n";
    }

    // A run's output replaces the last one's: earlier SVGs and CSVs in the output folders, and earlier
    // similarity plots and withinInRange.xlsx, are removed first.
    std::set<fs::path> bases;  // each selected folder's (or image folder's) output folder
    for (const Job& job : jobs) bases.insert(outDir / *job.outDir.lexically_relative(outDir).begin());
    for (const fs::path& base : bases) {
        std::error_code ec;
        std::vector<fs::path> stale;
        for (const auto& e : fs::recursive_directory_iterator(base, fs::directory_options::skip_permission_denied, ec)) {
            const std::string name = e.path().filename().string();
            const bool svg = name.size() > 11 && name.compare(name.size() - 11, 11, "_traces.svg") == 0;
            if (e.is_regular_file() && (svg || e.path().extension() == ".csv")) stale.push_back(e.path());
        }
        for (const fs::path& p : stale) fs::remove(p, ec);
    }
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(outDir, ec))
        if ((e.path().filename().string().rfind("similarityResult", 0) == 0 && e.path().extension() == ".svg") ||
            e.path().filename() == "withinInRange.xlsx")
            fs::remove(e.path(), ec);

    // Every image with a neck gets an SVG. The measurements of all images without a neck go to
    // analysisResult.xlsx, a sheet per folder; images with a neck or without both traces to its
    // "Special cases" sheet.
    std::vector<std::vector<std::string>> special = {{"Folder", "Prefix", "Image", "Suffix", "Comment"}};
    std::vector<std::string> specialLast;  // the last special case's name parts
    std::string specialFolder;
    std::map<std::string, std::vector<std::vector<std::string>>> results;  // folder -> its rows
    std::map<std::string, std::vector<std::string>> resultLast;           // folder -> its last row's name parts
    std::vector<Sample> samples;                                          // for similarityResult_<type>.svg
    int failures = 0;
    for (const Job& job : jobs) {
        const fs::path& file = job.image;
        try {
            const Gray img = readTiff(file);
            const Analysis analysis = analyzeImage(img);
            const std::vector<Trace>& traces = analysis.traces;
            const char* style = analysis.style == ImageStyle::Rims ? "rimmed" : "edged";
            const auto listed = [&](const std::string& comment) {  // the folder only where it changes
                const bool newFolder = job.folder != specialFolder;
                if (newFolder) specialLast.clear(), specialFolder = job.folder;
                addImageRow(special, specialLast, file.filename().string(), {comment}, {newFolder ? job.folder : ""});
            };
            if (traces.size() < 2) {
                const std::string comment = !analysis.error.empty() ? analysis.error : traces.empty() ? "no trace"
                                            : traces[0].vertical ? "only vertical trace" : "only horizontal trace";
                listed(comment);
                std::printf("%s: %s\n", file.string().c_str(), comment.c_str());
                continue;
            }
            const Trace& h = traces[0].vertical ? traces[1] : traces[0];
            const Trace& v = traces[0].vertical ? traces[0] : traces[1];
            std::printf("%s: %zu trace(s), %s", file.string().c_str(), traces.size(), style);
            if (h.connected) {
                listed("neck width: " + number(h.neckWidth) + " px");
            } else {
                std::vector<std::vector<std::string>>& rows = results[job.folder];
                if (rows.empty())
                    rows.push_back(kResultHeader);
                const std::string top = *fs::path(job.folder).begin();  // the selected folder
                samples.push_back({top.substr(0, top.find('_')), job.folder, file.filename().string(),
                                   h.width() * h.length() * depth(h), v.width() * v.length() * depth(v),
                                   {number(h.width()), number(h.length()), number(depth(h)), number(v.width()), number(v.length()), number(depth(v))}, file});
                addImageRow(rows, resultLast[job.folder], file.filename().string(), samples.back().measurements);
            }
            if (h.connected) {
                fs::create_directories(job.outDir);
                const fs::path out = job.outDir / (file.stem().string() + "_traces.svg");
                writeSvg(out, file.filename().string() + " (" + style + ")", img, traces);
                std::printf(" -> %s", out.string().c_str());
            }
            std::printf("\n");
            for (const Trace& t : traces) {
                std::printf("  %d (%s): x %.0f-%.0f  y %.0f-%.0f  length %4.0f  width %3.0f px  depth %.1f\n", t.id, t.vertical ? "vertical" : "horizontal", t.x0(), t.x1(), t.y0(), t.y1(), t.length(), t.width(), depth(t));
                if (t.connected) std::printf("  connection: neck into trace 2, narrowest width %.0f px\n", t.neckWidth);
            }
        } catch (const std::exception& e) {
            std::fprintf(stderr, "%s: %s\n", file.string().c_str(), e.what());
            ++failures;
        }
    }
    if (!results.empty() || special.size() > 1) {
        std::vector<Sheet> sheets;
        for (auto& [folder, rows] : results) sheets.push_back({fs::path(folder).filename().string(), std::move(rows), {}});
        sheets.push_back({"Special cases", std::move(special), {}});
        fs::create_directories(outDir);
        writeXlsx(outDir / "analysisResult.xlsx", sheets);
        std::printf("Measurements of %zu folder(s) and the special cases -> %s\n", sheets.size() - 1, (outDir / "analysisResult.xlsx").string().c_str());
    }
    std::map<std::string, std::vector<Sample>> byType;
    for (const Sample& s : samples) byType[s.type].push_back(s);
    for (const auto& [type, list] : byType) {
        const fs::path out = outDir / ("similarityResult_" + type + ".svg");
        fs::create_directories(outDir);
        writeSimilarity(out, type, list);
        std::printf("%s Analysis Result -> %s\n", type.c_str(), out.string().c_str());
    }

    // withinInRange.xlsx: the images within the standard deviation range of their type, i.e. within mean +- SD
    // on both similarity factors (the plot's shaded band), in analysisResult.xlsx's format: a sheet per folder.
    // Each shown prefix links to the image's raw data folder, each image name to the raw image.
    std::map<std::string, Sheet> inRange;
    std::map<std::string, std::vector<std::string>> inRangeLast;
    for (const auto& [type, list] : byType) {
        std::vector<double> hts, vts;
        for (const Sample& s : list) hts.push_back(s.ht), vts.push_back(s.vt);
        const auto [mh, sh] = meanSd(hts);
        const auto [mv, sv] = meanSd(vts);
        for (const Sample& s : list) {
            if (std::abs(s.ht - mh) > sh || std::abs(s.vt - mv) > sv) continue;
            Sheet& sheet = inRange[s.folder];
            if (sheet.rows.empty()) sheet.rows.push_back(kResultHeader);
            addImageRow(sheet.rows, inRangeLast[s.folder], s.image, s.measurements);
            const size_t r = sheet.rows.size() - 1;
            if (!sheet.rows[r][0].empty()) sheet.links[{r, 0}] = fileUrl(s.file.parent_path()) + "/";
            sheet.links[{r, 1}] = fileUrl(s.file);
        }
    }
    if (!inRange.empty()) {
        std::vector<Sheet> sheets;
        size_t count = 0;
        for (auto& [folder, sheet] : inRange) {
            count += sheet.rows.size() - 1;
            sheet.name = fs::path(folder).filename().string();
            sheets.push_back(std::move(sheet));
        }
        writeXlsx(outDir / "withinInRange.xlsx", sheets);
        std::printf("%zu image(s) within the standard deviation range -> %s\n", count, (outDir / "withinInRange.xlsx").string().c_str());
    }
    return failures ? 1 : 0;
}
