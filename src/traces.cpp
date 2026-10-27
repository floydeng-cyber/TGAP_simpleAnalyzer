#include "traces.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <tuple>

namespace {

constexpr double kRoughFraction = 0.5; // rough extent: where the mean-contrast profile is above half its peak
constexpr double kBodyFraction = 0.5;  // the body along the trace: mean contrast above this fraction of its peak
constexpr double kBrightSigmas = 4;    // a bright pixel (trace edge) is this many noise sigmas above the background
constexpr double kEdgeFraction = 0.03; // a side edge row/column has at least this fraction of bright pixels
constexpr int kEdgeSearch = 20;        // px around the rough extent searched for the side edges
constexpr int kMaxGap = 20;            // px of weak body tolerated along a trace (faint traces are patchy)
constexpr int kCapGap = 10;            // px without bright pixels tolerated beyond the body (end caps)
constexpr double kNeckReach = 1.5;     // a neck spans at most this many trace widths ...
constexpr int kNeckSlack = 20;         // ... plus this many px
constexpr double kNeckDip = 1.5;       // an end cap's edges (see halfWidth) are at least this many sigmas of local contrast
constexpr double kWallLow = 1.0;       // a wall: local contrast at least this (sigmas), connected to ...
constexpr double kWallHigh = 2.5;      // ... at least kWallPixels pixels of at least this
constexpr int kWallPixels = 60;
constexpr double kNeckDark = 4;        // inside a trace: 13x13 mean this many sigmas below the background
constexpr int kNeckGap = 7;            // px; the background floods in through openings at least 2 kNeckGap + 1 px wide
constexpr int kNeckJump = 3;           // px a neck's outline moves inward at most from one slice to the next (more: a leak)
constexpr double kNeckStart = 0.9;     // the trace narrows into a neck below this fraction of its width
constexpr double kNeckRows = 3;        // neck box rows hold an edge: local contrast at least this (sigmas) ...
constexpr int kNeckEdgePixels = 4;     // ... in at least this many pixels
constexpr double kNeckPath = -1.35;    // a neck's edges keep the widest path above this local contrast (sigmas)
constexpr int kMargin = 3;             // px added around the bright edges
constexpr double kMinStrength = 2.0;   // rough peak contrast over the background's below which there is no trace
constexpr int kRidgeInset = 2;         // px inside the ridges that count as the trace's dark interior

// Rimmed images (after processing: thick bright rims, the trace's inside about as dark as the background).
constexpr int kRimRun = 12;            // px; bright runs along rows at least this long (median): a rimmed image
constexpr double kRimSigmas = 4;       // a rim pixel: 5x5 mean this many noise sigmas above the background, and
constexpr double kRimLevel = 0.4;      // this fraction of the way from the background to the rims' level (98th
                                       // percentile), above the dimmer film of debris some sites carry
constexpr double kRimMinExtent = 0.2;  // a trace's rim spans at least this fraction of the image; shorter: a particle
constexpr int kRimBlobPixels = 20;     // a particle (or dot) holds at least this many rim pixels; fewer: noise
constexpr int kDotCount = 3;           // no trace but at least this many particles: a dot array (manufacture error)
constexpr double kRimTall = 0.3;       // the vertical trace's rims fill at least this fraction of a column
constexpr int kRimEdge = 8;            // px next to the vertical trace not searched for the waist (where the rims merge)

// Mean of `v` (an image-sized array) per column (columns = true) over rows [from, to), or per row
// over columns [from, to).
std::vector<double> profile(const std::vector<float>& v, int w, int h, bool columns, int from, int to) {
    std::vector<double> p(columns ? w : h, 0);
    for (int i = 0; i < static_cast<int>(p.size()); ++i) {
        for (int k = from; k < to; ++k) p[i] += v[columns ? size_t(k) * w + i : size_t(i) * w + k];
        p[i] /= std::max(1, to - from);
    }
    return p;
}

// 5 px running mean.
std::vector<double> smooth(const std::vector<double>& p) {
    std::vector<double> s(p.size());
    for (int i = 0; i < static_cast<int>(p.size()); ++i) {
        const int lo = std::max(0, i - 2), hi = std::min(static_cast<int>(p.size()) - 1, i + 2);
        for (int k = lo; k <= hi; ++k) s[i] += p[k] / (hi - lo + 1);
    }
    return s;
}

// (2r+1)x(2r+1) running mean (separable, shrinking at the border).
std::vector<float> blur(const Gray& img, int r) {
    std::vector<float> out = img.px;
    for (int pass = 0; pass < 2; ++pass) {
        const std::vector<float> in = out;
        const int n = pass ? img.h : img.w, lines = pass ? img.w : img.h;
        for (int l = 0; l < lines; ++l)
            for (int i = 0; i < n; ++i) {
                double sum = 0;
                const int lo = std::max(0, i - r), hi = std::min(n - 1, i + r);
                for (int k = lo; k <= hi; ++k) sum += in[pass ? size_t(k) * img.w + l : size_t(l) * img.w + k];
                out[pass ? size_t(i) * img.w + l : size_t(l) * img.w + i] = static_cast<float>(sum / (hi - lo + 1));
            }
    }
    return out;
}

// Number of set pixels in the (2r+1)x(2r+1) square around each pixel (clipped at the border).
std::vector<int> boxCount(const std::vector<uint8_t>& m, int w, int h, int r) {
    std::vector<int> sum(size_t(w + 1) * (h + 1), 0), out(m.size());
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
            sum[size_t(y + 1) * (w + 1) + x + 1] = m[size_t(y) * w + x] + sum[size_t(y) * (w + 1) + x + 1] + sum[size_t(y + 1) * (w + 1) + x] - sum[size_t(y) * (w + 1) + x];
    for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
            const int x0 = std::max(0, x - r), x1 = std::min(w, x + r + 1), y0 = std::max(0, y - r), y1 = std::min(h, y + r + 1);
            out[size_t(y) * w + x] = sum[size_t(y1) * (w + 1) + x1] - sum[size_t(y0) * (w + 1) + x1] - sum[size_t(y1) * (w + 1) + x0] + sum[size_t(y0) * (w + 1) + x0];
        }
    return out;
}

double median(std::vector<double> v) {
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
}

// The two ridges of a trace in intensity profile p within its box [lo, hi]: the brightest point on each
// side of the darkest point (the trace's inside), at most `reach` px from the box's end. Found on a 5 px
// mean (so noise cannot win), then moved to the profile's own peak within 2 px.
std::pair<int, int> ridges(const std::vector<double>& p, int lo, int hi, int reach) {
    const std::vector<double> s = smooth(p);
    const auto first = s.begin() + lo, last = s.begin() + std::max(lo, hi) + 1;
    const auto trough = std::min_element(first, last);
    auto peak = [&](int i) {
        int best = i;
        for (int k = std::max(lo, i - 2); k <= std::min(hi, i + 2); ++k)
            if (p[k] > p[best]) best = k;
        return best;
    };
    return {peak(static_cast<int>(std::max_element(first, std::min(trough, first + reach) + 1) - s.begin())),
            peak(static_cast<int>(std::max_element(std::max(trough, last - 1 - reach), last) - s.begin()))};
}

// Background (median) and noise (median deviation, returned in sigma) of a blurred image, re-estimated
// from pixels within 3 sigmas until stable, so traces do not inflate them.
double background(const std::vector<float>& s, double& sigma) {
    double bg = 0, previous = 0;
    sigma = HUGE_VAL;
    for (int i = 0; i < 100 && (i == 0 || std::abs(sigma - previous) > 1e-3 * sigma); ++i) {
        std::vector<double> kept, dev;
        for (float v : s)
            if (std::abs(v - bg) < 3 * sigma) kept.push_back(v);
        bg = median(kept);
        for (double v : kept) dev.push_back(std::abs(v - bg));
        previous = sigma, sigma = std::max(1e-3, 1.4826 * median(dev));
    }
    return bg;
}

class Finder {
public:
    explicit Finder(const Gray& img) : img_(img), w_(img.w), h_(img.h) {
        // Contrast (5x5 mean, for finding traces) and bright pixels (3x3 mean, which keeps thin edges
        // and end caps), each in noise sigmas of its own blur.
        const std::vector<float> s5 = blur(img, 2), s3 = blur(img, 1), s13 = blur(img, 6);
        double sigma5, sigma3;
        const double bg5 = background(s5, sigma5), bg3 = background(s3, sigma3);
        bg_ = bg5, sigma_ = sigma5;
        contrast_.resize(s5.size());
        bright_.resize(s3.size());
        wide_.resize(s13.size());
        for (size_t i = 0; i < s5.size(); ++i) {
            contrast_[i] = static_cast<float>(std::abs(s5[i] - bg5) / sigma5);
            bright_[i] = s3[i] - bg3 > kBrightSigmas * sigma3 ? 1.f : 0.f;
            wide_[i] = static_cast<float>((s13[i] - bg5) / sigma5);
        }
        // Local contrast: the 5x5 mean minus the 21x21 mean, in its own noise sigmas. Slow shading (the
        // dark haze around thick traces and necks) cancels out; edges stand out.
        const std::vector<float> s21 = blur(img, 10);
        std::vector<float> local(s5.size());
        for (size_t i = 0; i < s5.size(); ++i) local[i] = s5[i] - s21[i];
        double sigmaL;
        const double bgL = background(local, sigmaL);
        local_.resize(local.size());
        for (size_t i = 0; i < local.size(); ++i) local_[i] = static_cast<float>((local[i] - bgL) / sigmaL);
    }

    std::vector<Trace> run() {
        std::vector<Trace> found;
        // 1-2. The vertical trace: column profiles of the whole image, then row profiles over its columns.
        Trace v;
        v.id = 2, v.vertical = true;
        const bool hasV = find(v, 0, w_, 0, h_);
        if (hasV) capLength_ = capLength(v), vert_ = v, hasV_ = true;
        // 3-4. The horizontal trace, left of the vertical trace's box: row profiles, then column profiles
        //      over its rows. It starts at the image's left edge.
        Trace t;
        t.id = 1;
        const int split = hasV ? static_cast<int>(v.acrossOut0) : w_;
        if (split > 16 && find(t, 0, split, 0, h_)) found.push_back(t);
        if (hasV) found.push_back(v);
        return found;
    }

private:
    // Half the distance between a trace's two edges at slice s: the brightest local-contrast points on
    // each side of position `mid` within `reach` px (averaged over 5 slices), or -1 if either is under
    // kNeckDip sigmas.
    double halfWidth(bool vert, int s, int mid, int reach) const {
        const int along = vert ? h_ : w_, across = vert ? w_ : h_;
        auto value = [&](int k) {
            double sum = 0;
            const int c = std::clamp(mid + k, 0, across - 1);
            for (int ds = -2; ds <= 2; ++ds) {
                const int ss = std::clamp(s + ds, 0, along - 1);
                sum += local_[vert ? size_t(ss) * w_ + c : size_t(c) * w_ + ss];
            }
            return sum / 5;
        };
        int up = -reach, down = reach;
        for (int k = -reach; k <= 0; ++k)
            if (value(k) > value(up)) up = k;
        for (int k = 0; k <= reach; ++k)
            if (value(k) > value(down)) down = k;
        return std::min(value(up), value(down)) < kNeckDip ? -1 : (down - up) / 2.0;
    }

    // A free end cap is a half-ellipse: half-width W/2 at its start, 0 at its tip, so at distance d back
    // from the tip, g = sqrt(1 - (2h/W)^2) = 1 - d / a, a being the cap's length. Measured on the
    // vertical trace's two ends (same structure, no neck there), median over the slices where the edges
    // curve in; half the trace's width if that fails.
    double capLength(const Trace& v) const {
        const double W = v.width();
        const int mid = static_cast<int>(std::lround((v.across0 + v.across1) / 2)), reach = static_cast<int>(W / 2) + kEdgeSearch;
        std::vector<double> lengths;
        for (int end = 0; end < 2; ++end) {
            const double tip = end ? v.along1 : v.along0;
            const int dir = end ? -1 : 1;  // from the tip into the trace
            for (int d = 2; d <= static_cast<int>(2 * W); ++d) {
                const double h = halfWidth(true, static_cast<int>(tip) + dir * d, mid, reach);
                if (h < 0) continue;
                if (h >= 0.95 * W / 2) break;  // the cap has ended: full width
                const double g = std::sqrt(std::max(0.0, 1 - (2 * h / W) * (2 * h / W)));
                if (g > 0.15 && g < 0.95) lengths.push_back(d / (1 - g));
            }
        }
        return lengths.size() >= 3 ? median(lengths) : W / 2;
    }

    // The widest path in the local contrast over columns [x0, x1] and rows [y0, y1], from the first column
    // to the last (8-connected): of all such paths, the one whose dimmest pixel is brightest; returns that
    // pixel's value. A continuous bright structure across the region scores high; a stretch of background
    // anywhere across it caps the score near the noise level.
    double widestPath(int x0, int x1, int y0, int y1) const {
        const int W = x1 - x0 + 1, H = y1 - y0 + 1;
        if (W < 2 || H < 1) return HUGE_VAL;
        auto value = [&](int i) { return static_cast<double>(local_[size_t(y0 + i / W) * w_ + x0 + i % W]); };
        std::vector<double> best(size_t(W) * H, -HUGE_VAL);
        std::priority_queue<std::pair<double, int>> queue;
        for (int y = 0; y < H; ++y) best[size_t(y) * W] = value(y * W), queue.push({best[size_t(y) * W], y * W});
        while (!queue.empty()) {
            const auto [v, i] = queue.top();
            queue.pop();
            if (v < best[i]) continue;
            if (i % W == W - 1) return v;
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    const int x = i % W + dx, y = i / W + dy;
                    if (x < 0 || y < 0 || x >= W || y >= H) continue;
                    const int j = y * W + x;
                    const double nv = std::min(v, value(j));
                    if (nv > best[j]) best[j] = nv, queue.push({nv, j});
                }
        }
        return -HUGE_VAL;
    }

    // The background around the horizontal trace t: pixels that are neither wall nor inside a trace (13x13
    // mean more than kNeckDark sigmas below the background), flooded in from the image's border through
    // openings at least 2 kNeckGap + 1 px wide. A wall is a connected stretch of local contrast of at least
    // kWallLow sigmas holding at least kWallPixels pixels of kWallHigh: so a faint stretch of a long wall
    // still counts, while a noise speck does not. The background thus fills every bay right up to the
    // traces' outer walls but cannot squeeze through a faint wall into a neck. Each trace is closed where
    // it leaves the image. Returns 1 for the background, 2 for a wall, 0 for anything else.
    std::vector<uint8_t> backgroundMask(const Trace& t) const {
        std::vector<uint8_t> open(local_.size()), core(local_.size(), 0), out(local_.size(), 0);
        for (size_t i = 0; i < open.size(); ++i) open[i] = wide_[i] > -kNeckDark;
        std::vector<uint8_t> seen(local_.size(), 0), wall(local_.size(), 0);
        for (size_t i0 = 0; i0 < seen.size(); ++i0) {
            if (seen[i0] || local_[i0] < kWallLow) continue;
            std::vector<size_t> part{i0};
            seen[i0] = 1;
            int strong = 0;
            for (size_t q = 0; q < part.size(); ++q) {
                strong += local_[part[q]] >= kWallHigh;
                const int x = static_cast<int>(part[q] % w_), y = static_cast<int>(part[q] / w_);
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        const int xx = x + dx, yy = y + dy;
                        const size_t j = size_t(yy) * w_ + xx;
                        if (xx >= 0 && yy >= 0 && xx < w_ && yy < h_ && !seen[j] && local_[j] >= kWallLow) seen[j] = 1, part.push_back(j);
                    }
            }
            if (strong >= kWallPixels)
                for (size_t j : part) open[j] = 0, wall[j] = 1;
        }
        auto close = [&](double x0, double x1, double y0, double y1) {
            for (int y = std::max(0, int(y0)); y <= std::min(h_ - 1, int(y1)); ++y)
                for (int x = std::max(0, int(x0)); x <= std::min(w_ - 1, int(x1)); ++x) open[size_t(y) * w_ + x] = 0;
        };
        close(0, kNeckGap, t.acrossOut0, t.acrossOut1);
        if (hasV_) close(vert_.acrossOut0, vert_.acrossOut1, 0, kNeckGap), close(vert_.acrossOut0, vert_.acrossOut1, h_ - 1 - kNeckGap, h_ - 1);
        // Pixels whose whole (2 kNeckGap + 1)-square is open, flooded from the border; then that square
        // around each reached pixel, where open.
        std::vector<uint8_t> closed(open.size());
        for (size_t i = 0; i < open.size(); ++i) closed[i] = !open[i];
        const std::vector<int> blocked = boxCount(closed, w_, h_, kNeckGap);
        std::vector<size_t> stack;
        for (int x = 0; x < w_; ++x) stack.push_back(size_t(x)), stack.push_back(size_t(h_ - 1) * w_ + x);
        for (int y = 0; y < h_; ++y) stack.push_back(size_t(y) * w_), stack.push_back(size_t(y) * w_ + w_ - 1);
        while (!stack.empty()) {
            const size_t i = stack.back();
            stack.pop_back();
            if (core[i] || blocked[i]) continue;
            core[i] = 1;
            const int x = static_cast<int>(i % w_), y = static_cast<int>(i / w_);
            if (x > 0) stack.push_back(i - 1);
            if (x < w_ - 1) stack.push_back(i + 1);
            if (y > 0) stack.push_back(i - w_);
            if (y < h_ - 1) stack.push_back(i + w_);
        }
        const std::vector<int> nearCore = boxCount(core, w_, h_, kNeckGap);
        for (size_t i = 0; i < out.size(); ++i) out[i] = open[i] && nearCore[i] > 0 ? 1 : wall[i] ? 2 : 0;
        return out;
    }

    // Profile along the trace's direction or across it, over rows/columns [from, to) of the other axis.
    std::vector<double> profileOf(const std::vector<float>& v, const Trace& t, bool acrossIt, int from, int to) const {
        return profile(v, w_, h_, t.vertical == acrossIt, from, to);
    }

    // Finds trace t (its orientation and id set) within columns [x0, x1) and rows [y0, y1).
    bool find(Trace& t, int x0, int x1, int y0, int y1) const {
        const bool vert = t.vertical;
        const int a0 = vert ? y0 : x0, a1 = vert ? y1 : x1;  // along the trace
        const int c0 = vert ? x0 : y0, c1 = vert ? x1 : y1;  // across it
        // Rough extent across: where the mean contrast (over the whole region along the trace) is above
        // half its peak, around the peak.
        const std::vector<double> rough = smooth(profileOf(contrast_, t, true, a0, a1));
        const int peak = static_cast<int>(std::max_element(rough.begin() + c0, rough.begin() + c1) - rough.begin());
        const double noise = median(std::vector<double>(rough.begin() + c0, rough.begin() + c1));
        t.contrast = rough[peak] / std::max(1e-3, noise);
        if (t.contrast < kMinStrength) return false;
        int lo = peak, hi = peak;
        while (lo > c0 && rough[lo - 1] >= noise + kRoughFraction * (rough[peak] - noise)) --lo;
        while (hi < c1 - 1 && rough[hi + 1] >= noise + kRoughFraction * (rough[peak] - noise)) ++hi;

        const int half = (hi - lo) / 2, mid = (lo + hi) / 2;
        const int wideLo = std::max(c0, lo - kEdgeSearch), wideHi = std::min(c1 - 1, hi + kEdgeSearch);
        auto at = [&](const std::vector<float>& m, int s, int c) { return m[vert ? size_t(s) * w_ + c : size_t(c) * w_ + s]; };

        // Along: the body first, where the mean contrast over the band stays above kBodyFraction of its
        // typical level (gaps of kMaxGap px tolerated, so it carries on where an edge breaks,
        // e.g. at a neck); then on while slices have any bright pixel (within the widened band, gaps of
        // kCapGap px), which takes in the end caps. The horizontal trace always starts at the image's
        // left edge.
        std::vector<double> body(std::max(a1, vert ? h_ : w_), 0);
        for (int s = a0; s < a1; ++s) {
            for (int c = mid - half; c <= mid + half; ++c) body[s] += at(contrast_, s, c);
            body[s] /= 2 * half + 1;
        }
        body = smooth(body);
        const int bodyPeak = static_cast<int>(std::max_element(body.begin() + a0, body.begin() + a1) - body.begin());
        std::vector<double> strongSlices;
        for (int s = a0; s < a1; ++s)
            if (body[s] >= noise + 0.5 * (body[bodyPeak] - noise)) strongSlices.push_back(body[s]);
        const double bodyLevel = noise + kBodyFraction * (median(strongSlices) - noise);
        const std::vector<double> alongCount = profileOf(bright_, t, false, wideLo, wideHi + 1);
        // An end cap: the mean over the middle tenth of the band (and 3 slices) is a
        // noise sigma brighter than the background. The cap's bright arc crosses there; the trace's inside,
        // a neck's gray middle and the dark shading past a thick trace's tip are all darker. (A bright speck
        // inside a patchy trace can pass too.)
        const int centre = std::max(1, (hi - lo) / 20);
        auto cap = [&](int s) {
            double sum = 0;
            int count = 0;
            for (int k = std::max(a0, s - 1); k <= std::min(a1 - 1, s + 1); ++k)
                for (int c = mid - centre; c <= mid + centre; ++c) sum += at(img_.px, k, c), ++count;
            return sum / count > bg_ + sigma_;
        };
        bool capped = false;
        auto walk = [&](const std::vector<double>& p, double level, int from, int dir, int maxGap, bool stopAtCap) {
            int end = from;
            for (int s = from, gap = 0; s >= a0 && s < a1 && gap <= maxGap; s += dir) {
                if (stopAtCap && s != from && cap(s)) return capped = true, s;
                p[s] >= level ? (end = s, gap = 0) : ++gap;
            }
            return end;
        };
        // Only the horizontal trace's right end stops at a cap: that is where a cap must be told apart from
        // a neck or the dark shading in the gap. Elsewhere a bright speck inside a patchy trace could pass
        // for a cap and cut the trace short.
        int s0 = walk(alongCount, 1e-9, walk(body, bodyLevel, bodyPeak, -1, kMaxGap, false), -1, kCapGap, false);
        int s1 = walk(alongCount, 1e-9, walk(body, bodyLevel, bodyPeak, 1, kMaxGap, !vert), 1, kCapGap, false);
        if (!vert) {
            s0 = a0;
            // A neck joins the two traces if the trace ends close to the vertical trace's box (within
            // kNeckReach trace widths + kNeckSlack px) and a continuous bright structure (the neck's edges)
            // runs within the trace's rows from inside its body to the vertical trace's edge: the widest path
            // in the local contrast stays above kNeckPath. Across a gap, the background between the tip and
            // the vertical trace breaks every such path.
            const int band = hi - lo + 1;
            t.connected = a1 < w_ && a1 - s1 <= kNeckReach * band + kNeckSlack &&
                          widestPath(std::max(a0, s1 - band - 10), std::min(w_ - 1, a1 + 3), std::max(0, lo - 4),
                                     std::min(h_ - 1, hi + 4)) >= kNeckPath;
            if (t.connected) s1 = a1 - 1;
        }

        // Across: the outermost rows/columns near the rough edges that are part of a bright side edge:
        // at least kEdgeFraction of the trace's length (and 2 px) bright. A tilted edge spreads over a few
        // rows/columns, each holding only part of it.
        const std::vector<double> acrossCount = profileOf(bright_, t, true, s0, s1 + 1);
        const double edgeNeed = std::max(kEdgeFraction, 2.0 / (s1 - s0 + 1));
        int e0 = lo, e1 = hi;
        for (int c = std::max(c0, lo - kEdgeSearch); c <= std::min(c1 - 1, lo + kEdgeSearch); ++c)
            if (acrossCount[c] >= edgeNeed) { e0 = std::min(e0, c); break; }
        for (int c = std::min(c1 - 1, hi + kEdgeSearch); c >= std::max(c0, hi - kEdgeSearch); --c)
            if (acrossCount[c] >= edgeNeed) { e1 = std::max(e1, c); break; }

        const int limitA = vert ? h_ : x1, limitC = vert ? x1 : h_;
        t.acrossOut0 = std::max(c0, e0 - kMargin), t.acrossOut1 = std::min(limitC - 1, e1 + kMargin);
        t.alongOut0 = std::max(a0, s0 - kMargin), t.alongOut1 = std::min(limitA - 1, s1 + kMargin);

        // A slightly tilted trace drifts sideways along its length, which would smear its two side edges
        // together when averaged along it. So the profiles below follow the trace's centre line: a
        // straight line fitted to each slice's contrast centroid.
        double n = 0, sa = 0, sc = 0, saa = 0, sac = 0;
        for (int s = s0; s <= s1; ++s) {
            double weight = 0, moment = 0;
            for (int c = e0; c <= e1; ++c) {
                const double v = contrast_[vert ? size_t(s) * w_ + c : size_t(c) * w_ + s];
                weight += v, moment += v * c;
            }
            if (weight <= 0) continue;
            const double centre = moment / weight;
            n += 1, sa += s, sc += centre, saa += double(s) * s, sac += s * centre;
        }
        const double slope = n > 1 && n * saa - sa * sa > 0 ? (n * sac - sa * sc) / (n * saa - sa * sa) : 0;
        const int midA = (s0 + s1) / 2;
        auto pixel = [&](int s, int c) {  // image value at slice s, position c of the middle slice's frame
            const int cc = std::clamp(c + static_cast<int>(std::lround(slope * (s - midA))), 0, (vert ? w_ : h_) - 1);
            return static_cast<double>(img_.px[vert ? size_t(s) * w_ + cc : size_t(cc) * w_ + s]);
        };
        auto straightProfile = [&](bool acrossIt, int from, int to) {  // like profileOf, along the centre line
            std::vector<double> p(acrossIt ? (vert ? w_ : h_) : (vert ? h_ : w_), 0);
            for (int i = 0; i < static_cast<int>(p.size()); ++i) {
                for (int k = from; k < to; ++k) p[i] += acrossIt ? pixel(k, i) : pixel(i, k);
                p[i] /= std::max(1, to - from);
            }
            return p;
        };

        // Ridges (width and length) in the intensity profiles across the trace (averaged along it) and
        // along it (averaged over its dark interior).
        std::vector<double> across = straightProfile(true, s0, s1 + 1);
        std::tie(t.across0, t.across1) = ridges(across, e0, e1, e1 - e0);
        int from = static_cast<int>(t.across0) + kRidgeInset, to = static_cast<int>(t.across1) + 1 - kRidgeInset;
        if (to <= from) from = static_cast<int>(t.across0), to = static_cast<int>(t.across1) + 1;
        const std::vector<double> along = straightProfile(false, from, to);
        // An end cap lies within about a trace width of the box's end; further in, a bright spot of a
        // faint trace could outshine it.
        std::tie(t.along0, t.along1) = ridges(along, s0, s1, std::max(10, e1 - e0));
        if (!vert) t.along0 = 0;  // the horizontal trace runs on past the image's left edge

        if (t.connected) {
            // The neck's outline: what is left once the background is flooded in (see backgroundMask). In
            // each slice it runs outward from the trace's centre line to where the background starts. Its
            // ridges lie as far inside it as the trace's own ridges lie inside its outline where it has its full
            // width (the median over the first 20 slices): a neck's faint walls hold no ridge of their own to
            // find, and bright texture inside it and the trace's own cap, curving in to the centre, are no
            // outline. The width is ridge to ridge, like the traces'; the waist is the narrowest slice near
            // the end.
            const int mid = static_cast<int>(std::lround((t.across0 + t.across1) / 2));
            const int reach = static_cast<int>(t.width() / 2) + kEdgeSearch;
            const double normal = t.width();
            const int first = std::max(a0, a1 - 1 - static_cast<int>(2 * normal) - 20), last = a1 - 7;  // clear of the vertical trace's edge
            const std::vector<uint8_t> background = backgroundMask(t);
            const int span = static_cast<int>(1.5 * normal);
            std::vector<int> top(last - first + 1, -1), bottom(last - first + 1, -1);  // -1: no outline (open to the side)
            std::vector<double> inset0, inset1;
            for (int i = 0; i <= last - first; ++i) {
                const int s = first + i, c = static_cast<int>(std::lround(mid + slope * (s - midA)));
                auto inside = [&](int y) { return y >= 0 && y < h_ && std::abs(y - c) <= span && background[size_t(y) * w_ + s] != 1; };
                if (!inside(c)) continue;
                int y0 = c, y1 = c;
                while (inside(y0 - 1)) --y0;
                while (inside(y1 + 1)) ++y1;
                if (std::abs(y0 - 1 - c) > span || std::abs(y1 + 1 - c) > span) continue;
                top[i] = y0, bottom[i] = y1;
                if (i < 20) inset0.push_back(t.across0 + slope * (s - midA) - y0), inset1.push_back(y1 - t.across1 - slope * (s - midA));
            }
            // A wall curves smoothly; where the background has squeezed in through a faint stretch of it, the
            // outline jumps inward. Such slices (more than kNeckJump px further in than the last good slice,
            // plus 1 px per slice between them) are bridged by a straight line between the good slices on
            // either side; at the end, where there is no good slice to bridge to, they are left out.
            for (int side = 0; side < 2; ++side) {
                std::vector<int>& edge = side ? bottom : top;
                const int in = side ? -1 : 1;  // the inward direction
                std::vector<int> good;
                for (int i = 0; i <= last - first; ++i)
                    if (edge[i] >= 0 && (good.empty() || (edge[i] - edge[good.back()]) * in <= kNeckJump + (i - good.back() - 1))) good.push_back(i);
                std::vector<int> bridged(edge.size(), -1);
                for (size_t g = 0; g < good.size(); ++g) {
                    bridged[good[g]] = edge[good[g]];
                    if (g + 1 < good.size())
                        for (int i = good[g] + 1; i < good[g + 1]; ++i)
                            bridged[i] = static_cast<int>(std::lround(edge[good[g]] + double(edge[good[g + 1]] - edge[good[g]]) * (i - good[g]) / (good[g + 1] - good[g])));
                }
                edge = bridged;
            }
            for (int i = 0; i <= last - first; ++i)
                if (top[i] < 0 || bottom[i] < 0) top[i] = bottom[i] = -1;
            const int d0 = inset0.empty() ? 0 : static_cast<int>(std::lround(median(inset0))), d1 = inset1.empty() ? 0 : static_cast<int>(std::lround(median(inset1)));
            for (int i = 0; i <= last - first; ++i)
                if (top[i] >= 0) top[i] += d0, bottom[i] = std::max(top[i], bottom[i] - d1);
            // Where the background passes right through the neck (on both sides of the centre line), the neck
            // is a single thread too thin for an outline: its width is where its local contrast (averaged over
            // 7 slices) is at least halfway from its surroundings (the higher of the lowest points within
            // kEdgeSearch px each side) to its peak, the brightest point within kEdgeSearch px of the centre.
            auto value = [&](int s, int c) {
                double sum = 0;
                for (int ds = -3; ds <= 3; ++ds) sum += local_[size_t(std::clamp(c, 0, h_ - 1)) * w_ + std::clamp(s + ds, a0, a1 - 1)];
                return sum / 7;
            };
            std::vector<int> los(top), his(bottom);
            for (int i = 0; i <= last - first; ++i) {
                const int s = first + i, c = static_cast<int>(std::lround(mid + slope * (s - midA)));
                if (background[size_t(c) * w_ + s] != 1) continue;
                int p = c;
                for (int y = c - kEdgeSearch; y <= c + kEdgeSearch; ++y)
                    if (value(s, y) > value(s, p)) p = y;
                double above = value(s, p), below = above;
                for (int k = 1; k <= kEdgeSearch; ++k) above = std::min(above, value(s, p - k)), below = std::min(below, value(s, p + k));
                const double level = (value(s, p) + std::max(above, below)) / 2;
                int lo = p, hi = p;
                while (lo > p - kEdgeSearch && value(s, lo - 1) >= level) --lo;
                while (hi < p + kEdgeSearch && value(s, hi + 1) >= level) ++hi;
                los[i] = lo, his[i] = hi;
            }
            // Each slice's width is the median over 7 neighbouring slices, so a single slice cannot be the
            // narrowest by itself.
            int narrowest = -1, lo = 0, hi = 0;
            for (int i = 0; i <= last - first; ++i) {
                std::vector<std::pair<int, int>> near;  // (width, index) of the slices within 3
                for (int k = std::max(0, i - 3); k <= std::min(last - first, i + 3); ++k)
                    near.push_back({los[k] < 0 ? h_ : his[k] - los[k], k});
                std::nth_element(near.begin(), near.begin() + near.size() / 2, near.end());
                const auto [median, k] = near[near.size() / 2];
                if (narrowest < 0 || median < t.neckWidth) narrowest = first + k, t.neckWidth = median, lo = los[k], hi = his[k];
            }
            if (narrowest >= 0) {
                t.neckAt = narrowest, t.neckFrom = lo, t.neckTo = hi;
                // Where the trace would end without the neck. Walking back from the waist, the outline's ridges
                // curve in from where the trace is last at kNeckStart of its width (`start`) to the inflection,
                // where they turn outward into the neck (the steepest narrowing). Each slice there gives the
                // tip of a free cap through it (see capLength): x + a (1 - g). Box 1 ends at their median.
                const double W = normal, a = capLength_ > 0 ? capLength_ : W / 2;
                std::vector<double> half(a1, -1);
                for (int x = first; x <= last; ++x) half[x] = top[x - first] < 0 ? -1 : (bottom[x - first] - top[x - first]) / 2.0;  // the outline's ridges
                auto smoothHalf = [&](int x) {  // median of the measured half-widths within 2 slices
                    std::vector<double> near;
                    for (int k = std::max(first, x - 2); k <= std::min(a1 - 7, x + 2); ++k)
                        if (half[k] >= 0) near.push_back(half[k]);
                    return near.empty() ? -1.0 : median(near);
                };
                int start = narrowest;
                while (start > first && (smoothHalf(start) < 0 || smoothHalf(start) < kNeckStart * W / 2)) --start;
                int inflection = start;
                double steepest = 0;
                for (int x = start + 2; x <= narrowest - 2; ++x) {
                    const double h0 = smoothHalf(x - 2), h1 = smoothHalf(x + 2);
                    if (h0 >= 0 && h1 >= 0 && (h1 - h0) / 4 < steepest) steepest = (h1 - h0) / 4, inflection = x;
                }
                std::vector<double> tips;
                for (int x = start; x <= inflection; ++x) {
                    const double h = smoothHalf(x);
                    if (h < 0) continue;
                    const double g = std::sqrt(std::max(0.0, 1 - (2 * h / W) * (2 * h / W)));
                    tips.push_back(x + a * (1 - g));
                }
                const int tip = tips.empty() ? narrowest : std::clamp(static_cast<int>(std::lround(median(tips))), start, a1 - 1);
                s1 = tip;
                t.along1 = tip, t.alongOut1 = std::min(a1 - 1, tip + kMargin);
                across = straightProfile(true, s0, start + 1);  // width and cross-section from the full-width part
                std::tie(t.across0, t.across1) = ridges(across, e0, e1, e1 - e0);

                // The neck box: from box 1 to the vertical trace's box, over the rows holding the neck's
                // structure: from the topmost to the bottommost row, within `reach` of the trace's centre, with
                // an edge in those columns (including the arms flaring into the vertical trace): at least
                // kNeckEdgePixels pixels of local contrast kNeckRows sigmas or more. (The neck's dark middle
                // holds none; a single noise blob holds fewer.) It also holds the waist's two edges.
                t.neckX0 = t.alongOut1 + 1, t.neckX1 = a1 - 1;
                if (t.neckX1 > t.neckX0) {
                    auto edgeRow = [&](int y) {
                        int count = 0;
                        for (int x = static_cast<int>(t.neckX0); x <= t.neckX1; ++x) count += local_[size_t(y) * w_ + x] >= kNeckRows;
                        return count >= kNeckEdgePixels;
                    };
                    int y0 = mid, y1 = mid;
                    for (int y = std::max(0, mid - reach); y < mid; ++y)
                        if (edgeRow(y)) { y0 = y; break; }
                    for (int y = std::min(h_ - 1, mid + reach); y > mid; --y)
                        if (edgeRow(y)) { y1 = y; break; }
                    t.neckY0 = std::max(0, std::min(y0, lo) - kMargin), t.neckY1 = std::min(h_ - 1, std::max(y1, hi) + kMargin);  // and the waist
                    // Its cross-section at the waist: intensity at each of its rows, averaged over 7 slices; its
                    // edges are the waist's ridges.
                    const int x0 = std::max(a0, narrowest - 3), x1 = std::min(a1 - 1, narrowest + 3);
                    t.neckFrom0 = t.neckY0, t.neckEdge0 = lo, t.neckEdge1 = hi;
                    for (int y = static_cast<int>(t.neckY0); y <= t.neckY1; ++y) {
                        double sum = 0;
                        for (int x = x0; x <= x1; ++x) sum += img_.px[size_t(y) * w_ + x];
                        t.neckProfile.push_back(sum / (x1 - x0 + 1));
                    }
                }
            }
        }

        // Cross-section between the box edges (in the middle slice's frame).
        t.crossFrom = t.acrossOut0;
        t.crossSection.assign(across.begin() + static_cast<int>(t.acrossOut0), across.begin() + static_cast<int>(t.acrossOut1) + 1);
        t.plotFirst = 0;
        t.plotLast = static_cast<int>(t.crossSection.size()) - 1;
        return true;
    }

    const Gray& img_;
    int w_, h_;
    std::vector<float> contrast_;  // |5x5 mean - background| in noise sigmas
    std::vector<float> bright_;    // 1 where the 3x3 mean is kBrightSigmas above the background
    std::vector<float> local_;     // local contrast (5x5 mean - 21x21 mean) in noise sigmas
    std::vector<float> wide_;      // 13x13 mean - background, in noise sigmas of the 5x5 mean (signed)
    double bg_ = 0, sigma_ = 1;    // background level and noise of the 5x5 mean
    Trace vert_;                   // the vertical trace, once found
    bool hasV_ = false;
    double capLength_ = 0;         // length of a free end cap, from the vertical trace (0: unknown)
};

// The style of an image: rimmed if the bright pixels (5x5 mean kRimSigmas above the background) form runs along
// the rows of at least kRimRun px (median); thin bright edges give runs of a few px.
ImageStyle styleOf(const Gray& img) {
    const std::vector<float> s5 = blur(img, 2);
    double sigma;
    const double bg = background(s5, sigma);
    std::vector<int> runs;
    for (int y = 0; y < img.h; ++y)
        for (int x = 0, run = 0; x <= img.w; ++x) {
            if (x < img.w && s5[size_t(y) * img.w + x] > bg + kRimSigmas * sigma) ++run;
            else if (run) runs.push_back(run), run = 0;
        }
    if (runs.empty()) return ImageStyle::Edges;
    std::nth_element(runs.begin(), runs.begin() + runs.size() / 2, runs.end());
    return runs[runs.size() / 2] >= kRimRun ? ImageStyle::Rims : ImageStyle::Edges;
}

// Finds the two traces in a rimmed image, from each trace's inside outward, so that debris around a trace
// (with rims of its own) is never taken for it:
// 1. Rim pixels (5x5 mean kRimSigmas above the background and kRimLevel of the way to the rims' level), in
//    connected pieces; pieces shorter than kRimMinExtent of the image are particles and are dropped (no trace
//    at all but kDotCount or more particles: a dot array, a manufacture error).
// 2. The vertical trace: its two rims are the two fullest columns of rim pixels (at least kRimTall of the image);
//    its centre line runs between them. From the centre line, each row is walked outward to the first rim and
//    through it: the rim's outer edge is the trace's outline there. Its ends, likewise up and down the centre.
// 3. The horizontal trace: left of the vertical one, its two rims are the two fullest rows of rim pixels; from
//    its centre line each column is walked up and down the same way, from the image's left edge on, until no rim
//    is found within two trace widths (past its end cap). It is connected to the vertical trace if its outline
//    runs on into it; then its box ends at the narrowest column (the waist), and the rest is the neck.
// 4. Width, length and depth from the intensity profiles across and along each trace, between its ridges (the
//    rims' peaks), as for edged images.
class RimFinder {
public:
    explicit RimFinder(const Gray& img) : img_(img), w_(img.w), h_(img.h) {}

    Analysis run() {
        Analysis result;
        result.style = ImageStyle::Rims;
        const std::vector<float> s5 = blur(img_, 2);
        double sigma;
        const double bg = background(s5, sigma);
        std::vector<float> sorted(s5);
        std::nth_element(sorted.begin(), sorted.begin() + static_cast<long>(0.98 * (sorted.size() - 1)), sorted.end());
        const double rimLevel = sorted[static_cast<size_t>(0.98 * (sorted.size() - 1))];  // the rims' bright level
        const double threshold = std::max(bg + kRimSigmas * sigma, bg + kRimLevel * (rimLevel - bg));
        std::vector<uint8_t> raw(s5.size());
        for (size_t i = 0; i < s5.size(); ++i) raw[i] = s5[i] > threshold;

        // 1. Rim pieces: traces' rims are kept, particles dropped and counted.
        rim_.assign(raw.size(), 0);
        std::vector<uint8_t> seen(raw.size(), 0);
        int particles = 0;
        const int minExtent = static_cast<int>(kRimMinExtent * std::min(w_, h_));
        for (size_t i0 = 0; i0 < raw.size(); ++i0) {
            if (!raw[i0] || seen[i0]) continue;
            std::vector<size_t> part{i0};
            seen[i0] = 1;
            int x0 = w_, x1 = -1, y0 = h_, y1 = -1;
            for (size_t q = 0; q < part.size(); ++q) {
                const int x = static_cast<int>(part[q] % w_), y = static_cast<int>(part[q] / w_);
                x0 = std::min(x0, x), x1 = std::max(x1, x), y0 = std::min(y0, y), y1 = std::max(y1, y);
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        const int xx = x + dx, yy = y + dy;
                        const size_t j = size_t(yy) * w_ + xx;
                        if (xx >= 0 && yy >= 0 && xx < w_ && yy < h_ && raw[j] && !seen[j]) seen[j] = 1, part.push_back(j);
                    }
            }
            if (std::max(x1 - x0, y1 - y0) + 1 >= minExtent)
                for (size_t j : part) rim_[j] = 1;
            else if (static_cast<int>(part.size()) >= kRimBlobPixels) ++particles;
        }
        if (std::find(rim_.begin(), rim_.end(), 1) == rim_.end()) {
            if (particles >= kDotCount) result.error = "manufacture error (dot array)";
            return result;
        }

        // 2. The vertical trace.
        Trace v;
        v.id = 2, v.vertical = true;
        int vLeft = w_;  // its outline's left edge (the horizontal trace lies left of it)
        int cx = -1, vw = 0;
        if (centreBetweenRims(true, 0, w_, cx, vw)) {
            // Its ends: up and down the centre line from its middle, through the inside to the cap rims.
            int top = -1, bottom = -1;
            {
                int y = h_ / 2;
                while (y > 0 && !rimAt(cx, y)) --y;
                while (y > 0 && rimAt(cx, y - 1)) --y;
                top = y;
                y = h_ / 2;
                while (y < h_ - 1 && !rimAt(cx, y)) ++y;
                while (y < h_ - 1 && rimAt(cx, y + 1)) ++y;
                bottom = y;
            }
            // Its sides: the outline per row, left out where the rim walked through is over 1.5 times as thick
            // as that side's rim usually is (+3 px): there it runs on into the horizontal trace's rim, a neck or
            // debris. Of the rest, the 10th percentile from outside (a straight trace's edge, not a speck on it).
            struct Row { int l, r, il, ir; };
            std::vector<Row> rows;
            std::vector<double> thickL, thickR;
            for (int y = top; y <= bottom; ++y) {
                Row row{};
                if (!walk(true, y, cx, vw, row.l, row.r, &row.il, &row.ir)) continue;
                rows.push_back(row), thickL.push_back(row.il - row.l), thickR.push_back(row.r - row.ir);
            }
            int left = w_, right = -1;
            if (!rows.empty()) {
                const double usualL = median(thickL), usualR = median(thickR);
                std::vector<double> ls, rs;
                for (const Row& row : rows) {
                    if (row.il - row.l <= 1.5 * usualL + 3) ls.push_back(row.l);
                    if (row.r - row.ir <= 1.5 * usualR + 3) rs.push_back(-row.r);
                }
                auto tenth = [](std::vector<double> v) {
                    std::nth_element(v.begin(), v.begin() + v.size() / 10, v.end());
                    return v[v.size() / 10];
                };
                if (!ls.empty() && !rs.empty()) left = static_cast<int>(tenth(ls)), right = static_cast<int>(-tenth(rs));
            }
            if (right > left && bottom > top) {
                vLeft = left;
                v.acrossOut0 = std::max(0, left - kMargin), v.acrossOut1 = std::min(w_ - 1, right + kMargin);
                v.alongOut0 = std::max(0, top - kMargin), v.alongOut1 = std::min(h_ - 1, bottom + kMargin);
                const int W = right - left + 1, from = std::min(top + W, (top + bottom) / 2), to = std::max(bottom + 1 - W, from + 1);
                ridgesOf(v, profile(img_.px, w_, h_, true, from, to), left, right, top, bottom, false);
            }
        }
        const bool hasV = vLeft < w_;

        // 3. The horizontal trace.
        Trace t;
        t.id = 1;
        int cy = -1, hw = 0;
        const int split = hasV ? vLeft : w_;
        if (split > 16 && centreBetweenRims(false, 0, split, cy, hw)) {
            std::vector<int> tops, bottoms, innerTops, innerBottoms;  // the outline per column from the left edge on
            for (int x = 0; x < split; ++x) {
                int a, b, ia, ib;
                if (!walk(false, x, cy, hw, a, b, &ia, &ib)) break;
                tops.push_back(a), bottoms.push_back(b), innerTops.push_back(ia), innerBottoms.push_back(ib);
            }
            const int reached = static_cast<int>(tops.size()) - 1;  // the last column with the trace's outline
            if (reached >= minExtent / 2) {
                const int W0 = median(std::vector<double>(bottoms.begin(), bottoms.begin() + std::max(1, reached / 2))) -
                               median(std::vector<double>(tops.begin(), tops.begin() + std::max(1, reached / 2))) + 1;  // full width
                t.connected = hasV && reached >= split - 2;  // its outline runs on into the vertical trace's
                int end = reached, waist = -1;
                if (t.connected) {
                    // The waist: the lowest outline (median over 7 columns) within two widths of the vertical trace.
                    const int first = std::max(0, reached - 2 * W0), last = std::max(first, reached - kRimEdge);
                    double lowest = HUGE_VAL;
                    for (int x = first; x <= last; ++x) {
                        std::vector<double> near;
                        for (int k = std::max(first, x - 3); k <= std::min(last, x + 3); ++k) near.push_back(bottoms[k] - tops[k]);
                        const double m = median(near);
                        if (m < lowest - 1e-9) lowest = m, waist = x;
                    }
                    if (waist >= 0) end = waist;
                }
                // The box's rows: the outline up to its end, leaving out columns where the rim walked through is
                // over 1.5 times as thick as usual (+3 px): debris stuck on it.
                std::vector<double> thickT, thickB;
                for (int x = 0; x <= end; ++x) thickT.push_back(innerTops[x] - tops[x]), thickB.push_back(bottoms[x] - innerBottoms[x]);
                const double usualT = median(thickT), usualB = median(thickB);
                int y0 = h_, y1 = -1;
                for (int x = 0; x <= end; ++x) {
                    if (innerTops[x] - tops[x] <= 1.5 * usualT + 3) y0 = std::min(y0, tops[x]);
                    if (bottoms[x] - innerBottoms[x] <= 1.5 * usualB + 3) y1 = std::max(y1, bottoms[x]);
                }
                if (y1 < y0) y0 = tops[0], y1 = bottoms[0];
                t.alongOut0 = 0, t.alongOut1 = std::min(split - 1, end + (t.connected ? 0 : kMargin));
                t.acrossOut0 = std::max(0, y0 - kMargin), t.acrossOut1 = std::min(h_ - 1, y1 + kMargin);
                const int to = std::max(1, std::min(end - W0, end / 2 + 1));  // the full-width part, clear of the cap
                ridgesOf(t, profile(img_.px, w_, h_, false, 0, to), y0, y1, 0, end, true);
                if (t.connected && waist >= 0) {
                    // The neck: its width at the waist, ridge to ridge (the outline less the trace's own rim
                    // insets, measured on its full-width part); its box from the waist to the vertical trace.
                    const double inset0 = t.across0 - median(std::vector<double>(tops.begin(), tops.begin() + std::max(1, reached / 2)));
                    const double inset1 = median(std::vector<double>(bottoms.begin(), bottoms.begin() + std::max(1, reached / 2))) - t.across1;
                    t.neckAt = waist, t.neckFrom = tops[waist] + inset0, t.neckTo = std::max(t.neckFrom, bottoms[waist] - inset1);
                    t.neckWidth = t.neckTo - t.neckFrom;
                    t.along1 = waist;
                    t.neckX0 = waist + 1, t.neckX1 = (hasV ? v.acrossOut0 : w_) - 1;
                    int ny0 = static_cast<int>(t.neckFrom), ny1 = static_cast<int>(t.neckTo);
                    for (int x = waist; x <= reached; ++x) ny0 = std::min(ny0, tops[x]), ny1 = std::max(ny1, bottoms[x]);
                    t.neckY0 = std::max(0, ny0 - kMargin), t.neckY1 = std::min(h_ - 1, ny1 + kMargin);
                    t.neckFrom0 = t.neckY0, t.neckEdge0 = t.neckFrom, t.neckEdge1 = t.neckTo;
                    for (int y = static_cast<int>(t.neckY0); y <= t.neckY1; ++y) {
                        double sum = 0;
                        int n = 0;
                        for (int x = std::max(0, waist - 3); x <= std::min(w_ - 1, waist + 3); ++x) sum += img_.px[size_t(y) * w_ + x], ++n;
                        t.neckProfile.push_back(sum / n);
                    }
                }
                result.traces.push_back(t);
            }
        }
        if (hasV) result.traces.push_back(v);
        return result;
    }

private:
    bool rimAt(int x, int y) const { return x >= 0 && y >= 0 && x < w_ && y < h_ && rim_[size_t(y) * w_ + x]; }

    // The centre line between a trace's two rims: the two fullest columns (vertical) or rows (horizontal, over
    // columns [from, to)) of rim pixels, at least kRimTall of the image (vertical) or kRimMinExtent / 2
    // (horizontal) full and separated by the trace's inside. If no inside separates two (a thin trace, its rims
    // merged), the middle of the fullest band. Returns false if no column (row) is full enough.
    bool centreBetweenRims(bool vertical, int from, int to, int& centre, int& width) const {
        const int n = vertical ? w_ : h_;
        std::vector<int> count(n, 0);
        for (int i = 0; i < n; ++i)
            for (int k = vertical ? 0 : from; k < (vertical ? h_ : to); ++k) count[i] += vertical ? rimAt(i, k) : rimAt(k, i);
        const int need = static_cast<int>(vertical ? kRimTall * h_ : kRimMinExtent / 2 * std::min(w_, h_));
        const int p1 = static_cast<int>(std::max_element(count.begin(), count.end()) - count.begin());
        if (count[p1] < need) return false;
        // The other rim: the fullest position at least 4 px away whose run of rim-holding positions does not
        // reach p1 (i.e. separated from it by the trace's inside).
        int p2 = -1;
        for (int i = 0; i < n; ++i) {
            if (std::abs(i - p1) < 4 || count[i] < need) continue;
            bool apart = false;
            for (int k = std::min(i, p1); k <= std::max(i, p1) && !apart; ++k) apart = count[k] < count[p1] / 4;
            if (apart && (p2 < 0 || count[i] > count[p2])) p2 = i;
        }
        if (p2 < 0) {  // the rims too close to tell apart (a thin trace): its centre in the fullest band
            int lo = p1, hi = p1;
            while (lo > 0 && 2 * count[lo - 1] >= count[p1]) --lo;
            while (hi < n - 1 && 2 * count[hi + 1] >= count[p1]) ++hi;
            centre = (lo + hi) / 2, width = hi - lo + 1;
            return true;
        }
        centre = (p1 + p2) / 2, width = std::abs(p2 - p1);
        return true;
    }

    // From the centre (column c of row i for the vertical trace, row c of column i for the horizontal one) out
    // to each side: through the inside to the first rim, then through the rim; lo and hi are the rims' outer
    // edges (innerLo and innerHi, if given, their inner edges). If the centre itself is rim (a cap, or rims
    // merged), the rim run around it. False if either side has no rim within 2 widths of the centre.
    bool walk(bool vertical, int i, int c, int width, int& lo, int& hi, int* innerLo = nullptr, int* innerHi = nullptr) const {
        auto at = [&](int k) { return vertical ? rimAt(k, i) : rimAt(i, k); };
        const int n = vertical ? w_ : h_, reach = 2 * std::max(width, 4);
        lo = hi = c;
        if (!at(c)) {
            while (lo > 0 && c - lo < reach && !at(lo)) --lo;
            while (hi < n - 1 && hi - c < reach && !at(hi)) ++hi;
            if (!at(lo) || !at(hi)) return false;
        }
        if (innerLo) *innerLo = lo, *innerHi = hi;  // the rims' inner edges
        while (lo > 0 && at(lo - 1)) --lo;
        while (hi < n - 1 && at(hi + 1)) ++hi;
        return true;
    }

    // Width (ridges of the across profile within the outline [c0, c1]) and length (ridges of the profile along
    // the trace's inside, within [a0, a1]); the cross-section between the box's sides.
    void ridgesOf(Trace& t, const std::vector<double>& across, int c0, int c1, int a0, int a1, bool horizontal) {
        const int W = c1 - c0 + 1;
        std::tie(t.across0, t.across1) = ridges(across, c0, c1, W);
        int from = static_cast<int>(t.across0) + kRidgeInset, to = static_cast<int>(t.across1) + 1 - kRidgeInset;
        if (to <= from) from = static_cast<int>(t.across0), to = static_cast<int>(t.across1) + 1;
        const std::vector<double> along = profile(img_.px, w_, h_, horizontal, from, to);
        std::tie(t.along0, t.along1) = ridges(along, a0, a1, std::max(10, W));
        if (horizontal) t.along0 = 0;  // it runs on past the image's left edge
        t.crossFrom = t.acrossOut0;
        t.crossSection.assign(across.begin() + static_cast<int>(t.acrossOut0), across.begin() + static_cast<int>(t.acrossOut1) + 1);
        t.plotFirst = 0, t.plotLast = static_cast<int>(t.crossSection.size()) - 1;
    }

    const Gray& img_;
    int w_, h_;
    std::vector<uint8_t> rim_;  // the traces' rim pixels (particles dropped)
};

}  // namespace

Analysis analyzeImage(const Gray& img) {
    if (img.w < 16 || img.h < 16) return {};
    if (styleOf(img) == ImageStyle::Rims) return RimFinder(img).run();
    Analysis result;
    result.traces = Finder(img).run();
    return result;
}
