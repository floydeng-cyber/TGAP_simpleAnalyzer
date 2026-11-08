"""Step 6: summary tables, figures and headline numbers from OUT/metrics.csv.

Writes to OUT:
  summary_by_design.csv   medians per design width, and gaps per design gap
  hitachi_check.csv       our post-etch CDs next to Hitachi's, per site
  figures/fig_*.png       CD vs design, edge movement, gap vs design, BinDev heatmap, Hitachi check
  stats.json              headline numbers
"""

import csv
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

import config  # noqa: E402

WIDTHS, GAPS = config.WIDTHS, config.GAPS
# reference chart palette (light): categorical slots 1-3 (valid for any pairing), text and surface tokens
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
SEQ = ["#f4f8fd", "#d6e6f8", "#a9cbf0", "#6fa6e4", "#3b7fd2", "#1f5fae", "#13427c"]


def num(x) -> float:
    try:
        v = float(x)
        return v if np.isfinite(v) else np.nan
    except (TypeError, ValueError):
        return np.nan


def col(rows, k):
    return [num(r.get(k)) for r in rows]


def nanmed(v) -> float:
    v = [x for x in v if np.isfinite(x)]
    return float(np.median(v)) if v else float("nan")


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_visible(False)


def series(ax, xs, groups, color, label, offset=0.0):
    """Median with IQR bars per x, joined by a 2px line, with a direct label at the right end."""
    stats = []
    for g in groups:
        v = np.asarray([x for x in g if np.isfinite(x)])
        stats.append(np.percentile(v, [25, 50, 75]) if len(v) else [np.nan] * 3)
    q1, m, q3 = np.array(stats).T
    x = np.asarray(xs, float) + offset
    ax.errorbar(x, m, yerr=[m - q1, q3 - m], color=color, linewidth=2, marker="o", markersize=7,
                markeredgecolor=SURFACE, markeredgewidth=2, capsize=0, elinewidth=1.5, label=label)
    if np.isfinite(m[-1]):
        ax.annotate(label, (x[-1], m[-1]), xytext=(8, 0), textcoords="offset points", va="center", color=INK, fontsize=8)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(config.out("figures", name), dpi=160)
    plt.close(fig)


def main() -> dict:
    rows = list(csv.DictReader(open(config.OUT / "metrics.csv")))
    ok = [r for r in rows if r["pre_ok"] == "True" and r["post_ok"] == "True"]
    sep = [r for r in ok if r["pre_connected"] == "False" and r["post_connected"] == "False"]
    by_w = lambda k, rs=ok: [col([r for r in rs if int(r["design_width_nm"]) == w], k) for w in WIDTHS]  # noqa: E731
    by_g = lambda k, rs: [col([r for r in rs if int(r["design_gap_nm"]) == g], k) for g in GAPS]  # noqa: E731
    half = lambda k: [[v / 2 for v in g] for g in by_w(k)]  # noqa: E731

    # ---------- tables
    keys_w = ["pre_v_cd", "pre_h_cd", "post_v_inner_cd", "post_h_inner_cd", "post_v_cd", "post_h_cd", "hitachi_width_nm",
              "d_v_cd_bottom", "d_h_cd_bottom", "d_v_cd_top", "d_h_cd_top", "d_tip_extension",
              "pre_v_lwr3s", "post_v_lwr3s", "pre_h_lwr3s", "post_h_lwr3s", "pre_tip_radius", "post_tip_radius",
              "cmp_dev_bindev", "cmp_iou", "cmp_hausdorff"]
    keys_g = ["pre_gap", "post_gap", "d_gap", "hitachi_gap_nm"]
    table_w, table_g = [], []
    for w in WIDTHS:
        rs = [r for r in ok if int(r["design_width_nm"]) == w]
        table_w.append({"design_width_nm": w, "n": len(rs), **{k: nanmed(col(rs, k)) for k in keys_w}})
    for g in GAPS:
        rs = [r for r in ok if int(r["design_gap_nm"]) == g]
        s = [r for r in rs if r in sep]
        table_g.append({"design_gap_nm": g, "n_separated": len(s), **{k: nanmed(col(s, k)) for k in keys_g},
                        "n_neck_pre": sum(r["pre_connected"] == "True" for r in rs),
                        "n_neck_post": sum(r["post_connected"] == "True" for r in rs)})
    with open(config.out("summary_by_design.csv"), "w", newline="") as f:
        w = csv.writer(f)
        for table in (table_w, table_g):
            w.writerow(table[0].keys())
            w.writerows([[f"{v:.1f}" if isinstance(v, float) else v for v in r.values()] for r in table])
            w.writerow([])

    # ---------- Hitachi check
    hrows = [r for r in rows if r["post_ok"] == "True" and np.isfinite(num(r["hitachi_width_nm"]))]
    with open(config.out("hitachi_check.csv"), "w", newline="") as f:
        keys = ["chip", "col", "row", "post_image", "hitachi_width_nm", "post_v_cd", "post_h_cd", "hitachi_gap_nm", "post_gap"]
        w = csv.writer(f)
        w.writerow(keys)
        w.writerows([[r.get(k, "") for k in keys] for r in hrows])
    hit = np.array(col(hrows, "hitachi_width_nm"))
    wd = np.array(col(hrows, "design_width_nm"))
    hitachi = {}
    for t in ("v", "h"):
        ours = np.array(col(hrows, f"post_{t}_cd"))
        s = np.isfinite(ours)
        d = ours[s] - hit[s]
        resid = lambda y: y - np.polyval(np.polyfit(wd[s], y, 1), wd[s])  # noqa: E731
        hitachi[t] = {"n": int(s.sum()), "offset": float(np.median(d)),
                      "robust_sd": float(1.4826 * np.median(np.abs(d - np.median(d)))),
                      "r": float(np.corrcoef(ours[s], hit[s])[0, 1]),
                      "residual_r": float(np.corrcoef(resid(ours[s]), resid(hit[s]))[0, 1])}
    go, gh = np.array(col(hrows, "post_gap")), np.array(col(hrows, "hitachi_gap_nm"))
    s = np.isfinite(go) & np.isfinite(gh) & (go > 0)
    dg = go[s] - gh[s]
    hitachi["gap"] = {"n": int(s.sum()), "offset": float(np.median(dg)),
                      "robust_sd": float(1.4826 * np.median(np.abs(dg - np.median(dg)))),
                      "r": float(np.corrcoef(go[s], gh[s])[0, 1])}

    # ---------- figures
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE})

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, t, name in ((axes[0], "v", "Vertical trace"), (axes[1], "h", "Horizontal trace")):
        style(ax, f"{name}: CD vs design width", "Design width (nm)", "CD (nm), median and IQR")
        ax.plot([130, 420], [130, 420], color=INK2, linewidth=1, linestyle=(0, (4, 3)))
        ax.annotate("CD = design", (410, 410), xytext=(-4, 6), textcoords="offset points", ha="right", color=INK2, fontsize=8)
        series(ax, WIDTHS, by_w(f"post_{t}_cd"), C2, "Post-etch bottom", -3)
        series(ax, WIDTHS, by_w(f"pre_{t}_cd"), C1, "Pre-etch (resist)", 0)
        series(ax, WIDTHS, by_w(f"post_{t}_inner_cd"), C3, "Post-etch top", 3)
        ax.set_xticks(WIDTHS)
        ax.set_xlim(130, 470)
    axes[0].legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK)
    save(fig, "fig_cd_vs_design.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    style(axes[0], "Bottom (footprint): edge movement per side", "Design width (nm)", "Outward movement (nm), median and IQR")
    series(axes[0], WIDTHS, half("d_v_cd_bottom"), C1, "Vertical trace side", -4)
    series(axes[0], WIDTHS, half("d_h_cd_bottom"), C2, "Horizontal trace side", 0)
    series(axes[0], WIDTHS, by_w("d_tip_extension"), C3, "Tip (toward V)", 4)
    style(axes[1], "Top: edge movement per side", "Design width (nm)", "")
    series(axes[1], WIDTHS, half("d_v_cd_top"), C1, "Vertical trace side", -3)
    series(axes[1], WIDTHS, half("d_h_cd_top"), C2, "Horizontal trace side", 3)
    for ax in axes:
        ax.axhline(0, color=INK2, linewidth=1)
        ax.set_xticks(WIDTHS)
        ax.set_xlim(130, 490)
    axes[0].legend(loc="center right", frameon=False, fontsize=8, labelcolor=INK)
    save(fig, "fig_edge_shift.png")

    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    style(ax, "Tip-to-line gap vs design gap", "Design gap (nm)", "Gap (nm), median and IQR")
    ax.plot([130, 420], [130, 420], color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    series(ax, GAPS, by_g("pre_gap", sep), C1, "Pre-etch", -3)
    series(ax, GAPS, by_g("post_gap", sep), C2, "Post-etch (bottom)", 3)
    ax.set_xticks(GAPS)
    ax.set_xlim(130, 480)
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK)
    save(fig, "fig_gap_vs_design.png")

    heat = np.array([[nanmed([num(r["cmp_dev_bindev"]) for r in ok if int(r["design_gap_nm"]) == g and int(r["design_width_nm"]) == w])
                      for w in WIDTHS] for g in GAPS])
    fig, ax = plt.subplots(figsize=(6.4, 5))
    im = ax.imshow(heat, cmap=LinearSegmentedColormap.from_list("seq", SEQ), origin="lower")
    for i in range(6):
        for j in range(6):
            if np.isfinite(heat[i, j]):
                ax.text(j, i, f"{heat[i, j]:.3f}", ha="center", va="center", fontsize=8,
                        color=INK if heat[i, j] < np.nanpercentile(heat, 60) else "#ffffff")
    ax.set_xticks(range(6), WIDTHS)
    ax.set_yticks(range(6), GAPS)
    ax.set_xlabel("Design width (nm)", color=INK2, fontsize=9)
    ax.set_ylabel("Design gap (nm)", color=INK2, fontsize=9)
    ax.set_title("BinDev deviation, pre vs aligned post (median over chips)", color=INK, fontsize=10, loc="left")
    ax.tick_params(colors=INK2, labelsize=8, length=0)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.outline.set_visible(False)
    cb.ax.tick_params(colors=INK2, labelsize=8, length=0)
    save(fig, "fig_dev_heatmap.png")

    fig, ax = plt.subplots(figsize=(6.2, 5))
    style(ax, "Post-etch footprint CD vs Hitachi CD-SEM", "Hitachi line width (nm)", "Our footprint CD (nm)")
    lo, hi = np.nanmin(hit) - 20, np.nanmax(hit) + 20
    ax.plot([lo, hi], [lo, hi], color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate("equal", (hi, hi), xytext=(-4, 4), textcoords="offset points", ha="right", color=INK2, fontsize=8)
    for t, c, name in (("h", C1, "Horizontal trace"), ("v", C2, "Vertical trace")):
        ax.scatter(hit, col(hrows, f"post_{t}_cd"), s=22, color=c, edgecolors=SURFACE, linewidths=1,
                   label=f"{name}: offset {hitachi[t]['offset']:+.0f} nm, r={hitachi[t]['r']:.3f}")
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK)
    save(fig, "fig_hitachi.png")

    # ---------- headline numbers
    dev = col(ok, "cmp_dev_bindev")
    stats = {
        "n_sites": len(rows), "n_both_ok": len(ok),
        "not_measured": [{"kind": k, "image": r[f"{k}_image"], "why": r[f"{k}_why"]}
                         for r in rows for k in ("pre", "post") if r[f"{k}_ok"] != "True"],
        "necks_pre": sum(r["pre_connected"] == "True" for r in ok), "necks_post": sum(r["post_connected"] == "True" for r in ok),
        "iou_vs_classical": {k: nanmed(col(rows, f"{k}_iou_vs_classical")) for k in ("pre", "post")},
        "by_width": table_w, "by_gap": table_g, "hitachi": hitachi,
        "dev_bindev": {"median": nanmed(dev), "q1": float(np.nanpercentile(dev, 25)), "q3": float(np.nanpercentile(dev, 75))},
        "lwr_hitachi_median": nanmed(col(rows, "hitachi_lwr_3s_nm")),
        "neck_width_pre": [nanmed(col(ok, "pre_neck_width")), float(np.nanmin(col(ok, "pre_neck_width"))), float(np.nanmax(col(ok, "pre_neck_width")))],
        "neck_width_post": [nanmed(col(ok, "post_neck_width")), float(np.nanmin(col(ok, "post_neck_width"))), float(np.nanmax(col(ok, "post_neck_width")))],
    }
    json.dump(stats, open(config.out("stats.json"), "w"), indent=1)
    v, h = hitachi["v"], hitachi["h"]
    print(f"summarize: {len(ok)}/{len(rows)} sites measurable; Hitachi V r={v['r']:.3f} offset {v['offset']:+.1f} nm "
          f"(SD {v['robust_sd']:.1f}), H r={h['r']:.3f} offset {h['offset']:+.1f} nm; BinDev median {stats['dev_bindev']['median']:.3f}")
    return stats
