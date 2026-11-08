"""Step 4: segments every pre/post pair with the U-Net, measures both masks and compares them.

Writes to OUT:
  work/masks/{pre,post}/<image>.png  U-Net footprint masks
  metrics.csv                      one row per site: pre features, post features, pre->post changes
  check/chip_<chip>.png            per site: pre + outline (green) | post + outline (magenta) | both outlines
  check/seg_*.png                  classical (green) vs U-Net (magenta) masks: on the images whose classical
                                   mask failed QC, and on the held-out chip
"""

import csv

import cv2
import numpy as np
import tifffile

import config
import unet
from labels import read_qc, trace_widths
from measure import align, compare, inner_cd, measure, shift
from pairs import read_pairs

GREEN, MAGENTA = (0, 255, 0), (255, 0, 255)


def outlines(img, masks_colors):
    rgb = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
    for m, color in masks_colors:
        cnts, _ = cv2.findContours(m.astype(np.uint8), cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(rgb, cnts, -1, color, 1)
    return rgb


def labelled(tile, text):
    tile = cv2.copyMakeBorder(tile, 16, 2, 2, 2, cv2.BORDER_CONSTANT, value=(255, 255, 255))
    cv2.putText(tile, text, (3, 12), cv2.FONT_HERSHEY_SIMPLEX, 0.38, (0, 0, 0), 1, cv2.LINE_AA)
    return tile


def grid(tiles, ncols):
    rows = [np.hstack(tiles[i:i + ncols]) for i in range(0, len(tiles), ncols)]
    width = max(r.shape[1] for r in rows)
    return np.vstack([np.pad(r, ((0, 0), (0, width - r.shape[1]), (0, 0)), constant_values=255) for r in rows])


def save_rgb(path, rgb):
    cv2.imwrite(str(config.out(path)), cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR))


def site_record(model, r: dict) -> tuple[dict, dict, dict]:
    """Masks, measurements and comparison of one pre/post pair."""
    rec = {k: r[k] for k in ("chip", "col", "row", "design_width_nm", "design_gap_nm", "pre_image", "post_image",
                             "hitachi_width_nm", "hitachi_gap_nm", "hitachi_lwr_3s_nm")}
    imgs, masks, ms = {}, {}, {}
    for kind in ("pre", "post"):
        name = r[f"{kind}_image"]
        img = tifffile.imread(config.image(kind, name))
        mask = unet.predict(model, img)
        cv2.imwrite(str(config.out("work/masks", kind, f"{name}.png")), mask.astype(np.uint8) * 255)
        label = cv2.imread(str(config.OUT / "work/labels" / kind / f"{name}.png"), cv2.IMREAD_GRAYSCALE) > 0
        rec[f"{kind}_iou_vs_classical"] = float((mask & label).sum() / max((mask | label).sum(), 1))
        m = measure(mask)
        if kind == "post":
            m.update(inner_cd(img, mask, m))
        rec.update({f"{kind}_{k}": v for k, v in m.items() if k != "v_rows"})
        imgs[kind], masks[kind], ms[kind] = img, mask, m
    pre_m, post_m = ms["pre"], ms["post"]
    if pre_m["ok"] and post_m["ok"]:
        rec.update({f"cmp_{k}": v for k, v in compare(masks["pre"], masks["post"], pre_m, post_m).items()})
        for t in ("v", "h"):
            rec[f"d_{t}_cd_bottom"] = post_m[f"{t}_cd"] - pre_m[f"{t}_cd"]
            if f"{t}_inner_cd" in post_m:
                rec[f"d_{t}_cd_top"] = post_m[f"{t}_inner_cd"] - pre_m[f"{t}_cd"]
        rec["d_gap"] = post_m["gap"] - pre_m["gap"]
        # tip extension: how far the tip moved toward V, net of V's own growth (alignment-free)
        if not pre_m["connected"] and not post_m["connected"]:
            rec["d_tip_extension"] = -(post_m["gap"] - pre_m["gap"]) - (post_m["v_cd"] - pre_m["v_cd"]) / 2
    return rec, imgs, masks


def qa_sheets(pairs: list[dict]) -> None:
    """Classical vs U-Net masks where the classical mask failed QC, and on the held-out chip."""
    in_pairs = {("pre", r["pre_image"]) for r in pairs} | {("post", r["post_image"]) for r in pairs}
    rejected = [(r["kind"], r["image"]) for r in read_qc() if r["ok"] != "True" and (r["kind"], r["image"]) in in_pairs]
    holdout = [(k, r[f"{k}_image"]) for r in pairs if r["chip"] == config.HOLDOUT_CHIP for k in ("pre", "post")][::3]
    for tag, items in (("classical_rejects", rejected), ("holdout_chip", holdout)):
        tiles, passed = [], 0
        for kind, name in items:
            img = tifffile.imread(config.image(kind, name))
            label = cv2.imread(str(config.OUT / "work/labels" / kind / f"{name}.png"), cv2.IMREAD_GRAYSCALE) > 0
            mask = cv2.imread(str(config.OUT / "work/masks" / kind / f"{name}.png"), cv2.IMREAD_GRAYSCALE) > 0
            ok = trace_widths(mask)["ok"]
            passed += ok
            iou = (label & mask).sum() / max((label | mask).sum(), 1)
            tile = cv2.resize(outlines(img, [(label, GREEN), (mask, MAGENTA)]), (220, 220), interpolation=cv2.INTER_AREA)
            tiles.append(labelled(tile, f"{kind} {name[:9]} QC {'pass' if ok else 'FAIL'} IoU {iou:.2f}"))
        if tiles:
            save_rgb(f"check/seg_{tag}.png", grid(tiles, 6))
        print(f"analyze: {tag}: {len(items)} images, U-Net mask passes QC on {passed}")


def main() -> None:
    model, ck = unet.load()
    pairs = read_pairs()
    rows, tiles = [], {}
    for i, r in enumerate(pairs):
        rec, imgs, masks = site_record(model, r)
        rows.append(rec)
        if "cmp_shift_dy_px" in rec:
            dy, dx = rec["cmp_shift_dy_px"], rec["cmp_shift_dx_px"]
        else:
            dy, dx, _ = align(masks["pre"], masks["post"])
        three = [outlines(imgs["pre"], [(masks["pre"], GREEN)]), outlines(imgs["post"], [(masks["post"], MAGENTA)]),
                 outlines(imgs["pre"], [(masks["pre"], GREEN), (shift(masks["post"], dy, dx), MAGENTA)])]
        tile = np.hstack([cv2.resize(t, (200, 200), interpolation=cv2.INTER_AREA) for t in three])
        text = f"W{r['design_width_nm']} G{r['design_gap_nm']} {r['pre_image'][:9]}|{r['post_image'][:9]}"
        tiles.setdefault(r["chip"], []).append(((int(r["row"]), int(r["col"])), labelled(tile, text)))
        if (i + 1) % 60 == 0:
            print(f"analyze: {i + 1}/{len(pairs)} sites", flush=True)

    lead = ["chip", "col", "row", "design_width_nm", "design_gap_nm", "pre_image", "post_image"]
    cols = lead + sorted({k for r in rows for k in r} - set(lead))
    with open(config.out("metrics.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    for chip, items in tiles.items():
        save_rgb(f"check/chip_{chip.replace(',', '_')}.png", grid([t for _, t in sorted(items, key=lambda x: x[0])], 6))
    print(f"analyze: U-Net from step {ck['step']}; {len(rows)} sites; measurable pre "
          f"{sum(r['pre_ok'] for r in rows)}, post {sum(r['post_ok'] for r in rows)}")
    qa_sheets(pairs)
