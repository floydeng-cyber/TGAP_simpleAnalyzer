"""Step 2: labels for the U-Net. Classical segmentation of every image, then a QC of each mask.

Segmentation gives the structure's footprint: everything inside the outer edge of its bright edge band,
by the same rule for both kinds of image.
  - pre-etch (resist): a dark inside with a thin bright edge, on a lighter background;
  - post-etch (Si):    a thick bright rim around an inside about as dark as the background.
The outer edge is where the smoothed intensity crosses halfway from the background to the edge band's
bright level (half maximum). Holes are filled after mirror-padding the image on every side, so that a
trace leaving the image (the horizontal trace always enters from the left edge) still closes into a ring.

A mask passes QC when
  - it holds a horizontal trace entering from the left edge and a vertical trace to its right (two
    pieces, or one where a neck joins them);
  - both traces have steady widths: along the middle of each trace, the width's spread and largest
    jump stay under relative limits (MAX_CV, MAX_JUMP), with a floor in px for thin traces.

Writes OUT/work/labels/{pre,post}/<image>.png (0/255) and OUT/work/labels_qc.csv.
"""

import csv
from concurrent.futures import ProcessPoolExecutor

import cv2
import numpy as np
import tifffile
from scipy import ndimage as ndi

import config
from pairs import read_sites

# ---------------------------------------------------------------- classical segmentation

PAD = 512


def robust_background(img: np.ndarray) -> tuple[float, float]:
    """Median and robust sigma of the background, re-estimated from pixels within 3 sigma until stable."""
    v = img.astype(np.float32).ravel()
    med, sig = float(np.median(v)), float(np.std(v))
    for _ in range(10):
        sel = v[np.abs(v - med) < 3 * sig]
        new_med = float(np.median(sel))
        new_sig = 1.4826 * float(np.median(np.abs(sel - new_med)))
        if abs(new_med - med) < 0.1 and abs(new_sig - sig) < 0.1:
            break
        med, sig = new_med, max(new_sig, 1e-3)
    return med, sig


MAX_HOLE_THICKNESS = 140  # px; a trace's inside is thinner than this, the background between traces is not


def fill_mirrored(mask: np.ndarray) -> np.ndarray:
    """Fills the insides of traces. Where a neck joins the traces, the mirrored copies can also close off
    the background in a corner; such holes are thick (a large circle fits inside) and are left unfilled.
    A trace's inside stays thin even where, through a neck, it runs on into the other trace's inside."""
    padded = np.pad(mask, PAD, mode="symmetric")
    filled = ndi.binary_fill_holes(padded)[PAD:-PAD, PAD:-PAD]
    holes = filled & ~mask
    n, lab = cv2.connectedComponents(holes.astype(np.uint8), connectivity=4)
    if n <= 1:
        return mask
    # thickness of each hole: twice its largest inscribed radius (image borders count as open, not walls)
    dist = cv2.distanceTransform(np.pad(holes, 1, constant_values=1).astype(np.uint8), cv2.DIST_L2, 5)[1:-1, 1:-1]
    thick = ndi.maximum(dist, lab, index=np.arange(n)) * 2
    ok = np.asarray(thick) < MAX_HOLE_THICKNESS
    ok[0] = False
    return mask | ok[lab]


def keep_large(mask: np.ndarray, min_area: int) -> np.ndarray:
    n, lab, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8), connectivity=8)
    keep = np.zeros(n, bool)
    keep[1:] = stats[1:, cv2.CC_STAT_AREA] >= min_area
    return keep[lab]


def segment(img: np.ndarray, kind: str) -> np.ndarray:
    """Footprint mask (bool) of the structures in a 512x512 uint8 CD-SEM image; kind is "pre" or "post"."""
    f = img.astype(np.float32)
    bg, sig = robust_background(img)
    if kind == "post":
        sm = cv2.GaussianBlur(f, (0, 0), 2.0)
        rim = np.percentile(sm, 98)
        edge = sm > bg + 0.5 * (rim - bg)
        edge = keep_large(edge, 2000)  # drop particles
        edge = cv2.morphologyEx(edge.astype(np.uint8), cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8)) > 0
        mask = fill_mirrored(edge)
    else:
        den8 = cv2.fastNlMeansDenoising(img, None, h=25, templateWindowSize=7, searchWindowSize=21)
        # flatten the slow background shading (a dark haze around traces, gradients across the field),
        # estimated from background pixels only: anything 3 sigma from the background, widened, is left out
        d = den8.astype(np.float32)
        bg0, sig0 = robust_background(den8)
        bgpx = (np.abs(d - bg0) < 3 * sig0).astype(np.uint8)
        bgpx = cv2.erode(bgpx, np.ones((9, 9), np.uint8)).astype(np.float32)
        num = cv2.GaussianBlur(d * bgpx, (0, 0), 30)
        den_w = cv2.GaussianBlur(bgpx, (0, 0), 30)
        shading = num / np.maximum(den_w, 1e-3)
        den = d - shading + 128
        bg, sig = robust_background(np.clip(den, 0, 255).astype(np.uint8))
        sm = cv2.GaussianBlur(den, (0, 0), 1.0)
        top = np.percentile(sm, 99.5)
        # The two edges of a trace differ a lot in brightness (the lower edge of a horizontal trace is
        # often dim), so the bar is set by the noise, not by the brightest edge: 4 sigma above the
        # background, i.e. where the edge rises out of the noise.
        edge = sm > bg + max(4 * sig, 0.15 * (top - bg))
        dark_sm = cv2.GaussianBlur(den, (0, 0), 2.0)
        dark_level = np.percentile(dark_sm, 1)
        inside = dark_sm < bg - 0.5 * (bg - dark_level)
        inside = keep_large(inside, 150)
        # bright edge pieces count only where they touch a dark inside (keeps out bright specks)
        near = cv2.dilate(inside.astype(np.uint8), np.ones((7, 7), np.uint8)) > 0
        n, lab = cv2.connectedComponents(edge.astype(np.uint8), connectivity=8)
        touching = np.zeros(n, bool)
        touching[np.unique(lab[near & edge])] = True
        touching[0] = False
        body = inside | touching[lab]
        body = cv2.morphologyEx(body.astype(np.uint8), cv2.MORPH_CLOSE, np.ones((5, 5), np.uint8)) > 0
        mask = fill_mirrored(body)
    mask = keep_large(mask, 1500)
    mask = ndi.binary_fill_holes(mask)
    mask = cv2.medianBlur(mask.astype(np.uint8) * 255, 5) > 0
    return mask


# ---------------------------------------------------------------- QC of a mask

MAX_CV = 0.12
MAX_JUMP = 0.35


def trace_widths(mask: np.ndarray) -> dict:
    """Splits a TGAP mask into its vertical and horizontal traces and measures their width profiles."""
    h, w = mask.shape
    col_fill = mask.sum(axis=0)
    # vertical trace: the columns that are filled over the longest stretch
    vcols = np.flatnonzero(col_fill > 0.5 * col_fill.max())
    if len(vcols) == 0:
        return {"ok": False, "why": "no vertical trace"}
    # take the widest contiguous run of such columns, right of the image middle third
    runs = np.split(vcols, np.flatnonzero(np.diff(vcols) > 1) + 1)
    vrun = max(runs, key=len)
    v0, v1 = int(vrun[0]), int(vrun[-1])
    vrows = np.flatnonzero(mask[:, (v0 + v1) // 2])
    if len(vrows) < 0.3 * h:
        return {"ok": False, "why": "vertical trace too short"}
    top, bot = int(vrows[0]), int(vrows[-1])
    # vertical width per row, in the middle 60% of its length, away from the horizontal trace
    left_rows = np.flatnonzero(mask[:, 2])
    if len(left_rows) == 0:
        return {"ok": False, "why": "no horizontal trace at left edge"}
    hr0, hr1 = int(left_rows[0]), int(left_rows[-1])
    lo, hi = top + int(0.2 * (bot - top)), bot - int(0.2 * (bot - top))
    vw = []
    for y in range(lo, hi):
        if hr0 - 40 <= y <= hr1 + 40:
            continue  # junction / neck region
        xs = np.flatnonzero(mask[y, max(0, v0 - 60):min(w, v1 + 60)])
        if len(xs):
            vw.append(xs[-1] - xs[0] + 1)
    # horizontal width per column, from the left edge to 40 px before the vertical trace's left side
    hw = []
    x_end = None
    for x in range(0, v0):
        ys = np.flatnonzero(mask[max(0, hr0 - 60):min(h, hr1 + 60), x])
        if len(ys) == 0:
            x_end = x
            break
        hw.append(ys[-1] - ys[0] + 1)
    connected = x_end is None
    # stop well before the tip, or, where a neck joins the traces, before the trace flares into it
    stop = (x_end - 40) if x_end is not None else int(0.6 * v0)
    hw = hw[5:max(6, stop)]
    out = {"ok": True, "why": "", "connected": connected}
    for key, arr in (("v", vw), ("h", hw)):
        arr = np.asarray(arr, float)
        if len(arr) < 20:
            return {"ok": False, "why": f"{key} trace profile too short"}
        med = float(np.median(arr))
        out[f"{key}_width_px"] = med
        out[f"{key}_cv"] = float(arr.std() / med)
        out[f"{key}_jump"] = float(np.max(np.abs(arr - med)) / med)
    n_parts = cv2.connectedComponents(mask.astype(np.uint8))[0] - 1
    out["n_parts"] = n_parts
    if n_parts not in (1, 2) or (n_parts == 1) != connected:
        out.update(ok=False, why=f"{n_parts} pieces")
    # relative limits, with a floor in px so that thin traces (about 10 px) are not failed for 1-2 px of jitter
    elif any(out[f"{k}_cv"] * out[f"{k}_width_px"] > max(MAX_CV * out[f"{k}_width_px"], 2.0) for k in "vh"):
        out.update(ok=False, why="unsteady width")
    elif any(out[f"{k}_jump"] * out[f"{k}_width_px"] > max(MAX_JUMP * out[f"{k}_width_px"], 6.0) for k in "vh"):
        out.update(ok=False, why="width jump")
    return out


def _one(args):
    kind, name = args
    img = tifffile.imread(config.image(kind, name))
    mask = segment(img, kind)
    cv2.imwrite(str(config.out("work/labels", kind, f"{name}.png")), mask.astype(np.uint8) * 255)
    return {"kind": kind, "image": name, **trace_widths(mask)}


def main() -> None:
    jobs = [(s["kind"], s["image"]) for s in read_sites()]
    with ProcessPoolExecutor(max_workers=6) as ex:
        rows = list(ex.map(_one, jobs, chunksize=8))
    cols = ["kind", "image", "ok", "why", "connected", "n_parts", "v_width_px", "v_cv", "v_jump", "h_width_px", "h_cv", "h_jump"]
    with open(config.out("work/labels_qc.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    for kind in ("pre", "post"):
        k = [r for r in rows if r["kind"] == kind]
        fails = {}
        for r in k:
            if not r["ok"]:
                fails[r["why"]] = fails.get(r["why"], 0) + 1
        print(f"labels: {kind} {sum(r['ok'] for r in k)}/{len(k)} pass QC; failures {fails}")


def read_qc() -> list[dict]:
    return list(csv.DictReader(open(config.OUT / "work/labels_qc.csv")))
