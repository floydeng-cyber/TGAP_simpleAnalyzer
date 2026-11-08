"""Used by step 4: features and dimensions of a TGAP footprint mask, and the pre/post shape comparison.

A TGAP site: a horizontal trace (H) entering from the left edge, a vertical trace (V) to its right, and a
designed tip-to-line gap between them, which a neck may bridge. All lengths in nm.

Per mask:
  - V and H: CD (median width), LWR (3 sigma of the width), LER of each edge (3 sigma of its residual
    about a straight-line fit), on the trace's straight part;
  - V length; H tip position; tip radius (circle fitted to the tip outline);
  - gap: from H's tip to V's facing edge along H's centre line (0 where a neck joins them);
  - neck width: narrowest H width before V where H and V are one piece;
  - post-etch only: top CD, between the inner edges of the bright rims (inner_cd).
Per pre/post pair (compare):
  - Dev (BinDev): fraction of differing pixels after BinDev's XNOR-maximizing shift;
  - IoU, XOR fraction, Hausdorff distance and area gained/lost, with the traces' centre lines aligned.
"""

import numpy as np
import cv2
from scipy import ndimage as ndi

from config import NM_PER_PX

MAX_SHIFT = 48  # px searched in each direction when aligning post onto pre


def _edges_rows(mask, rows, c0, c1):
    """Left and right edge x of the mask on each row, within columns [c0, c1)."""
    out = []
    for y in rows:
        xs = np.flatnonzero(mask[y, c0:c1])
        if len(xs):
            out.append((y, c0 + xs[0], c0 + xs[-1] + 1))
    return np.array(out, float).reshape(-1, 3)


def _line_stats(pos, e0, e1):
    """CD, LWR and LER of a straight trace from its two edges e0, e1 sampled at positions pos."""
    w = e1 - e0
    res = []
    for e in (e0, e1):
        fit = np.polyfit(pos, e, 1)
        res.append(e - np.polyval(fit, pos))
    return {
        "cd": float(np.median(w)) * NM_PER_PX,
        "lwr3s": float(3 * np.std(w)) * NM_PER_PX,
        "ler0_3s": float(3 * np.std(res[0])) * NM_PER_PX,
        "ler1_3s": float(3 * np.std(res[1])) * NM_PER_PX,
        "n": int(len(w)),
    }


def _circle_fit(pts):
    """Algebraic least-squares circle through points (x, y): radius."""
    x, y = pts[:, 0], pts[:, 1]
    A = np.c_[2 * x, 2 * y, np.ones_like(x)]
    b = x**2 + y**2
    (cx, cy, c), *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.sqrt(c + cx**2 + cy**2))


def measure(mask: np.ndarray) -> dict:
    """Features and dimensions of one TGAP footprint mask (bool, 512x512)."""
    mask = mask.astype(bool)
    h, w = mask.shape
    out = {"ok": False}
    col_fill = mask.sum(axis=0)
    vcols = np.flatnonzero(col_fill > 0.5 * col_fill.max())
    if len(vcols) == 0:
        return {**out, "why": "no vertical trace"}
    runs = np.split(vcols, np.flatnonzero(np.diff(vcols) > 1) + 1)
    vrun = max(runs, key=len)
    vc = (vrun[0] + vrun[-1]) // 2
    vrows = np.flatnonzero(mask[:, vc])
    left_rows = np.flatnonzero(mask[:, 1])
    if len(vrows) < 0.3 * h or len(left_rows) == 0:
        return {**out, "why": "missing trace"}
    top, bot = int(vrows[0]), int(vrows[-1]) + 1
    hr0, hr1 = int(left_rows[0]), int(left_rows[-1]) + 1

    # --- vertical trace: rows of its straight part, away from its ends and from the junction
    lo, hi = top + int(0.12 * (bot - top)), bot - int(0.12 * (bot - top))
    # junction rows: where something (H's tip, its flare, a neck, a bump) sits just left of V, plus 15 rows
    probe = mask[:, max(0, int(vrun[0]) - 6):max(1, int(vrun[0]) - 2)].any(axis=1)
    near_h = np.zeros(h, bool)
    near_h[max(0, hr0 - 15):min(h, hr1 + 15)] = True  # H's own rows, even if its tip ends far away
    junction = ndi.binary_dilation(probe & ndi.binary_dilation(near_h, iterations=80), iterations=15) | near_h
    rows = [y for y in range(lo, hi) if not junction[y]]
    ve = _edges_rows(mask, rows, max(0, vrun[0] - 80), min(w, vrun[-1] + 80))
    if len(ve) < 30:
        return {**out, "why": "vertical trace too short"}
    v = _line_stats(ve[:, 0], ve[:, 1], ve[:, 2])
    v_left = float(np.median(ve[:, 1]))
    v_right = float(np.median(ve[:, 2]))

    # --- horizontal trace: joined to V by a neck if they are one connected piece
    _, lab = cv2.connectedComponents(mask.astype(np.uint8), connectivity=8)
    h_id = lab[int(left_rows[len(left_rows) // 2]), 1]
    v_id = lab[(top + bot) // 2, vc]
    connected = h_id == v_id
    hmask = lab == h_id
    if not connected:
        hmask_cols = np.flatnonzero(hmask.any(axis=0))
        x_tip = int(hmask_cols[-1]) + 1
        # another piece on H's rows between this tip and V: H is broken, so its tip and gap are not real
        between = np.zeros_like(mask)
        between[max(0, hr0 - 5):min(h, hr1 + 5), x_tip:max(x_tip, int(vrun[0]) - 2)] = True
        others = np.unique(lab[between & mask & ~hmask])
        others = [i for i in others if i not in (0, v_id) and (lab == i).sum() >= 50]
        if others:
            return {**out, "why": "broken horizontal trace"}
    else:
        x_tip = None
    # columns from the left edge to the tip (or to V, through the neck)
    band0, band1 = max(0, hr0 - 60), min(h, hr1 + 60)
    cols_h, tops_h, bots_h = [], [], []
    for x in range(0, int(v_left) if connected else x_tip):
        src = mask if connected else hmask
        ys = np.flatnonzero(src[band0:band1, x])
        if len(ys) == 0:
            break
        cols_h.append(x)
        tops_h.append(band0 + ys[0])
        bots_h.append(band0 + ys[-1] + 1)
    cols_h, tops_h, bots_h = map(np.asarray, (cols_h, tops_h, bots_h))
    widths_h = bots_h - tops_h
    w_ref = float(np.median(widths_h[5:max(6, len(widths_h) // 2)]))
    if connected:
        straight_end = int(0.6 * v_left)
    else:
        straight_end = x_tip - int(max(1.0 * w_ref, 20))
    sel = (cols_h >= 5) & (cols_h < straight_end)
    if sel.sum() < 30:
        return {**out, "why": "horizontal trace too short"}
    hs = _line_stats(cols_h[sel].astype(float), tops_h[sel].astype(float), bots_h[sel].astype(float))

    out.update(ok=True, why="", connected=connected)
    for k, val in v.items():
        out[f"v_{k}"] = val
    for k, val in hs.items():
        out[f"h_{k}"] = val
    out["v_length"] = (bot - top) * NM_PER_PX
    out["v_rows"] = [int(y) for y in ve[:, 0]]  # the straight-part rows the V CD was measured on
    out["v_left_x"] = v_left
    out["v_right_x"] = v_right
    out["h_center_y"] = float(np.median((tops_h[sel] + bots_h[sel]) / 2))
    out["area"] = float(mask.sum()) * NM_PER_PX**2 / 1e6  # um^2

    if connected:
        # neck: narrowest H width over the last part before V (the waist)
        part = (cols_h >= straight_end) & (cols_h < v_left - 2)
        if part.any():
            k = np.argmin(widths_h[part])
            out["neck_width"] = float(widths_h[part][k]) * NM_PER_PX
            out["neck_x"] = float(cols_h[part][k])
        out["gap"] = 0.0
        out["h_tip_x"] = float("nan")
        out["tip_radius"] = float("nan")
    else:
        # gap along H's centre line: from the tip to V's left edge at that height
        yc = int(round(out["h_center_y"]))
        xs = np.flatnonzero(hmask[yc])
        tip_x = xs[-1] + 1 if len(xs) else x_tip
        # V's facing edge at that height (V's own piece only)
        vxs = np.flatnonzero((lab[yc] == v_id)[int(tip_x):])
        v_at = tip_x + vxs[0] if len(vxs) else v_left
        out["h_tip_x"] = float(tip_x)
        out["gap"] = float(v_at - tip_x) * NM_PER_PX
        # tip radius: outline points within one half-width of the tip
        cnts, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        pts = np.concatenate([c[:, 0, :] for c in cnts]).astype(float)
        near = pts[(pts[:, 0] > tip_x - 0.5 * w_ref) & (pts[:, 0] <= tip_x + 1)
                   & (np.abs(pts[:, 1] - out["h_center_y"]) < w_ref)]
        out["tip_radius"] = _circle_fit(near) * NM_PER_PX if len(near) >= 8 else float("nan")
    return out


def align(pre: np.ndarray, post: np.ndarray) -> tuple[int, int, float]:
    """Shift (dy, dx) of post onto pre, within MAX_SHIFT, maximizing pixel agreement over pre's central window.

    The window excludes a MAX_SHIFT margin, so every shift compares the same number of pixels; post is
    extended by its edge pixels, exactly as shift() does. Among equally good shifts, the smallest is taken.
    Returns the shift and its agreement (fraction of the window's pixels that match).
    """
    m = MAX_SHIFT
    window = pre[m:-m, m:-m].astype(np.float32) * 2 - 1
    search = np.pad(post, m, mode="edge").astype(np.float32)[m:3 * m + window.shape[0], m:3 * m + window.shape[1]] * 2 - 1
    score = cv2.matchTemplate(search, window, cv2.TM_CCORR)  # sum of +-1 products: matches - mismatches
    dys, dxs = np.nonzero(score >= score.max() - 0.5)
    k = np.argmin((m - dys) ** 2 + (m - dxs) ** 2)  # score index i means shift m - i
    agreement = (score[dys[k], dxs[k]] / window.size + 1) / 2
    return int(m - dys[k]), int(m - dxs[k]), float(agreement)


def shift(mask: np.ndarray, dy: int, dx: int) -> np.ndarray:
    """Shifts a mask, extending edge pixels into the uncovered border (traces run on past the field)."""
    m = np.pad(mask, MAX_SHIFT, mode="edge")
    return m[MAX_SHIFT - dy:MAX_SHIFT - dy + mask.shape[0], MAX_SHIFT - dx:MAX_SHIFT - dx + mask.shape[1]]


def signed_distance(mask: np.ndarray) -> np.ndarray:
    """Signed distance (px) to the mask boundary: positive outside, negative inside."""
    return ndi.distance_transform_edt(~mask) - ndi.distance_transform_edt(mask)


def compare(pre: np.ndarray, post: np.ndarray, pre_m: dict | None = None, post_m: dict | None = None) -> dict:
    """Shape difference of a pre/post pair of footprint masks.

    Two alignments:
      - BinDev's: the shift maximizing pixel agreement (XNOR), used for dev_bindev as defined in the paper;
      - centre lines: V's centre line onto V's, H's onto H's (etch assumed symmetric about a trace's centre).
        Used for everything else. Pixel agreement alone is unreliable here: the post-etch footprint is about
        twice as wide, so many shifts score almost the same and edge shifts would trade off left against right.
    """
    pre, post = pre.astype(bool), post.astype(bool)
    m = MAX_SHIFT
    by, bx, _ = align(pre, post)
    xor_b = pre[m:-m, m:-m] ^ shift(post, by, bx)[m:-m, m:-m]
    if pre_m and post_m and pre_m.get("ok") and post_m.get("ok"):
        dx = int(round((pre_m["v_left_x"] + pre_m["v_right_x"]) / 2 - (post_m["v_left_x"] + post_m["v_right_x"]) / 2))
        dy = int(round(pre_m["h_center_y"] - post_m["h_center_y"]))
        dx, dy = int(np.clip(dx, -m, m)), int(np.clip(dy, -m, m))
        how = "centre lines"
    else:
        dy, dx, how = by, bx, "xnor"
    post_a = shift(post, dy, dx)
    # compare away from a margin, where the shifted post mask was extended
    P, Q = pre[m:-m, m:-m], post_a[m:-m, m:-m]
    xor = P ^ Q
    out = {
        "bindev_shift_dy_px": by, "bindev_shift_dx_px": bx,
        "dev_bindev": float(xor_b.mean()),
        "align": how, "shift_dy_px": dy, "shift_dx_px": dx,
        "xor_frac_centre_aligned": float(xor.mean()),
        "iou": float((P & Q).sum() / max((P | Q).sum(), 1)),
        "area_gain": float((Q & ~P).sum()) * NM_PER_PX**2 / 1e6,
        "area_loss": float((P & ~Q).sum()) * NM_PER_PX**2 / 1e6,
    }
    # Hausdorff distance between the two outlines, within the compared field
    def outline(mask):
        cnts, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        pts = np.concatenate([c[:, 0, :] for c in cnts])
        return pts[(pts[:, 0] >= m) & (pts[:, 0] < 512 - m) & (pts[:, 1] >= m) & (pts[:, 1] < 512 - m)]

    p_pts, q_pts = outline(pre), outline(post_a)
    if len(p_pts) and len(q_pts):
        d1 = np.abs(signed_distance(post_a)[p_pts[:, 1], p_pts[:, 0]]).max()
        d2 = np.abs(signed_distance(pre)[q_pts[:, 1], q_pts[:, 0]]).max()
        out["hausdorff"] = float(max(d1, d2)) * NM_PER_PX
    return out


def _inner_edges(profile: np.ndarray) -> tuple[float, float] | None:
    """Inner edges of a bright-rim profile across a post-etch trace (outer edge to outer edge): on each
    side, where the intensity falls halfway from the rim's peak to the dark inside's floor."""
    n = len(profile)
    if n < 12:
        return None
    half = n // 2
    pl = int(np.argmax(profile[:half]))
    pr = half + int(np.argmax(profile[half:]))
    if pr - pl < 4:
        return None
    floor_i = pl + int(np.argmin(profile[pl:pr + 1]))
    floor = profile[floor_i]
    lvl_l = (profile[pl] + floor) / 2
    lvl_r = (profile[pr] + floor) / 2
    # walk inward from each peak to the half level (linear interpolation)
    xl = None
    for i in range(pl, floor_i):
        if profile[i + 1] <= lvl_l:
            xl = i + (profile[i] - lvl_l) / max(profile[i] - profile[i + 1], 1e-6)
            break
    xr = None
    for i in range(pr, floor_i, -1):
        if profile[i - 1] <= lvl_r:
            xr = i - (profile[i] - lvl_r) / max(profile[i] - profile[i - 1], 1e-6)
            break
    if xl is None or xr is None or xr <= xl:
        return None
    return xl, xr


def inner_cd(img: np.ndarray, mask: np.ndarray, m: dict) -> dict:
    """Post-etch top CD: the width of the dark inside between the rims' inner edges (top of the sloped
    sidewalls), on the same straight parts of V and H that measure() used. Profiles are averaged over 5
    rows/columns."""
    if not m.get("ok"):
        return {}
    sm = cv2.GaussianBlur(img.astype(np.float32), (0, 0), 1.5)
    mask = mask.astype(bool)
    out = {}
    h = mask.shape[0]
    # V: the straight-part rows measure() used for the outer CD
    v_rows = m["v_rows"]
    ws = []
    for y in v_rows[::3]:
        y0, y1 = max(0, y - 2), min(h, y + 3)
        xs = np.flatnonzero(mask[y, int(m["v_left_x"]) - 10:int(m["v_right_x"]) + 10])
        if len(xs) < 12:
            continue
        a = int(m["v_left_x"]) - 10 + xs[0]
        b = int(m["v_left_x"]) - 10 + xs[-1] + 1
        e = _inner_edges(sm[y0:y1, a:b].mean(axis=0))
        if e:
            ws.append(e[1] - e[0])
    if len(ws) >= 10:
        out["v_inner_cd"] = float(np.median(ws)) * NM_PER_PX
    # H: columns of its straight part
    yc = m["h_center_y"]
    hw = m["h_cd"] / NM_PER_PX
    x_end = (m["h_tip_x"] - 1.0 * hw - 20) if not m["connected"] else 0.6 * m["v_left_x"]
    ws = []
    for x in range(5, int(x_end), 3):
        r0, r1 = int(yc - hw), int(yc + hw)
        col = mask[max(0, r0):min(h, r1), x]
        ys = np.flatnonzero(col)
        if len(ys) < 12:
            continue
        a, b = max(0, r0) + ys[0], max(0, r0) + ys[-1] + 1
        e = _inner_edges(sm[a:b, max(0, x - 2):x + 3].mean(axis=1))
        if e:
            ws.append(e[1] - e[0])
    if len(ws) >= 10:
        out["h_inner_cd"] = float(np.median(ws)) * NM_PER_PX
    return out
