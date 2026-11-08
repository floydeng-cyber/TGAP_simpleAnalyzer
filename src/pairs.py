"""Step 2: pairs each pre-etch site with its post-etch image, with design values and Hitachi results.

Pairs are matched by site (chip, col, row), not by file name: in these two runs the same file name
names different sites. Design width/gap follow config.design(); the rule is checked against every
layout mask in the paired set. Hitachi values come from the post-etch .msr file.

Writes OUT/pairs.csv.
"""

import csv
import re
import sys
from collections import defaultdict

import numpy as np
from PIL import Image

import config
from extract import read_sites


def _runs(line: np.ndarray) -> list[tuple[int, int]]:
    d = np.diff(np.concatenate([[0], line.astype(np.int8), [0]]))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def mask_design(path) -> tuple[int, int] | None:
    """Design (vertical bar width, gap) of a TGAP layout mask, rounded to 50 nm; None if bars touch."""
    m = np.array(Image.open(path).convert("L")) > 127
    left = _runs(m[:, 2])
    if len(left) != 1:
        return None
    r0, r1 = left[0]
    row = _runs(m[(r0 + r1) // 2])
    if len(row) < 2:
        return None
    gap = row[1][0] - row[0][1]
    far = 20 if (r0 + r1) // 2 > 256 else 491
    v = [r for r in _runs(m[far]) if r[1] > row[1][0] - 5]
    if len(v) != 1:
        return None
    snap = lambda px: int(round(px * config.NM_PER_PX / 50) * 50)  # noqa: E731
    return snap(v[0][1] - v[0][0]), snap(gap)


def read_msr(path) -> dict[str, dict]:
    """Post-etch image name -> Hitachi line width (Mean'), LWR/LER (3 sigma) and gap.

    Each site has a width point (on its -01MS image) and a gap point (on its -02MS image) at the same chip
    coordinate. A result header's last field is the value in 0.1 nm (1883.0 -> 188.3 nm in Results.xlsx).
    """
    text = path.read_text(errors="replace")
    width_at, gap_at = {}, {}
    for b in re.split(r"\n\s*\$result ", text)[1:]:
        head = b.split("\n", 1)[0].split()
        p_no, key, value = int(head[2]), (head[3], head[4]), head[-1]
        kind = re.search(r"Measurement\s*:\s*(\w+)", b)
        kind = kind.group(1) if kind else ""
        info = dict(re.findall(r":\s*(Mean'|3Sigma(?:\([LR]\.E\.\))?)\s*:\s*([\d.]+) nm", b))
        if kind == "Width":
            width_at[key] = {
                "hitachi_p_no": p_no,
                "hitachi_width_nm": float(info["Mean'"]) if "Mean'" in info else "",
                "hitachi_lwr_3s_nm": float(info["3Sigma"]) if "3Sigma" in info else "",
                "hitachi_ler_left_3s_nm": float(info["3Sigma(L.E.)"]) if "3Sigma(L.E.)" in info else "",
                "hitachi_ler_right_3s_nm": float(info["3Sigma(R.E.)"]) if "3Sigma(R.E.)" in info else "",
            }
        elif kind == "Gap":
            tif = next((i for i in re.findall(r"&mp_image_name \d+ (\S+)", b) if i.endswith("-02MS.tif")), None)
            if tif:
                gap_at[key] = (tif[:-4], "" if value == "*" else float(value) / 10)
    return {name: {**width_at.get(key, {}), "hitachi_gap_nm": gap} for key, (name, gap) in gap_at.items()}


def main() -> None:
    sites = read_sites()
    post_at = defaultdict(list)
    for s in sites:
        if s["kind"] == "post":
            post_at[(s["chip"], s["col"], s["row"])].append(s["image"])
    masks = {p.stem: p for p in config.DESIGN_MASKS.glob("*/*.tif")}
    hitachi = read_msr(config.POST_MSR)

    rows, checked, disagree = [], 0, 0
    for s in sorted((s for s in sites if s["kind"] == "pre"), key=lambda s: (s["chip"], int(s["row"]), int(s["col"]))):
        col, row = int(s["col"]), int(s["row"])
        w, g = config.design(col, row)
        if s["image"] in masks and (d := mask_design(masks[s["image"]])) is not None:
            checked += 1
            disagree += d != (w, g)
        match = post_at.get((s["chip"], s["col"], s["row"]), [])
        post = match[0] if len(match) == 1 else ""
        rows.append({"chip": s["chip"], "col": col, "row": row, "design_width_nm": w, "design_gap_nm": g,
                     "pre_image": s["image"], "post_image": post, **hitachi.get(post, {})})
    if disagree:
        print(f"WARNING: {disagree} layout masks disagree with the design rule", file=sys.stderr)
    cols = ["chip", "col", "row", "design_width_nm", "design_gap_nm", "pre_image", "post_image", "hitachi_p_no",
            "hitachi_width_nm", "hitachi_gap_nm", "hitachi_lwr_3s_nm", "hitachi_ler_left_3s_nm", "hitachi_ler_right_3s_nm"]
    with open(config.out("pairs.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=cols)
        wr.writeheader()
        wr.writerows(rows)
    print(f"pairs: {sum(bool(r['post_image']) for r in rows)}/{len(rows)} paired; design rule checked on "
          f"{checked} layout masks, {disagree} disagree; Hitachi values for "
          f"{sum(r.get('hitachi_width_nm', '') != '' for r in rows)} sites")


def read_pairs() -> list[dict]:
    return list(csv.DictReader(open(config.OUT / "pairs.csv")))
