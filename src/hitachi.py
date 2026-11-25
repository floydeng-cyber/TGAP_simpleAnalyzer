"""Hitachi CD-SEM's own critical dimensions: query by structure, wafer, chip and design cell.

Standalone: the Python standard library only, and data/hitachi_cd.csv (one row per Hitachi site, extracted from
Hitachi's .msr files). It does not need the images, the U-Net or the rest of the analyzer, and its values do not
change when the analyzer's own measurements are re-run.

What is in the data
  - pre-etch: every structure (TGAP, DTGAP, TTGAP, UL, ETE, TETE, FLG, TCL), wafers 1-9, 5 chips each
    ("-4,0", "0,-5", "0,0", "0,5", "4,0"), 6 x 6 design cells per chip;
  - post-etch: TGAP only, 17 chips per wafer (45 on wafer 1).
  - grid_col 1-6: the design's line width (narrowest at 1); grid_row 1-6: its gap (smallest at 1).
  - width, gap in nm; x_um, y_um: the site's position in Hitachi's coordinates.

How to read the values
  - Gap 0: Hitachi found the gap closed (bridged / neck), not a measured distance of zero. gaps(...) leaves these
    out unless bridged=True.
  - None: Hitachi did not measure the site (1-13 % of sites per structure; most on TCL).
  - FLG and TCL have no Gap.
  - TCL: Hitachi measures the middle of the three lines, not the outermost one.
  - UL row 1 and TETE on wafer 2: Hitachi counts the grey haze in front of a line end as line.
  - TTGAP grid col 1, rows 5-6: Hitachi's Gap falls while the design gap grows; treat with care.

Python:
    from hitachi import width, gap, d_width, d_gap, sites, widths, gaps, change
    width("TGAP", 1, "0,0", 4, 4)                              # 255.51 (pre); kind="post" -> 533.47
    d_width("TGAP", 1, "0,0", 4, 4)                            # 277.96
    sites("TGAP", wafer=1, chip="0,0", row=4, col=4)          # list of matching sites (dicts)
    widths("DTGAP", wafer=3, kind="pre")                       # list of widths, unmeasured left out
    gaps("TGAP", wafer=1, kind="post", bridged=False)          # list of gaps, Gap 0 and unmeasured left out
    change("TGAP", wafer=1)                                    # TGAP pre and post per cell, on the 5 shared chips

Command line:
    python src/hitachi.py TGAP 1 --chip 0,0 --row 4 --col 4
    python src/hitachi.py ETE 2 --kind pre --fields width gap
    python src/hitachi.py TGAP 1 --change --json
"""

from __future__ import annotations  # "float | None" annotations on Python 3.7-3.9 (macOS system python3)

import argparse
import csv
import json
import sys
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "hitachi_cd.csv"
STRUCTURES = ["TGAP", "DTGAP", "TTGAP", "UL", "ETE", "TETE", "FLG", "TCL"]
INTS = ("wafer", "grid_col", "grid_row")
FLOATS = ("x_um", "y_um", "width", "gap")


def _num(v: str):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else round(x, 2)  # NaN -> None


@lru_cache(maxsize=4)
def _load(path: str) -> tuple:
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for k in INTS:
                r[k] = int(r[k])
            for k in FLOATS:
                r[k] = _num(r[k])
            rows.append(r)
    return tuple(rows)


def _chip(chip) -> str | None:
    """'0,0', '0, 0', (0, 0) -> '0,0'."""
    if chip is None:
        return None
    if isinstance(chip, (tuple, list)):
        return f"{int(chip[0])},{int(chip[1])}"
    return ",".join(str(int(p)) for p in str(chip).replace("(", "").replace(")", "").split(","))


def sites(structure: str, wafer: int = None, kind: str = None, chip=None, row: int = None, col: int = None,
          path=None) -> list[dict]:
    """Every Hitachi site matching the filters (None: any). kind: "pre" or "post" (post exists for TGAP only).
    Each site is a dict: structure, wafer, kind, chip, grid_col, grid_row, x_um, y_um, width, gap."""
    s = structure.upper()
    if s not in STRUCTURES:
        raise ValueError(f"unknown structure {structure!r}; one of {', '.join(STRUCTURES)}")
    if kind is not None and kind not in ("pre", "post"):
        raise ValueError("kind is 'pre' or 'post'")
    c = _chip(chip)
    return [dict(r) for r in _load(str(path or DATA))
            if r["structure"] == s and (wafer is None or r["wafer"] == int(wafer))
            and (kind is None or r["kind"] == kind) and (c is None or r["chip"] == c)
            and (row is None or r["grid_row"] == int(row)) and (col is None or r["grid_col"] == int(col))]


def widths(structure: str, wafer: int = None, kind: str = "pre", chip=None, row: int = None, col: int = None,
           path=None) -> list[float]:
    """Hitachi's widths (nm) of the matching sites; unmeasured sites left out."""
    return [r["width"] for r in sites(structure, wafer, kind, chip, row, col, path) if r["width"] is not None]


def gaps(structure: str, wafer: int = None, kind: str = "pre", chip=None, row: int = None, col: int = None,
         bridged: bool = False, path=None) -> list[float]:
    """Hitachi's gaps (nm) of the matching sites; unmeasured sites left out, and Gap 0 (bridged) too unless
    bridged=True."""
    return [r["gap"] for r in sites(structure, wafer, kind, chip, row, col, path)
            if r["gap"] is not None and (bridged or r["gap"] > 0)]


def change(structure: str = "TGAP", wafer: int = None, chip=None, row: int = None, col: int = None,
           path=None) -> list[dict]:
    """Pre- and post-etch values per cell where Hitachi measured both (TGAP: the 5 shared chips), with
    d_width = post - pre and d_gap = post - pre (None where either is missing, or a gap is 0)."""
    pre = {(r["wafer"], r["chip"], r["grid_col"], r["grid_row"]): r
           for r in sites(structure, wafer, "pre", chip, row, col, path)}
    out = []
    for r in sites(structure, wafer, "post", chip, row, col, path):
        p = pre.get((r["wafer"], r["chip"], r["grid_col"], r["grid_row"]))
        if p is None:
            continue
        dw = round(r["width"] - p["width"], 2) if r["width"] is not None and p["width"] is not None else None
        dg = round(r["gap"] - p["gap"], 2) if r["gap"] and p["gap"] else None
        out.append({"structure": r["structure"], "wafer": r["wafer"], "chip": r["chip"],
                    "grid_col": r["grid_col"], "grid_row": r["grid_row"],
                    "pre_width": p["width"], "post_width": r["width"], "d_width": dw,
                    "pre_gap": p["gap"], "post_gap": r["gap"], "d_gap": dg})
    return out


def _one(structure, wafer, chip, row, col, key, kind=None, path=None):
    rows = (change(structure, wafer, chip, row, col, path) if kind is None
            else sites(structure, wafer, kind, chip, row, col, path))
    return rows[0][key] if len(rows) == 1 else None


def width(structure: str, wafer: int, chip, row: int, col: int, kind: str = "pre", path=None) -> float | None:
    """One site's width (nm), or None if not measured."""
    return _one(structure, wafer, chip, row, col, "width", kind, path)


def gap(structure: str, wafer: int, chip, row: int, col: int, kind: str = "pre", path=None) -> float | None:
    """One site's gap (nm), or None if not measured; 0 = bridged."""
    return _one(structure, wafer, chip, row, col, "gap", kind, path)


def d_width(structure: str, wafer: int, chip, row: int, col: int, path=None) -> float | None:
    """One cell's post - pre width (nm; TGAP), or None."""
    return _one(structure, wafer, chip, row, col, "d_width", None, path)


def d_gap(structure: str, wafer: int, chip, row: int, col: int, path=None) -> float | None:
    """One cell's post - pre gap (nm; TGAP), or None."""
    return _one(structure, wafer, chip, row, col, "d_gap", None, path)


def main() -> int:
    ap = argparse.ArgumentParser(description="Query Hitachi CD-SEM critical dimensions (data/hitachi_cd.csv).")
    ap.add_argument("structure", help=", ".join(STRUCTURES))
    ap.add_argument("wafer", type=int, nargs="?", help="wafer 1-9 (default: all)")
    ap.add_argument("--kind", choices=["pre", "post"], help="pre- or post-etch (default: both)")
    ap.add_argument("--chip", help='chip, e.g. "0,0" or "-4,0"')
    ap.add_argument("--row", "-r", type=int, help="grid row 1-6 (design gap)")
    ap.add_argument("--col", "-c", type=int, help="grid column 1-6 (design width)")
    ap.add_argument("--change", action="store_true", help="pre and post per cell, and their difference (TGAP)")
    ap.add_argument("--fields", "-f", nargs="+", help="only these fields (identifiers are always kept)")
    ap.add_argument("--json", "-j", action="store_true", help="JSON instead of CSV")
    ap.add_argument("--data", help=f"the CSV (default: {DATA})")
    a = ap.parse_args()
    try:
        rows = (change if a.change else sites)(a.structure, a.wafer, **({} if a.change else {"kind": a.kind}),
                                               chip=a.chip, row=a.row, col=a.col, path=a.data)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 2
    if not rows:
        print("no matching sites", file=sys.stderr)
        return 1
    if a.fields:
        ids = ["structure", "wafer", "kind", "chip", "grid_col", "grid_row"]
        rows = [{k: r[k] for k in ids + a.fields if k in r} for r in rows]
    if a.json:
        print(json.dumps(rows, indent=1))
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows({k: ("" if v is None else v) for k, v in r.items()} for r in rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
