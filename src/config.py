"""Paths and constants. The analyzer only reads the data folders; everything it writes goes to OUT.

Like floAnalyzer, outputs go to analysisOutput/ in the folder holding etchAnalyzer, wherever that is.
"""

from pathlib import Path

ANALYZER = Path(__file__).resolve().parent.parent
PROJECT = ANALYZER.parent
OUT = PROJECT / "analysisOutput" / "etch_W1"  # replaced by run.py -o

# --- inputs (read only): Wafer 1, TGAP
IMAGES = {  # image folders (unzipped): *-02MS.tif, each with its *-02MS.tif_cnd/cond.txt
    "pre": PROJECT / "02_HitachiRun02_9waferProcess_preEtch/TGAP/TGAP_Part1/20200121_195959_01-TGAP_Jan20_W1",
    "post": PROJECT / "03_HitachiRun03_9waferProcess-postETCH/TGAP/20200317_122851_01-TGAP_Mar20_W1P1",
}
POST_MSR = PROJECT / "03_HitachiRun03_9waferProcess-postETCH/TGAP/MSR file/20200317_122851_01-TGAP_Mar20_W1P1.msr"
DESIGN_MASKS = PROJECT / "02_HitachiRun02_9waferProcess_preEtch/TGAP/postprocesseddata_paired_Ziwang/PreEtch_9wafer_TGAP_ZiwangPairs/A"

# --- geometry
NM_PER_PX = 1999.0 / 512  # Field_Size 1.999 um over 512 px
CELL_PITCH_UM = 50  # test cells on a chip: 6 x 6 grid
WIDTHS = [150, 200, 250, 300, 350, 400]
GAPS = [150, 200, 250, 300, 350, 400]


def design(col: int, row: int) -> tuple[int, int]:
    """Design width and gap (nm) of a grid cell; checked against the layout masks by pairs.py."""
    return 150 + 50 * col, 150 + 50 * row


# --- U-Net
HOLDOUT_CHIP = "0,5"  # held out of training, to measure agreement on unseen sites


def image(kind: str, name: str) -> Path:
    """Path of an input image; kind is "pre" or "post", name e.g. "E04_M0002-02MS"."""
    return IMAGES[kind] / f"{name}.tif"


def out(*parts: str) -> Path:
    p = OUT.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
