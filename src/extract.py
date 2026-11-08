"""Step 1: extracts the -02MS.tif images of both runs and records each image's site.

Writes OUT/work/images/{pre,post}/*.tif and OUT/work/sites.csv: kind, image, chip, chip coordinate (um)
and the cell (col, row) on the chip's 6x6 grid, counted from the run's smallest coordinate.
"""

import csv
import re
import zipfile

import config


def extract_run(kind: str, zip_path, folder: str) -> list[dict]:
    z = zipfile.ZipFile(zip_path)
    names = set(z.namelist())
    pat = re.compile(re.escape(folder) + r"/([ES]\d\d_M\d{4}-02MS)\.tif_cnd/cond\.txt$")
    sites = []
    for member in sorted(names):
        m = pat.search(member)
        tif = member.replace(".tif_cnd/cond.txt", ".tif")
        if not m or tif not in names:
            continue
        text = z.read(member).decode(errors="replace")
        chip = re.search(r"Chip_number\s+(\S+)", text).group(1)
        cx, cy = map(int, re.search(r"Chip_coordinate\s+(\d+) um, (\d+) um", text).groups())
        config.out("work/images", kind, f"{m.group(1)}.tif").write_bytes(z.read(tif))
        sites.append({"kind": kind, "image": m.group(1), "chip": chip, "cx": cx, "cy": cy})
    x0, y0 = min(s["cx"] for s in sites), min(s["cy"] for s in sites)
    for s in sites:
        s["col"] = round((s["cx"] - x0) / config.CELL_PITCH_UM)
        s["row"] = round((s["cy"] - y0) / config.CELL_PITCH_UM)
    return sites


def main() -> None:
    sites = extract_run("pre", config.PRE_ZIP, config.PRE_FOLDER) + extract_run("post", config.POST_ZIP, config.POST_FOLDER)
    with open(config.out("work/sites.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["kind", "image", "chip", "cx", "cy", "col", "row"])
        w.writeheader()
        w.writerows(sites)
    for kind in ("pre", "post"):
        print(f"extract: {kind} {sum(s['kind'] == kind for s in sites)} images")


def read_sites() -> list[dict]:
    return list(csv.DictReader(open(config.OUT / "work/sites.csv")))
