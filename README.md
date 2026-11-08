# etchAnalyzer

Compares pre-etch and post-etch CD-SEM images of Wafer 1 (TGAP)
- Segment:  u-net
- Measure:  critical dimensions, gap/neck
- Quantify: change through etch
- Check:    post-etch against Hitachi's measurements

## Where to put input data

Put data folders next to `etchAnalyzer`, under in the same folder:

```
<folder>/
├── etchAnalyzer/
├── 02_HitachiRun02_9waferProcess_preEtch/TGAP/
│   ├── TGAP_Part1/20200121_195959_01-TGAP_Jan20_W1/    (pre-etch images, unzipped)
│   └── postprocesseddata_paired_Ziwang/...             (layout masks)
└── 03_HitachiRun03_9waferProcess-postETCH/TGAP/
    ├── 20200317_122851_01-TGAP_Mar20_W1P1/              (post-etch images, unzipped)
    └── MSR file/20200317_122851_01-TGAP_Mar20_W1P1.msr  (Hitachi results)
```
Images are read in place, never copied or changed.
To use other locations, edit the paths at top of `src/config.py`.

## Run
Once, from `<folder>`:

```bash
python3 -m venv etchAnalyzer/.venv
```

```bash
etchAnalyzer/.venv/bin/pip install -r etchAnalyzer/requirements.txt
```

Then:

```bash
etchAnalyzer/.venv/bin/python -I etchAnalyzer/src/run.py
```

runs the program automatically.

## What happens

Everything is written to `<folder>/analysisOutput/etch_W1/`. Each step prints one or two lines.

| Step | Time | Writes |
|---|---|---|
| `pairs` | seconds | `work/sites.csv`: chip and cell of each image; `work/pairs.csv`: each pre-etch site with its post-etch image, design width/gap, Hitachi values |
| `labels` | ~1 min | `work/labels/`: classical masks used to train the U-Net; `work/labels_qc.csv`: which masks are clean |
| `train` | ~40 min (CPU) | `work/model/unet.pt`, `work/model/train_log.csv`. Skipped while `work/model/unet.pt` exists. |
| `analyze` | ~7 min | `work/masks/`, `metrics.csv`, `check/chip_*.png`, `check/seg_*.png` |
| `summarize` | seconds | `summary_by_design.csv`, `figures/fig_*.png`; prints the agreement with Hitachi |

The analysis includes:
- **`metrics.csv`:** one row per site, with pre-etch, post-etch and change values.
- **`check/chip_*.png`:** each image with its outlines, for checking by eye.
- **`check/seg_*.png`:** the U-Net masks against the classical ones.

The summary includes:
- **`summary_by_design.csv`:** medians per design width and gap.
- **`figures/`:** CD vs design, edge movement, gap, BinDev heatmap, and our post-etch CD vs Hitachi's.

Columns in `metrics.csv`:
- **`pre_*` / `post_*`:** one image's features.
- **`d_*`:** change through etch.
- **`cmp_*`:** shape comparison.
- **`hitachi_*`:** Hitachi's values.

All lengths are in nm (3.904 nm per pixel).
