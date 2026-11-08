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
│   └── ...
└── 03_HitachiRun03_9waferProcess-postETCH/TGAP/
    └── ..
```
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
| `extract` | seconds | `work/images/` (the Wafer 1 `.tif` files), `work/sites.csv` (chip and cell of each image) |
| `pairs` | seconds | `pairs.csv`: each pre-etch site with its post-etch image, design width/gap, Hitachi values |
| `labels` | ~1 min | `work/labels/`: classical masks used to train the U-Net; `work/labels_qc.csv`: which masks are clean |
| `train` | ~40 min (CPU) | `model/unet.pt`, `model/train_log.csv`. Skipped while `model/unet.pt` exists. |
| `analyze` | ~7 min | `masks/`, `metrics.csv`, `overlays/chip_*.png`, `qa/seg_*.png` |
| `summarize` | seconds | `summary_by_design.csv`, `hitachi_check.csv`, `figures/fig_*.png`, `stats.json` |

The analysis includes:
- **`metrics.csv`:** one row per site, with pre-etch, post-etch and change values.
- **`overlays/chip_*.png`:** each image with its outlines, for checking by eye.
- **`qa/seg_*.png`:** the U-Net masks against the classical ones.

The summary includes:
- **`summary_by_design.csv`:** medians per design width and gap.
- **`hitachi_check.csv`:** our post-etch values next to Hitachi's.

Columns in `metrics.csv`:
- **`pre_*` / `post_*`:** one image's features.
- **`d_*`:** change through etch.
- **`cmp_*`:** shape comparison.
- **`hitachi_*`:** Hitachi's values.

All lengths are in nm (3.904 nm per pixel).
