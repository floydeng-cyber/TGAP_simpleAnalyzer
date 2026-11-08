"""etchAnalyzer: Part 1, pre- vs post-etch CD-SEM comparison of Wafer 1 (TGAP).

Usage (from the folder holding etchAnalyzer, or anywhere):
  etchAnalyzer/.venv/bin/python -I etchAnalyzer/src/run.py                 # every step whose output is missing
  etchAnalyzer/.venv/bin/python -I etchAnalyzer/src/run.py analyze summarize  # just these steps
  etchAnalyzer/.venv/bin/python -I etchAnalyzer/src/run.py --force         # every step, again
  ... run.py -o DIR                                                        # write to DIR instead

Steps, in order: extract, pairs, labels, train, analyze, summarize.
Training takes about 40 minutes on the CPU and is skipped while OUT/model/unet.pt exists.
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # -I leaves the script's folder off sys.path

import config  # noqa: E402

STEPS = {  # step: (module, output that marks it done)
    "extract": ("extract", "work/sites.csv"),
    "pairs": ("pairs", "pairs.csv"),
    "labels": ("labels", "work/labels_qc.csv"),
    "train": ("train", "model/unet.pt"),
    "analyze": ("analyze", "metrics.csv"),
    "summarize": ("summarize", "stats.json"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("steps", nargs="*", help=f"steps to run: {', '.join(STEPS)} (default: the missing ones, in order)")
    ap.add_argument("--force", action="store_true", help="run the steps even if their outputs exist")
    ap.add_argument("-o", "--out", type=Path, help=f"output folder (default: {config.OUT})")
    args = ap.parse_args()
    if unknown := [s for s in args.steps if s not in STEPS]:
        ap.error(f"unknown step(s): {', '.join(unknown)}")
    if args.out:
        config.OUT = args.out.resolve()
    print(f"output: {config.OUT}")
    chosen = args.steps or list(STEPS)
    for step in STEPS:
        if step not in chosen:
            continue
        module, done = STEPS[step]
        if not args.steps and not args.force and (config.OUT / done).exists():
            print(f"{step}: done already ({done})")
            continue
        t0 = time.time()
        __import__(module).main()
        print(f"{step}: {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
