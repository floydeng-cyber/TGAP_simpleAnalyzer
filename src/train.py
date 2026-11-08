"""Step 4: trains the U-Net on the classical masks that passed QC, pre- and post-etch images together.

Chip config.HOLDOUT_CHIP is held out to measure agreement with clean labels on unseen sites. Pre-etch
images are few (about 160 vs 1600), so each batch is half pre-etch, half post-etch. Training runs on the
CPU (about 1 s per step): this torch build's MPS backend fails on convolution backward on macOS 14.

Writes OUT/model/unet.pt and OUT/model/train_log.csv.
"""

import csv
import time

import cv2
import numpy as np
import tifffile
import torch
import torch.nn.functional as F

import config
from extract import read_sites
from labels import read_qc
from unet import MODEL_PATH, UNet, normalize

CROP = 192
BATCH = 8
BASE = 16
SEED = 0


def load_set() -> dict:
    chip_of = {(s["kind"], s["image"]): s["chip"] for s in read_sites()}
    data = {"train": {"pre": [], "post": []}, "val": {"pre": [], "post": []}}
    for r in read_qc():
        if r["ok"] != "True":
            continue
        kind, name = r["kind"], r["image"]
        img = tifffile.imread(config.OUT / "work/images" / kind / f"{name}.tif")
        mask = cv2.imread(str(config.OUT / "work/labels" / kind / f"{name}.png"), cv2.IMREAD_GRAYSCALE) > 0
        split = "val" if chip_of[(kind, name)] == config.HOLDOUT_CHIP else "train"
        data[split][kind].append((img, mask))
    return data


def augment(img, mask, rng):
    """Random crop, rotation by 90 deg, flip, contrast/brightness and extra noise."""
    y, x = rng.integers(0, 512 - CROP + 1, size=2)
    img, mask = img[y:y + CROP, x:x + CROP], mask[y:y + CROP, x:x + CROP]
    k = rng.integers(4)
    img, mask = np.rot90(img, k), np.rot90(mask, k)
    if rng.random() < 0.5:
        img, mask = img[:, ::-1], mask[:, ::-1]
    img = img.astype(np.float32) * rng.uniform(0.8, 1.2) + rng.uniform(-20, 20)
    img = img + rng.normal(0, rng.uniform(0, 8), img.shape)
    return np.ascontiguousarray(img, dtype=np.float32), np.ascontiguousarray(mask)


def loss_fn(logits, y):
    """Binary cross-entropy + Dice."""
    p = torch.sigmoid(logits)
    inter = (p * y).sum(dim=(1, 2, 3))
    dice = 1 - (2 * inter + 1) / (p.sum(dim=(1, 2, 3)) + y.sum(dim=(1, 2, 3)) + 1)
    return F.binary_cross_entropy_with_logits(logits, y) + dice.mean()


@torch.no_grad()
def evaluate(model, items) -> float:
    """Mean IoU with the labels on full 512x512 images."""
    model.eval()
    ious = []
    for img, m in items:
        p = (model(normalize(torch.from_numpy(img.astype(np.float32))[None, None]))[0, 0] > 0).numpy()
        ious.append((p & m).sum() / max((p | m).sum(), 1))
    model.train()
    return float(np.mean(ious)) if ious else float("nan")


def main(steps: int = 2400) -> None:
    torch.manual_seed(SEED)
    torch.set_num_threads(8)
    rng = np.random.default_rng(SEED)
    data = load_set()
    print("train: " + ", ".join(f"{s} {k} {len(v)}" for s in data for k, v in data[s].items()), flush=True)
    model = UNet(base=BASE)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=1e-3, total_steps=steps, pct_start=0.1)
    log = csv.writer(open(config.out("model/train_log.csv"), "w", newline=""))
    log.writerow(["step", "loss", "val_iou_pre", "val_iou_post", "seconds"])
    t0, run = time.time(), 0.0
    for step in range(1, steps + 1):
        pre = [data["train"]["pre"][i] for i in rng.integers(len(data["train"]["pre"]), size=BATCH // 2)]
        post = [data["train"]["post"][i] for i in rng.integers(len(data["train"]["post"]), size=BATCH - BATCH // 2)]
        imgs, masks = zip(*(augment(img, m, rng) for img, m in pre + post))
        x = normalize(torch.from_numpy(np.stack(imgs))[:, None])
        y = torch.from_numpy(np.stack(masks).astype(np.float32))[:, None]
        loss = loss_fn(model(x), y)
        opt.zero_grad()
        loss.backward()
        opt.step()
        sched.step()
        run = 0.98 * run + 0.02 * loss.item() if step > 1 else loss.item()
        if step % 300 == 0 or step == steps:
            vp, vq = evaluate(model, data["val"]["pre"]), evaluate(model, data["val"]["post"][:40])
            log.writerow([step, f"{run:.4f}", f"{vp:.4f}", f"{vq:.4f}", f"{time.time() - t0:.0f}"])
            print(f"train: step {step}: loss {run:.4f}, held-out IoU pre {vp:.4f} post {vq:.4f} ({time.time() - t0:.0f} s)", flush=True)
            torch.save({"state_dict": model.state_dict(), "base": BASE, "holdout_chip": config.HOLDOUT_CHIP, "step": step},
                       config.out(MODEL_PATH))
