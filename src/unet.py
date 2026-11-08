"""The U-Net (Ronneberger et al. 2015): one-channel CD-SEM image -> footprint mask. Model, loading, prediction."""

import cv2
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy import ndimage as ndi

import config


class DoubleConv(nn.Sequential):
    def __init__(self, cin, cout):
        super().__init__(
            nn.Conv2d(cin, cout, 3, padding=1, padding_mode="reflect", bias=False),
            nn.BatchNorm2d(cout),
            nn.ReLU(inplace=True),
            nn.Conv2d(cout, cout, 3, padding=1, padding_mode="reflect", bias=False),
            nn.BatchNorm2d(cout),
            nn.ReLU(inplace=True),
        )


class UNet(nn.Module):
    def __init__(self, base: int = 32, depth: int = 4):
        super().__init__()
        chs = [base * 2**i for i in range(depth + 1)]
        self.inc = DoubleConv(1, chs[0])
        self.downs = nn.ModuleList(DoubleConv(chs[i], chs[i + 1]) for i in range(depth))
        # bilinear upsampling + 1x1 conv; transposed-conv backward crashes on the MPS backend of macOS 14
        self.ups = nn.ModuleList(
            nn.Sequential(nn.Upsample(scale_factor=2, mode="bilinear", align_corners=False), nn.Conv2d(chs[i + 1], chs[i], 1))
            for i in reversed(range(depth))
        )
        self.dec = nn.ModuleList(DoubleConv(chs[i] * 2, chs[i]) for i in reversed(range(depth)))
        self.out = nn.Conv2d(chs[0], 1, 1)

    def forward(self, x):
        skips = [self.inc(x)]
        for down in self.downs:
            skips.append(down(F.max_pool2d(skips[-1], 2)))
        x = skips.pop()
        for up, dec in zip(self.ups, self.dec):
            x = up(x)
            x = dec(torch.cat([skips.pop(), x], dim=1))
        return self.out(x)


def normalize(img: torch.Tensor) -> torch.Tensor:
    """Per-image standardization (images of both kinds differ in brightness and contrast)."""
    m = img.mean(dim=(-2, -1), keepdim=True)
    s = img.std(dim=(-2, -1), keepdim=True).clamp_min(1e-3)
    return (img - m) / s


MODEL_PATH = "model/unet.pt"  # under config.OUT


def load() -> tuple[UNet, dict]:
    ck = torch.load(config.OUT / MODEL_PATH, map_location="cpu")
    model = UNet(base=ck["base"])
    model.load_state_dict(ck["state_dict"])
    model.eval()
    return model, ck


@torch.no_grad()
def predict(model: UNet, img: np.ndarray) -> np.ndarray:
    """Footprint mask: mean probability over the 4 flips > 0.5, specks under 300 px dropped (thin
    W=150 nm traces can be under 2000 px), holes under 500 px filled."""
    x = normalize(torch.from_numpy(img.astype(np.float32))[None, None])
    probs = []
    for dims in ([], [2], [3], [2, 3]):
        xi = torch.flip(x, dims) if dims else x
        p = torch.sigmoid(model(xi))
        probs.append(torch.flip(p, dims) if dims else p)
    mask = torch.stack(probs).mean(0)[0, 0].numpy() > 0.5
    n, lab, stats, _ = cv2.connectedComponentsWithStats(mask.astype(np.uint8), connectivity=8)
    keep = np.zeros(n, bool)
    keep[1:] = stats[1:, cv2.CC_STAT_AREA] >= 300
    mask = keep[lab]
    holes = ndi.binary_fill_holes(mask) & ~mask
    hn, hlab, hstats, _ = cv2.connectedComponentsWithStats(holes.astype(np.uint8), connectivity=4)
    small = np.zeros(hn, bool)
    small[1:] = hstats[1:, cv2.CC_STAT_AREA] < 500
    return mask | small[hlab]
