"""Per-angle ring sampling + 2-point linear calibration.

For each printed ring we sample mean RGB at every angular bin within
the ring's radial band (using the inset 80 % of the band to avoid
anti-aliased edges and the angular tick marks on the outside of the
black ring). Per-angle slope / intercept then maps the printed black
ring to 0 and the printed gray ring to its known mid-tone linear
value. Output is unclipped: agar can read darker than the printed
black ink and calibrates to slightly negative values; clipping there
would discard real signal.
"""
from __future__ import annotations

import numpy as np

from .config import CALIB_MIN_GAP, RING_TRUTH
from .geometry import Geometry, _angle_grid, _radius_grid


def _smooth_circular(arr: np.ndarray, window_deg: int) -> np.ndarray:
    n = arr.shape[0]
    w = max(3, int(round(window_deg * n / 360.0)))
    if w % 2 == 0:
        w += 1
    pad = w // 2
    padded = np.concatenate([arr[-pad:], arr, arr[:pad]], axis=0)
    kernel = np.ones(w, dtype=np.float32) / w
    out = np.empty_like(arr)
    for ch in range(arr.shape[1]):
        out[:, ch] = np.convolve(padded[:, ch], kernel, mode='valid')
    return out


def sample_ring_per_angle(img: np.ndarray, geom: Geometry,
                          ring_name: str, n_angles: int = 360) -> np.ndarray:
    """Mean RGB at each angle bin within the ring's radial band.

    Uses the inset 80 % of the radial band (10 % shaved off each
    edge) to avoid anti-aliased edges and the angular tick marks
    printed on the outside of the black ring. Per-bin mean is
    computed via `np.bincount`; the subsequent circular smoothing
    damps any per-bin outlier influence.
    """
    h, w = img.shape[:2]
    theta = _angle_grid(geom.card_cx, geom.card_cy, h, w)
    rgrid = _radius_grid(geom.card_cx, geom.card_cy, h, w)
    r_in = getattr(geom, f'r_{ring_name}_in')
    r_out = getattr(geom, f'r_{ring_name}_out')
    band = (r_out - r_in) * 0.10
    in_ring = (rgrid >= (r_in + band)) & (rgrid <= (r_out - band))
    if not in_ring.any():
        raise RuntimeError(f'ring {ring_name!r} band has no pixels (geometry off)')

    pixel_theta = theta[in_ring]
    pixel_rgb = img[in_ring]
    step = 360.0 / n_angles
    bin_idx = np.minimum((pixel_theta / step).astype(np.int32), n_angles - 1)

    counts = np.bincount(bin_idx, minlength=n_angles).astype(np.float32)
    counts_safe = np.maximum(counts, 1.0)
    out = np.empty((n_angles, 3), dtype=np.float32)
    for ch in range(3):
        sums = np.bincount(bin_idx,
                           weights=pixel_rgb[:, ch].astype(np.float32),
                           minlength=n_angles)
        out[:, ch] = sums / counts_safe
    out[counts == 0] = 0.0
    return _smooth_circular(out, window_deg=15)


def build_calibration(black_per_angle: np.ndarray,
                      gray_per_angle: np.ndarray) -> np.ndarray:
    """Per-channel, per-angle 2-point linear fit. Returns shape
    (n, 3, 2): [..., 0]=slope, [..., 1]=intercept such that
        slope * raw_black + intercept = 0
        slope * raw_gray  + intercept = gray_truth
    """
    gray_truth = RING_TRUTH['gray'][None, :]
    gap = gray_per_angle - black_per_angle
    if np.any(gap < CALIB_MIN_GAP):
        raise RuntimeError(
            f'calibration ill-conditioned: gray-black gap < {CALIB_MIN_GAP} '
            f'(min gap = {float(gap.min()):.4f})'
        )
    slope = gray_truth / gap
    intercept = -slope * black_per_angle
    return np.stack([slope, intercept], axis=-1)


def apply_calibration(img: np.ndarray, geom: Geometry,
                      calibration: np.ndarray) -> np.ndarray:
    """Per-pixel linear correction looked up at the pixel's polar angle.

    Output is NOT clipped: the printed "black" ink has some non-zero
    reflectance, so agar (often darker than printed black) calibrates
    to slightly negative values. Clipping there would discard real
    signal. Saving to PNG clips for display only.
    """
    h, w = img.shape[:2]
    theta = _angle_grid(geom.card_cx, geom.card_cy, h, w)
    n_angles = calibration.shape[0]
    step = 360.0 / n_angles
    bin_idx = np.minimum((theta / step).astype(np.int32), n_angles - 1)
    out = np.empty_like(img)
    for ch in range(3):
        s = calibration[bin_idx, ch, 0]
        b = calibration[bin_idx, ch, 1]
        out[..., ch] = s * img[..., ch] + b
    return out
