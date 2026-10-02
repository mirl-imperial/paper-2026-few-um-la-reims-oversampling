"""Step-edge spatial-resolution fit: an ideal step edge blurred by a
Gaussian point-spread function, fit to an intensity profile sampled across
the edge with an error function (Zhang & Bergholm, Int. J. Comput. Vis.
24, 219-250, 1997), reported as a 16-84% (+/-1 sigma) rise-distance
resolution metric per the convention of Kompauer et al. (Nat. Methods, 2017).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf, erfinv

# erf edge-blur model (Zhang & Bergholm, Int. J. Comput. Vis. 24, 219-250,
# 1997), per the 16-84% (+/-1 sigma) convention of Kompauer et al.
# (Nat. Methods, 2017).
LOW_FRAC, HIGH_FRAC = 0.16, 0.84


def edge_model(x, low, high, x0, sigma):
    """Error-function edge model: an ideal step blurred by a Gaussian
    point-spread function of width `sigma`."""
    return low + (high - low) * 0.5 * (1 + erf((x - x0) / (np.sqrt(2) * sigma)))


def fit_edge(distance_um, intensity):
    lo_guess, hi_guess = intensity.min(), intensity.max()
    x0_guess = distance_um[np.argmin(np.abs(intensity - 0.5 * (lo_guess + hi_guess)))]
    sigma_guess = (distance_um.max() - distance_um.min()) / 8.0
    if sigma_guess <= 0:
        sigma_guess = 1.0
    p0 = [lo_guess, hi_guess, x0_guess, sigma_guess]
    bounds = (
        [-np.inf, -np.inf, distance_um.min(), 1e-3],
        [np.inf, np.inf, distance_um.max(), max(distance_um.max() - distance_um.min(), 1.0)],
    )
    popt, pcov = curve_fit(edge_model, distance_um, intensity, p0=p0, bounds=bounds,
                           maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr


def resolution_16_84(sigma):
    """16-84% (+/-1 sigma) rise distance for the erf edge model, given its sigma."""
    z = erfinv(HIGH_FRAC - LOW_FRAC)  # erf(z) = 0.68 for 16-84%
    return 2 * np.sqrt(2) * sigma * z
