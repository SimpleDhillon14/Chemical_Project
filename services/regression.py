"""
regression.py
-------------
Thin wrapper around scipy.stats.linregress used by every CRE method so that
slope / intercept / R2 are computed the same way everywhere, and so the
resulting straight line (for plotting) is generated consistently.
"""
import numpy as np
from scipy import stats


def linear_fit(x, y):
    """
    Fit y = slope*x + intercept by least squares.

    Returns a dict with slope, intercept, r2, r_value, p_value, std_err,
    and ready-to-plot line coordinates (line_x, line_y) spanning the data.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if len(x) < 2:
        raise ValueError("At least two points are required for linear regression.")

    result = stats.linregress(x, y)
    slope = float(result.slope)
    intercept = float(result.intercept)
    r2 = float(result.rvalue ** 2)

    x_min, x_max = float(np.min(x)), float(np.max(x))
    if x_min == x_max:
        line_x = [x_min, x_min]
    else:
        line_x = list(np.linspace(x_min, x_max, 50))
    line_y = [slope * xv + intercept for xv in line_x]

    return {
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "r_value": float(result.rvalue),
        "p_value": float(result.pvalue),
        "std_err": float(result.stderr),
        "line_x": line_x,
        "line_y": line_y,
        "equation": f"y = {slope:.6g}x + {intercept:.6g}",
    }


def round_sig(value, sig=4):
    """Round to a given number of significant figures for clean display."""
    if value == 0 or value is None:
        return 0.0
    import math
    d = sig - int(math.floor(math.log10(abs(value)))) - 1
    return round(value, d)
