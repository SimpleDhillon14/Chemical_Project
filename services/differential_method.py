"""
differential_method.py
-----------------------
Implements the DIFFERENTIAL METHOD OF ANALYSIS exactly as presented in
CRE_Diff.pdf (Chapter 3, "Differential Method of Analysis of Data"):

  1. Plot C_A vs t, draw a smooth curve through the data.
  2. Determine the slope dC_A/dt at each concentration -> these slopes are
     the rates of reaction, -r_A, at each composition (Fig. 3.17).
  3. To test an nth-order rate form  -r_A = k*C_A^n, take logarithms:

         log10(-dC_A/dt) = log10(k) + n*log10(C_A)
                y                intercept      slope   x

     and plot log10(-r_A) vs log10(C_A) (Fig. 3.18). The SLOPE of the best
     straight line is the reaction order n, and the INTERCEPT is log10(k).

Because we do not have a hand-drawn smooth curve, dC_A/dt is estimated
numerically (central differences for interior points, forward/backward
differences at the endpoints), which is the standard numerical stand-in for
"draw a smooth curve and take its slope" -- this mirrors Example 3.2 in the
material, where the tabulated slope column is exactly -dC_A/dt at each
reported concentration.

The module also supports the Michaelis-Menten-type rate form taught in the
same chapter (Eq. 57):

        -r_A = k1*C_A / (1 + k2*C_A)

with its two linearizations (Eqs. 61 and 62, Fig. 3.19):

    Eq. 61:  1/(-r_A) = 1/(k1*C_A) + k2/k1     -> plot 1/(-r_A) vs 1/C_A
             slope = 1/k1,  intercept = k2/k1

    Eq. 62:  (-r_A) = k1/k2 - (1/k2)*[(-r_A)/C_A]   -> plot -r_A vs (-r_A)/C_A
             slope = -1/k2,  intercept = k1/k2
"""
import numpy as np
from .validation import validate_time_concentration, ValidationError
from .regression import linear_fit


def numerical_rate(t, C):
    """
    Estimate -dC_A/dt at every data point.
    Central difference for interior points, forward/backward at the ends,
    following standard numerical differentiation practice for unevenly or
    evenly spaced batch-reactor data (the numerical equivalent of drawing a
    smooth curve and reading its slope, as instructed in CRE_Diff.pdf).
    """
    t = np.asarray(t, dtype=float)
    C = np.asarray(C, dtype=float)
    n = len(t)
    dCdt = np.zeros(n)

    # forward difference for first point
    dCdt[0] = (C[1] - C[0]) / (t[1] - t[0])
    # backward difference for last point
    dCdt[-1] = (C[-1] - C[-2]) / (t[-1] - t[-2])
    # central difference for interior points (non-uniform spacing safe)
    for i in range(1, n - 1):
        h1 = t[i] - t[i - 1]
        h2 = t[i + 1] - t[i]
        dCdt[i] = (C[i + 1] * h1 ** 2 - C[i - 1] * h2 ** 2 + C[i] * (h2 ** 2 - h1 ** 2)) / (h1 * h2 * (h1 + h2))

    neg_rate = -dCdt  # -r_A = -dC_A/dt
    return neg_rate


def analyze(times, concentrations):
    t, C = validate_time_concentration(times, concentrations, min_points=3)
    t = np.array(t)
    C = np.array(C)

    neg_rate = numerical_rate(t, C)

    rate_table = [{"time": float(t[i]), "concentration": float(C[i]), "neg_dCdt": float(neg_rate[i])}
                  for i in range(len(t))]

    # keep only physically valid points (rate must be > 0 for log)
    valid_idx = [i for i in range(len(t)) if neg_rate[i] > 0]
    if len(valid_idx) < 3:
        raise ValidationError(
            "Fewer than 3 points have a positive computed rate (-dC_A/dt > 0). "
            "This usually means the concentration data is not monotonically "
            "decreasing, which is required to compute a rate of disappearance."
        )

    logC = np.log10(C[valid_idx])
    logR = np.log10(neg_rate[valid_idx])

    fit = linear_fit(logC, logR)
    n_order = fit["slope"]
    k = 10 ** fit["intercept"]

    transformed_table = [
        {"time": float(t[i]), "concentration": float(C[i]), "neg_dCdt": float(neg_rate[i]),
         "log_C": float(np.log10(C[i])), "log_rate": float(np.log10(neg_rate[i])) if neg_rate[i] > 0 else None}
        for i in range(len(t))
    ]

    rounded_hint = None
    for candidate in [0, 0.5, 1, 1.5, 2, 2.5, 3]:
        if abs(round(n_order, 2) - candidate) <= 0.1:
            rounded_hint = candidate
            break

    steps = [
        {"title": "Step 1: Plot C_A vs t and estimate slopes",
         "text": "The experimental C_A vs t data was used to numerically estimate the slope dC_A/dt at each point (central difference for interior points, forward/backward difference at the endpoints), following Fig. 3.17 of the course material."},
        {"title": "Step 2: Rate of reaction",
         "text": "The rate of reaction at each composition is -r_A = -dC_A/dt (see rate table)."},
        {"title": "Step 3: Assume an nth-order rate form",
         "text": "-r_A = k * C_A^n"},
        {"title": "Step 4: Linearize by taking logarithms",
         "text": "log10(-dC_A/dt) = log10(k) + n * log10(C_A)"},
        {"title": "Step 5: Graph coordinates",
         "text": "X-axis = log10(C_A), Y-axis = log10(-dC_A/dt), exactly as in Fig. 3.18."},
        {"title": "Step 6: Linear regression",
         "text": f"slope = {fit['slope']:.6g}, intercept = {fit['intercept']:.6g}, R^2 = {fit['r2']:.5f}."},
        {"title": "Step 7: Interpretation",
         "text": "Slope = n (reaction order). Intercept = log10(k)."},
        {"title": "Step 8: Reaction order",
         "text": f"n = {n_order:.4g}" + (f"  (approximately {rounded_hint:g} order)" if rounded_hint is not None else "")},
        {"title": "Step 9: Rate constant",
         "text": f"k = 10^({fit['intercept']:.6g}) = {k:.6g}"},
        {"title": "Step 10: Final rate law",
         "text": f"-r_A = {k:.6g} * C_A^{n_order:.4g}"},
    ]

    # --- optional Michaelis-Menten style check (Eqs. 61/62) ---------------
    mm_result = None
    try:
        inv_C = 1.0 / C[valid_idx]
        inv_rate = 1.0 / neg_rate[valid_idx]
        fit61 = linear_fit(inv_C, inv_rate)
        k1_from61 = 1.0 / fit61["slope"] if fit61["slope"] != 0 else None
        k2_from61 = fit61["intercept"] * k1_from61 if k1_from61 else None

        rate_over_C = neg_rate[valid_idx] / C[valid_idx]
        fit62 = linear_fit(rate_over_C, neg_rate[valid_idx])
        k2_from62 = -1.0 / fit62["slope"] if fit62["slope"] != 0 else None
        k1_from62 = fit62["intercept"] * k2_from62 if k2_from62 else None

        if k1_from61 and k1_from61 > 0 and k2_from61 and k2_from61 > 0:
            mm_result = {
                "eq61": {"slope": fit61["slope"], "intercept": fit61["intercept"], "r2": fit61["r2"],
                         "k1": k1_from61, "k2": k2_from61,
                         "x_label": "1/C_A", "y_label": "1/(-r_A)",
                         "points_x": inv_C.tolist(), "points_y": inv_rate.tolist(),
                         "line_x": fit61["line_x"], "line_y": fit61["line_y"]},
                "eq62": {"slope": fit62["slope"], "intercept": fit62["intercept"], "r2": fit62["r2"],
                         "k1": k1_from62, "k2": k2_from62,
                         "x_label": "(-r_A)/C_A", "y_label": "-r_A",
                         "points_x": rate_over_C.tolist(), "points_y": neg_rate[valid_idx].tolist(),
                         "line_x": fit62["line_x"], "line_y": fit62["line_y"]},
            }
    except Exception:
        mm_result = None

    return {
        "method": "Differential Method",
        "input_data": [{"time": float(t[i]), "concentration": float(C[i])} for i in range(len(t))],
        "original_equation": "-r_A = -dC_A/dt = k * C_A^n",
        "rate_table": rate_table,
        "transformed_table": transformed_table,
        "regression": fit,
        "n_order": n_order,
        "rounded_order_hint": rounded_hint,
        "k": k,
        "rate_law_str": f"-r_A = {k:.6g} * C_A^{n_order:.4g}",
        "x_label": "log10(C_A)",
        "y_label": "log10(-dC_A/dt)",
        "graph": {
            "x_axis": "log10(C_A)",
            "y_axis": "log10(-dC_A/dt)",
            "points_x": logC.tolist(),
            "points_y": logR.tolist(),
            "line_x": fit["line_x"],
            "line_y": fit["line_y"],
        },
        "steps": steps,
        "michaelis_menten_check": mm_result,
        "prediction": {"order": n_order, "k": k, "CA0": float(C[0])},
    }


def predict_concentration(order, k, CA0, t):
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    if abs(order - 1.0) < 1e-9:
        C = CA0 * np.exp(-k * t)
    else:
        base = CA0 ** (1.0 - order) + (order - 1.0) * k * t
        if base <= 0:
            return None
        C = base ** (1.0 / (1.0 - order))
    if C <= 0 or np.isnan(C) or np.isinf(C):
        return None
    return float(C)


def predict_time(order, k, CA0, C):
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    if C > CA0:
        raise ValidationError("Target concentration cannot exceed the initial concentration C_A0.")
    if abs(order - 1.0) < 1e-9:
        t = -np.log(C / CA0) / k
    else:
        t = (C ** (1.0 - order) - CA0 ** (1.0 - order)) / ((order - 1.0) * k)
    if t < 0 or np.isnan(t) or np.isinf(t):
        return None
    return float(t)


def predict_rate(order, k, C):
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    return float(k * (C ** order))
