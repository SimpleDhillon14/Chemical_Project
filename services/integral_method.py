"""
integral_method.py
-------------------
Implements the INTEGRAL METHOD of analysis exactly as taught for a single
reactant nth-order irreversible reaction in a constant-volume batch reactor:

    -rA = -dCA/dt = k * CA^n

Integrated forms (the "guess an order, integrate, plot, check linearity"
procedure):

    n = 1 :   -ln(CA/CA0) = k*t            -> plot  -ln(CA/CA0)  vs  t
              slope = k,  intercept = 0

    n != 1:   CA^(1-n) - CA0^(1-n) = (n-1)*k*t
              -> plot  CA^(1-n)  vs  t
              slope = (n-1)*k ,  intercept = CA0^(1-n)

    n = 0 :   CA0 - CA = k*t                -> plot CA vs t (special case of above)
              slope = -k, intercept = CA0

The method searches over reaction order n (0 to 3, step 0.05) linearizes the
data for each n, fits a straight line, and reports R^2 for every order so the
best (including fractional) order can be identified -- this mirrors the
textbook instruction to "try zero, first, second... order plots and see
which gives a straight line."
"""
import numpy as np
from .validation import validate_time_concentration, ValidationError
from .regression import linear_fit


def _linearize_for_order(t, C, n):
    """Return (x, y) arrays for the integral-method straight-line test of order n."""
    t = np.asarray(t, dtype=float)
    C = np.asarray(C, dtype=float)
    if abs(n - 1.0) < 1e-9:
        y = -np.log(C / C[0])
        x = t
        y_label = "-ln(C_A / C_A0)"
    else:
        y = np.power(C, 1.0 - n)
        x = t
        y_label = f"C_A^(1-{n:g})"
    return x, y, y_label


def analyze(times, concentrations, candidate_orders=None):
    t, C = validate_time_concentration(times, concentrations, min_points=3)
    CA0 = C[0]

    if candidate_orders is None:
        candidate_orders = [0.0, 1.0, 2.0, 3.0]

    # ---- 1. Standard integer-order comparison table -----------------------
    comparison = []
    for n in candidate_orders:
        try:
            x, y, y_label = _linearize_for_order(t, C, n)
            fit = linear_fit(x, y)
            if abs(n - 1.0) < 1e-9:
                k = fit["slope"]
            else:
                k = fit["slope"] / (n - 1.0)
            physically_valid = k > 0
            comparison.append({
                "order": n,
                "y_label": y_label,
                "slope": fit["slope"],
                "intercept": fit["intercept"],
                "r2": fit["r2"],
                "k": k,
                "physically_valid": physically_valid,
            })
        except Exception:
            continue

    # ---- 2. Fine grid search for the best-fit (possibly fractional) order -
    grid = np.arange(0.0, 3.01, 0.05)
    best = None
    grid_results = []
    for n in grid:
        try:
            x, y, _ = _linearize_for_order(t, C, float(n))
            fit = linear_fit(x, y)
            if abs(n - 1.0) < 1e-9:
                k = fit["slope"]
            else:
                k = fit["slope"] / (n - 1.0)
            if k <= 0:
                continue
            grid_results.append({"order": float(n), "r2": fit["r2"], "k": k,
                                  "slope": fit["slope"], "intercept": fit["intercept"]})
            if best is None or fit["r2"] > best["r2"]:
                best = {"order": float(n), "r2": fit["r2"], "k": k,
                        "slope": fit["slope"], "intercept": fit["intercept"]}
        except Exception:
            continue

    if best is None:
        raise ValidationError(
            "No physically valid (k > 0) straight-line fit could be found for any "
            "tested reaction order. Please check the experimental data."
        )

    best_n = round(best["order"], 2)
    rounded_hint = None
    for candidate in [0, 0.5, 1, 1.5, 2, 2.5, 3]:
        if abs(best_n - candidate) <= 0.1:
            rounded_hint = candidate
            break

    # ---- 3. Best-order transformed data table (for display) --------------
    x_best, y_best, y_label_best = _linearize_for_order(t, C, best_n)
    best_fit = linear_fit(x_best, y_best)
    k_best = best_fit["slope"] if abs(best_n - 1.0) < 1e-9 else best_fit["slope"] / (best_n - 1.0)

    transformed_table = []
    for i in range(len(t)):
        transformed_table.append({
            "time": t[i],
            "concentration": C[i],
            "x": float(x_best[i]),
            "y": float(y_best[i]),
        })

    if abs(best_n - 1.0) < 1e-9:
        rate_equation = f"-dC_A/dt = k*C_A,  k = {k_best:.6g} (1/time)"
        rate_law_str = f"-r_A = {k_best:.6g} * C_A"
    else:
        rate_equation = f"-dC_A/dt = k*C_A^{best_n:g},  k = {k_best:.6g}"
        rate_law_str = f"-r_A = {k_best:.6g} * C_A^{best_n:g}"

    steps = [
        {"title": "Step 1: Original kinetic equation",
         "text": "Assume a single-reactant nth-order irreversible reaction: -r_A = -dC_A/dt = k*C_A^n"},
        {"title": "Step 2: Integrated equation",
         "text": ("For n = 1: -ln(C_A/C_A0) = k t.  "
                   "For n != 1: C_A^(1-n) - C_A0^(1-n) = (n-1) k t.")},
        {"title": "Step 3: Data transformation",
         "text": f"Using the best-fit order n = {best_n:g}, the data was transformed to y = {y_label_best} and plotted against x = t."},
        {"title": "Step 4: Graph coordinates",
         "text": f"X-axis = t (time), Y-axis = {y_label_best}."},
        {"title": "Step 5: Linear regression",
         "text": f"A least-squares straight line was fit: slope = {best_fit['slope']:.6g}, intercept = {best_fit['intercept']:.6g}, R^2 = {best_fit['r2']:.5f}."},
        {"title": "Step 6: Slope and intercept interpretation",
         "text": ("For n = 1: slope = k. "
                   "For n != 1: slope = (n-1) k, and intercept = C_A0^(1-n) "
                   f"(k = {k_best:.6g}).")},
        {"title": "Step 7: Reaction order",
         "text": f"By comparing R^2 across candidate orders (grid search 0 to 3), the best-fitting order is n = {best_n:g}."
                  + (f" This is approximately a {rounded_hint:g}-order reaction." if rounded_hint is not None else "")},
        {"title": "Step 8: Rate constant, k",
         "text": f"k = {k_best:.6g} (units depend on order and concentration/time units used)."},
        {"title": "Step 9: Final rate equation",
         "text": rate_law_str},
    ]

    return {
        "method": "Integral Method",
        "input_data": [{"time": t[i], "concentration": C[i]} for i in range(len(t))],
        "CA0": CA0,
        "original_equation": "-r_A = -dC_A/dt = k * C_A^n",
        "comparison_table": comparison,
        "grid_search": grid_results,
        "best_order": best_n,
        "rounded_order_hint": rounded_hint,
        "k": k_best,
        "regression": best_fit,
        "y_label": y_label_best,
        "x_label": "t (time)",
        "transformed_table": transformed_table,
        "rate_equation": rate_equation,
        "rate_law_str": rate_law_str,
        "steps": steps,
        "graph": {
            "x_axis": "t (time)",
            "y_axis": y_label_best,
            "points_x": list(x_best.tolist() if hasattr(x_best, "tolist") else x_best),
            "points_y": list(y_best.tolist() if hasattr(y_best, "tolist") else y_best),
            "line_x": best_fit["line_x"],
            "line_y": best_fit["line_y"],
        },
        "prediction": {
            "order": best_n,
            "k": k_best,
            "CA0": CA0,
        },
    }


def predict_concentration(order, k, CA0, t):
    """C_A(t) from the integrated nth-order law. Returns None if non-physical."""
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    if abs(order - 1.0) < 1e-9:
        C = CA0 * np.exp(-k * t)
    else:
        base = CA0 ** (1.0 - order) + (order - 1.0) * k * t
        if base <= 0:
            return None  # reaction would have gone to completion before time t (n>1) or invalid
        C = base ** (1.0 / (1.0 - order))
    if C <= 0 or np.isnan(C) or np.isinf(C):
        return None
    return float(C)


def predict_time(order, k, CA0, C):
    """t at which C_A = C is reached, from the integrated nth-order law."""
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    if C > CA0:
        raise ValidationError("Target concentration cannot exceed the initial concentration C_A0 "
                               "for a reactant being consumed.")
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
