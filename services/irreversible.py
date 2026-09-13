"""
irreversible.py
----------------
Implements the IRREVERSIBLE SERIES REACTION analysis exactly as derived in
the CRE_Reactions.pdf handwritten notes:

    A --(k1)--> R --(k2)--> S      (both steps irreversible, first order)

    -r_A = -dC_A/dt = k1*C_A            ->  C_A = C_A0 * e^(-k1 t)
                                             ->  -ln(C_A/C_A0) = k1 t   (LINEAR in t)

    r_R = k1*C_A - k2*C_R = dC_R/dt

    C_R = C_A0 * k1 * [ e^(-k1 t)/(k2-k1) + e^(-k2 t)/(k1-k2) ]

    t_max = ln(k2/k1) / (k2 - k1)
    C_R,max / C_A0 = (k1/k2)^(k2/(k2-k1))

Step 1 (mandatory) uses the C_A vs t data to get k1 from the straight line
    X-axis = t,  Y-axis = -ln(C_A/C_A0),  slope = k1.

Step 2 (optional, only if the user also supplies C_R vs t data) fits the
non-linear C_R(t) expression above by least squares to estimate k2, exactly
using the equation reproduced from the course notes (no alternative model is
substituted).
"""
import numpy as np
from scipy.optimize import curve_fit
from .validation import validate_time_concentration, ValidationError
from .regression import linear_fit


def _cr_model(t, CA0, k1, k2):
    if abs(k2 - k1) < 1e-9:
        # degenerate case k1 == k2
        return CA0 * k1 * t * np.exp(-k1 * t)
    return CA0 * k1 * (np.exp(-k1 * t) / (k2 - k1) + np.exp(-k2 * t) / (k1 - k2))


def analyze(times, concentrations_A, times_R=None, concentrations_R=None):
    t, CA = validate_time_concentration(times, concentrations_A, min_points=3)
    t = np.array(t)
    CA = np.array(CA)
    CA0 = CA[0]

    y = -np.log(CA / CA0)
    x = t
    fit = linear_fit(x, y)
    k1 = fit["slope"]
    if k1 <= 0:
        raise ValidationError("The regression produced a non-positive k1. C_A must be "
                               "monotonically decreasing for this first-order-in-A series reaction.")

    transformed_table = [{"time": float(t[i]), "CA": float(CA[i]), "y": float(y[i])} for i in range(len(t))]

    steps = [
        {"title": "Step 1: Original kinetic scheme",
         "text": "A --k1--> R --k2--> S  (irreversible, first order in each step)"},
        {"title": "Step 2: Rate equation for A",
         "text": "-r_A = -dC_A/dt = k1*C_A  ->  C_A = C_A0 * e^(-k1 t)"},
        {"title": "Step 3: Linearization",
         "text": "-ln(C_A/C_A0) = k1 * t"},
        {"title": "Step 4: Graph coordinates",
         "text": "X-axis = t (time), Y-axis = -ln(C_A/C_A0)."},
        {"title": "Step 5: Linear regression",
         "text": f"slope = {fit['slope']:.6g}, intercept = {fit['intercept']:.6g}, R^2 = {fit['r2']:.5f}."},
        {"title": "Step 6: Rate constant k1",
         "text": f"k1 = slope = {k1:.6g}"},
    ]

    result = {
        "method": "Irreversible Series Reaction (A -> R -> S)",
        "input_data": [{"time": float(t[i]), "concentration": float(CA[i])} for i in range(len(t))],
        "original_equation": "A --k1--> R --k2--> S,  -r_A = k1*C_A",
        "CA0": float(CA0),
        "transformed_table": transformed_table,
        "regression": fit,
        "k1": k1,
        "x_label": "t (time)",
        "y_label": "-ln(C_A / C_A0)",
        "graph": {
            "x_axis": "t (time)",
            "y_axis": "-ln(C_A / C_A0)",
            "points_x": x.tolist(),
            "points_y": y.tolist(),
            "line_x": fit["line_x"],
            "line_y": fit["line_y"],
        },
        "steps": steps,
        "k2": None,
        "t_max": None,
        "CR_max_over_CA0": None,
        "prediction": {"k1": k1, "k2": None, "CA0": float(CA0)},
    }

    # ---- optional: fit k2 from C_R(t) data --------------------------------
    if times_R and concentrations_R and len(times_R) >= 3:
        try:
            tR, CR = validate_time_concentration(times_R, concentrations_R, min_points=3,
                                                   require_monotonic_time=True)
            tR = np.array(tR)
            CR = np.array(CR)
            popt, _ = curve_fit(lambda tt, k2: _cr_model(tt, CA0, k1, k2), tR, CR,
                                 p0=[max(k1 * 2, 1e-3)], maxfev=10000)
            k2 = float(popt[0])
            if k2 <= 0:
                raise ValidationError("Fitted k2 was non-positive; C_R data may not match the series model.")

            pred_CR = _cr_model(tR, CA0, k1, k2)
            ss_res = float(np.sum((CR - pred_CR) ** 2))
            ss_tot = float(np.sum((CR - np.mean(CR)) ** 2))
            r2_k2 = 1 - ss_res / ss_tot if ss_tot > 0 else None

            if abs(k2 - k1) > 1e-9:
                t_max = np.log(k2 / k1) / (k2 - k1)
                cr_max_ratio = (k1 / k2) ** (k2 / (k2 - k1))
            else:
                t_max = 1.0 / k1
                cr_max_ratio = 1.0 / np.e

            result["k2"] = k2
            result["t_max"] = float(t_max)
            result["CR_max_over_CA0"] = float(cr_max_ratio)
            result["cr_fit_r2"] = r2_k2
            result["cr_data"] = [{"time": float(tR[i]), "CR_experimental": float(CR[i]),
                                   "CR_model": float(pred_CR[i])} for i in range(len(tR))]
            result["prediction"]["k2"] = k2
            result["steps"].append({
                "title": "Step 7: Fitting k2 from C_R(t) data",
                "text": ("C_R(t) = C_A0*k1*[ e^(-k1 t)/(k2-k1) + e^(-k2 t)/(k1-k2) ] was fit by "
                         f"nonlinear least squares using the known k1 to obtain k2 = {k2:.6g} "
                         f"(fit R^2 = {r2_k2:.5f}).")
            })
            result["steps"].append({
                "title": "Step 8: Time and concentration of maximum R",
                "text": (f"t_max = ln(k2/k1)/(k2-k1) = {t_max:.6g}; "
                         f"C_R,max/C_A0 = (k1/k2)^(k2/(k2-k1)) = {cr_max_ratio:.6g}")
            })
        except ValidationError:
            raise
        except Exception:
            pass  # k2 estimation is optional; silently skip if fit fails

    return result


def predict_concentration_A(k1, CA0, t):
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    return float(CA0 * np.exp(-k1 * t))


def predict_time_A(k1, CA0, C):
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    if C > CA0:
        raise ValidationError("Target concentration cannot exceed C_A0.")
    return float(-np.log(C / CA0) / k1)


def predict_rate_A(k1, C):
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    return float(k1 * C)


def predict_concentration_R(k1, k2, CA0, t):
    if k2 is None:
        raise ValidationError("k2 has not been determined. Provide C_R vs t data to estimate k2.")
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    return float(_cr_model(np.array([t]), CA0, k1, k2)[0])
