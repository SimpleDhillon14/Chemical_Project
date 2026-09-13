"""
reversible.py
--------------
Implements the REVERSIBLE REACTION analysis exactly as derived in the
CRE_Reactions.pdf handwritten notes:

    A <=(k1, k2)=> R

At equilibrium the net rate is zero, giving the equilibrium constant

    K_C = k1 / k2 = C_Re / C_Ae

Let  M = C_R0 / C_A0  and  X_Ae = 1 - C_Ae / C_A0  (equilibrium conversion,
obtained from the user-supplied equilibrium concentration C_Ae). Integrating
the first-order-in-each-direction rate law gives the straight-line test:

    -ln(1 - X_A/X_Ae) = [(M+1)/(M+X_Ae)] * k1 * t = -ln[(C_A - C_Ae)/(C_A0 - C_Ae)]

    Plot:  X-axis = t
           Y-axis = -ln[(C_A - C_Ae)/(C_A0 - C_Ae)]

    slope = [(M+1)/(M+X_Ae)] * k1   ->   k1 = slope * (M + X_Ae) / (M + 1)
    k2 = k1 / K_C
"""
import numpy as np
from .validation import validate_time_concentration, validate_positive, validate_non_negative, ValidationError
from .regression import linear_fit


def analyze(times, concentrations, CR0, CAe):
    t, CA = validate_time_concentration(times, concentrations, min_points=3)
    t = np.array(t)
    CA = np.array(CA)
    CA0 = CA[0]
    validate_non_negative(CR0, "Initial concentration of R (C_R0)")
    validate_positive(CAe, "Equilibrium concentration of A (C_Ae)")

    if CAe >= CA0:
        raise ValidationError("Equilibrium concentration C_Ae must be less than the initial "
                               "concentration C_A0 (A is being consumed as it approaches equilibrium).")

    M = CR0 / CA0
    C0_total = CA0 + CR0
    CRe = C0_total - CAe
    XAe = 1.0 - CAe / CA0
    KC = CRe / CAe

    denom = CA0 - CAe
    y_arg = (CA - CAe) / denom
    if np.any(y_arg <= 0):
        raise ValidationError("Some concentration values are at or below the equilibrium "
                               "concentration C_Ae, which is not valid for this transformation. "
                               "Check that C_Ae represents the true equilibrium (long-time) value.")

    y = -np.log(y_arg)
    x = t

    fit = linear_fit(x, y)
    k1 = fit["slope"] * (M + XAe) / (M + 1.0)
    if k1 <= 0:
        raise ValidationError("The regression produced a non-positive forward rate constant k1. "
                               "Check the trend of the data and the supplied C_Ae.")
    k2 = k1 / KC

    transformed_table = [
        {"time": float(t[i]), "CA": float(CA[i]), "y": float(y[i])}
        for i in range(len(t))
    ]

    steps = [
        {"title": "Step 1: Original kinetic equation",
         "text": "A <=> R (forward k1, reverse k2). At equilibrium: K_C = k1/k2 = C_Re/C_Ae."},
        {"title": "Step 2: Equilibrium data",
         "text": f"C_A0 = {CA0:.6g}, C_R0 = {CR0:.6g}, M = C_R0/C_A0 = {M:.6g}, "
                  f"C_Ae = {CAe:.6g}, C_Re = {CRe:.6g}, X_Ae = {XAe:.6g}, K_C = {KC:.6g}."},
        {"title": "Step 3: Integrated linear form",
         "text": "-ln[(C_A - C_Ae)/(C_A0 - C_Ae)] = [(M+1)/(M+X_Ae)] * k1 * t"},
        {"title": "Step 4: Data transformation",
         "text": "y = -ln[(C_A - C_Ae)/(C_A0 - C_Ae)] computed for every time point."},
        {"title": "Step 5: Graph coordinates",
         "text": "X-axis = t (time), Y-axis = -ln[(C_A - C_Ae)/(C_A0 - C_Ae)]."},
        {"title": "Step 6: Linear regression",
         "text": f"slope = {fit['slope']:.6g}, intercept = {fit['intercept']:.6g}, R^2 = {fit['r2']:.5f}."},
        {"title": "Step 7: Interpretation",
         "text": "slope = [(M+1)/(M+X_Ae)] * k1  ->  k1 = slope * (M + X_Ae) / (M + 1). Then k2 = k1 / K_C."},
        {"title": "Step 8: Rate constants",
         "text": f"k1 (forward) = {k1:.6g},  k2 (reverse) = {k2:.6g},  K_C = {KC:.6g}."},
        {"title": "Step 9: Final rate law",
         "text": f"-r_A = k1*C_A - k2*C_R = {k1:.6g}*C_A - {k2:.6g}*C_R"},
    ]

    return {
        "method": "Reversible Reaction",
        "input_data": [{"time": float(t[i]), "concentration": float(CA[i])} for i in range(len(t))],
        "original_equation": "A <=> R,   K_C = k1/k2 = C_Re/C_Ae",
        "CA0": float(CA0), "CR0": float(CR0), "CAe": float(CAe), "CRe": float(CRe),
        "M": float(M), "XAe": float(XAe), "KC": float(KC),
        "transformed_table": transformed_table,
        "regression": fit,
        "k1": k1, "k2": k2,
        "rate_law_str": f"-r_A = {k1:.6g}*C_A - {k2:.6g}*C_R",
        "x_label": "t (time)",
        "y_label": "-ln[(C_A - C_Ae) / (C_A0 - C_Ae)]",
        "graph": {
            "x_axis": "t (time)",
            "y_axis": "-ln[(C_A - C_Ae) / (C_A0 - C_Ae)]",
            "points_x": x.tolist(),
            "points_y": y.tolist(),
            "line_x": fit["line_x"],
            "line_y": fit["line_y"],
        },
        "steps": steps,
        "prediction": {"k1": k1, "k2": k2, "CA0": float(CA0), "CAe": float(CAe)},
    }


def predict_concentration(k1, k2, CA0, CAe, t):
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    kobs = k1 + k2  # since slope*(M+1)/(M+XAe) reduces to (k1+k2) at the CA-only level
    C = CAe + (CA0 - CAe) * np.exp(-kobs * t)
    if C <= 0:
        return None
    return float(C)


def predict_time(k1, k2, CA0, CAe, C):
    if C <= CAe:
        raise ValidationError("Concentration cannot go below the equilibrium concentration C_Ae.")
    if C > CA0:
        raise ValidationError("Target concentration cannot exceed C_A0.")
    kobs = k1 + k2
    t = -np.log((C - CAe) / (CA0 - CAe)) / kobs
    if t < 0:
        return None
    return float(t)


def predict_rate(k1, k2, CA, CR):
    if CA <= 0 or CR < 0:
        raise ValidationError("Concentrations must be non-negative (C_A > 0).")
    return float(k1 * CA - k2 * CR)
