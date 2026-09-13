"""
autocatalytic.py
-----------------
Implements the AUTOCATALYTIC REACTION analysis exactly as derived in the
CRE_Reactions.pdf handwritten notes:

    A + R -> R + R          (product R catalyzes its own formation)

Mole balance (rate expressed as second order, first order in each of A and R):

    -r_A = -dC_A/dt = k * C_A * C_R

Total concentration is constant because one mole of A is converted to one
mole of R:

    C0 = C_A + C_R          (C0 = C_A0 + C_R0, constant throughout)

Let  M = C_R0 / C_A0   (R can never be zero, i.e. M > 0, otherwise the
reaction never starts -- exactly as noted: "R can't be zero").

With  X_A = 1 - C_A/C_A0  (conversion of A), integrating the rate law gives
the linear form used for graphical analysis:

    ln[ (M + X_A) / (M (1 - X_A)) ]  =  C_A0 (M + 1) k t

    Plot:  X-axis = t
           Y-axis = ln[ (M + X_A) / (M (1 - X_A)) ]

    slope = C_A0 (M + 1) k   ->   k = slope / [C_A0 (M + 1)]
    intercept (theoretical) = 0
"""
import numpy as np
from .validation import validate_time_concentration, validate_positive, ValidationError
from .regression import linear_fit


def analyze(times, concentrations, CR0):
    t, CA = validate_time_concentration(times, concentrations, min_points=3)
    t = np.array(t)
    CA = np.array(CA)
    CA0 = CA[0]
    validate_positive(CR0, "Initial concentration of R (C_R0)")

    M = CR0 / CA0
    if M <= 0:
        raise ValidationError("C_R0 must be greater than zero: an autocatalytic reaction "
                               "cannot start with zero product present (R can't be zero).")

    C0_total = CA0 + CR0
    XA = 1.0 - CA / CA0

    if np.any(XA >= 1.0) or np.any(XA < 0):
        raise ValidationError("Computed conversion X_A must lie in [0, 1). Check that "
                               "concentration is decreasing and does not exceed C_A0.")

    denom = M * (1.0 - XA)
    if np.any(denom <= 0) or np.any((M + XA) <= 0):
        raise ValidationError("Invalid autocatalytic transformation: (M + X_A) and "
                               "M(1 - X_A) must both be positive.")

    y = np.log((M + XA) / denom)
    x = t

    fit = linear_fit(x, y)
    k = fit["slope"] / (CA0 * (M + 1.0))
    if k <= 0:
        raise ValidationError("The regression produced a non-positive rate constant k. "
                               "The autocatalytic model may not fit this data; check the "
                               "sign/trend of the concentration data.")

    transformed_table = [
        {"time": float(t[i]), "CA": float(CA[i]), "CR": float(C0_total - CA[i]),
         "XA": float(XA[i]), "y": float(y[i])}
        for i in range(len(t))
    ]

    steps = [
        {"title": "Step 1: Original kinetic equation",
         "text": "A + R -> R + R,   -r_A = -dC_A/dt = k * C_A * C_R"},
        {"title": "Step 2: Constant total concentration",
         "text": f"C0 = C_A + C_R = C_A0 + C_R0 = {C0_total:.6g} (constant). M = C_R0/C_A0 = {M:.6g}."},
        {"title": "Step 3: Integrated linear form",
         "text": "ln[ (M + X_A) / (M (1 - X_A)) ] = C_A0 (M + 1) k t"},
        {"title": "Step 4: Data transformation",
         "text": "X_A = 1 - C_A/C_A0 was computed for every point, then y = ln[(M + X_A)/(M(1 - X_A))]."},
        {"title": "Step 5: Graph coordinates",
         "text": "X-axis = t (time), Y-axis = ln[(M + X_A)/(M(1 - X_A))]."},
        {"title": "Step 6: Linear regression",
         "text": f"slope = {fit['slope']:.6g}, intercept = {fit['intercept']:.6g}, R^2 = {fit['r2']:.5f}."},
        {"title": "Step 7: Interpretation",
         "text": "slope = C_A0 (M + 1) k  ->  k = slope / [C_A0 (M + 1)]"},
        {"title": "Step 8: Rate constant",
         "text": f"k = {k:.6g}"},
        {"title": "Step 9: Final rate law",
         "text": f"-r_A = {k:.6g} * C_A * C_R"},
    ]

    return {
        "method": "Autocatalytic Reaction",
        "input_data": [{"time": float(t[i]), "concentration": float(CA[i])} for i in range(len(t))],
        "original_equation": "A + R -> R + R,   -r_A = -dC_A/dt = k*C_A*C_R",
        "CA0": float(CA0), "CR0": float(CR0), "M": float(M), "C0_total": float(C0_total),
        "transformed_table": transformed_table,
        "regression": fit,
        "k": k,
        "rate_law_str": f"-r_A = {k:.6g} * C_A * C_R",
        "x_label": "t (time)",
        "y_label": "ln[(M + X_A) / (M(1 - X_A))]",
        "graph": {
            "x_axis": "t (time)",
            "y_axis": "ln[(M + X_A) / (M(1 - X_A))]",
            "points_x": x.tolist(),
            "points_y": y.tolist(),
            "line_x": fit["line_x"],
            "line_y": fit["line_y"],
        },
        "steps": steps,
        "prediction": {"k": k, "CA0": float(CA0), "CR0": float(CR0), "M": float(M), "C0_total": float(C0_total)},
    }


def predict_concentration(k, CA0, M, C0_total, t):
    """Solve ln[(M+XA)/(M(1-XA))] = CA0(M+1)kt for XA, then CA = CA0(1-XA)."""
    if t < 0:
        raise ValidationError("Time cannot be negative.")
    rhs = CA0 * (M + 1.0) * k * t
    # ln[(M+XA)/(M(1-XA))] = rhs  ->  (M+XA) = M(1-XA) * exp(rhs)
    E = np.exp(rhs)
    # M + XA = M*E - M*E*XA  ->  XA(1 + M*E) = M*E - M  -> XA = M(E-1)/(1+M*E)
    XA = M * (E - 1.0) / (1.0 + M * E)
    if XA < 0 or XA >= 1.0 or np.isnan(XA):
        return None
    CA = CA0 * (1.0 - XA)
    if CA <= 0:
        return None
    return float(CA)


def predict_time(k, CA0, M, C0_total, C):
    if C <= 0:
        raise ValidationError("Concentration must be positive.")
    if C > CA0:
        raise ValidationError("Target concentration cannot exceed C_A0.")
    XA = 1.0 - C / CA0
    if XA >= 1.0 or XA < 0:
        return None
    denom = M * (1.0 - XA)
    if denom <= 0 or (M + XA) <= 0:
        return None
    y = np.log((M + XA) / denom)
    t = y / (CA0 * (M + 1.0) * k)
    if t < 0 or np.isnan(t):
        return None
    return float(t)


def predict_rate(k, C, CR):
    if C <= 0 or CR < 0:
        raise ValidationError("Concentrations must be non-negative (C_A > 0).")
    return float(k * C * CR)
