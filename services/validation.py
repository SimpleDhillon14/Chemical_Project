"""
validation.py
--------------
Shared input validation for the Chemical Reaction Kinetics Analyzer.

Every analysis method (Integral, Differential, Autocatalytic, Reversible,
Irreversible) receives raw (time, concentration) pairs typed by the user or
uploaded via CSV. This module makes sure the data is numeric, physically
sane, and large enough for a meaningful regression BEFORE it reaches the
math in services/*.py, so the Flask routes can return a clean JSON error
instead of raising an exception.
"""
import math


class ValidationError(Exception):
    """Raised when experimental data fails a physical/numerical sanity check."""
    pass


def validate_time_concentration(times, concentrations, min_points=3, require_monotonic_time=True):
    """
    Validate parallel lists of time / concentration values.

    Returns (times, concentrations) as clean Python float lists.
    Raises ValidationError with a human-readable message otherwise.
    """
    if times is None or concentrations is None:
        raise ValidationError("Time and Concentration data are required.")

    if len(times) != len(concentrations):
        raise ValidationError(
            f"Time column has {len(times)} values but Concentration column has "
            f"{len(concentrations)} values. They must match row for row."
        )

    if len(times) < min_points:
        raise ValidationError(
            f"At least {min_points} data points are required for a meaningful "
            f"regression. Only {len(times)} were provided."
        )

    clean_t, clean_c = [], []
    for i, (t_raw, c_raw) in enumerate(zip(times, concentrations), start=1):
        if t_raw is None or c_raw is None or str(t_raw).strip() == "" or str(c_raw).strip() == "":
            raise ValidationError(f"Row {i} has an empty Time or Concentration value.")
        try:
            t = float(t_raw)
            c = float(c_raw)
        except (TypeError, ValueError):
            raise ValidationError(f"Row {i} contains a non-numeric value ('{t_raw}', '{c_raw}').")

        if math.isnan(t) or math.isnan(c) or math.isinf(t) or math.isinf(c):
            raise ValidationError(f"Row {i} contains NaN or infinite values.")

        if t < 0:
            raise ValidationError(f"Row {i}: Time cannot be negative ({t}).")

        if c <= 0:
            raise ValidationError(
                f"Row {i}: Concentration must be strictly greater than zero ({c}). "
                f"Zero or negative concentrations are not physically valid and cannot "
                f"be used in logarithmic transformations."
            )

        clean_t.append(t)
        clean_c.append(c)

    # sort by time (defensive: user may paste rows out of order)
    order = sorted(range(len(clean_t)), key=lambda i: clean_t[i])
    clean_t = [clean_t[i] for i in order]
    clean_c = [clean_c[i] for i in order]

    if require_monotonic_time and len(set(clean_t)) != len(clean_t):
        raise ValidationError("Duplicate time values were found. Each time point must be unique.")

    return clean_t, clean_c


def safe_log(x, label="value"):
    if x <= 0:
        raise ValidationError(f"Cannot take the logarithm of a non-positive {label} ({x}).")
    return math.log(x)


def safe_divide(numerator, denominator, label="expression"):
    if denominator == 0:
        raise ValidationError(f"Division by zero encountered while computing {label}.")
    return numerator / denominator


def validate_positive(value, name):
    if value is None or value <= 0:
        raise ValidationError(f"{name} must be a positive number.")
    return value


def validate_non_negative(value, name):
    if value is None or value < 0:
        raise ValidationError(f"{name} cannot be negative.")
    return value
