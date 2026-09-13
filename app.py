"""
app.py
------
Flask application for the Chemical Reaction Kinetics Analyzer.

Routes:
    GET  /                          Home page with method selector
    GET  /integral                  Integral Method page
    GET  /differential              Differential Method page
    GET  /autocatalytic             Autocatalytic Reaction page
    GET  /reversible                Reversible Reaction page
    GET  /irreversible              Irreversible (series) Reaction page
    GET  /which-method               Integral vs Differential method advisor

    POST /analyze/integral          JSON in -> JSON results out
    POST /analyze/differential
    POST /analyze/autocatalytic
    POST /analyze/reversible
    POST /analyze/irreversible

    POST /predict/<method>          JSON in (kinetic params + query) -> prediction

    POST /upload_csv                multipart CSV -> parsed {time:[], concentration:[]}
"""
import io
import os
import traceback

import pandas as pd
from flask import Flask, render_template, request, jsonify

from services.validation import ValidationError
from services import integral_method, differential_method, autocatalytic, reversible, irreversible

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024  # 5 MB upload cap


# --------------------------------------------------------------------------
# Page routes
# --------------------------------------------------------------------------
@app.route("/")
def index():
    return render_template("index.html")


@app.route("/integral")
def integral_page():
    return render_template("integral.html")


@app.route("/differential")
def differential_page():
    return render_template("differential.html")


@app.route("/autocatalytic")
def autocatalytic_page():
    return render_template("autocatalytic.html")


@app.route("/reversible")
def reversible_page():
    return render_template("reversible.html")


@app.route("/irreversible")
def irreversible_page():
    return render_template("irreversible.html")


@app.route("/which-method")
def which_method_page():
    return render_template("which_method.html")


# --------------------------------------------------------------------------
# CSV upload helper
# --------------------------------------------------------------------------
@app.route("/upload_csv", methods=["POST"])
def upload_csv():
    try:
        if "file" not in request.files:
            return jsonify({"error": "No file uploaded."}), 400
        f = request.files["file"]
        if f.filename == "":
            return jsonify({"error": "No file selected."}), 400

        raw = f.read().decode("utf-8-sig", errors="replace")
        df = pd.read_csv(io.StringIO(raw))
        df.columns = [c.strip().lower() for c in df.columns]

        time_col = next((c for c in df.columns if "time" in c or c == "t"), None)
        conc_col = next((c for c in df.columns if "conc" in c or c in ("c", "ca")), None)
        extra_col = next((c for c in df.columns if "cr" in c or "product" in c), None)

        if time_col is None or conc_col is None:
            return jsonify({"error": "CSV must contain a 'Time' column and a 'Concentration' column."}), 400

        payload = {
            "time": df[time_col].tolist(),
            "concentration": df[conc_col].tolist(),
        }
        if extra_col is not None:
            payload["concentration_r"] = df[extra_col].tolist()

        return jsonify(payload)
    except Exception as exc:
        return jsonify({"error": f"Could not parse CSV file: {exc}"}), 400


# --------------------------------------------------------------------------
# Analysis API routes
# --------------------------------------------------------------------------
def _get_lists(data, *keys):
    return [data.get(k) for k in keys]


@app.route("/analyze/integral", methods=["POST"])
def analyze_integral():
    try:
        data = request.get_json(force=True)
        result = integral_method.analyze(data.get("time"), data.get("concentration"))
        return jsonify({"ok": True, "result": result})
    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during analysis. Please check your data."}), 500


@app.route("/analyze/differential", methods=["POST"])
def analyze_differential():
    try:
        data = request.get_json(force=True)
        result = differential_method.analyze(data.get("time"), data.get("concentration"))
        return jsonify({"ok": True, "result": result})
    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during analysis. Please check your data."}), 500


@app.route("/analyze/autocatalytic", methods=["POST"])
def analyze_autocatalytic():
    try:
        data = request.get_json(force=True)
        CR0 = float(data.get("CR0"))
        result = autocatalytic.analyze(data.get("time"), data.get("concentration"), CR0)
        return jsonify({"ok": True, "result": result})
    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except (TypeError, ValueError):
        return jsonify({"ok": False, "error": "C_R0 must be a valid positive number."}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during analysis. Please check your data."}), 500


@app.route("/analyze/reversible", methods=["POST"])
def analyze_reversible():
    try:
        data = request.get_json(force=True)
        CR0 = float(data.get("CR0"))
        CAe = float(data.get("CAe"))
        result = reversible.analyze(data.get("time"), data.get("concentration"), CR0, CAe)
        return jsonify({"ok": True, "result": result})
    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except (TypeError, ValueError):
        return jsonify({"ok": False, "error": "C_R0 and C_Ae must be valid numbers."}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during analysis. Please check your data."}), 500


@app.route("/analyze/irreversible", methods=["POST"])
def analyze_irreversible():
    try:
        data = request.get_json(force=True)
        times_r = data.get("time_r") or None
        conc_r = data.get("concentration_r") or None
        result = irreversible.analyze(data.get("time"), data.get("concentration"), times_r, conc_r)
        return jsonify({"ok": True, "result": result})
    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during analysis. Please check your data."}), 500


# --------------------------------------------------------------------------
# Prediction API routes
# --------------------------------------------------------------------------
@app.route("/predict/<method>", methods=["POST"])
def predict(method):
    try:
        data = request.get_json(force=True)
        mode = data.get("mode")  # "concentration" | "time" | "rate"
        params = data.get("params", {})

        if method == "integral" or method == "differential":
            mod = integral_method if method == "integral" else differential_method
            order = float(params["order"])
            k = float(params["k"])
            CA0 = float(params["CA0"])
            if mode == "concentration":
                t = float(data["value"])
                val = mod.predict_concentration(order, k, CA0, t)
                if val is None:
                    return jsonify({"ok": False, "error": "The requested time leads to a non-physical "
                                                            "(negative or undefined) concentration for this model."})
                return jsonify({"ok": True, "value": val, "label": "Concentration"})
            elif mode == "time":
                C = float(data["value"])
                val = mod.predict_time(order, k, CA0, C)
                if val is None:
                    return jsonify({"ok": False, "error": "The requested concentration cannot be reached "
                                                            "using this kinetic model within the valid physical domain."})
                return jsonify({"ok": True, "value": val, "label": "Time"})
            elif mode == "rate":
                C = float(data["value"])
                val = mod.predict_rate(order, k, C)
                return jsonify({"ok": True, "value": val, "label": "Rate"})

        elif method == "autocatalytic":
            k = float(params["k"]); CA0 = float(params["CA0"])
            M = float(params["M"]); C0_total = float(params["C0_total"])
            if mode == "concentration":
                t = float(data["value"])
                val = autocatalytic.predict_concentration(k, CA0, M, C0_total, t)
                if val is None:
                    return jsonify({"ok": False, "error": "The requested time leads to a non-physical result."})
                return jsonify({"ok": True, "value": val, "label": "Concentration"})
            elif mode == "time":
                C = float(data["value"])
                val = autocatalytic.predict_time(k, CA0, M, C0_total, C)
                if val is None:
                    return jsonify({"ok": False, "error": "The requested concentration cannot be reached "
                                                            "using this kinetic model."})
                return jsonify({"ok": True, "value": val, "label": "Time"})
            elif mode == "rate":
                C = float(data["value"]); CR = C0_total - C
                val = autocatalytic.predict_rate(k, C, CR)
                return jsonify({"ok": True, "value": val, "label": "Rate"})

        elif method == "reversible":
            k1 = float(params["k1"]); k2 = float(params["k2"])
            CA0 = float(params["CA0"]); CAe = float(params["CAe"])
            if mode == "concentration":
                t = float(data["value"])
                val = reversible.predict_concentration(k1, k2, CA0, CAe, t)
                if val is None:
                    return jsonify({"ok": False, "error": "The requested time leads to a non-physical result."})
                return jsonify({"ok": True, "value": val, "label": "Concentration"})
            elif mode == "time":
                C = float(data["value"])
                val = reversible.predict_time(k1, k2, CA0, CAe, C)
                if val is None:
                    return jsonify({"ok": False, "error": "This concentration cannot be reached; it may be "
                                                            "below the equilibrium concentration C_Ae."})
                return jsonify({"ok": True, "value": val, "label": "Time"})
            elif mode == "rate":
                CA = float(data["value"])
                CR = (CA0 + params.get("CR0", 0)) - CA if params.get("CR0") is not None else None
                val = reversible.predict_rate(k1, k2, CA, CR if CR is not None else 0.0)
                return jsonify({"ok": True, "value": val, "label": "Rate"})

        elif method == "irreversible":
            k1 = float(params["k1"]); CA0 = float(params["CA0"])
            k2 = params.get("k2")
            k2 = float(k2) if k2 not in (None, "", "null") else None
            species = params.get("species", "A")
            if mode == "concentration":
                t = float(data["value"])
                if species == "R" and k2 is not None:
                    val = irreversible.predict_concentration_R(k1, k2, CA0, t)
                else:
                    val = irreversible.predict_concentration_A(k1, CA0, t)
                return jsonify({"ok": True, "value": val, "label": "Concentration"})
            elif mode == "time":
                C = float(data["value"])
                val = irreversible.predict_time_A(k1, CA0, C)
                if val is None:
                    return jsonify({"ok": False, "error": "This concentration cannot be reached using this model."})
                return jsonify({"ok": True, "value": val, "label": "Time"})
            elif mode == "rate":
                C = float(data["value"])
                val = irreversible.predict_rate_A(k1, C)
                return jsonify({"ok": True, "value": val, "label": "Rate"})

        return jsonify({"ok": False, "error": "Unknown method or mode."}), 400

    except ValidationError as ve:
        return jsonify({"ok": False, "error": str(ve)}), 400
    except (KeyError, TypeError, ValueError) as e:
        return jsonify({"ok": False, "error": f"Invalid input for prediction: {e}"}), 400
    except Exception:
        app.logger.error(traceback.format_exc())
        return jsonify({"ok": False, "error": "Unexpected error during prediction."}), 500


# --------------------------------------------------------------------------
@app.errorhandler(404)
def not_found(e):
    return render_template("base.html", content_404=True), 404


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=False)
