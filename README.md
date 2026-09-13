# Chemical Reaction Kinetics Analyzer

An interactive Flask web application for analyzing experimental concentration–time
data using the graphical and mathematical methods taught in Chemical Reaction
Engineering (CRE): the **Integral Method**, the **Differential Method**, and the
specialized linearizations for **Autocatalytic**, **Reversible**, and
**Irreversible (series)** reactions.

The application is built for transparency: every equation, transformation,
table, graph, slope, intercept, R², and final kinetic parameter is shown
step-by-step, so a professor or student can see exactly how each answer was
derived — not just the final number.

---

## Features

- Five CRE analysis methods, each following the exact equations and graphical
  coordinates from the course material (see *CRE Methodology Note* below).
- Manual, editable Time/Concentration data table (add/remove/edit rows) **or**
  CSV upload.
- Full input validation (numeric-only, positive concentrations, no
  division-by-zero or invalid logarithms, minimum data points) with clear
  error messages instead of crashes.
- Step-by-step explanation panel for every analysis.
- Interactive Plotly graphs: experimental points, best-fit line, axis labels,
  slope/intercept/R² annotations.
- Order-comparison table (Integral Method) covering zero, first, second,
  third, and fractional orders via a systematic R² grid search.
- Prediction calculator: concentration at a given time, time to reach a given
  concentration, and instantaneous rate — all validated for physical
  meaningfulness (no negative concentrations or times).
- "Quiet Luxury" color palette (cream / gold / dusty rose / taupe) applied
  throughout the UI.

## Supported Kinetic Methods

| # | Method | Equation | Straight-line test |
|---|--------|----------|---------------------|
| 1 | Integral Method | -r_A = kC_A^n | n=1: -ln(C_A/C_A0) vs t; n≠1: C_A^(1-n) vs t |
| 2 | Differential Method | -r_A = kC_A^n | log10(-dC_A/dt) vs log10(C_A) |
| 3 | Autocatalytic (A+R→R+R) | -r_A = kC_AC_R | ln[(M+X_A)/(M(1-X_A))] vs t |
| 4 | Reversible (A⇌R) | -r_A = k1C_A - k2C_R | -ln[(C_A-C_Ae)/(C_A0-C_Ae)] vs t |
| 5 | Irreversible series (A→R→S) | -r_A = k1C_A | -ln(C_A/C_A0) vs t (+ optional nonlinear fit for k2) |

## Project Structure

```
chemical-kinetics-analyzer/
├── app.py
├── requirements.txt
├── Procfile
├── runtime.txt
├── README.md
├── .gitignore
├── services/
│   ├── __init__.py
│   ├── validation.py
│   ├── regression.py
│   ├── integral_method.py
│   ├── differential_method.py
│   ├── autocatalytic.py
│   ├── reversible.py
│   └── irreversible.py
├── templates/
│   ├── base.html
│   ├── index.html
│   ├── integral.html
│   ├── differential.html
│   ├── autocatalytic.html
│   ├── reversible.html
│   ├── irreversible.html
│   └── results.html
└── static/
    ├── css/style.css
    └── js/main.js
```

## Installation (local)

```bash
git clone <your-repo-url>
cd chemical-kinetics-analyzer
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## Running Locally

```bash
python app.py
```

Then open http://127.0.0.1:5000 in your browser.

For a production-style local run:

```bash
gunicorn app:app --bind 0.0.0.0:5000
```

## Input Format

CSV upload expects a header row containing a time-like column (`time` or `t`)
and a concentration-like column (`concentration`, `c`, or `ca`). An optional
third column containing `cr` or `product` in its name is read as C_R data
(used by the Irreversible-reaction page).

```csv
Time,Concentration
0,10
20,8
40,6
60,5
120,3
180,2
300,1
```

## Deployment (Render)

See the step-by-step deployment guide provided alongside this project
(Phase 7). In short:

- **Environment:** Python
- **Build Command:** `pip install -r requirements.txt`
- **Start Command:** `gunicorn app:app --bind 0.0.0.0:$PORT` (already in `Procfile`)

## Technologies Used

- **Backend:** Python, Flask
- **Analysis:** NumPy, Pandas, SciPy (`scipy.stats.linregress`, `scipy.optimize.curve_fit`)
- **Graphs:** Plotly.js (CDN)
- **Frontend:** HTML, CSS (custom "Quiet Luxury" palette), Bootstrap 5, vanilla JavaScript
- **Deployment:** Gunicorn + Render

## CRE Methodology Note

All equations and graphical coordinates implemented here were taken directly
from the two supplied course-material excerpts:

- `CRE_Diff.pdf` — Differential Method of Analysis (Fig. 3.17, 3.18, 3.19,
  Example 3.2), including the Michaelis–Menten-type linearizations (Eqs. 61–62).
- `CRE_Reactions.pdf` — handwritten notes covering first/second/nth/zero-order
  integral analysis, autocatalytic reactions, reversible reactions, and
  irreversible series reactions (A→R→S), including the t_max and C_R,max
  relationships.

Where the reversible-reaction derivation in the notes used a general `M`
(ratio of initial R to initial A) framework, that same general form is
implemented here (rather than assuming C_R0 = 0), so it remains consistent
with the notes even when a non-zero initial product concentration is given.

## Limitations and Assumptions

- The Differential Method estimates dC_A/dt numerically (central/forward/backward
  differences) rather than by hand-drawing a smooth curve, since the app has
  no way to "draw by eye." This is the standard numerical stand-in for that
  step and reproduces the same log–log linearization taught in the notes.
- The Irreversible-reaction k2 estimate requires C_R vs t data in addition to
  C_A vs t data; without it, only k1 is reported.
- Fractional-order results should be interpreted alongside the R² comparison
  table — a very close R² between two orders (e.g., 1.4 and 1.5) usually
  means either is a reasonable engineering approximation.
- This tool is intended for coursework and educational use; always sanity-check
  fitted parameters against known chemistry for the system being studied.
