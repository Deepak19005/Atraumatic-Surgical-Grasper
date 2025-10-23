# Atraumatic Surgical Grasper — Pressure Measurement System

**Author:** L. G. Deepak (2022MEB1322)  
**Project:** Atraumatic Surgical Grasper with Integrated Pressure Measurement System.  
Reference: Design and Development report. :contentReference[oaicite:1]{index=1}

---

## Short overview
This repository contains models, simulation inputs and calibration data for an **integrated pressure measurement system** in an atraumatic surgical grasper. The measurement system infers applied pressure at the grasper tips by measuring strain in limb-4 (where a strain gauge is mounted) and the angular position `θ` of the mechanism (obtained from an LVDT). The pipeline:

1. Derive analytical relation between tension `T`, applied pressure `P` and configuration angle `θ`.  
2. Simulate the assembly in ANSYS to obtain `strain` in limb-4 for a grid of `(θ, P)` points.  
3. Save the simulated dataset `strain(θ, P)` to CSV files (calibration table).  
4. During practical use measure limb-4 position by LVDT → infer `θ`; measure strain → look up/interpolate the data table to recover `P`.  

---

## Key mathematical relation (from design notes)
From the derivation used in the design:

\[
\boxed{\,T(\theta,P) \;=\; 7.17\times10^{-5}\; P \;\cot\theta \,}
\]

**Important unit convention:**  
- `P` must be in **Pascals (Pa = N/m²)**.  
- `θ` must be in **radians** when evaluating `cot(θ)`.  
- `T` will be in **Newtons (N)** with the above units.

> This relation gives the axial tension `T` in limb-4 as a function of geometry and applied pressure. Use this to sanity-check ANSYS results and to convert between pressure and internal tension.

---

## ANSYS simulation plan (data-generation)
To create the calibration table:

1. **Choose θ-grid**: simulate for **≥ 10** distinct values of `θ` spanning the expected operating range (e.g., 10°, 20°, … or as needed). Convert degrees → radians for calculations.
2. **Choose pressure levels**: for each `θ`, run ANSYS for multiple pressure magnitudes `P` (e.g., 5–12 levels between minimum and maximum expected pressures).
3. **Record outputs**: for each `(θ_i, P_j)` save:
   - global metadata: `theta_deg, theta_rad, P_Pa`
   - simulation outputs at the strain-gauge location: `strain` (unitless or microstrain as configured)
   - optionally `Tension_T_N` (if you compute T from reaction forces) for cross-validation.
4. **File format**: store as CSV (`calibration_grid.csv`) with header:


5. **Number of points**: recommended grid: `N_theta >= 10` × `N_pressure >= 8` → `>= 80` rows.

---

## Data storage / naming convention
- `ANSYS_exports/` — raw export of ANSYS result files (for traceability).  
- `data/calibration_grid.csv` — processed CSV containing the grid of `(theta, P, strain)`.  
- `scripts/generate_lookup.py` — scripts to read CSV and build interpolation objects.  
- `scripts/calibrate_from_measurement.py` — example runtime code to infer `P` from measured `theta` and `strain`.

---

## Calibration / runtime workflow (practical use)
1. **Measure LVDT displacement** (linear displacement of limb-4 tip).  
2. **Convert LVDT value → θ** using the geometry mapping (either pre-computed function or a small lookup table). Example: `theta = f_lvd(t_disp)` (implement this from your mechanism geometry or a simple polynomial fit).  
3. **Read strain** from strain gauge (ε_meas).  
4. **Recover pressure P**:
   - Extract the slice of the calibration table at the nearest `θ` (or interpolate between two `θ` rows) to get `strain_vs_P` curve.
   - Invert the curve (interpolate) to find `P` corresponding to `ε_meas`.
   - If the measured θ lies between grid θ values, first interpolate strain values between the two nearest theta grids, then invert.
5. **If result lies between discrete pressure points** use linear interpolation (or higher order if desired).

---

## Interpolation & inversion — example Python snippet

> This example assumes a regular grid: `theta_unique` and `pressure_unique` and a matrix `strain_matrix` shaped (N_theta, N_pressure). It shows how to compute `P` from measured `(theta, strain)`.

```python
# scripts/calibrate_from_measurement.py
import numpy as np
from scipy.interpolate import interp1d

# --- load calibration CSV (precomputed) ---
# CSV expected columns: theta_deg, theta_rad, P_Pa, strain
data = np.genfromtxt('data/calibration_grid.csv', delimiter=',', names=True)

# build unique axes
theta_unique = np.unique(data['theta_rad'])
P_unique = np.unique(data['P_Pa'])

# reshape strain into (N_theta, N_pressure)
Ntheta = len(theta_unique)
Np = len(P_unique)
strain_grid = np.zeros((Ntheta, Np))
for i, th in enumerate(theta_unique):
    sel = data['theta_rad'] == th
    row = data[sel]
    # ensure pressures align with P_unique ordering
    # sort row by P_Pa
    order = np.argsort(row['P_Pa'])
    strain_grid[i, :] = row['strain'][order]

def estimate_pressure(theta_meas_rad, strain_meas):
    # 1) find two nearest theta indices for interpolation
    idx = np.searchsorted(theta_unique, theta_meas_rad)
    if idx == 0:
        i_lo, i_hi = 0, 0
    elif idx >= Ntheta:
        i_lo, i_hi = Ntheta-1, Ntheta-1
    else:
        i_lo, i_hi = idx-1, idx

    # 2) for each theta row, build 1D interp function strain->P
    def invert_in_row(i):
        # strain increases (or decreases) monotonically with P in each row.
        s_row = strain_grid[i, :]
        p_row = P_unique
        # if monotonic decreasing, flip arrays for interp
        if s_row[0] > s_row[-1]:
            s_row = s_row[::-1]
            p_row = p_row[::-1]
        inv = interp1d(s_row, p_row, bounds_error=False, fill_value='extrapolate')
        return inv(strain_meas)

    P_lo = invert_in_row(i_lo)
    P_hi = invert_in_row(i_hi)

    # 3) linear interpolation along theta
    if i_lo == i_hi:
        return float(P_lo)
    t_lo, t_hi = theta_unique[i_lo], theta_unique[i_hi]
    w = (theta_meas_rad - t_lo) / (t_hi - t_lo)
    P_est = (1-w)*P_lo + w*P_hi
    return float(P_est)

# Example usage:
theta_meas_deg = 25.0
theta_meas = np.deg2rad(theta_meas_deg)
strain_meas = 1500e-6  # example: 1500 microstrain (if your CSV is in microstrain convert accordingly)
P_est = estimate_pressure(theta_meas, strain_meas)
print("Estimated pressure (Pa):", P_est)
