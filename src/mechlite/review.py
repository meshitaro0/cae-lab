"""Phase 6 review metrics, not automatic physical acceptance criteria.

The energy ledger below is the exercise's convention, not an LS-DYNA reader.
Momentum and forces refer to one common direction and consistent units.
"""

import numpy as np

from ._validation import finite_array, scalar


def checked_history(data):
    """Validate t and drive/support/P/K/U/D/H/W, all arrays (n,), n >= 2.

    t strictly increases. P is momentum, K/U kinetic/strain energy, D physical
    dissipation, H artificial energy, W external work. D/H/W start at zero;
    K/U/P may start nonzero. No solver-specific accounting is inferred.
    """
    t = finite_array(data["t"])
    if t.ndim != 1 or len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError("t must be a strictly increasing vector of length >= 2")
    out = {"t": t}
    for key in ("drive", "support", "P", "K", "U", "D", "H", "W"):
        out[key] = finite_array(data[key], t.shape)
    if any(abs(out[key][0]) > 1e-12 for key in ("D", "H", "W")):
        raise ValueError("D, H, W must start at zero")
    return out


def audit_history(data, energy_ref, force_ref):
    """Return normalized energy/impulse residuals and their absolute maxima.

    r_E = (K-K0+U-U0+D+H-W)/energy_ref;
    r_P = (integral(drive+support)dt-P+P0)/(force_ref*(t_end-t_start)).
    Trapezoidal integration supports nonuniform times. eta_K=max(K)/energy_ref
    and eta_H=max(abs(H))/energy_ref. Both references are positive scalars.
    Inputs are preserved; a small residual alone does not establish validity.
    """
    data = checked_history(data)
    energy_ref = scalar(energy_ref, positive=True)
    force_ref = scalar(force_ref, positive=True)
    r_E = (data["K"] - data["K"][0] + data["U"] - data["U"][0]
           + data["D"] + data["H"] - data["W"]) / energy_ref
    force = data["drive"] + data["support"]
    impulse = np.concatenate(([0.0], np.cumsum((force[1:] + force[:-1]) / 2 * np.diff(data["t"]))))
    r_P = (impulse - data["P"] + data["P"][0]) / (force_ref * (data["t"][-1] - data["t"][0]))
    return dict(r_E=r_E, r_P=r_P, max_abs_r_E=float(np.max(np.abs(r_E))),
                max_abs_r_P=float(np.max(np.abs(r_P))), eta_K=float(np.max(data["K"]) / energy_ref),
                eta_H=float(np.max(np.abs(data["H"])) / energy_ref))


def mesh_changes(values, floor):
    """Return successive absolute/relative changes for coarse-to-fine values.

    relative[i] = abs(v[i+1]-v[i])/max(abs(v[i+1]), floor).
    values is a finite vector of length >= 2; floor > 0 has the same units.
    Changes are sensitivities, not error estimates against an exact solution.
    """
    values = finite_array(values)
    floor = scalar(floor, positive=True)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("Expected at least two values")
    absolute = np.abs(np.diff(values))
    return dict(absolute=absolute, relative=absolute / np.maximum(np.abs(values[1:]), floor))
