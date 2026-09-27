"""Elementary linear-elastic member formulas, using N, mm and MPa."""

import numpy as np

from ._validation import finite_array, scalar


def axial_response(N_N, lengths_mm, areas_mm2, E_mpa):
    """Return stress [MPa], strain, elongation [mm] per serial segment.

    Uniform axial force, centroidal loading and small strain are assumed.
    lengths/areas are nonempty positive vectors of equal shape; E is scalar.
    """
    N_N, E_mpa = scalar(N_N), scalar(E_mpa, positive=True)
    lengths = finite_array(lengths_mm)
    areas = finite_array(areas_mm2, lengths.shape)
    if lengths.ndim != 1 or lengths.size == 0 or np.any(lengths <= 0) or np.any(areas <= 0):
        raise ValueError("lengths and areas must be nonempty positive vectors")
    sigma = N_N / areas
    strain = sigma / E_mpa
    return sigma, strain, strain * lengths


def cantilever_tip_load(P_N, L_mm, b_mm, h_mm, E_mpa):
    """Return Iz [mm^4], |root stress| [MPa], |tip deflection| [mm].

    Rectangular Euler-Bernoulli cantilever, tip point load, small deflection;
    shear deformation and root stress concentration are excluded.
    """
    P = abs(scalar(P_N))
    L, b, h, E = [scalar(v, positive=True) for v in (L_mm, b_mm, h_mm, E_mpa)]
    Iz = b * h**3 / 12
    return Iz, P * L * h / (2 * Iz), P * L**3 / (3 * E * Iz)


def euler_buckling(E_mpa, area_mm2, inertias_mm4, L_mm, K):
    """Return radii [mm], slenderness, Pcr [N], mean critical stress [MPa].

    Ideal slender elastic column; K is a positive effective length factor.
    Does not predict imperfection sensitivity, plastic or local buckling.
    inertias may be a positive scalar or array; outputs retain its shape.
    """
    E, area, L, K = [scalar(v, positive=True) for v in (E_mpa, area_mm2, L_mm, K)]
    inertias = finite_array(inertias_mm4)
    if inertias.size == 0 or np.any(inertias <= 0):
        raise ValueError("inertias must be positive and nonempty")
    radii = np.sqrt(inertias / area)
    slenderness = K * L / radii
    Pcr = np.pi**2 * E * area / slenderness**2
    return radii, slenderness, Pcr, Pcr / area
