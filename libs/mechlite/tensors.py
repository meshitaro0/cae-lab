"""Symmetric 3D tensors; Voigt order is xx, yy, zz, yz, xz, xy.

Strain vectors use engineering shear (2*epsilon_ij); stress vectors do not.
All functions operate on one tensor and leave input arrays unchanged.
"""

import numpy as np

from ._validation import finite_array, symmetric_tensor


def stress_to_voigt(sigma):
    """Convert a symmetric stress tensor (3,3) to a vector (6,)."""
    a = symmetric_tensor(sigma)
    return np.array([a[0, 0], a[1, 1], a[2, 2], a[1, 2], a[0, 2], a[0, 1]])


def strain_to_voigt(eps):
    """Convert a small strain tensor (3,3) to engineering strain (6,)."""
    v = stress_to_voigt(eps)
    v[3:] *= 2
    return v


def voigt_to_stress(v):
    """Convert stress components (6,) without dividing shear by two."""
    v = finite_array(v, (6,))
    return np.array([[v[0], v[5], v[4]], [v[5], v[1], v[3]], [v[4], v[3], v[2]]])


def voigt_to_strain(v):
    """Convert engineering strain (6,) to a symmetric tensor (3,3)."""
    v = finite_array(v, (6,)).copy()
    v[3:] /= 2
    return voigt_to_stress(v)


def invariants(sigma):
    """Return (mean stress, deviator, J2, von Mises).

    Mean stress is positive in tension (not pressure positive in compression).
    J2 has stress-squared units; mean and Mises have the input stress units.
    """
    sigma = symmetric_tensor(sigma)
    mean = np.trace(sigma) / 3
    dev = sigma - mean * np.eye(3)
    j2 = np.sum(dev * dev) / 2
    return float(mean), dev, float(j2), float(np.sqrt(3 * j2))


def mises(stress):
    """Return von Mises stress for a symmetric (3,3) tensor."""
    return invariants(stress)[3]
