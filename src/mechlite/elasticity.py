"""Small-strain, 3D isotropic linear elasticity (no plane-stress reduction)."""

from dataclasses import dataclass

import numpy as np

from ._validation import scalar, symmetric_tensor


@dataclass(frozen=True)
class IsotropicElastic:
    """Young's modulus E > 0 (stress units), and -1 < nu < 0.5."""

    E: float
    nu: float

    def __post_init__(self):
        object.__setattr__(self, "E", scalar(self.E, positive=True))
        object.__setattr__(self, "nu", scalar(self.nu))
        if not -1 < self.nu < 0.5:
            raise ValueError("nu must satisfy -1 < nu < 0.5")

    @property
    def G(self):
        """Shear modulus in stress units."""
        return self.E / (2 * (1 + self.nu))

    @property
    def K(self):
        """Bulk modulus in stress units."""
        return self.E / (3 * (1 - 2 * self.nu))

    def stress(self, eps):
        """Map elastic strain (3,3), dimensionless, to stress (3,3)."""
        eps = symmetric_tensor(eps)
        tr = np.trace(eps)
        return self.K * tr * np.eye(3) + 2 * self.G * (eps - tr / 3 * np.eye(3))

    def matrix(self):
        """Return D (6,6): engineering strain -> stress, xx,yy,zz,yz,xz,xy."""
        lam = self.K - 2 * self.G / 3
        D = np.zeros((6, 6))
        D[:3, :3] = lam
        D[:3, :3] += 2 * self.G * np.eye(3)
        D[3:, 3:] = self.G * np.eye(3)
        return D


def elastic_stress(eps, E, nu):
    """Notebook-compatible functional form of IsotropicElastic.stress."""
    return IsotropicElastic(E, nu).stress(eps)


def elasticity_matrix_3d(E, nu):
    """Return 3D isotropic engineering-Voigt stiffness (6,6)."""
    return IsotropicElastic(E, nu).matrix()


def stress_from_strain(eps, E, nu):
    """Phase 1 entry point, with corrected shear stress conversion."""
    return elastic_stress(eps, E, nu)
