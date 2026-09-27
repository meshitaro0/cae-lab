"""Small-strain 3D J2 radial return with linear isotropic hardening.

No finite rotation, plane-stress enforcement, kinematic hardening or consistent
tangent is supplied. All stress-like quantities must use one consistent unit.
The caller owns the committed state; update never commits a trial implicitly.
"""

from dataclasses import dataclass, field

import numpy as np

from ._validation import finite_array, scalar, symmetric_tensor, positive_integer
from .elasticity import IsotropicElastic
from .tensors import invariants


@dataclass(frozen=True)
class J2State:
    """Stress, plastic strain (3,3), and nonnegative equivalent strain alpha.

    Arrays are copied and read-only to prevent accidental trial-state mutation.
    The caller must supply a physically consistent, previously committed state.
    """

    stress: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    plastic_strain: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    alpha: float = 0.0

    def __post_init__(self):
        for name in ("stress", "plastic_strain"):
            a = symmetric_tensor(getattr(self, name)).copy()
            a.setflags(write=False)
            object.__setattr__(self, name, a)
        object.__setattr__(self, "alpha", scalar(self.alpha, nonnegative=True))


@dataclass(frozen=True)
class J2Material:
    """E, nu, initial yield stress > 0, hardening H >= 0.

    yield_tolerance is an absolute stress tolerance (default 1e-7, as in the
    MPa notebook); adapt it when using a different stress unit.
    """

    E: float
    nu: float
    sigma_y0: float
    H: float
    yield_tolerance: float = 1e-7

    def __post_init__(self):
        elastic = IsotropicElastic(self.E, self.nu)
        object.__setattr__(self, "E", elastic.E)
        object.__setattr__(self, "nu", elastic.nu)
        for name, positive in (("sigma_y0", True), ("H", False), ("yield_tolerance", False)):
            object.__setattr__(self, name, scalar(getattr(self, name), positive=positive, nonnegative=True))

    def update(self, deps, state):
        """Return (new J2State, dalpha) for a total strain increment (3,3)."""
        if not isinstance(state, J2State):
            raise TypeError("state must be a J2State")
        elastic = IsotropicElastic(self.E, self.nu)
        trial = state.stress + elastic.stress(deps)
        _, dev, _, q = invariants(trial)
        f = q - self.sigma_y0 - self.H * state.alpha
        if f <= self.yield_tolerance:
            return J2State(trial, state.plastic_strain, state.alpha), 0.0
        dalpha = f / (3 * elastic.G + self.H)
        dep = 1.5 * dalpha * dev / q
        return J2State(trial - 2 * elastic.G * dep,
                       state.plastic_strain + dep, state.alpha + dalpha), float(dalpha)


def j2_update(deps, sigma_n, ep_n, alpha_n, E, nu, sigma_y0, H):
    """Phase 5 signature; return (stress, plastic strain, alpha, dalpha).

    Returned arrays are independent writable copies, including elastic steps.
    """
    new, da = J2Material(E, nu, sigma_y0, H).update(deps, J2State(sigma_n, ep_n, alpha_n))
    return new.stress.copy(), new.plastic_strain.copy(), new.alpha, da


def accumulate_alpha(dep_history):
    """Return [0, cumulative sqrt(2/3 dep:dep)] for history (m,3,3).

    Increments are symmetric plastic strain tensors, not engineering Voigt
    vectors. Shape (0,3,3) returns [0]. J2 plasticity assumes isochoric flow.
    """
    dep = finite_array(dep_history)
    if dep.ndim != 3 or dep.shape[1:] != (3, 3):
        raise ValueError("Expected shape (m,3,3)")
    for value in dep:
        symmetric_tensor(value)
    norms = np.sqrt(2 / 3 * np.sum(dep * dep, axis=(1, 2)))
    return np.concatenate(([0.0], np.cumsum(norms)))


def integrate_strain_path(strain_points, material, initial_state=None, subdivisions=1):
    """Integrate a piecewise-linear total strain path (m,3,3), m >= 1.

    initial_state belongs to strain_points[0]; default is an unstressed state.
    Every segment uses the same positive integer number of subdivisions.
    Return stress/ep (n,3,3), alpha (n,), dalpha (n-1,), including initial state.
    Unlike notebook drivers, material parameters and tensor paths are explicit.
    """
    points = finite_array(strain_points)
    if points.ndim != 3 or points.shape[1:] != (3, 3) or len(points) < 1:
        raise ValueError("Expected a nonempty path (m,3,3)")
    for point in points:
        symmetric_tensor(point)
    nd = positive_integer(subdivisions)
    if not isinstance(material, J2Material):
        raise TypeError("material must be a J2Material")
    state = J2State() if initial_state is None else initial_state
    if not isinstance(state, J2State):
        raise TypeError("initial_state must be a J2State")
    n = (len(points) - 1) * nd + 1
    stress, ep = np.empty((2, n, 3, 3))
    alpha, da = np.empty(n), np.empty(n - 1)
    stress[0], ep[0], alpha[0] = state.stress, state.plastic_strain, state.alpha
    step = 0
    for start, end in zip(points[:-1], points[1:]):
        deps = (end - start) / nd
        for _ in range(nd):
            state, da[step] = material.update(deps, state)
            step += 1
            stress[step], ep[step], alpha[step] = state.stress, state.plastic_strain, state.alpha
    return stress, ep, alpha, da
