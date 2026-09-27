"""Two-node axial bar FEM: linear displacement, small strain, dense solver.

Use consistent units (e.g. mm, N, MPa). No bending, buckling or dynamics.
Connectivity must orient every element from smaller to larger x.
"""

import numpy as np

from ._validation import finite_array, scalar, positive_integer


def bar_element(xa, xb, EA_of_x, p_of_x, n_gauss=2):
    """Return ke (2,2) [force/length], fe (2,) [force].

    EA_of_x and p_of_x accept one scalar coordinate, returning scalar axial
    rigidity [force] > 0 and distributed load [force/length], respectively.
    Both callbacks are evaluated once per quadrature point.
    """
    xa, xb = scalar(xa), scalar(xb)
    if xa >= xb:
        raise ValueError("Expected xa < xb")
    le = scalar(xb - xa, positive=True)
    xi, weights = np.polynomial.legendre.leggauss(positive_integer(n_gauss))
    points = (1 - xi) * xa / 2 + (1 + xi) * xb / 2
    EA = np.array([scalar(EA_of_x(float(x)), positive=True) for x in points])
    load = np.array([scalar(p_of_x(float(x))) for x in points])
    N = np.column_stack(((1 - xi) / 2, (1 + xi) / 2))
    B = np.array([-1.0, 1.0]) / le
    return np.outer(B, B) * np.dot(weights, EA) * le / 2, N.T @ (weights * load) * le / 2


def _mesh(x, conn, E, A):
    x = finite_array(x)
    conn = np.asarray(conn)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("Expected at least two nodal coordinates")
    if conn.ndim != 2 or conn.shape[1] != 2 or len(conn) == 0:
        raise ValueError("Expected connectivity (ne,2), ne > 0")
    if not np.issubdtype(conn.dtype, np.integer):
        raise TypeError("Connectivity must contain integer node indices")
    if np.any(conn < 0) or np.any(conn >= len(x)):
        raise IndexError("Connectivity index out of range")
    lengths = x[conn[:, 1]] - x[conn[:, 0]]
    if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0):
        raise ValueError("Each element requires finite xa < xb")
    E, A = finite_array(E, (len(conn),)), finite_array(A, (len(conn),))
    if np.any(E <= 0) or np.any(A <= 0):
        raise ValueError("E and A must be positive")
    return x, conn, E, A


def assemble_bar(x, conn, E, A, p, point_forces):
    """Assemble K (n,n), f (n,) with element-constant E/A/p (ne,).

    x and point_forces have shape (n,); conn has shape (ne,2).
    Inputs are not modified; integer loads are converted to floating point.
    """
    x, conn, E, A = _mesh(x, conn, E, A)
    p = finite_array(p, (len(conn),))
    f = finite_array(point_forces, x.shape).copy()
    K = np.zeros((len(x), len(x)))
    for i, c in enumerate(conn):
        ke, fe = bar_element(*x[c], lambda _: E[i] * A[i], lambda _: p[i])
        K[np.ix_(c, c)] += ke
        f[c] += fe
    return K, f


def solve_dirichlet(K, f, prescribed):
    """Return displacement d (n,), residual K@d-f (n,) using original K/f.

    prescribed maps integer DOF indices to finite displacement values.
    Constrained residuals are reactions; free residuals should approach zero.
    Singular free systems raise numpy.linalg.LinAlgError.
    """
    f = finite_array(f)
    if f.ndim != 1 or f.size == 0:
        raise ValueError("Expected nonempty load vector")
    K = finite_array(K, (len(f), len(f)))
    if not prescribed:
        raise ValueError("At least one prescribed DOF is required")
    d = np.zeros(len(f))
    for index, value in prescribed.items():
        if isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)):
            raise TypeError("Prescribed DOF indices must be integers")
        if not 0 <= index < len(f):
            raise IndexError("Prescribed DOF out of range")
        d[index] = scalar(value)
    fixed = np.array(list(prescribed), dtype=int)
    free = np.setdiff1d(np.arange(len(f)), fixed)
    if free.size:
        d[free] = np.linalg.solve(K[np.ix_(free, free)], f[free] - K[np.ix_(free, fixed)] @ d[fixed])
    return d, K @ d - f


def recover_bar(x, conn, E, A, d):
    """Return element strain, stress, axial force, each (ne,)."""
    x, conn, E, A = _mesh(x, conn, E, A)
    d = finite_array(d, x.shape)
    strain = (d[conn[:, 1]] - d[conn[:, 0]]) / (x[conn[:, 1]] - x[conn[:, 0]])
    stress = E * strain
    return strain, stress, A * stress


def uniform_bar_mesh(L, ne):
    """Return x (ne+1,), conn (ne,2) on [0,L], L > 0.

    Loads are specified explicitly by node index in assemble_bar, avoiding the
    notebook helper's silent snapping of a requested load position to a node.
    """
    L, ne = scalar(L, positive=True), positive_integer(ne)
    return np.linspace(0, L, ne + 1), np.column_stack((np.arange(ne), np.arange(1, ne + 1)))
