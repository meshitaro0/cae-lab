"""Shared validation for real, finite numerical inputs."""

import numpy as np


def finite_array(value, shape=None):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError("Expected real values")
    a = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(a)) or (shape is not None and a.shape != shape):
        raise ValueError(f"Expected finite values with shape {shape}")
    return a


def scalar(value, *, positive=False, nonnegative=False):
    value = float(finite_array(value, ()))
    if (positive and value <= 0) or (nonnegative and value < 0):
        raise ValueError("Scalar is outside the permitted range")
    return value


def symmetric_tensor(value):
    a = finite_array(value, (3, 3))
    if not np.allclose(a, a.T, atol=1e-12, rtol=0):
        raise ValueError("Expected a symmetric 3x3 tensor")
    return a


def positive_integer(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError("Expected a positive integer")
    return int(value)
