"""Independent physical checks plus regression against submitted functions.

Run from the repository root: python -m unittest discover -s tests -v
"""

import ast
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "libs"))

from mechlite import IsotropicElastic, J2Material, J2State
from mechlite import bar, elasticity, members, plasticity, review, tensors


def notebook_functions(phase):
    """Load only function definitions, without executing plots or experiments."""
    notebook = json.loads((ROOT / "notebooks" / f"Phase{phase}.ipynb").read_text(encoding="utf-8"))
    env = {"np": np, "I": np.eye(3), "E": 210000.0, "nu": 0.3, "sigma_y0": 250.0, "H": 1000.0}
    nodes = []
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            nodes.extend(n for n in ast.parse("".join(cell["source"])).body if isinstance(n, ast.FunctionDef))
    exec(compile(ast.Module(body=nodes, type_ignores=[]), f"Phase{phase}.ipynb", "exec"), env)
    return env


class TensorTests(unittest.TestCase):
    def test_voigt_roundtrip_and_work_conjugacy(self):
        eps = np.array([[.001, .002, -.003], [.002, -.001, .004], [-.003, .004, .002]])
        sig = np.array([[200., 30., -40.], [30., -20., 50.], [-40., 50., 10.]])
        assert_allclose(tensors.voigt_to_strain(tensors.strain_to_voigt(eps)), eps)
        assert_allclose(tensors.voigt_to_stress(tensors.stress_to_voigt(sig)), sig)
        self.assertAlmostEqual(np.sum(sig * eps), tensors.stress_to_voigt(sig) @ tensors.strain_to_voigt(eps))

    def test_pure_shear_regression_fix(self):
        eps = np.array([[0., .002, 0.], [.002, 0., 0.], [0., 0., 0.]])
        mat = IsotropicElastic(210000, .3)
        expected = 2 * (210000 / 2.6) * eps
        assert_allclose(elasticity.stress_from_strain(eps, mat.E, mat.nu), expected)
        assert_allclose(tensors.voigt_to_stress(mat.matrix() @ tensors.strain_to_voigt(eps)), expected)
        # The old Phase 1 function incorrectly halves stress shear components.
        assert_allclose(notebook_functions(1)["stress_from_strain"](eps, mat.E, mat.nu), expected / 2)

    def test_rotation_and_mises(self):
        sig = np.array([[200., 70., 0.], [70., -50., 40.], [0., 40., 10.]])
        Q, _ = np.linalg.qr(np.random.default_rng(42).normal(size=(3, 3)))
        vals = np.linalg.eigvalsh(sig)
        expected = np.sqrt(((vals[0]-vals[1])**2 + (vals[1]-vals[2])**2 + (vals[2]-vals[0])**2) / 2)
        self.assertAlmostEqual(tensors.mises(Q.T @ sig @ Q), expected)
        self.assertEqual(tensors.mises(100 * np.eye(3)), 0.)
        mat = IsotropicElastic(100000, .2)
        assert_allclose(mat.stress(Q.T @ sig @ Q * 1e-6), Q.T @ mat.stress(sig * 1e-6) @ Q, atol=1e-10)

    def test_invalid_inputs(self):
        for E, nu in ((0, .3), (np.nan, .3), (1, .5), (1, -1)):
            with self.assertRaises(ValueError):
                IsotropicElastic(E, nu)
        for a in (np.zeros((2, 2)), np.full((3, 3), np.nan), np.triu(np.ones((3, 3))), np.eye(3)*1j):
            with self.assertRaises(ValueError):
                tensors.mises(a)


class MemberTests(unittest.TestCase):
    def test_notebook_regression(self):
        old = notebook_functions(3)
        for name, args in (
            ("axial_response", (30000, [300, 200], [200, 100], 210000)),
            ("cantilever_tip_load", (100, 1000, 20, 40, 210000)),
            ("euler_buckling", (210000, 800, [26666, 106666], 1000, 1.0)),
        ):
            for actual, expected in zip(getattr(members, name)(*args), old[name](*args)):
                assert_allclose(actual, expected)

    def test_physical_scalings_and_validation(self):
        a = members.cantilever_tip_load(100, 1000, 20, 40, 210000)
        b = members.cantilever_tip_load(-100, 2000, 20, 40, 210000)
        assert_allclose(b, [a[0], 2*a[1], 8*a[2]])
        p = members.euler_buckling(210000, 800, [26666], 1000, 1)[2]
        assert_allclose(members.euler_buckling(210000, 800, [26666], 2000, 1)[2], p / 4)
        with self.assertRaises(ValueError):
            members.axial_response(1, [1], [0], 1)
        with self.assertRaises(ValueError):
            members.cantilever_tip_load(1, 1, 1, 1, np.nan)


class BarTests(unittest.TestCase):
    def test_element_quadrature_callbacks(self):
        for n in (2, 3, 4):
            calls = []
            def rigidity(x):
                self.assertIsInstance(x, float)
                calls.append(x)
                return 60.0
            ke, fe = bar.bar_element(0, 2, rigidity, lambda x: 3*x, n)
            assert_allclose(ke, [[30, -30], [-30, 30]])
            assert_allclose(fe, [2, 4])
            assert_allclose(ke @ [1, 1], 0)
            self.assertEqual(len(calls), n)

    def test_assembly_reactions_and_nonzero_displacement(self):
        x, conn = bar.uniform_bar_mesh(2, 2)
        loads = np.array([0, 0, 3])
        K, f = bar.assemble_bar(x, conn, [100, 100], [2, 2], [2, 2], loads)
        assert_array_equal(loads, [0, 0, 3])
        d, r = bar.solve_dirichlet(K, f, {0: .1})
        assert_allclose(d, [.1, .13, .15])
        assert_allclose(r, [-7, 0, 0], atol=1e-12)
        strain, stress, force = bar.recover_bar(x, conn, [100, 100], [2, 2], d)
        assert_allclose(strain, [.03, .02])
        assert_allclose(stress, [3, 2])
        assert_allclose(force, [6, 4])
        all_d, all_r = bar.solve_dirichlet(K, f, {0: .1, 1: .1, 2: .1})
        assert_allclose(all_d, .1)
        assert_allclose(all_r, -f)

    def test_notebook_regression_and_order(self):
        old = notebook_functions(4)
        x, conn = bar.uniform_bar_mesh(100, 4)
        E, A, p, loads = np.full(4, 210000.), np.full(4, 100.), np.full(4, 3.), np.zeros(5)
        args = (x, conn, E, A, p, loads)
        K, f = bar.assemble_bar(*args)
        old_K, old_f = old["assemble_bar"](*args)
        assert_allclose(K, old_K)
        assert_allclose(f, old_f)
        K2, f2 = bar.assemble_bar(x, conn[::-1], E[::-1], A[::-1], p[::-1], loads)
        assert_allclose(K2, K)
        assert_allclose(f2, f)
        d, r = bar.solve_dirichlet(K, f, {0: .2})
        for a, b in zip((d, r), old["solve_dirichlet"](K, f, {0: .2})):
            assert_allclose(a, b, atol=1e-10)
        for a, b in zip(bar.recover_bar(x, conn, E, A, d), old["recover_bar"](x, conn, E, A, d)):
            assert_allclose(a, b)

    def test_field_convergence(self):
        errors = []
        xi, weights = np.polynomial.legendre.leggauss(3)
        for ne in (2, 4, 8):
            x, conn = bar.uniform_bar_mesh(1, ne)
            ones = np.ones(ne)
            K, f = bar.assemble_bar(x, conn, ones, ones, ones, np.zeros(ne+1))
            d, _ = bar.solve_dirichlet(K, f, {0: 0})
            _, sig, _ = bar.recover_bar(x, conn, ones, ones, d)
            eu, es = 0., 0.
            for i, (a, b) in enumerate(conn):
                xx = (1-xi)*x[a]/2+(1+xi)*x[b]/2
                uh = (1-xi)*d[a]/2+(1+xi)*d[b]/2
                eu += np.dot(weights, (uh-(xx-xx**2/2))**2) * (x[b]-x[a])/2
                es += np.dot(weights, (sig[i]-(1-xx))**2) * (x[b]-x[a])/2
            errors.append(np.sqrt([eu, es]))
        assert_allclose(np.array(errors[:-1])/errors[1:], [[4, 2], [4, 2]], rtol=1e-10)

    def test_invalid_and_singular(self):
        for xa, xb, EA, p, ng in ((0, 0, 1, 1, 2), (0, np.inf, 1, 1, 2), (0, 1, np.nan, 1, 2), (0, 1, 1, np.inf, 2), (0, 1, 1, 1, 1.5)):
            with self.assertRaises(ValueError):
                bar.bar_element(xa, xb, lambda _: EA, lambda _: p, ng)
        with self.assertRaises(TypeError):
            bar.solve_dirichlet(np.eye(2), [0, 0], {0.5: 0})
        with self.assertRaises(np.linalg.LinAlgError):
            bar.solve_dirichlet(np.zeros((2, 2)), [0, 0], {0: 0})
        with self.assertRaises(TypeError):
            bar.assemble_bar([0, 1], [[0., 1.]], [1], [1], [0], [0, 0])


class PlasticityTests(unittest.TestCase):
    def test_reference_and_notebook_regression(self):
        old = notebook_functions(5)
        for E, nu, sy, H in ((210000, .3, 250, 1000), (100000, .2, 100, 0)):
            sig, ep, alpha = np.zeros((3, 3)), np.zeros((3, 3)), 0.
            for inc in (.003, -.001, -.005, .002):
                deps = np.diag([inc, 0, 0])
                expected = old["j2_update"](deps, sig, ep, alpha, E, nu, sy, H)
                actual = plasticity.j2_update(deps, sig, ep, alpha, E, nu, sy, H)
                for a, b in zip(actual, expected):
                    assert_allclose(a, b, atol=1e-10)
                sig, ep, alpha, da = actual
                self.assertGreaterEqual(da, 0)
                self.assertLessEqual(tensors.mises(sig)-sy-H*alpha, 1e-7)
        sig, ep, alpha, da = plasticity.j2_update(np.diag([.003, 0, 0]), np.zeros((3,3)), np.zeros((3,3)), 0, 210000, .3, 250, 1000)
        assert_allclose(np.diag(sig), [692.309516282, 441.345241859, 441.345241859], atol=1e-7, rtol=0)
        self.assertAlmostEqual(alpha, .000964274423, places=11)

    def test_elastic_hydrostatic_and_state_ownership(self):
        initial = np.zeros((3, 3))
        state = J2State(initial, initial)
        initial[0, 0] = 100
        self.assertEqual(state.stress[0, 0], 0)
        mat = J2Material(210000, .3, 250, 1000)
        new, da = mat.update(.001*np.eye(3), state)
        assert_allclose(new.stress, 525*np.eye(3))
        self.assertEqual(da, 0)
        self.assertFalse(np.shares_memory(state.plastic_strain, new.plastic_strain))
        with self.assertRaises(ValueError):
            state.stress[0, 0] = 10
        out = plasticity.j2_update(np.zeros((3,3)), np.zeros((3,3)), initial, 0, 210000, .3, 250, 0)
        self.assertFalse(np.shares_memory(out[1], initial))

    def test_path_reversal_rotation_and_history(self):
        mat = J2Material(210000, .3, 250, 1000)
        direction = np.diag([1., -.5, -.5])
        path = np.array([0., .004, -.004, 0.])[:, None, None] * direction
        stress, ep, alpha, da = plasticity.integrate_strain_path(path, mat, subdivisions=40)
        self.assertEqual(stress.shape, (121, 3, 3))
        self.assertTrue(np.all(np.diff(alpha) >= 0))
        assert_allclose(plasticity.accumulate_alpha(np.diff(ep, axis=0)), alpha, atol=1e-12)
        assert_allclose(np.trace(ep, axis1=1, axis2=2), 0, atol=1e-12)
        eps = np.concatenate([np.linspace(a, b, 40, endpoint=False) for a, b in zip(path[:-1], path[1:])] + [path[-1:]])
        assert_allclose(np.array([IsotropicElastic(mat.E, mat.nu).stress(e-e_p) for e,e_p in zip(eps, ep)]), stress, atol=1e-9)
        q = np.array([tensors.mises(s) for s in stress])
        assert_allclose((q-mat.sigma_y0-mat.H*alpha)[1:][da > 0], 0, atol=1e-9)
        self.assertGreater(np.linalg.norm(stress[-1]), 0)
        old = notebook_functions(5)["history_driver"](np.array([0., .004, -.004, 0.]), np.zeros((3,3)), np.zeros((3,3)), 0, 40)
        for a,b in zip((stress,ep,alpha,da),old):
            assert_allclose(a,b,atol=1e-9)
        Q, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(3,3)))
        rotated = plasticity.integrate_strain_path(Q.T @ path @ Q, mat, subdivisions=40)
        assert_allclose(rotated[0], Q.T @ stress @ Q, atol=1e-9)
        assert_allclose(rotated[2], alpha, atol=1e-12)

    def test_empty_and_invalid_history(self):
        assert_array_equal(plasticity.accumulate_alpha(np.empty((0,3,3))), [0])
        mat = J2Material(210000, .3, 250, 0)
        out = plasticity.integrate_strain_path(np.zeros((1,3,3)), mat)
        self.assertEqual(out[3].shape, (0,))
        for nd in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                plasticity.integrate_strain_path(np.zeros((2,3,3)), mat, subdivisions=nd)
        with self.assertRaises(ValueError):
            J2State(alpha=-1)
        with self.assertRaises(ValueError):
            J2Material(210000, .3, 250, -1)


class ReviewTests(unittest.TestCase):
    def history(self):
        t = np.array([0., .1, .4, 1.])
        # F=2t integrates exactly by trapezoids to t^2.
        return dict(t=t, drive=2*t, support=0*t, P=3+t*t, K=47+t, U=29+2*t,
                    D=t, H=.1*t, W=4.1*t)

    def test_nonuniform_time_offsets_scaling_and_input_preservation(self):
        data = self.history()
        original = {k: v.copy() for k,v in data.items()}
        out = review.audit_history(data, 10, 2)
        assert_allclose(out["r_E"], 0, atol=1e-14)
        assert_allclose(out["r_P"], 0, atol=1e-14)
        for k in data:
            assert_array_equal(data[k], original[k])
        data["P"][-1] += 1
        data["W"][-1] += 2
        out = review.audit_history(data, 10, 2)
        self.assertAlmostEqual(out["r_E"][-1], -.2)
        self.assertAlmostEqual(out["r_P"][-1], -.5)
        for key in ("r_E", "r_P"):
            assert_allclose(review.audit_history(data, 20, 4)[key], out[key]/2)
        old = notebook_functions(6)["audit_history"](data, 10, 2)
        for key in out:
            assert_allclose(out[key], old[key])

    def test_mesh_changes_and_bad_history(self):
        out = review.mesh_changes([1, 0, -1], .5)
        assert_allclose(out["absolute"], [1, 1])
        assert_allclose(out["relative"], [2, 1])
        data = self.history()
        data["t"][1] = 0
        with self.assertRaises(ValueError):
            review.audit_history(data, 10, 2)
        with self.assertRaises(ValueError):
            review.mesh_changes([1, 2], 0)


if __name__ == "__main__":
    unittest.main()
