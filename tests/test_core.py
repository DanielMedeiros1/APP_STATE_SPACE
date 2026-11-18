"""Testes do núcleo matemático (pytest). Execução: pytest -q"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
import pytest
import sympy as sp

from statespace import core
from statespace.core import DesignError
from statespace.parsing import InputError, parse_number, parse_poles, parse_matrix_text
from statespace.render import Fmt

F = Fmt(4, True)
A1 = sp.Matrix([[-4, 3, 0], [-1, -4, -1], [0, 0, -2]]); B1 = sp.Matrix([0, 0, 3]); C1 = sp.Matrix([[1, 0, 0]])
A2 = sp.Matrix([[-5, 1, 0], [2, 7, -1], [0, -3, 6]]); B2 = sp.Matrix([0, 0, 1]); C2 = sp.Matrix([[1, -2, 3]])
s = sp.Symbol("s")


def poly_of(poles, v=s):
    return sp.Poly(sp.expand(sp.prod([v - p for p in poles])), v).all_coeffs()


def test_parsing_fracoes_e_rejeicao():
    assert parse_number("-25/3") == sp.Rational(-25, 3)
    assert parse_number("0.25") == sp.Rational(1, 4)
    for bad in ["abc", "__import__('os')", "1/0", "2**3", ""]:
        with pytest.raises(InputError):
            parse_number(bad)
    with pytest.raises(InputError):
        parse_matrix_text("1 2; 3", 2, 2)


def test_polos_complexos_e_repetidos():
    p = parse_poles("-1±0.1j; -5; -5")
    assert len(p) == 4 and p.count(-5) == 2
    assert sp.im(p[0]) == sp.Rational(1, 10)


def test_polinomio_estabilidade():
    r = core.analyze_stability(A1, False, F)
    assert r.data["coeffs"] == [1, 10, 35, 38] and r.data["stability"]["status"] == "assintoticamente estável"
    r2 = core.analyze_stability(A2, False, F)
    assert r2.data["coeffs"] == [1, -8, -28, 207] and r2.data["stability"]["status"] == "instável"


def test_estabilidade_fronteira_jordan():
    assert core.stability(sp.Matrix([[0, 1], [0, 0]]), False)["status"] == "instável"
    assert core.stability(sp.zeros(2, 2), False)["status"] == "marginalmente estável"
    assert core.stability(sp.Matrix([[1, 1], [0, 1]]), True)["status"] == "instável"
    assert core.stability(sp.eye(2), True)["status"] == "marginalmente estável"
    assert core.stability(sp.Matrix([[0, 1], [-1, 0]]), False)["status"] == "marginalmente estável"


def test_controlabilidade_observabilidade():
    assert core.controllability(A1, B1, F).data["det"] == -81
    assert core.observability(A1, C1, F).data["det"] == -9
    assert core.controllability(A2, B2, F).data["det"] == -1
    nc = core.controllability(sp.diag(1, 2), sp.Matrix([1, 0]), F)
    assert not nc.data["controllable"] and nc.data["bad_modes"]
    assert core.design_feedback(sp.diag(1, 2), sp.Matrix([1, 0]), [-1, -2], F).data["K"] is None


@pytest.mark.parametrize("form", ["controlável", "observável"])
def test_forma_canonica_reconstroi_tf(form):
    r = core.tf_to_canonical([3, 4, -2], [1, 3, 7, 5], form, F)
    assert r.data["ok"]
    if form == "controlável":
        assert r.data["A"] == sp.Matrix([[0, 1, 0], [0, 0, 1], [-5, -7, -3]]) and list(r.data["C"]) == [-2, 4, 3]
    r2 = core.tf_to_canonical([3], [1, 2, -4], "observável", F)
    assert r2.data["A"] == sp.Matrix([[0, 4], [1, -2]]) and list(r2.data["B"]) == [3, 0] and list(r2.data["C"]) == [0, 1]


def test_tf_propria_e_impropria():
    r = core.tf_to_canonical([2, 1, 1], [1, 3, 2], "controlável", F)
    assert r.data["D"][0, 0] == 2 and r.data["ok"]
    with pytest.raises(DesignError):
        core.tf_to_canonical([1, 0, 0, 1], [1, 1], "controlável", F)


def test_ackermann_realimentacao_e_observador():
    poles = [-10, -2 + sp.I, -2 - sp.I]
    r = core.design_feedback(A1, B1, poles, F)
    assert r.data["K"] == sp.Matrix([[sp.Rational(8, 3), -sp.Rational(22, 3), -sp.Rational(4, 3)]])
    assert [sp.simplify(c) for c in core.char_coeffs(A1 + B1 * r.data["K"])] == poly_of(poles)
    r = core.design_feedback(A2, B2, poles, F)
    assert core.char_coeffs(A2 + B2 * r.data["K"]) == poly_of(poles)
    lo = core.design_observer(A1, C1, [-10, -1 + sp.I, -1 - sp.I], F)
    assert list(lo.data["L"]) == [2, -sp.Rational(25, 3), -sp.Rational(16, 3)]
    assert core.char_coeffs(A1 - lo.data["L"] * C1) == poly_of([-10, -1 + sp.I, -1 - sp.I])


def test_polos_repetidos_e_validacoes():
    t = core.tf_to_canonical([3, 4, -2], [1, 3, 7, 5], "controlável", F)
    k = core.design_feedback(t.data["A"], t.data["B"], [-4, -4, -5], F)
    assert list(k.data["K"]) == [-75, -49, -10]
    assert core.char_coeffs(t.data["A"] + t.data["B"] * k.data["K"]) == [1, 13, 56, 80]
    l = core.design_observer(t.data["A"], t.data["C"], [-1, -1, -2], F)
    assert core.char_coeffs(t.data["A"] - l.data["L"] * t.data["C"]) == [1, 4, 5, 2]
    with pytest.raises(DesignError):
        core.design_feedback(A1, B1, [-1, -2], F)
    with pytest.raises(DesignError):
        core.design_feedback(A1, B1, [-1, -2 + sp.I, -3], F)  # sem par conjugado


def test_seguidores():
    D0 = sp.Matrix([[0]])
    t = core.tf_to_canonical([3, 4, -2], [1, 3, 7, 5], "controlável", F)
    Aa_expected = sp.Matrix([[0, -2, 4, 3], [0, 0, 1, 0], [0, 0, 0, 1], [0, -5, -7, -3]])
    poles = parse_poles("-1±0.1j; -5; -12")
    r = core.design_tracker(t.data["A"], t.data["B"], t.data["C"], D0, poles, F)
    assert r.data["Aa"] == Aa_expected and r.data["ok"]
    assert core.controllability(Aa_expected, r.data["Ba"], F).data["det"] == -2
    assert core.char_coeffs(r.data["Acl"]) == poly_of(poles)
    assert abs(float(r.data["k1"]) - 30.3) < 1e-9
    r = core.design_tracker(sp.Matrix([[0, 1], [0, -3]]), sp.Matrix([0, 1]), sp.Matrix([[1, 2]]), D0, [-1, -1, -1], F)
    assert list(r.data["Ka"]) == [-1, -1, 0]
    r = core.design_tracker(sp.Matrix([[0, 4], [1, -2]]), sp.Matrix([3, 0]), sp.Matrix([[0, 1]]), D0, [-2, -2, -3], F)
    assert r.data["ok"] and core.controllability(r.data["Aa"], r.data["Ba"], F).data["det"] == 27
    # com D != 0
    r = core.design_tracker(sp.Matrix([[0, 1], [0, -3]]), sp.Matrix([0, 1]), sp.Matrix([[1, 2]]), sp.Matrix([[1]]), [-1, -2, -3], F)
    assert r.data["ok"]
    with pytest.raises(DesignError):
        core.design_tracker(A1, B1, C1, D0, [-1, -2, -3], F)  # precisa de n+1 polos


def test_seguidor_simulacao_converge_para_referencia():
    A = sp.Matrix([[0, 1], [0, -3]]); B = sp.Matrix([0, 1]); C = sp.Matrix([[1, 2]]); D = sp.Matrix([[0]])
    r = core.design_tracker(A, B, C, D, [-1, -1, -1], F)
    Acl = core.to_np(r.data["Acl"]); N = core.to_np(r.data["N"])
    Cz = np.array([[1, 2, 0]]); t, Z, Y = core.simulate_ct(Acl, N, np.zeros(3), 2.0, 20, 400, Cz, [0])
    assert abs(Y[-1, 0] - 2.0) < 1e-3


def test_discretizacao_A_singular():
    d = core.discretize(sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([0, 1]), sp.Integer(1), F)
    assert d.data["G"] == sp.Matrix([[1, 1], [0, 1]]) and d.data["H"] == sp.Matrix([sp.Rational(1, 2), 1])
    rendered = " ".join(value for step in d.steps for kind, value in step.items if kind == "tex")
    assert r"G(T)=e^{AT}" in rendered
    assert r"H(T)=\int_0^T e^{A\tau}B" in rendered
    assert r"A^{2}" in rendered and r"e^{A\tau}B" in rendered
    assert r"\right]_{0}^{1}" in rendered
    d2 = core.discretize(sp.Matrix([[0, 1], [0, -3]]), sp.Matrix([0, 1]), sp.Rational(1, 2), F)
    G = core.to_np(d2.data["G"])
    assert abs(G[1, 1] - np.exp(-1.5)) < 1e-10 and abs(G[0, 1] - (1 - np.exp(-1.5)) / 3) < 1e-10


def test_projeto_discreto_atividade3():
    G = sp.Matrix([[1, 1], [0, 1]]); H = sp.Matrix([sp.Rational(1, 2), 1]); C = sp.Matrix([[1, 1]])
    k = core.design_feedback(G, H, parse_poles("0.25±0.5j"), F, True)
    assert list(k.data["K"]) == [-sp.Rational(13, 16), -sp.Rational(35, 32)]
    l = core.design_observer(G, C, [-sp.Rational(1, 2), sp.Rational(1, 2)], F, True)
    assert list(l.data["L"]) == [sp.Rational(5, 4), sp.Rational(3, 4)]
    e = core.estimated_state_feedback(G, H, C, sp.Matrix([[0]]), k.data["K"], l.data["L"], F, True)
    assert e.data["separation_ok"]
    Gq1 = sp.Matrix([[1, sp.Rational(1, 10)], [sp.Rational(1, 2), sp.Rational(1, 10)]])
    assert core.stability(Gq1, True)["status"] == "instável"


def test_ss_para_tf_e_cancelamento():
    r = core.ss_to_tf(A1, B1, C1, sp.Matrix([[0]]), F)
    assert r.data["den"] == [1, 10, 35, 38] or not r.data["cancelled"]
