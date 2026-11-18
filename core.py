"""Núcleo matemático: análise, projeto e simulação em espaço de estados.

Convenções da disciplina:
  contínuo : x' = Ax + Bu,      y = Cx + Du
  discreto : x[k+1] = Gx[k] + Hu[k],  y[k] = Cx[k] + Du[k]
  realimentação u = Kx  ->  A_MF = A + BK  (G_MF = G + HK)
  Ackermann : K = -e_n^T U^{-1} q_c(A) ;  L = q_o(A) V^{-1} e_n   (autovalores de A - LC)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from scipy.linalg import expm, svdvals

from .render import Fmt


class DesignError(ValueError):
    """Erro matemático/de projeto com mensagem compreensível."""


@dataclass
class Step:
    title: str
    items: list


@dataclass
class Result:
    steps: list = field(default_factory=list)
    data: dict = field(default_factory=dict)
    conclusion: str = ""
    warnings: list = field(default_factory=list)


def TEX(s: str):
    return ("tex", s)


def TXT(s: str):
    return ("txt", s)


# ----------------------------------------------------------------- utilidades
def to_np(M) -> np.ndarray:
    M = sp.Matrix(M)
    a = np.array([[complex(sp.N(M[i, j])) for j in range(M.cols)] for i in range(M.rows)], dtype=complex)
    if np.allclose(a.imag, 0, atol=1e-14):
        return a.real.astype(float)
    return a


def from_np(a, digits: int = 12) -> sp.Matrix:
    """Converte um array real em matriz racional (12 algarismos significativos)."""
    a = np.asarray(a, dtype=float)
    return sp.Matrix(
        a.shape[0], a.shape[1],
        lambda i, j: sp.Rational(f"{(0.0 if abs(a[i, j]) < 1e-13 else a[i, j]):.{digits}g}"),
    )


def names(discrete: bool):
    return ("G", "H", "z") if discrete else ("A", "B", "s")


def powers(F: sp.Matrix, n: int):
    out = [sp.eye(F.rows)]
    for _ in range(n):
        out.append(F * out[-1])
    return out


def rank_info(M) -> dict:
    M = sp.Matrix(M)
    r_exact = M.rank()
    a = to_np(M).astype(complex)
    sv = svdvals(a)
    tol = max(a.shape) * np.finfo(float).eps * (sv.max() if sv.size else 0.0)
    r_num = int((sv > tol).sum())
    cond = float(sv.max() / sv.min()) if sv.size and sv.min() > 0 else float("inf")
    return dict(rank=r_exact, rank_num=r_num, tol=float(tol), cond=cond)


def _rank_warnings(info: dict, label: str) -> list:
    w = []
    if info["rank"] != info["rank_num"]:
        w.append(f"{label}: o posto exato ({info['rank']}) difere do posto numérico ({info['rank_num']}); "
                 "a matriz é numericamente mal condicionada.")
    elif info["cond"] > 1e8 and np.isfinite(info["cond"]):
        w.append(f"{label}: número de condição ≈ {info['cond']:.2e}; o resultado numérico pode ser sensível.")
    return w


def char_coeffs(F) -> list:
    x = sp.Symbol("x")
    return list(sp.Matrix(F).charpoly(x).all_coeffs())


def eig_info(F) -> list:
    """Autovalores com multiplicidade algébrica (exatos quando possível)."""
    F = sp.Matrix(F)
    n = F.rows
    x = sp.Symbol("x")
    groups = None
    try:
        r = sp.roots(sp.Poly(char_coeffs(F), x), cubics=False, quartics=False)
        if sum(r.values()) == n:
            groups = []
            for v, m in r.items():
                c = complex(sp.N(v, 30))
                if abs(c.imag) < 1e-12:
                    c = complex(c.real, 0.0)
                groups.append(dict(value=c, mult=int(m), exact=v))
    except Exception:
        groups = None
    if groups is None:
        lam = np.linalg.eigvals(to_np(F).astype(complex))
        used = [False] * len(lam)
        groups = []
        for i, l in enumerate(lam):
            if used[i]:
                continue
            idx = [j for j in range(len(lam)) if not used[j] and abs(lam[j] - l) < 1e-5 * max(1.0, abs(l))]
            for j in idx:
                used[j] = True
            val = complex(np.mean([lam[j] for j in idx]))
            if abs(val.imag) < 1e-10:
                val = complex(val.real, 0.0)
            groups.append(dict(value=val, mult=len(idx), exact=None))
    groups.sort(key=lambda g: (round(g["value"].real, 9), round(g["value"].imag, 9)))
    return groups


def flat_poles(info: list) -> list:
    out = []
    for g in info:
        out += [g["value"]] * g["mult"]
    return out


def _eig_tex(g: dict, fmt: Fmt) -> str:
    if fmt.exact and g["exact"] is not None:
        t = fmt.num(g["exact"])
        if len(t) <= 70:
            return t
    return Fmt(max(fmt.dec, 4), False).num(g["value"])


def _eig_list_tex(info: list, fmt: Fmt, var: str = "\\lambda") -> str:
    parts, k = [], 1
    for g in info:
        for _ in range(g["mult"]):
            parts.append(f"{var}_{{{k}}} = {_eig_tex(g, fmt)}")
            k += 1
    return r",\quad ".join(parts)


def _sI_minus(F, var: str, fmt: Fmt) -> str:
    F = sp.Matrix(F)
    rows = []
    for i in range(F.rows):
        cells = []
        for j in range(F.cols):
            cells.append(fmt.poly([1, -F[i, j]], var) if i == j else fmt.num(-F[i, j]))
        rows.append(" & ".join(cells))
    return r"\begin{bmatrix}" + r" \\ ".join(rows) + r"\end{bmatrix}"


# ------------------------------------------------- polinômio, polos, estabilidade
def stability(F, discrete: bool, info: list | None = None) -> dict:
    F = sp.Matrix(F)
    n = F.rows
    info = info or eig_info(F)
    Fn = to_np(F).astype(complex)
    tol = 1e-6
    unstable, boundary = False, []
    for g in info:
        l = g["value"]
        d = (abs(l) - 1.0) if discrete else l.real
        if d < -tol:
            continue
        if d > tol:
            unstable = True
            continue
        sv = svdvals(Fn - l * np.eye(n))
        r = int((sv > 1e-6 * max(1.0, sv.max())).sum())
        geo = n - r
        boundary.append(dict(value=l, mult=g["mult"], geo=geo, ok=geo >= g["mult"]))
    if unstable or any(not b["ok"] for b in boundary):
        status = "instável"
    elif boundary:
        status = "marginalmente estável"
    else:
        status = "assintoticamente estável"
    return dict(status=status, boundary=boundary, info=info)


def analyze_stability(F, discrete: bool, fmt: Fmt) -> Result:
    F = sp.Matrix(F)
    nf, _, var = names(discrete)
    coeffs = char_coeffs(F)
    info = eig_info(F)
    st = stability(F, discrete, info)
    steps = []
    steps.append(Step("Fórmula geral", [
        TEX(rf"p({var}) = \det({var}I - {nf})"),
        TXT("Os polos do sistema são as raízes de $p$ (autovalores de $" + nf + "$)."),
    ]))
    steps.append(Step("Substituição dos dados", [
        TEX(rf"{var}I-{nf} = " + _sI_minus(F, var, fmt)),
    ]))
    det_items = [TEX(rf"p({var}) = " + fmt.poly(coeffs, var))]
    if fmt.exact and info and all(g["exact"] is not None for g in info):
        facs = []
        for g in info:
            if g["value"].imag >= 0 or g["mult"] > 1:
                pass
        det_items.append(TXT("Polinômio característico monônico, obtido do determinante expandido."))
    steps.append(Step("Polinômio característico (expandido)", det_items))
    steps.append(Step("Polos (autovalores)", [
        TEX(_eig_list_tex(info, fmt)),
    ]))
    if discrete:
        crit = [TXT("**Critério discreto:** estabilidade assintótica se todos os polos têm $|z|<1$. "
                    "Módulos: " + ", ".join(f"{abs(g['value']):.4f}" for g in info))]
    else:
        crit = [TXT("**Critério contínuo:** estabilidade assintótica se todos os polos têm $\\mathrm{Re}(s)<0$. "
                    "Partes reais: " + ", ".join(f"{g['value'].real:.4f}" for g in info))]
    for b in st["boundary"]:
        crit.append(TXT(
            f"Polo na fronteira {Fmt(4, False).num(b['value'])}: multiplicidade algébrica {b['mult']}, "
            f"geométrica {b['geo']} → " + ("blocos de Jordan unitários (permitido na estabilidade marginal)."
                                         if b["ok"] else "bloco de Jordan não trivial (instável).")))
    steps.append(Step("Critério de estabilidade", crit))
    concl = f"O sistema é **{st['status']}**."
    return Result(steps, dict(coeffs=coeffs, info=info, poles=flat_poles(info), stability=st), concl)


# ------------------------------------------------------------- controlabilidade
def pbh_modes(F, B, mode: str = "ctrl") -> list:
    """Autovalores modos não controláveis ('ctrl') ou não observáveis ('obs') — teste PBH."""
    F, X = sp.Matrix(F), sp.Matrix(B)
    n = F.rows
    Fn, Xn = to_np(F).astype(complex), to_np(X).astype(complex)
    bad = []
    for g in eig_info(F):
        l = g["value"]
        M = np.hstack([l * np.eye(n) - Fn, Xn]) if mode == "ctrl" else np.vstack([l * np.eye(n) - Fn, Xn])
        sv = svdvals(M)
        if int((sv > 1e-9 * max(1.0, sv.max())).sum()) < n:
            bad.append(l)
    return bad


def controllability(F, B, fmt: Fmt, discrete: bool = False) -> Result:
    F, B = sp.Matrix(F), sp.Matrix(B)
    n = F.rows
    nf, nb, _ = names(discrete)
    nm = "W_c" if discrete else "U"
    steps = [Step("Fórmula geral", [
        TEX(rf"{nm} = \begin{{bmatrix}} {nb} & {nf}{nb} & {nf}^2{nb} & \cdots & {nf}^{{{n - 1}}}{nb} \end{{bmatrix}}"),
        TXT("O sistema é completamente controlável se $\\mathrm{posto}(" + nm + ")=n$."),
    ]), ]
    cols, items = [B], [TEX(rf"{nb} = " + fmt.mat(B))]
    for k in range(1, n):
        prev = nb if k == 1 else rf"{nf}^{{{k - 1}}}{nb}"
        cols.append(F * cols[-1])
        items.append(TEX(rf"{nf}^{{{k}}}{nb} = {nf}\,{prev} = " + fmt.mat(F) + fmt.mat(cols[-2]) + " = " + fmt.mat(cols[-1])))
    steps.append(Step("Produtos intermediários", items))
    U = sp.Matrix.hstack(*cols)
    ri = rank_info(U)
    its = [TEX(rf"{nm} = " + fmt.mat(U))]
    det = None
    if U.rows == U.cols:
        det = U.det()
        its.append(TEX(rf"\det({nm}) = " + fmt.num(det)))
    its.append(TEX(rf"\mathrm{{posto}}({nm}) = {ri['rank']}"))
    its.append(TXT(f"Tolerância numérica usada no posto (SVD): {ri['tol']:.2e}; número de condição: "
                   + ("infinito (singular)" if not np.isfinite(ri["cond"]) else f"{ri['cond']:.3e}") + "."))
    steps.append(Step("Matriz de controlabilidade, determinante e posto", its))
    ok = ri["rank"] == n
    warns = _rank_warnings(ri, nm)
    bad = [] if ok else pbh_modes(F, B, "ctrl")
    if ok:
        concl = f"posto({nm}) = n = {n}: o sistema é **completamente controlável**; os polos de malha fechada podem ser alocados arbitrariamente (pares complexos conjugados)."
    else:
        mt = ", ".join(Fmt(4, False).num(b) for b in bad) or "—"
        concl = (f"posto({nm}) = {ri['rank']} < n = {n}: o sistema **não é completamente controlável**. "
                 f"Modos não controláveis (teste PBH): {mt}. A alocação arbitrária de todos os polos não é possível; "
                 "somente os modos controláveis podem ser movidos, e a estabilização exige que os modos "
                 "não controláveis já sejam estáveis (estabilizabilidade).")
    return Result(steps, dict(U=U, rank=ri["rank"], det=det, controllable=ok, info=ri, bad_modes=bad), concl, warns)


def observability(F, C, fmt: Fmt, discrete: bool = False) -> Result:
    F, C = sp.Matrix(F), sp.Matrix(C)
    n = F.rows
    nf, _, _ = names(discrete)
    nm = "W_o" if discrete else "V"
    col = [r"C", r"C" + nf] + [rf"C{nf}^{{{k}}}" for k in range(2, n)]
    steps = [Step("Fórmula geral", [
        TEX(rf"{nm} = \begin{{bmatrix}} " + r" \\ ".join(col[:n]) + r" \end{bmatrix}"),
        TXT("O sistema é completamente observável se $\\mathrm{posto}(" + nm + ")=n$."),
    ])]
    rows, items = [C], [TEX(r"C = " + fmt.mat(C))]
    for k in range(1, n):
        prev = "C" if k == 1 else rf"C{nf}^{{{k - 1}}}"
        rows.append(rows[-1] * F)
        items.append(TEX(rf"C{nf}^{{{k}}} = {prev}\,{nf} = " + fmt.mat(rows[-2]) + fmt.mat(F) + " = " + fmt.mat(rows[-1])))
    steps.append(Step("Produtos intermediários", items))
    V = sp.Matrix.vstack(*rows)
    ri = rank_info(V)
    its = [TEX(rf"{nm} = " + fmt.mat(V))]
    det = None
    if V.rows == V.cols:
        det = V.det()
        its.append(TEX(rf"\det({nm}) = " + fmt.num(det)))
    its.append(TEX(rf"\mathrm{{posto}}({nm}) = {ri['rank']}"))
    its.append(TXT(f"Tolerância numérica usada no posto (SVD): {ri['tol']:.2e}; número de condição: "
                   + ("infinito (singular)" if not np.isfinite(ri["cond"]) else f"{ri['cond']:.3e}") + "."))
    steps.append(Step("Matriz de observabilidade, determinante e posto", its))
    ok = ri["rank"] == n
    bad = [] if ok else pbh_modes(F, C, "obs")
    if ok:
        concl = f"posto({nm}) = n = {n}: o sistema é **completamente observável**; os polos do observador podem ser alocados arbitrariamente."
    else:
        mt = ", ".join(Fmt(4, False).num(b) for b in bad) or "—"
        concl = (f"posto({nm}) = {ri['rank']} < n = {n}: o sistema **não é completamente observável**. "
                 f"Modos não observáveis (PBH): {mt}. Não é possível alocar arbitrariamente todos os polos do observador; "
                 "a detectabilidade exige que os modos não observáveis sejam estáveis.")
    return Result(steps, dict(V=V, rank=ri["rank"], det=det, observable=ok, info=ri, bad_modes=bad),
                  concl, _rank_warnings(ri, nm))


# --------------------------------------------------------- polinômio desejado
def desired_poly(poles, var: str, fmt: Fmt):
    n = len(poles)
    v = sp.Symbol(var)
    expr = sp.expand(sp.prod([v - p for p in poles]))
    coeffs = [sp.simplify(c) for c in sp.Poly(expr, v).all_coeffs()]
    if any(not sp.im(c) == 0 for c in coeffs):
        raise DesignError("Os polos complexos devem aparecer em pares conjugados (ganhos reais exigem polinômio real).")
    items = [
        TXT("Polos desejados: " + ", ".join(Fmt(6, False).num(p) if not fmt.exact else fmt.num(p) for p in poles) + "."),
        TEX(rf"q({var}) = " + fmt.factors(poles, var)),
        TEX(rf"q({var}) = " + fmt.poly(coeffs, var)),
    ]
    return coeffs, items


def matrix_poly(coeffs, F):
    n = len(coeffs) - 1
    P = powers(F, n)
    Q = sp.zeros(F.rows, F.rows)
    for i, c in enumerate(coeffs):
        Q += c * P[n - i]
    return Q, P


def _check_poly(Acl, coeffs, fmt, var, label):
    got = char_coeffs(Acl)
    same = all(sp.simplify(a - b) == 0 for a, b in zip(got, coeffs))
    info = eig_info(Acl)
    items = [TEX(rf"\det({var}I - {label}) = " + fmt.poly(got, var)),
             TEX(_eig_list_tex(info, fmt)),
             TXT("✅ O polinômio obtido coincide com o polinômio desejado." if same
                 else "⚠️ O polinômio obtido NÃO coincide com o desejado.")]
    return same, info, items


def design_feedback(F, B, poles, fmt: Fmt, discrete: bool = False) -> Result:
    """K = -e_n^T U^{-1} q_c(F)  (convenção u = Kx, F_MF = F + BK)."""
    F, B = sp.Matrix(F), sp.Matrix(B)
    n = F.rows
    nf, nb, var = names(discrete)
    nm = "W_c" if discrete else "U"
    if B.cols != 1:
        raise DesignError("O projeto por Ackermann requer sistema com uma única entrada (SISO).")
    if len(poles) != n:
        raise DesignError(f"São necessários exatamente {n} polos desejados; foram informados {len(poles)}.")
    steps = [Step("Fórmula geral (Ackermann, convenção $u=Kx$)", [
        TEX(rf"K = -\begin{{bmatrix}}0&\cdots&0&1\end{{bmatrix}}\,{nm}^{{-1}}\,q_c({nf})"),
        TEX(rf"{nf}_{{MF}} = {nf} + {nb}K"),
        TXT("Atenção: esta é a convenção $u=Kx$ (sinal já incluído em $K$); não use diretamente ganhos de rotinas que retornam $u=-Kx$."),
    ])]
    cr = controllability(F, B, fmt, discrete)
    steps += [Step("Controlabilidade (pré-requisito) — " + s.title, s.items) for s in cr.steps[1:]]
    if not cr.data["controllable"]:
        steps.append(Step("Conclusão da análise", [TXT(cr.conclusion)]))
        return Result(steps, dict(K=None, controllable=False), cr.conclusion, cr.warnings)
    coeffs, items = desired_poly(poles, var, fmt)
    steps.append(Step("Polinômio desejado", items))
    Q, P = matrix_poly(coeffs, F)
    its = [TEX(rf"q_c({nf}) = " + " + ".join(
        ([rf"{nf}^{{{n - i}}}" if n - i > 1 else (nf if n - i == 1 else "I")] if c == 1 and i == 0 else
         [rf"({fmt.num(c)})" + (rf"{nf}^{{{n - i}}}" if n - i > 1 else (nf if n - i == 1 else "I"))])[0]
        for i, c in enumerate(coeffs)))]
    for k in range(2, n + 1):
        its.append(TEX(rf"{nf}^{{{k}}} = " + fmt.mat(P[k])))
    its.append(TEX(rf"q_c({nf}) = " + fmt.mat(Q)))
    steps.append(Step(f"Substituição de ${var}$ por ${nf}$ e avaliação de $q_c({nf})$ (Cayley–Hamilton)", its))
    U = cr.data["U"]
    Uinv = U.inv()
    en = sp.Matrix([[0] * (n - 1) + [1]])
    row = en * Uinv
    K = -(row * Q)
    steps.append(Step(f"Inversa de ${nm}$", [TEX(rf"{nm}^{{-1}} = " + fmt.mat(Uinv))]))
    steps.append(Step("Cálculo do ganho $K$", [
        TEX(r"e_n^T " + nm + r"^{-1} = " + fmt.mat(row)),
        TEX(r"K = -(" + fmt.mat(row) + fmt.mat(Q) + r") = " + fmt.mat(K)),
    ]))
    Acl = F + B * K
    steps.append(Step(f"Matriz de malha fechada ${nf}+{nb}K$", [
        TEX(rf"{nf}+{nb}K = " + fmt.mat(F) + " + " + fmt.mat(B) + fmt.mat(K) + " = " + fmt.mat(Acl))]))
    same, info, its = _check_poly(Acl, coeffs, fmt, var, rf"({nf}+{nb}K)")
    steps.append(Step("Conferência dos polos obtidos", its))
    return Result(steps, dict(K=K, Acl=Acl, U=U, q=Q, coeffs=coeffs, poles_cl=flat_poles(info), ok=same, controllable=True),
                  f"K = {fmt.mat(K)} (convenção u = Kx); polos de malha fechada conferidos: " + ("✅ coincidem." if same else "⚠️ divergem."),
                  cr.warnings)


def design_observer(F, C, poles, fmt: Fmt, discrete: bool = False) -> Result:
    """L = q_o(F) V^{-1} e_n  (autovalores de F - LC)."""
    F, C = sp.Matrix(F), sp.Matrix(C)
    n = F.rows
    nf, _, var = names(discrete)
    nm = "W_o" if discrete else "V"
    if C.rows != 1:
        raise DesignError("O observador por Ackermann requer uma única saída (SISO).")
    if len(poles) != n:
        raise DesignError(f"São necessários exatamente {n} polos desejados para o observador; foram informados {len(poles)}.")
    steps = [Step("Fórmula geral (Ackermann para o observador)", [
        TEX(rf"L = q_o({nf})\,{nm}^{{-1}}\begin{{bmatrix}}0\\ \vdots \\ 1\end{{bmatrix}}"),
        TEX(r"\dot{\hat x} = A\hat x + Bu + L(y-\hat y),\quad \hat y = C\hat x + Du" if not discrete else
            r"\hat x_{k+1} = G\hat x_k + Hu_k + L(y_k-\hat y_k),\quad \hat y_k = C\hat x_k + Du_k"),
        TEX(rf"\dot{{\tilde x}} = ({nf}-LC)\tilde x" if not discrete else rf"\tilde x_{{k+1}} = ({nf}-LC)\tilde x_k"),
    ])]
    orr = observability(F, C, fmt, discrete)
    steps += [Step("Observabilidade (pré-requisito) — " + s.title, s.items) for s in orr.steps[1:]]
    if not orr.data["observable"]:
        steps.append(Step("Conclusão da análise", [TXT(orr.conclusion)]))
        return Result(steps, dict(L=None, observable=False), orr.conclusion, orr.warnings)
    coeffs, items = desired_poly(poles, var, fmt)
    steps.append(Step("Polinômio desejado do observador", items))
    Q, P = matrix_poly(coeffs, F)
    its = [TEX(rf"{nf}^{{{k}}} = " + fmt.mat(P[k])) for k in range(2, n + 1)]
    its.append(TEX(rf"q_o({nf}) = " + fmt.mat(Q)))
    steps.append(Step(f"Substituição de ${var}$ por ${nf}$ e avaliação de $q_o({nf})$", its))
    V = orr.data["V"]
    Vinv = V.inv()
    en = sp.Matrix([0] * (n - 1) + [1])
    col = Vinv * en
    L = Q * col
    steps.append(Step(f"Inversa de ${nm}$", [TEX(rf"{nm}^{{-1}} = " + fmt.mat(Vinv))]))
    steps.append(Step("Cálculo do ganho $L$", [
        TEX(rf"{nm}^{{-1}}e_n = " + fmt.mat(col)),
        TEX(r"L = " + fmt.mat(Q) + fmt.mat(col) + " = " + fmt.mat(L))]))
    Ao = F - L * C
    steps.append(Step(f"Matriz do erro ${nf}-LC$", [TEX(rf"{nf}-LC = " + fmt.mat(F) + " - " + fmt.mat(L) + fmt.mat(C) + " = " + fmt.mat(Ao))]))
    same, info, its = _check_poly(Ao, coeffs, fmt, var, rf"({nf}-LC)")
    steps.append(Step("Conferência dos polos do observador", its))
    return Result(steps, dict(L=L, Ao=Ao, V=V, q=Q, coeffs=coeffs, poles_obs=flat_poles(info), ok=same, observable=True),
                  f"L = {fmt.mat(L)}; autovalores de {nf}−LC conferidos: " + ("✅ coincidem." if same else "⚠️ divergem."),
                  orr.warnings)


def estimated_state_feedback(F, B, C, D, K, L, fmt: Fmt, discrete: bool = False) -> Result:
    """Interligação planta + observador + u = K x̂ e verificação do princípio da separação."""
    F, B, C, D, K, L = [sp.Matrix(M) for M in (F, B, C, D, K, L)]
    n = F.rows
    nf, nb, _ = names(discrete)
    Z = sp.zeros(n, n)
    Acomb = sp.Matrix.vstack(sp.Matrix.hstack(F, B * K), sp.Matrix.hstack(L * C, F + B * K - L * C))
    Asep = sp.Matrix.vstack(sp.Matrix.hstack(F + B * K, -B * K), sp.Matrix.hstack(Z, F - L * C))
    if discrete:
        eqs = [TEX(r"u_k = K\hat x_k"), TEX(r"x_{k+1} = Gx_k + HK\hat x_k"),
               TEX(r"\hat x_{k+1} = G\hat x_k + HK\hat x_k + LC(x_k-\hat x_k)")]
    else:
        eqs = [TEX(r"u = K\hat x"), TEX(r"\dot x = Ax + BK\hat x"),
               TEX(r"\dot{\hat x} = A\hat x + BK\hat x + L\,C(x-\hat x)")]
    steps = [Step("Equações da interligação", eqs + [TXT(
        r"Como $\hat y = C\hat x + Du$ e $y = Cx + Du$, o termo $D u$ se cancela em $y-\hat y = C(x-\hat x)$.")])]
    steps.append(Step("Sistema completo nos estados $[x;\\hat x]$", [
        TEX(r"\begin{bmatrix}x\\ \hat x\end{bmatrix}' = " + fmt.mat(Acomb) + r"\begin{bmatrix}x\\ \hat x\end{bmatrix}")]))
    steps.append(Step("Mudança de variável $\\tilde x = x-\\hat x$ (princípio da separação)", [
        TEX(rf"\begin{{bmatrix}}x\\ \tilde x\end{{bmatrix}}' = \begin{{bmatrix}} {nf}+{nb}K & -{nb}K \\ 0 & {nf}-LC \end{{bmatrix}}\begin{{bmatrix}}x\\ \tilde x\end{{bmatrix}}"),
        TEX(fmt.mat(Asep)),
        TXT("A matriz é bloco-triangular: os autovalores são a união dos autovalores de $" + nf + "+" + nb + "K$ (regulador) e de $" + nf + "-LC$ (observador)."),
    ]))
    if discrete:
        for step in steps:
            step.items = [(kind, value.replace(r"\begin{bmatrix}x\\ \hat x\end{bmatrix}'", r"\begin{bmatrix}x_{k+1}\\ \hat x_{k+1}\end{bmatrix}").replace(r"\begin{bmatrix}x\\ \tilde x\end{bmatrix}'", r"\begin{bmatrix}x_{k+1}\\ \tilde x_{k+1}\end{bmatrix}")) for kind,value in step.items]
    ev = np.sort_complex(np.linalg.eigvals(to_np(Acomb).astype(complex)))
    ex = np.sort_complex(np.array(flat_poles(eig_info(F + B * K)) + flat_poles(eig_info(F - L * C)), dtype=complex))
    variable=sp.Symbol("s")
    ok = sp.expand(Acomb.charpoly(variable).as_expr()-(F+B*K).charpoly(variable).as_expr()*(F-L*C).charpoly(variable).as_expr())==0
    steps.append(Step("Verificação numérica", [
        TXT("Autovalores do sistema completo: " + ", ".join(Fmt(4, False).num(e) for e in ev)),
        TXT("✅ Polinômio conferido exatamente: coincide com o produto dos polinômios do regulador e do observador." if ok else "⚠️ Divergência numérica; verifique mal condicionamento."),
    ]))
    return Result(steps, dict(Acomb=Acomb, Asep=Asep, eig=list(ev), separation_ok=ok),
                  "Princípio da separação verificado." if ok else "Não foi possível confirmar a separação numericamente.")


# ------------------------------------------------------------------ seguidor
def design_tracker(A, B, C, D, poles, fmt: Fmt) -> Result:
    """Seguidor de referência (degrau) com ação integral: u = k1 ξ + K2 x, ξ' = y - r."""
    A, B, C, D = [sp.Matrix(M) for M in (A, B, C, D)]
    n = A.rows
    if B.cols != 1 or C.rows != 1:
        raise DesignError("O seguidor requer sistema SISO.")
    if len(poles) != n + 1:
        raise DesignError(f"O sistema aumentado tem ordem {n + 1}; informe exatamente {n + 1} polos (foram {len(poles)}).")
    d = D[0, 0]
    Aa = sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(1, 1), C), sp.Matrix.hstack(sp.zeros(n, 1), A))
    Ba = sp.Matrix.vstack(sp.Matrix([[d]]), B)
    steps = [Step("Construção do sistema aumentado", [
        TXT("Erro $e=y-r$ e estado integral $\\xi=\\int (y-r)\\,dt$. Para o projeto considera-se $r=0$ e estado aumentado $x_a=[\\xi;\\,x]$."),
        TEX(r"A_a = \begin{bmatrix} 0 & C \\ 0 & A\end{bmatrix} = " + fmt.mat(Aa)),
        TEX(r"B_a = \begin{bmatrix} D \\ B\end{bmatrix} = " + fmt.mat(Ba)),
        TXT("Com $D=0$ (caso usual) $B_a=[0;\\,B]$. Se $D\\neq0$, $\\dot\\xi=Cx+Du$ e a primeira entrada de $B_a$ passa a ser $D$ (adaptação aplicada automaticamente)."
            + ("" if d == 0 else " **Aqui $D\\neq0$.**")),
    ])]
    des = design_feedback(Aa, Ba, poles, fmt, False)
    steps += des.steps
    if des.data.get("K") is None:
        return Result(steps, dict(Ka=None, controllable=False), des.conclusion, des.warnings)
    Ka = des.data["K"]
    k1, K2 = Ka[0, 0], Ka[0, 1:]
    steps.append(Step("Separação do ganho e implementação física", [
        TEX(r"K_a = [k_1\;\;K_2] = " + fmt.mat(Ka)),
        TEX(r"k_1 = " + fmt.num(k1) + r",\qquad K_2 = " + fmt.mat(K2)),
        TEX(r"\xi = \int (y-r)\,dt,\qquad u = k_1\xi + K_2 x"),
        TXT("O projeto usa $[\\xi;\\,x]$ (integrador primeiro); a simulação abaixo usa a ordem $[x;\\,\\xi]$, que apenas permuta linhas e colunas e preserva os autovalores."),
    ]))
    Acl = sp.Matrix.vstack(sp.Matrix.hstack(A + B * K2, B * k1), sp.Matrix.hstack(C + D * K2, sp.Matrix([[d * k1]])))
    N = sp.Matrix.vstack(sp.zeros(n, 1), sp.Matrix([[-1]]))
    steps.append(Step("Sistema para simulação", [
        TEX(r"\begin{bmatrix}\dot x\\ \dot\xi\end{bmatrix} = " + (
            r"\begin{bmatrix} A+BK_2 & Bk_1 \\ C & 0\end{bmatrix}" if d == 0 else
            r"\begin{bmatrix} A+BK_2 & Bk_1 \\ C+DK_2 & Dk_1\end{bmatrix}") +
            r"\begin{bmatrix}x\\ \xi\end{bmatrix} + \begin{bmatrix}0\\-1\end{bmatrix} r"),
        TEX(r"= " + fmt.mat(Acl) + r"\begin{bmatrix}x\\ \xi\end{bmatrix} + " + fmt.mat(N) + "r"),
        TXT("No regime permanente $\\dot\\xi=0\\Rightarrow y=r$ (erro nulo para degrau)."),
    ]))
    ev = np.sort_complex(np.linalg.eigvals(to_np(Acl).astype(complex)))
    ok_ev = all(sp.simplify(a-b)==0 for a,b in zip(char_coeffs(Acl), des.data["coeffs"]))
    steps.append(Step("Conferência", [
        TXT("Autovalores do sistema de simulação: " + ", ".join(Fmt(4, False).num(e) for e in ev)),
        TXT("✅ Coincidem com os polos desejados." if ok_ev else "⚠️ Divergência com os polos desejados."),
    ]))
    return Result(steps, dict(Ka=Ka, k1=k1, K2=K2, Aa=Aa, Ba=Ba, Acl=Acl, N=N, ok=ok_ev and des.data["ok"], controllable=True),
                  f"K_a = {fmt.mat(Ka)}; k₁ = {fmt.num(k1)}.", des.warnings)


# ------------------------------------------------------ função de transferência
def tf_to_canonical(num, den, form: str, fmt: Fmt, var: str = "s") -> Result:
    """Formas canônicas controlável/observável conforme as tabelas da disciplina."""
    num, den = [sp.sympify(c) for c in num], [sp.sympify(c) for c in den]
    while len(num) > 1 and num[0] == 0:
        num = num[1:]
    while den and den[0] == 0:
        den = den[1:]
    if not den:
        raise DesignError("O denominador não pode ser nulo.")
    n = len(den) - 1
    if n < 1:
        raise DesignError("O denominador deve ter grau ≥ 1.")
    if len(num) - 1 > n:
        raise DesignError("Função de transferência imprópria (grau do numerador > grau do denominador): "
                          "não é suportada pela realização em espaço de estados adotada.")
    a0 = den[0]
    den_m = [c / a0 for c in den]
    num_m = [c / a0 for c in num]
    num_p = [sp.Integer(0)] * (n + 1 - len(num_m)) + num_m
    D = num_p[0]
    beta = [num_p[i] - D * den_m[i] for i in range(1, n + 1)]   # beta[0]=β1 (s^{n-1}) ... beta[n-1]=βn
    alpha = den_m[1:]                                          # alpha[0]=α1 ... alpha[n-1]=αn
    steps = [Step("Normalização do denominador", [
        TEX(rf"G({var}) = \frac{{" + fmt.poly(num, var) + r"}{" + fmt.poly(den, var) + "}"),
        TEX(rf"G({var}) = \frac{{" + fmt.poly(num_m, var) + r"}{" + fmt.poly(den_m, var) + "}"
            + rf" \quad(\text{{denominador monônico: }}{var}^n+\alpha_1{var}^{{n-1}}+\dots+\alpha_n)"),
        TEX(r"\alpha = (" + ", ".join(fmt.num(a) for a in alpha) + r")"),
    ])]
    if D != 0:
        steps.append(Step("Termo direto $D$", [
            TXT("Função própria (grau do numerador = grau do denominador): divide-se para separar $D$."),
            TEX(r"D = " + fmt.num(D)),
            TEX(r"N_{sp}(" + var + r") = N(" + var + r") - D\,\mathrm{den}(" + var + r") = " + fmt.poly(beta, var)),
        ]))
    steps.append(Step("Coeficientes $\\beta$", [TEX(r"\beta = (" + ", ".join(fmt.num(b) for b in beta) + r")\;\;(\beta_1\text{ de }" + var + r"^{n-1}\text{ até }\beta_n\text{ constante})")]))
    A = sp.zeros(n, n)
    if form == "controlável":
        for i in range(n - 1):
            A[i, i + 1] = 1
        for j in range(n):
            A[n - 1, j] = -den_m[n - j]
        B = sp.zeros(n, 1)
        B[n - 1, 0] = 1
        C = sp.Matrix([[beta[n - 1 - j] for j in range(n)]])
        form_tex = [TEX(r"A=\begin{bmatrix}0&1&\cdots&0\\ \vdots&&\ddots&\vdots\\ 0&0&\cdots&1\\ -\alpha_n&-\alpha_{n-1}&\cdots&-\alpha_1\end{bmatrix},\;B=\begin{bmatrix}0\\ \vdots\\0\\1\end{bmatrix},\;C=[\beta_n\;\cdots\;\beta_1]")]
    else:
        for i in range(n - 1):
            A[i + 1, i] = 1
        for i in range(n):
            A[i, n - 1] = -den_m[n - i]
        B = sp.Matrix([beta[n - 1 - i] for i in range(n)])
        C = sp.zeros(1, n)
        C[0, n - 1] = 1
        form_tex = [TEX(r"A=\begin{bmatrix}0&\cdots&0&-\alpha_n\\ 1&\cdots&0&-\alpha_{n-1}\\ \vdots&\ddots&&\vdots\\ 0&\cdots&1&-\alpha_1\end{bmatrix},\;B=\begin{bmatrix}\beta_n\\ \vdots\\ \beta_2\\ \beta_1\end{bmatrix},\;C=[0\;\cdots\;0\;1]")]
    steps.append(Step(f"Forma canônica {form} (organização dos estados da disciplina)", form_tex + [
        TEX(r"A = " + fmt.mat(A)), TEX(r"B = " + fmt.mat(B)), TEX(r"C = " + fmt.mat(C)), TEX(r"D = " + fmt.num(D))]))
    v = sp.Symbol("v")
    Gr = sp.cancel(sp.together((C * (v * sp.eye(n) - A).adjugate() * B)[0, 0] / (v * sp.eye(n) - A).det() + D))
    Go = sp.cancel(sp.Poly(num_m, v).as_expr() / sp.Poly(den_m, v).as_expr())
    ok = sp.simplify(Gr - Go) == 0
    steps.append(Step("Conferência: $C(sI-A)^{-1}B+D$", [
        TEX(r"C(" + var + r"I-A)^{-1}B + D = " + sp.latex(Gr.subs(v, sp.Symbol(var)))),
        TXT("✅ Reconstrói exatamente a função de transferência dada." if ok else "⚠️ A reconstrução não coincide."),
    ]))
    return Result(steps, dict(A=A, B=B, C=C, D=sp.Matrix([[D]]), n=n, ok=ok), f"Realização {form} de ordem {n}.")


def ss_to_tf(A, B, C, D, fmt: Fmt, var: str = "s") -> Result:
    A, B, C, D = [sp.Matrix(M) for M in (A, B, C, D)]
    if B.cols != 1 or C.rows != 1:
        raise DesignError("A conversão para função de transferência é implementada para sistemas SISO.")
    n = A.rows
    v = sp.Symbol("v")
    den = sp.Poly(sp.Matrix(A).charpoly(v).as_expr(), v)
    adj = (v * sp.eye(n) - A).adjugate()
    numr = sp.expand((C * adj * B)[0, 0] + D[0, 0] * den.as_expr())
    g = sp.cancel(numr / den.as_expr())
    nn, dd = sp.fraction(g)
    cn, cd = sp.Poly(nn, v).all_coeffs(), sp.Poly(dd, v).all_coeffs()
    lead = cd[0]
    cn, cd = [c / lead for c in cn], [c / lead for c in cd]
    cancelled = len(cd) - 1 < n
    steps = [Step("Fórmula", [TEX(rf"G({var}) = C({var}I-A)^{{-1}}B + D = \frac{{C\,\mathrm{{adj}}({var}I-A)B}}{{\det({var}I-A)}}+D")]),
             Step("Resultado", [TEX(rf"G({var}) = \frac{{" + fmt.poly(cn, var) + "}{" + fmt.poly(cd, var) + "}")])]
    if cancelled:
        steps[-1].items.append(TXT("Houve cancelamento polo-zero: a realização não é mínima (perda de controlabilidade e/ou observabilidade)."))
    return Result(steps, dict(num=cn, den=cd, cancelled=cancelled), "")


# --------------------------------------------------------------- discretização
def discretize(A, B, T, fmt: Fmt, symbolic: bool = True) -> Result:
    """Discretização exata com entrada constante entre amostras (ZOH), sem usar A^{-1}."""
    A, B = sp.Matrix(A), sp.Matrix(B)
    T = sp.nsimplify(T) if not isinstance(T, sp.Basic) else T
    if T <= 0:
        raise DesignError("O período de amostragem T deve ser positivo.")
    n, m = A.rows, B.cols
    steps = [Step("Equações de discretização (conforme a atividade)", [
        TEX(r"G(T)=e^{AT}"),
        TEX(r"H(T)=\int_0^T e^{A\tau}B\,d\tau"),
        TXT("A matriz G é obtida pela exponencial de A; a matriz H é obtida integrando e^{Aτ}B no intervalo de 0 a T."),
    ])]
    P = powers(A, n)
    nilpotent_index = next((k for k in range(1, n + 1) if A ** k == sp.zeros(n, n)), None)
    if nilpotent_index is not None:
        # Para matrizes nilpotentes, a série da exponencial termina. A apresentação
        # desenvolve separadamente G(T) e a integral H(T), como no roteiro do professor.
        tau = sp.Symbol("tau", real=True)
        E_tau = sum((A ** k * tau ** k / sp.factorial(k) for k in range(nilpotent_index)),
                    sp.zeros(n, n))
        G = E_tau.subs(tau, T).applyfunc(sp.simplify)
        integrand = (E_tau * B).applyfunc(sp.simplify)
        primitive = integrand.applyfunc(lambda expr: sp.integrate(expr, tau))
        H = integrand.applyfunc(lambda expr: sp.simplify(sp.integrate(expr, (tau, 0, T))))

        power_items = [TEX(r"A=" + fmt.mat(A))]
        for k in range(2, nilpotent_index + 1):
            power_items.append(TEX(rf"A^{{{k}}}=" + fmt.mat(A ** k)))
        series_g = [r"I"] + [
            (r"AT" if k == 1 else rf"\frac{{A^{{{k}}}T^{{{k}}}}}{{{k}!}}")
            for k in range(1, nilpotent_index)
        ]
        series_tau = [r"I"] + [
            (r"A\tau" if k == 1 else rf"\frac{{A^{{{k}}}\tau^{{{k}}}}}{{{k}!}}")
            for k in range(1, nilpotent_index)
        ]
        steps.append(Step("Potências de A e cálculo de G(T)", power_items + [
            TEX(r"e^{AT}=" + "+".join(series_g)),
            TEX(r"G(T)=e^{AT}=" + "+".join(series_g) + r"=" + fmt.mat(G)),
            TEX(r"G(T)=I+AT=" + fmt.mat(sp.eye(n)) + "+" + fmt.num(T) + fmt.mat(A) + "=" + fmt.mat(G))
                if nilpotent_index == 2 else TEX(r"G(T)=" + fmt.mat(G)),
        ]))
        steps.append(Step("Cálculo de H(T) pela integral", [
            TEX(r"e^{A\tau}=" + "+".join(series_tau) + r"=" + fmt.mat(E_tau)),
            TEX(r"e^{A\tau}B=" + fmt.mat(E_tau) + fmt.mat(B) + "=" + fmt.mat(integrand)),
            TEX(r"H(T)=\int_0^T e^{A\tau}B\,d\tau=\int_0^T " + fmt.mat(integrand) + r"\,d\tau"),
            TEX(r"=\left[" + fmt.mat(primitive) + r"\right]_{0}^{" + fmt.num(T) + r"}=" + fmt.mat(H)),
        ]))
        exact = True
    elif symbolic and n <= 3 and sum(sp.roots(A.charpoly().as_expr(), cubics=False, quartics=False).values()) == n:
        tau = sp.Symbol("tau", real=True)
        Et = (A*tau).exp().applyfunc(sp.simplify)
        G = Et.subs(tau,T).applyfunc(sp.simplify)
        H = (Et*B).applyfunc(lambda e: sp.simplify(sp.integrate(e,(tau,0,T))))
        steps.append(Step("Exponencial e integral simbólicas", [
            TEX(r"e^{A\tau} = " + fmt.mat(Et)),
            TEX(r"G(T) = e^{AT} = " + fmt.mat(G)),
            TEX(r"H(T)=\int_0^T e^{A\tau}B\,d\tau = " + fmt.mat(H)),
        ]))
        exact = True
    else:
        big = sp.zeros(n + m, n + m)
        big[:n, :n] = A * T
        big[:n, n:] = B * T
        E = expm(to_np(big).astype(float))
        G, H = from_np(E[:n, :n]), from_np(E[:n, n:])
        dfmt = Fmt(max(fmt.dec, 6), False)
        steps.append(Step("Matriz em blocos $M=\\begin{bmatrix}A&B\\\\0&0\\end{bmatrix}T$", [TEX(r"M = " + dfmt.mat(big))]))
        steps.append(Step("Exponencial matricial e extração de $G$ e $H$", [
            TEX(r"e^{M} = " + dfmt.mat(sp.Matrix(from_np(E)))),
            TEX(r"G = " + dfmt.mat(G) + r",\qquad H = " + dfmt.mat(H)),
            TXT("Valores numéricos (SciPy `expm`, precisão dupla), mantidos com 12 algarismos significativos."),
        ]))
        exact = False
    extra = []
    if A.det() != 0:
        Hc = A.inv() * (G - sp.eye(n)) * B
        err = float(max(abs(complex(sp.N(e))) for e in (Hc - H)))
        extra.append(TXT(f"Verificação independente (A invertível): $H=A^{{-1}}(G-I)B$ difere em {err:.2e}."))
    else:
        extra.append(TXT("$A$ é singular: a identidade $H=A^{-1}(G-I)B$ não pode ser utilizada; mantém-se o cálculo pela integral original."))
    ev = np.linalg.eigvals(to_np(A).astype(complex))
    evd = np.linalg.eigvals(to_np(G).astype(complex))
    mapped = np.sort_complex(np.exp(ev * float(T)))
    extra.append(TXT("Polos: $z=e^{sT}$ → " + ", ".join(Fmt(4, False).num(e) for e in mapped)
                     + "; autovalores de $G$: " + ", ".join(Fmt(4, False).num(e) for e in np.sort_complex(evd))))
    steps.append(Step("Verificações", extra))
    return Result(steps, dict(G=G, H=H, exact=exact, T=T), f"Sistema discreto equivalente para T = {T} s.")


# ------------------------------------------------------------------ simulação
def simulate_ct(M, N, z0, w, t_end, npts, Cz, Dz):
    """Simula z' = Mz + N w (w constante) com a exponencial em blocos; saída = Cz z + Dz w."""
    M = np.asarray(M, dtype=float)
    N = np.asarray(N, dtype=float).reshape(-1)
    n = M.shape[0]
    t = np.linspace(0.0, float(t_end), int(npts))
    dt = t[1] - t[0]
    big = np.zeros((n + 1, n + 1))
    big[:n, :n] = M * dt
    big[:n, n] = N * dt
    E = expm(big)
    Phi, Gam = E[:n, :n], E[:n, n]
    Z = np.zeros((len(t), n))
    Z[0] = np.asarray(z0, dtype=float).reshape(-1)
    for k in range(len(t) - 1):
        Z[k + 1] = Phi @ Z[k] + Gam * w
    Y = Z @ np.asarray(Cz, dtype=float).T + w * np.asarray(Dz, dtype=float).reshape(1, -1)
    return t, Z, Y


def simulate_dt(M, N, z0, w, nsteps, T, Cz, Dz):
    M = np.asarray(M, dtype=float)
    N = np.asarray(N, dtype=float).reshape(-1)
    k = np.arange(int(nsteps) + 1)
    Z = np.zeros((len(k), M.shape[0]))
    Z[0] = np.asarray(z0, dtype=float).reshape(-1)
    for i in range(len(k) - 1):
        Z[i + 1] = M @ Z[i] + N * w
    Y = Z @ np.asarray(Cz, dtype=float).T + w * np.asarray(Dz, dtype=float).reshape(1, -1)
    return k * float(T), Z, Y
