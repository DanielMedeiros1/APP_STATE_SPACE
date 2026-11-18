"""Formatação (LaTeX / Markdown) independente da interface."""
from __future__ import annotations

import sympy as sp


class Fmt:
    """Formata números, matrizes e polinômios em LaTeX.

    exact=True  -> mantém frações/radicais (sympy);
    exact=False -> aproximação decimal com `decimals` casas (só na apresentação).
    """

    def __init__(self, decimals: int = 4, exact: bool = True):
        self.dec = int(decimals)
        self.exact = bool(exact)

    def _f(self, v: float) -> str:
        s = f"{v:.{self.dec}f}"
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        return "0" if s in ("-0", "") else s

    def num(self, x) -> str:
        x = sp.sympify(x)
        if self.exact:
            return sp.latex(x, imaginary_unit="j")
        try:
            v = complex(sp.N(x))
        except (TypeError, ValueError):
            return sp.latex(x, imaginary_unit="j")
        eps = 0.5 * 10 ** (-self.dec)
        if abs(v.imag) < eps:
            return self._f(v.real)
        re, im = self._f(v.real), self._f(abs(v.imag))
        sign = "-" if v.imag < 0 else "+"
        if re == "0":
            return ("-" if v.imag < 0 else "") + im + "j"
        return f"{re} {sign} {im}j"

    def mat(self, M) -> str:
        M = sp.Matrix(M)
        rows = [" & ".join(self.num(M[i, j]) for j in range(M.cols)) for i in range(M.rows)]
        return r"\begin{bmatrix}" + r" \\ ".join(rows) + r"\end{bmatrix}"

    def poly(self, coeffs, var: str = "s") -> str:
        """Polinômio a partir dos coeficientes em potências decrescentes."""
        n = len(coeffs) - 1
        out = ""
        for i, c in enumerate(coeffs):
            c = sp.sympify(c)
            p = n - i
            if c == 0:
                continue
            neg = bool(c.is_extended_negative)
            mag = -c if neg else c
            if mag == 1 and p > 0:
                cs = ""
            else:
                cs = self.num(mag)
                if not mag.is_real:
                    cs = r"\left(" + cs + r"\right)"
            pw = "" if p == 0 else (var if p == 1 else f"{var}^{{{p}}}")
            term = cs + pw
            if out == "":
                out = ("-" if neg else "") + term
            else:
                out += (" - " if neg else " + ") + term
        return out or "0"

    def factors(self, poles, var: str = "s") -> str:
        parts = []
        for p in poles:
            p = sp.sympify(p)
            if p.is_real:
                if p == 0:
                    parts.append(var)
                else:
                    parts.append(f"({var} {'-' if p > 0 else '+'} {self.num(abs(p))})")
            else:
                parts.append(r"\left(" + var + r"-\left(" + self.num(p) + r"\right)\right)")
        return "".join(parts)


def steps_to_markdown(title: str, steps, extra: str = "") -> str:
    """Converte etapas estruturadas em Markdown com LaTeX."""
    out = [f"# {title}", ""]
    if extra:
        out += [extra, ""]
    for st in steps:
        out += [f"## {st.title}", ""]
        for kind, content in st.items:
            if kind == "tex":
                out += ["$$", content, "$$", ""]
            else:
                out += [content, ""]
    return "\n".join(out)
