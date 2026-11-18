"""Leitura segura de números, matrizes e polos (sem eval)."""
from __future__ import annotations

import re

import sympy as sp


class InputError(ValueError):
    """Erro de preenchimento com mensagem compreensível."""


_NUM = r"(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"


def _clean(t: str) -> str:
    return (t or "").strip().replace("−", "-").replace("–", "-").replace("\u2212", "-")


def parse_number(tok: str) -> sp.Rational:
    """Inteiro, decimal (com ponto) ou fração: 3, -0.5, 1/3, -25/3."""
    t = _clean(tok).replace(" ", "")
    m = re.fullmatch(rf"([+-]?)({_NUM})(?:/({_NUM}))?", t)
    if not m:
        raise InputError(
            f'Valor inválido: "{tok}". Use inteiros, decimais com ponto (0.25) ou frações (1/3, -25/3).'
        )
    val = sp.Rational(m.group(2))
    if m.group(3):
        den = sp.Rational(m.group(3))
        if den == 0:
            raise InputError(f'Divisão por zero em "{tok}".')
        val = val / den
    return -val if m.group(1) == "-" else val


def split_tokens(text: str) -> list[str]:
    t = _clean(text)
    if ";" in t or "\n" in t:
        toks = re.split(r"[;\n]+", t)
    else:
        toks = re.split(r"[\s,]+", t)
    return [x.strip() for x in toks if x.strip()]


def parse_list(text: str) -> list[sp.Rational]:
    return [parse_number(t) for t in split_tokens(text)]


def parse_matrix_text(text: str, rows: int, cols: int, name: str = "matriz") -> sp.Matrix:
    """Linhas separadas por quebra de linha ou ';'; elementos por espaço ou vírgula."""
    lines = [l for l in re.split(r"[;\n]+", _clean(text)) if l.strip()]
    data = [[parse_number(x) for x in re.split(r"[\s,]+", l.strip()) if x] for l in lines]
    if len(data) != rows or any(len(r) != cols for r in data):
        got = f"{len(data)}x{max((len(r) for r in data), default=0)}"
        raise InputError(f"A matriz {name} deve ter dimensão {rows}x{cols}; foi informada {got}.")
    return sp.Matrix(data)


def cells_to_matrix(cells: list[list[str]], rows: int, cols: int, name: str) -> sp.Matrix:
    if len(cells) != rows or any(len(r) != cols for r in cells):
        raise InputError(f"A matriz {name} deve ter dimensão {rows}x{cols}.")
    out = []
    for i in range(rows):
        row = []
        for j in range(cols):
            c = cells[i][j]
            try:
                row.append(parse_number(str(c)))
            except InputError as e:
                raise InputError(f"Matriz {name}, linha {i + 1}, coluna {j + 1}: {e}")
        out.append(row)
    return sp.Matrix(out)


def text_to_cells(text: str, rows: int, cols: int, name: str) -> list[list[str]]:
    lines = [l for l in re.split(r"[;\n]+", _clean(text)) if l.strip()]
    cells = [[x for x in re.split(r"[\s,]+", l.strip()) if x] for l in lines]
    if len(cells) != rows or any(len(r) != cols for r in cells):
        raise InputError(f"A matriz {name} deve ter dimensão {rows}x{cols}.")
    return cells


def cells_to_text(cells: list[list[str]]) -> str:
    return "\n".join(" ".join(str(c) for c in r) for r in cells)


def fit_cells(cells, rows: int, cols: int) -> list[list[str]]:
    """Redimensiona preservando o conteúdo (preenche com 0)."""
    out = []
    for i in range(rows):
        row = []
        for j in range(cols):
            try:
                row.append(str(cells[i][j]))
            except (IndexError, TypeError):
                row.append("0")
        out.append(row)
    return out


def _parse_complex_token(tok: str) -> list:
    s = _clean(tok).replace(" ", "").replace("+-", "±").replace("-+", "±")
    if not s:
        return []
    if s[-1] in "ij":
        body = s[:-1]
        pos = -1
        for p in range(len(body) - 1, 0, -1):
            if body[p] in "+-±" and body[p - 1] not in "eE":
                pos = p
                break
        if pos == -1:
            re_part, im_tok = sp.Integer(0), body
        else:
            re_part, im_tok = parse_number(body[:pos]), body[pos:]
        sign = im_tok[0] if im_tok and im_tok[0] in "+-±" else "+"
        mag = im_tok.lstrip("+-±")
        im = parse_number(mag) if mag else sp.Integer(1)
        if sign == "±":
            return [re_part + sp.I * im, re_part - sp.I * im]
        return [re_part + sp.I * (-im if sign == "-" else im)]
    if "±" in s:
        raise InputError(f'Polo "{tok}": o símbolo ± exige parte imaginária (ex.: -1±0.1j).')
    return [parse_number(s)]


def parse_poles(text: str) -> list:
    """Polos reais ou complexos: '-2; -2+1j; -2-1j' ou '-1±0.1j -5'.
    Separe por ';' (ou espaço/vírgula, se não houver ';'). Use ponto decimal."""
    poles: list = []
    for tok in split_tokens(text):
        poles += _parse_complex_token(tok)
    return poles


def parse_vector(text: str, n: int, name: str) -> sp.Matrix:
    vals = parse_list(text)
    if len(vals) != n:
        raise InputError(f"{name} deve ter {n} elementos; foram informados {len(vals)}.")
    return sp.Matrix(vals)
