"""Entrada por linhas e espaços; vírgulas decimais aceitas sem ambiguidade."""
import re
import sympy as sp
from statespace.parsing import InputError, parse_number, parse_poles

def number(text):
    t=str(text).strip().replace(',','.')
    if len(t)>70: raise InputError('Número muito longo.')
    val=parse_number(t)
    if abs(val)>10**12: raise InputError('Use valores de módulo até 10¹².')
    return val

def matrix(text):
    text=text.strip().strip('[]').replace('−','-')
    if len(text)>4000: raise InputError('Entrada muito longa.')
    rows=[r.strip().strip('[]') for r in re.split(r';|\n',text) if r.strip()]
    data=[[number(x) for x in row.split()] for row in rows]
    if not data or not data[0] or any(len(x)!=len(data[0]) for x in data):
        raise InputError('Use o mesmo número de elementos por linha. Separe colunas por espaço e linhas por ; ou Enter.')
    return sp.Matrix(data)

def vector(text):
    m=matrix(text)
    if min(m.shape)!=1: raise InputError('Informe um vetor, com valores separados por espaço.')
    return list(m)

def poles(text):
    if len(text)>500: raise InputError('Entrada de polos muito longa.')
    return parse_poles(text.replace(',','.'))

def system(a,b,c,d):
    A,B,C,D=map(matrix,(a,b,c,d));n=A.rows
    if not 1<=n<=6 or A.cols!=n: raise InputError('A (ou G) deve ser quadrada, de ordem 1 a 6.')
    if B.shape==(1,n):B=B.T
    if C.shape==(n,1):C=C.T
    if B.shape!=(n,1) or C.shape!=(1,n) or D.shape!=(1,1):
        raise InputError(f'Para ordem {n}, B/H precisa de {n} elementos, C de {n} elementos e D de um único valor.')
    return A,B,C,D
