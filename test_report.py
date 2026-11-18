import pytest
import sympy as s
from statespace.render import Fmt
from statespace.examples import EXEMPLOS
from statespace import core
import inputs as inp
import report as rp
F=Fmt()

@pytest.mark.parametrize('name',list(EXEMPLOS))
def test_examples(name):
    ex=EXEMPLOS[name]
    if ex['entrada']=='Matrizes':A,B,C,D=inp.system(ex['A'],ex['B'],ex['C'],ex['D'])
    else:
        r=core.tf_to_canonical(inp.vector(ex['num']),inp.vector(ex['den']),ex['canon'],F)
        A,B,C,D=[r.data[k] for k in ('A','B','C','D')]
    discrete=ex['tipo']=='Discreto'
    if 'Discretização' in ex['modulo']:
        r=core.discretize(A,B,s.Rational(ex['T']),F);A,B=r.data['G'],r.data['H'];discrete=True
    for mode,M in [('ctrl',B),('obs',C)]:
        sheet,r=rp.rank_sheet(A,M,mode,discrete,F)
        assert r.data['rank']==A.rows
    for observer,M,pk in [(False,B,'poles_K'),(True,C,'poles_L')]:
        q,r=rp.gain(A,M,inp.poles(ex[pk]),observer,discrete,F)
        assert r.data['ok']
        assert any('boxed' in v for k,v in q.lines)
    if ex['poles_F']:
        q,r=rp.tracker(A,B,C,D,inp.poles(ex['poles_F']),discrete,F)
        assert r.data['ok']

def test_discrete_tracker_and_repeated_poles():
    G=s.Matrix([[1,1],[0,1]]);H=s.Matrix([s.Rational(1,2),1]);C=s.Matrix([[1,1]])
    q,r=rp.tracker(G,H,C,s.zeros(1,1),[s.Rational(1,2)]*3,True,F)
    assert core.char_coeffs(r.data['Acl'])==[1,-s.Rational(3,2),s.Rational(3,4),-s.Rational(1,8)]

def test_symbolic_exponential_non_nilpotent():
    A=s.Matrix([[0,1],[0,-2]]);B=s.Matrix([0,1])
    r=core.discretize(A,B,s.Integer(1),F)
    assert r.data['exact']
    assert s.simplify(r.data['G'][0,1]-(1-s.exp(-2))/2)==0
    assert s.simplify(r.data['H'][0]-(1+s.exp(-2))/4)==0

def test_decimal_comma_and_invalid():
    assert inp.matrix('1,5 0; -2 0,25')==s.Matrix([[s.Rational(3,2),0],[-2,s.Rational(1,4)]])
    assert len(inp.poles('-1±0,1j; -5; -12'))==4
    with pytest.raises(ValueError):inp.number('__import__("os")')

def test_export_is_escaped():
    html=rp.export_html([rp.Sheet('<script>').text('<img src=x>')],'<test>')
    assert '&lt;img src=x&gt;' in html and '&lt;script&gt;' in html
