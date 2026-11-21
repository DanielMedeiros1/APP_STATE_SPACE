import pytest
import sympy as s
from statespace.render import Fmt
from statespace.examples import EXEMPLOS, tasks_for_module
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
        if ex.get(pk):
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


def test_list4_question1_canonical_discretization_and_tracker():
    ex=EXEMPLOS["A4 · Q1 — (3s−2)/(s²+s), canônica controlável, discretização e seguidor"]
    conv=core.tf_to_canonical(inp.vector(ex['num']),inp.vector(ex['den']),ex['canon'],F)
    A,B,C,D=[conv.data[k] for k in ('A','B','C','D')]
    assert A==s.Matrix([[0,1],[0,-1]])
    assert B==s.Matrix([0,1])
    assert C==s.Matrix([[-2,3]])
    disc=core.discretize(A,B,s.Integer(1),F)
    G,H=disc.data['G'],disc.data['H']
    assert G==s.Matrix([[1,1-s.exp(-1)],[0,s.exp(-1)]])
    assert H==s.Matrix([s.exp(-1),1-s.exp(-1)])
    q,r=rp.tracker(G,H,C,D,inp.poles(ex['poles_F']),True,F)
    assert r.data['ok']
    assert r.data['coeffs']==[1,-s.Rational(9,10),s.Rational(19,100),s.Rational(29,1000)]
    assert abs(float(r.data['k1'].evalf())+0.2523252847)<1e-8
    assert max(abs(float(v.evalf())-e) for v,e in zip(r.data['Ka'], [0.50465057,0.43737330,0.46787944])) < 1e-7
    assert max(abs(float(v.evalf())-e) for v,e in zip(-r.data['K2'], [2.15757803,1.52977190])) < 1e-7
    assert r.data['K_A'] == r.data['K']
    assert r.data['Ka'] == -r.data['K_A']
    assert any(r'K_A=-e_n^TW_c^{-1}q_c(\widehat G)' == line for kind,line in q.lines if kind=='tex')
    assert any(r'\widehat K=-K_A=' in line for kind,line in q.lines if kind=='tex')
    assert any(line.startswith(r'\det(zI-(\widehat G+\widehat H K_A))=') for kind,line in q.lines if kind=='tex')
    assert len([line for kind,line in q.lines if kind=='tex'])>10
    assert any(r'[k_2\;k_1]=(\widehat K+[0\;\cdots\;0\;1])Q^{-1}' in line for kind,line in q.lines if kind=='tex')
    assert not any('Fórmula de transformação impressa' in line or 'Nota de convenção' in line or 'Lei de controle: u = Kx' in line for _,line in q.lines)


def test_list4_question2_discrete_tracker_is_possible():
    ex=EXEMPLOS["A4 · Q2 — sistema discreto de 3ª ordem, seguidor com polos ±0,5±0,5i"]
    G,H,C,D=inp.system(ex['A'],ex['B'],ex['C'],ex['D'])
    q,r=rp.tracker(G,H,C,D,inp.poles(ex['poles_F']),True,F)
    assert r.data['ok']
    assert r.data['coeffs']==[1,0,0,0,s.Rational(1,4)]
    assert r.data['k1']==s.Rational(5,48)
    assert r.data['K2']==s.Matrix([[-s.Rational(263,192),-s.Rational(667,384),-s.Rational(91,64)]])
    assert r.data['K_A'] == r.data['K']
    assert r.data['Ka'] == -r.data['K_A']
    assert any(r'[k_2\;k_1]=(\widehat K+[0\;\cdots\;0\;1])Q^{-1}' in line for kind,line in q.lines if kind=='tex')
    Q=(G-s.eye(3)).row_join(H).col_join((C*G).row_join(C*H))
    assert Q.det()==-12


def test_list4_question1_loads_tracker_module():
    ex=EXEMPLOS["A4 · Q1 — (3s−2)/(s²+s), canônica controlável, discretização e seguidor"]
    tasks=tasks_for_module(ex["modulo"])
    assert "Seguidor de degrau" in tasks
    assert "Realimentação K" not in tasks
    assert "Observador L" not in tasks
    assert ex["poles_F"] == "0.5+0.2j; 0.5-0.2j; -0.1"
