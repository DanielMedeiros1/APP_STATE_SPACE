"""Caderno de cálculo: fórmulas, matrizes e conferências, sem dependência de Streamlit."""
from dataclasses import dataclass, field
from html import escape
import json
import sympy as s
from statespace import core as c
from statespace.render import Fmt

@dataclass
class Sheet:
    title: str
    lines: list = field(default_factory=list)
    def tex(self, value): self.lines.append(('tex', value)); return self
    def text(self, value): self.lines.append(('txt', value)); return self


def system(A,B,C,D,f,discrete=False):
    a,b,_=c.names(discrete)
    q=Sheet('Sistema informado')
    q.tex(r'x_{k+1}=Gx_k+Hu_k,\qquad y_k=Cx_k+Du_k' if discrete else r'\dot x=Ax+Bu,\qquad y=Cx+Du')
    for name,m in [(a,A),(b,B),('C',C),('D',D)]: q.tex(name+'='+f.mat(m))
    return q


def stability(A,discrete,f):
    r=c.analyze_stability(A,discrete,f)
    q=Sheet('Polinômio característico, polos e estabilidade')
    for st in r.steps:
        for kind,line in st.items:
            if kind=='tex': q.tex(line)
    q.text(r.conclusion.replace('**',''))
    q.text('Critério: todos os módulos menores que 1.' if discrete else 'Critério: todas as partes reais estritamente negativas.')
    for b in r.data['stability']['boundary']:
        q.text(f"Polo na fronteira: multiplicidade algébrica {b['mult']}, geométrica {b['geo']}.")
    return q,r


def rank_sheet(A,M,mode,discrete,f):
    ctrl=mode=='ctrl'; n=A.rows; a,b,_=c.names(discrete)
    name=('W_c' if discrete else 'U') if ctrl else ('W_o' if discrete else 'V')
    r=c.controllability(A,M,f,discrete) if ctrl else c.observability(A,M,f,discrete)
    q=Sheet('Controlabilidade' if ctrl else 'Observabilidade')
    for k in range(n):
        nm=(b if k==0 else (a if k==1 else a+f'^{{{k}}}')+b) if ctrl else ('C' if k==0 else 'C'+(a if k==1 else a+f'^{{{k}}}'))
        val=A**k*M if ctrl else M*A**k
        q.tex(nm+'='+f.mat(val))
    blocks=([b]+[a+(f'^{{{k}}}' if k>1 else '')+b for k in range(1,n)]) if ctrl else (['C']+['C'+a+(f'^{{{k}}}' if k>1 else '') for k in range(1,n)])
    symbol=r'\begin{bmatrix}'+(' & ' if ctrl else r'\\').join(blocks)+r'\end{bmatrix}'
    mat=r.data['U' if ctrl else 'V']
    q.tex(name+'='+symbol).tex(name+'='+f.mat(mat))
    q.tex(r'\det('+name+')='+f.num(mat.det())+r',\qquad \operatorname{posto}('+name+')='+str(r.data['rank']))
    q.text(('Sistema controlável.' if ctrl else 'Sistema observável.') if r.data['rank']==n else ('Sistema não completamente controlável.' if ctrl else 'Sistema não completamente observável.'))
    q.text(f"Ordem n = {n}; critério: posto igual a n.")
    for warning in r.warnings: q.text(warning)
    return q,r


def gain(A,M,poles,observer,discrete,f):
    a,b,v=c.names(discrete); n=A.rows
    r=c.design_observer(A,M,poles,f,discrete) if observer else c.design_feedback(A,M,poles,f,discrete)
    q=Sheet('Projeto do observador de estados' if observer else 'Projeto da realimentação de estados')
    name=('W_o' if discrete else 'V') if observer else ('W_c' if discrete else 'U')
    label='L' if observer else 'K'; sub='o' if observer else 'c'
    q.tex(r'\mathcal{P}_d=\{'+r',\;'.join(f.num(x) for x in poles)+r'\}')
    q.tex(r'q_'+sub+'('+v+')='+f.factors(poles,v))
    if r.data.get(label) is None:
        q.text(r.conclusion.replace('**','')); return q,r
    coeff=r.data['coeffs']; Q=r.data['q']; U=r.data['V' if observer else 'U']
    q.tex(f'q_{sub}({v})='+f.poly(coeff,v))
    q.tex(name+'='+f.mat(U))
    for k in range(2,n+1): q.tex(a+f'^{{{k}}}='+f.mat(A**k))
    # A avaliação de q(F) é apresentada antes da inversa, na mesma sequência
    # algébrica usada na resolução manuscrita da Lista 3, questão 3.
    terms=[]
    for i,x in enumerate(coeff):
        power=n-i
        if x==0: continue
        term=('I' if power==0 else a if power==1 else a+f'^{{{power}}}')
        terms.append(('' if x==1 else '('+f.num(x)+')')+term)
    q.tex(f'q_{sub}({a})='+'+'.join(terms)).tex(f'q_{sub}({a})='+f.mat(Q))
    q.tex(name+'^{-1}='+f.mat(U.inv()))
    en=s.zeros(n,1);en[-1]=1
    if observer:
        row=U.inv()*en
        q.tex('L=q_o('+a+')'+name+r'^{-1}e_n')
        q.tex(r'e_n='+f.mat(en))
        q.tex(name+r'^{-1}e_n='+f.mat(row))
        q.tex('L='+f.mat(Q)+r'\,'+f.mat(row))
        if discrete:
            q.tex(r'\hat x_{k+1}=G\hat x_k+Hu_k+L(y_k-\hat y_k),\qquad \hat y_k=C\hat x_k+Du_k')
        else:
            q.tex(r'\dot{\hat x}=A\hat x+Bu+L(y-\hat y),\qquad \hat y=C\hat x+Du')
    else:
        row=en.T*U.inv()
        q.tex('K=-e_n^T'+name+r'^{-1}q_c('+a+')')
        q.tex('e_n^T='+f.mat(en.T))
        q.tex('e_n^T'+name+r'^{-1}='+f.mat(row))
        q.tex('K=-'+f.mat(row)+r'\,'+f.mat(Q))
    q.tex(r'\boxed{'+label+'='+f.mat(r.data[label])+'}')
    cl=r.data['Ao' if observer else 'Acl']; nm=a+'-LC' if observer else a+'+'+b+'K'
    q.tex(nm+'='+f.mat(cl))
    q.tex(r'\det('+v+'I-('+nm+'))='+f.poly(c.char_coeffs(cl),v))
    q.tex(c._eig_list_tex(c.eig_info(cl),f,'z' if discrete else r'\lambda'))
    q.text('Conferência exata: o polinômio obtido coincide com o desejado.' if r.data['ok'] else 'Atenção: falha na conferência do polinômio.')
    q.text('Lei de controle: u = Kx; o sinal já está incorporado em K.' if not observer else 'Dinâmica do erro: '+nm+'.')
    return q,r


def conversion(result,title):
    q=Sheet(title)
    for step in result.steps:
        for kind,line in step.items: q.lines.append((kind,line))
    return q


def tracker(A,B,C,D,poles,discrete,f):
    n=A.rows
    if not discrete:
        r=c.design_tracker(A,B,C,D,poles,f)
        aa=s.zeros(n+1);aa[0,1:]=C;aa[1:,1:]=A
        ba=s.Matrix.vstack(D,B)
        q=Sheet('Seguidor contínuo de referência — degrau')
        q.tex(r'e=y-r,\quad z=\dot x,\quad w=\dot u')
        q.text('Para D = 0, o projeto usa os estados [e; z] da apostila. A realização integral [ξ; x] tem as mesmas matrizes aumentadas. Se D ≠ 0, aplica-se diretamente a realização integral.')
        q.tex(r'A_a=\begin{bmatrix}0&C\\0&A\end{bmatrix}='+f.mat(aa))
        q.tex(r'B_a=\begin{bmatrix}D\\B\end{bmatrix}='+f.mat(ba))
        q.text('No cálculo de Ackermann abaixo, A e B representam as matrizes aumentadas Aₐ e Bₐ; o ganho K calculado é Kₐ.')
        design,_=gain(aa,ba,poles,False,False,f)
        q.lines+=design.lines
        if r.data.get('Ka') is not None:
            q.tex(r'K_a=[k_1\;K_2]='+f.mat(r.data['Ka']))
            q.tex('k_1='+f.num(r.data['k1'])+r',\qquad K_2='+f.mat(r.data['K2']))
            q.tex(r'\dot\xi=y-r,\qquad u=k_1\xi+K_2x')
            q.text('Implementação na ordem [x; ξ]:')
            q.tex(r'A_f='+f.mat(r.data['Acl'])).tex(r'B_r='+f.mat(r.data['N']))
            q.tex(r'\dot X=A_fX+B_rr,\qquad X=[x^T\;\xi]^T')
            q.text('Erro de regime nulo para referência constante quando a malha fechada é assintoticamente estável.')
        return q,r
    if D!=s.zeros(1,1): raise ValueError('O seguidor discreto da apostila requer D = 0.')
    # Apostila 7.18: states [x_e;u_e], input w, then transform Khat into [k2,k1].
    aa=A.row_join(B).col_join(s.zeros(1,n+1));ba=s.zeros(n+1,1);ba[-1]=1
    q=Sheet('Seguidor discreto de referência — degrau')
    q.tex(r'\widehat G=\begin{bmatrix}G&H\\0&0\end{bmatrix}='+f.mat(aa))
    q.tex(r'\widehat H=\begin{bmatrix}0\\1\end{bmatrix}='+f.mat(ba))
    q.text('No cálculo abaixo, G e H representam Ĝ e Ĥ; K corresponde ao ganho K̂ do sistema aumentado.')
    ds,rr=gain(aa,ba,poles,False,True,f);q.lines+=ds.lines
    if rr.data.get('K') is None:return q,rr
    kh=rr.data['K'];Q=(A-s.eye(n)).row_join(B).col_join((C*A).row_join(C*B))
    q.tex(r'Q=\begin{bmatrix}G-I&H\\CG&CH\end{bmatrix}='+f.mat(Q))
    if Q.det()==0: raise ValueError('A transformação do seguidor é singular. Não se aplica a alocação arbitrária por esse método.')
    gains=(s.Matrix([[0]*n+[1]])-kh)*Q.inv();k2=gains[:,:n];k1=gains[0,n]
    q.text('Conversão compatível com a convenção w = K̂[xₑ; uₑ]. Ela decorre de K̂ = [0 … 1] − [k₂ k₁]Q.')
    q.tex(r'[k_2\;k_1]=([0\;\cdots\;0\;1]-\widehat K)Q^{-1}')
    q.tex('Q^{-1}='+f.mat(Q.inv())).tex(r'[k_2\;k_1]='+f.mat(gains))
    q.tex('k_1='+f.num(k1)+r',\qquad k_2='+f.mat(k2))
    q.tex(r'v_k=v_{k-1}+r_k-y_k,\qquad u_k=-k_2x_k+k_1v_k')
    cl=A.row_join(B).col_join((k2-k2*A-k1*C*A).row_join(s.ones(1,1)-k2*B-k1*C*B))
    q.tex('G_f='+f.mat(cl))
    q.tex(r'\det(zI-G_f)='+f.poly(c.char_coeffs(cl),'z'))
    q.text('Conferência exata: '+('polinômio coincide com o desejado.' if c.char_coeffs(cl)==rr.data['coeffs'] else 'divergência encontrada.'))
    rr.data.update(Ka=kh,k1=k1,K2=-k2,Acl=cl,discrete_tracker=True)
    return q,rr


def export_html(sheets,subtitle):
    parts=[]
    for i,q in enumerate(sheets,1):
        body=''.join('<div class="math">\\['+escape(t)+'\\]</div>' if k=='tex' else '<p>'+escape(t.replace('**',''))+'</p>' for k,t in q.lines)
        parts.append('<section><h2>'+str(i).zfill(2)+' · '+escape(q.title)+'</h2>'+body+'</section>')
    return '''<!doctype html><html lang="pt-BR"><meta charset="utf-8"><title>Resolução — Espaço de estados</title>
<script>window.MathJax={tex:{inlineMath:[['$','$'],['\\\\(','\\\\)']]},svg:{fontCache:'global'}};</script>
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>
<style>body{font:16px system-ui;color:#172033;background:#eef2f7;margin:0}main{max-width:1050px;margin:35px auto}h1{font-size:30px}section{background:white;border:1px solid #dce2ed;border-radius:12px;padding:24px 32px;margin:20px 0;break-inside:auto}h2{font-size:21px;border-bottom:2px solid #2563eb;padding-bottom:14px;break-after:avoid}.math{overflow:auto;margin:20px 0;break-inside:avoid}p{line-height:1.6}button{padding:12px 20px;border:0;border-radius:8px;background:#2563eb;color:white;cursor:pointer}@media print{body{background:white}main{margin:0;max-width:none}button{display:none}section{border-radius:0;padding:14px;border:1px solid #ccc}.math{font-size:85%;overflow:visible}h1{font-size:23px}@page{size:A4;margin:12mm}}</style><main><button onclick="window.print()">Imprimir / salvar em PDF</button><h1>Espaço de estados · Resolução</h1><p>'''+escape(subtitle)+'''</p><p>Valores decimais são aproximações de apresentação; os cálculos mantêm a precisão das entradas. Conecte à internet ao abrir para carregar as equações.</p>'''+''.join(parts)+'</main></html>'
