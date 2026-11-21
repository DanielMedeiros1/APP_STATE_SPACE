"""Aplicativo Streamlit para resolução matricial, com modo de captura."""
import hashlib
import streamlit as st
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from statespace import core
from statespace.render import Fmt, steps_to_markdown
from statespace.examples import EXEMPLOS, tasks_for_module
import inputs as inp
import report as rp

st.set_page_config(page_title='Espaço de Estados | Daniel Medeiros',page_icon='🧮',layout='wide')
TASKS=['Estabilidade e polos','Controlabilidade','Observabilidade','Realimentação K','Observador L','Estados estimados','Seguidor de degrau']
DEFAULT=dict(A='-4 3 0; -1 -4 -1; 0 0 -2',B='0; 0; 3',C='1 0 0',D='0',num='3 4 -2',den='1 3 7 5',pk='-2±1j; -10',pl='-1±1j; -10',pf='-1±0,1j; -5; -12',T='1',tipo='Contínuo',entrada='Matrizes',canon='controlável',discretizar=False,task=TASKS[:5])
for key,val in DEFAULT.items():
    if key not in st.session_state: st.session_state[key]=val

def load_example():
    ex=EXEMPLOS[st.session_state.example]
    for key in ('A','B','C','D','num','den','canon','tipo','entrada','T'):
        st.session_state[key]=ex.get(key,DEFAULT[key])
    for dest,src in [('pk','poles_K'),('pl','poles_L'),('pf','poles_F')]: st.session_state[dest]=ex.get(src,'')
    mod=ex.get('modulo','')
    if st.session_state.example.startswith('A3 · Q3'):
        # A questão 3 pede a discretização, seguida do projeto K e do observador L.
        # Controlabilidade e observabilidade aparecem como pré-requisitos desses projetos;
        # estabilidade não é solicitada no enunciado.
        st.session_state.task=['Controlabilidade','Observabilidade','Estados estimados']
    else:
        st.session_state.task=tasks_for_module(mod)
    st.session_state.discretizar='Discretização' in mod
    st.session_state.pop('saved',None)

with st.sidebar:
    st.markdown('## ⚙️ Dados da questão')
    with st.expander('Carregar exemplo das listas'):
        st.selectbox('Exemplo',list(EXEMPLOS),key='example')
        st.button('Usar este exemplo',on_click=load_example,use_container_width=True)
    st.radio('Domínio',['Contínuo','Discreto'],key='tipo',horizontal=True)
    st.radio('Representação',['Matrizes','Função de transferência'],key='entrada')
    if st.session_state.entrada=='Matrizes':
        st.caption('Colunas: espaço · Linhas: Enter ou ; · Decimal: ponto ou vírgula. Vetores podem ficar em uma linha.')
        a,b=('G','H') if st.session_state.tipo=='Discreto' else ('A','B')
        st.text_area(a,key='A',height=110)
        st.text_area(b,key='B',height=80)
        st.text_input('C',key='C');st.text_input('D',key='D')
    else:
        st.caption('Coeficientes em ordem decrescente. Separe por espaço. Ex.: 1 3 7 5.')
        st.text_input('Numerador',key='num');st.text_input('Denominador',key='den')
        st.radio('Forma canônica',['controlável','observável'],key='canon',horizontal=True)
    if st.session_state.tipo=='Contínuo':st.checkbox('Discretizar antes de projetar',key='discretizar')
    if st.session_state.tipo=='Discreto' or st.session_state.discretizar:st.text_input('Período T (s)',key='T')
    st.multiselect('O que resolver',TASKS,key='task')
    if 'Realimentação K' in st.session_state.task or 'Estados estimados' in st.session_state.task:
        st.text_input('Polos para K',key='pk',help='Ex.: -2±1j; -10. No discreto, informe polos no plano z.')
    if 'Observador L' in st.session_state.task or 'Estados estimados' in st.session_state.task:st.text_input('Polos para L',key='pl')
    if 'Seguidor de degrau' in st.session_state.task:st.text_input('Polos do seguidor (n + 1)',key='pf')
    st.caption('Separe polos por ;. Aceita ±, i ou j. Ex.: -1±0,1j; -5; -12.')
    calculate=st.button('▶ Resolver questão',type='primary',use_container_width=True)
    st.divider()
    style=st.radio('Apresentação',['Decimais','Frações exatas'],horizontal=True)
    decimals=st.slider('Casas decimais',2,10,5)
    scale=st.slider('Tamanho das equações',75,140,105,5)
    capture=st.checkbox('Modo de captura',key='capture',help='Esconde a lateral e simplifica a página. Use o controle da barra lateral para voltar.')
    graphs=st.checkbox('Mostrar simulação e gráficos',value=False)
    symbolic=st.checkbox('Exponencial simbólica quando viável',value=True,help='Até ordem 3 com autovalores simples de expressar. Desative para obter uma aproximação rápida por exponencial em blocos.')

css='''<style>.block-container{max-width:1120px;padding-top:2rem;padding-bottom:3rem}h1{letter-spacing:-1px}h2{letter-spacing:-.4px}.badge{display:inline-block;padding:5px 12px;border-radius:30px;background:#eaf0ff;color:#2451a6;font-size:13px;font-weight:600}.stLatex{overflow-x:auto;overflow-y:hidden;padding:8px 0}.stLatex .katex-display{text-align:left!important}.stLatex .katex-display>.katex{text-align:left!important}div[data-testid="stVerticalBlockBorderWrapper"]{border-radius:12px}h3{font-size:1.18rem!important}@media print{[data-testid="stSidebar"],header,footer,button,[data-testid="stToolbar"]{display:none!important}.block-container{padding:0;max-width:none}.stLatex{overflow:visible!important}}</style>'''
st.markdown(css,unsafe_allow_html=True)
st.markdown('<style>.stLatex .katex{font-size:'+str(scale/100*1.1)+'em}</style>',unsafe_allow_html=True)
if 'question_label' in st.session_state: st.session_state.question_label=st.session_state.question_label
if capture:st.markdown('<style>[data-testid="stSidebar"]{display:none}header{background:transparent!important}.block-container{max-width:1040px}</style>',unsafe_allow_html=True)
if capture:
    st.button('← Voltar à edição',on_click=lambda: st.session_state.update(capture=False))
st.title('🧮 Espaço de Estados')
st.markdown('<span class="badge">DCA 3701 · Caderno de cálculos para a 3ª unidade</span>',unsafe_allow_html=True)
st.caption('Daniel Medeiros Soares Carneiro · u = Kx · malha fechada A + BK / G + HK · observador A − LC / G − LC')
if not capture:st.write('Informe os dados na lateral, escolha os itens da questão e clique em **Resolver questão**. A resolução fica aberta, pronta para conferir e capturar.')

params={key:st.session_state[key] for key in DEFAULT}
fingerprint=hashlib.sha256(repr(params).encode()).hexdigest()
if calculate:
    st.session_state.saved=dict(params=params,fingerprint=fingerprint)
if 'saved' not in st.session_state:
    st.info('← Preencha o sistema ou carregue um exemplo das listas. Depois clique em Resolver questão.')
    st.markdown('### O que aparece na resolução')
    st.write('Matrizes informadas; determinantes e postos; polinômios; potências matriciais; inversas; avaliação de q(A); ganhos K e L; matrizes de malha fechada e conferência dos polos. Cada item fica em um bloco para facilitar o print.')
    st.stop()
if fingerprint!=st.session_state.saved['fingerprint']:st.warning('Os dados foram alterados. Clique em Resolver questão para atualizar os resultados; abaixo permanece a última resolução.')
p=st.session_state.saved['params'];f=Fmt(decimals,style=='Frações exatas');sheets=[];errors=[];results={};discrete=p['tipo']=='Discreto'
try:
    if p['entrada']=='Matrizes':A,B,C,D=inp.system(p['A'],p['B'],p['C'],p['D'])
    else:
        num,den=inp.vector(p['num']),inp.vector(p['den'])
        if not 2<=len(den)<=7:raise ValueError('Use denominador de grau 1 a 6.')
        rr=core.tf_to_canonical(num,den,p['canon'],f,'z' if discrete else 's')
        sheets.append(rp.conversion(rr,'Conversão para a forma canônica '+p['canon']))
        A,B,C,D=[rr.data[k] for k in ('A','B','C','D')]
    sheets.insert(0,rp.system(A,B,C,D,f,discrete))
    if p['discretizar'] and not discrete:
        rr=core.discretize(A,B,inp.number(p['T']),f,symbolic=symbolic)
        sheets.append(rp.conversion(rr,'Discretização exata com entrada constante entre amostras'))
        A,B=rr.data['G'],rr.data['H'];discrete=True
        sheets.append(rp.system(A,B,C,D,f,True))
    active=p['task']
    def execute(label,fn):
        try:
            q,r=fn();sheets.append(q);results[label]=r
        except (ValueError,TypeError,ZeroDivisionError,NotImplementedError) as exc: errors.append(label+': '+str(exc))
    if 'Estabilidade e polos' in active:execute('stability',lambda:rp.stability(A,discrete,f))
    if 'Controlabilidade' in active:execute('ctrl',lambda:rp.rank_sheet(A,B,'ctrl',discrete,f))
    if 'Observabilidade' in active:execute('obs',lambda:rp.rank_sheet(A,C,'obs',discrete,f))
    if 'Realimentação K' in active or 'Estados estimados' in active:execute('K',lambda:rp.gain(A,B,inp.poles(p['pk']),False,discrete,f))
    if 'Observador L' in active or 'Estados estimados' in active:execute('L',lambda:rp.gain(A,C,inp.poles(p['pl']),True,discrete,f))
    if 'Estados estimados' in active:
        if results.get('K') and results['K'].data.get('K') is not None and results.get('L') and results['L'].data.get('L') is not None:
            rr=core.estimated_state_feedback(A,B,C,D,results['K'].data['K'],results['L'].data['L'],f,discrete)
            # The state labels must be discrete, not derivatives, for discrete systems.
            sheets.append(rp.conversion(rr,'Realimentação com estados estimados e princípio da separação'));results['estimated']=rr
        else:errors.append('Estados estimados: é necessário obter K e L antes de montar o sistema combinado.')
    if 'Seguidor de degrau' in active:execute('tracker',lambda:rp.tracker(A,B,C,D,inp.poles(p['pf']),discrete,f))
except (ValueError,TypeError,ZeroDivisionError,NotImplementedError) as exc:
    st.error(str(exc));st.stop()

if not capture:
    st.divider()
    left,right=st.columns([3,1])
    with left: name=st.text_input('Identificação no print (opcional)',key='question_label',placeholder='Ex.: Questão 2 — itens a e b')
    with right: st.caption('Apresentação');st.write(f'{"Discreto" if discrete else "Contínuo"} · ordem {A.rows}')
else:name=st.session_state.get('question_label','')
if name:st.subheader(name)
st.caption('Os valores exibidos são aproximações; os cálculos usam as entradas sem arredondamento intermediário.' if not f.exact else 'Resultados em frações e expressões exatas. Matrizes grandes podem exigir redução do tamanho das equações.')
if len(sheets)>2:
    focus=st.selectbox('Bloco para capturar',['Resolução completa']+[f'{i:02d} · {x.title}' for i,x in enumerate(sheets,1)],key='focus')
else:focus='Resolução completa'
chosen=sheets if focus=='Resolução completa' else [sheets[int(focus.split(' · ')[0])-1]]
for i,q in enumerate(sheets,1):
    if q not in chosen:continue
    with st.container(border=True):
        st.markdown(f'### {i:02d} · {q.title}')
        for kind,line in q.lines:
            if kind=='tex':
                # Long chains of equality are split vertically; no matrix is truncated.
                pieces=line.split(' = ') if len(line)>170 else [line]
                if len(pieces)>1:
                    for j,x in enumerate(pieces):st.latex(('' if j==0 else '= ')+x)
                else:st.latex(line)
            else:st.markdown(line)
for error in errors:st.error(error)

if not capture:
    st.divider()
    subtitle=(name+' · ' if name else '')+('Discreto' if discrete else 'Contínuo')+f' · ordem {A.rows} · '+style
    html=rp.export_html(sheets,subtitle)
    md='# Espaço de estados\n\n'+subtitle+'\n\n'+'\n'.join('## '+q.title+'\n\n'+'\n\n'.join('$$\n'+v+'\n$$' if k=='tex' else v for k,v in q.lines) for q in sheets)
    c1,c2=st.columns(2)
    c1.download_button('↓ Resolução para imprimir (HTML)',html,'resolucao_espaco_estados.html','text/html',use_container_width=True)
    c2.download_button('↓ Equações em Markdown / LaTeX',md,'resolucao_espaco_estados.md','text/markdown',use_container_width=True)
    st.caption('O HTML abre uma folha limpa com botão para imprimir ou salvar em PDF. Requer internet para carregar o renderizador de equações.')

if graphs and not capture:
    with st.expander('Simulação e plano de polos',expanded=True):
        which=st.selectbox('Sistema a simular',['Original']+(['Malha fechada com K'] if results.get('K') and results['K'].data.get('K') is not None else [])+(['Seguidor'] if results.get('tracker') and results['tracker'].data.get('Ka') is not None and not discrete else []))
        x0=st.text_input('Estado inicial x(0)',value=' '.join(['0']*A.rows))
        amp=st.number_input('Amplitude da entrada / referência',value=1.)
        duration=st.number_input('Número de amostras' if discrete else 'Tempo final (s)',min_value=1,value=30 if discrete else 15,max_value=500)
        try:
            M,N=core.to_np(A),core.to_np(B);cz,dz=core.to_np(C),core.to_np(D);z=np.array([float(x) for x in inp.vector(x0)])
            if len(z)!=A.rows:raise ValueError('O estado inicial deve ter n elementos.')
            if which=='Malha fechada com K':
                K=core.to_np(results['K'].data['K']);M=core.to_np(results['K'].data['Acl']);cz=cz+dz@K
            elif which=='Seguidor':
                tr=results['tracker'].data;M,N=core.to_np(tr['Acl']),core.to_np(tr['N']);z=np.r_[z,0.]
                cz=np.hstack([core.to_np(C)+core.to_np(D)@core.to_np(tr['K2']),core.to_np(D)*float(tr['k1'])]);dz=np.zeros((1,1))
            if discrete:t,X,Y=core.simulate_dt(M,N,z,amp,int(duration),float(inp.number(p['T'])),cz,dz)
            else:t,X,Y=core.simulate_ct(M,N,z,amp,float(duration),500,cz,dz)
            if not np.all(np.isfinite(Y)):raise ValueError('Resposta excedeu a faixa numérica.')
            fig,axs=plt.subplots(1,2,figsize=(12,4));axs[0].plot(t,np.asarray(Y).reshape(-1),'.-' if discrete else '-',label='Saída y')
            if which=='Seguidor':axs[0].axhline(amp,color='gray',ls='--',label='Referência')
            axs[0].set_xlabel('Tempo (s)');axs[0].grid(alpha=.2);axs[0].legend()
            ev=np.linalg.eigvals(M);axs[1].scatter(ev.real,ev.imag,marker='x',s=70)
            if discrete:axs[1].add_patch(plt.Circle((0,0),1,fill=False,color='gray',ls='--'))
            axs[1].axvline(0,color='gray',lw=.7);axs[1].axhline(0,color='gray',lw=.7);axs[1].set_xlabel('Parte real');axs[1].set_ylabel('Parte imaginária');axs[1].grid(alpha=.2)
            fig.tight_layout();st.pyplot(fig);plt.close(fig)
            import pandas as pd
            frame=pd.DataFrame(X,columns=[f'x{i+1}' for i in range(X.shape[1])]);frame.insert(0,'t',t);frame['y']=np.asarray(Y).reshape(-1)
            st.download_button('↓ Dados da simulação',frame.to_csv(index=False),'simulacao.csv','text/csv')
        except (ValueError,OverflowError) as exc:st.warning(str(exc))
