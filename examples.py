"""Exemplos pré-cadastrados (transcritos dos PDFs das atividades da 3ª unidade, DCA3701 2026.2).

Todos os resultados são calculados pelo aplicativo; aqui há apenas os dados de entrada.
Campos: tipo ('Contínuo'/'Discreto'), entrada ('Matrizes'/'Função de transferência'),
A,B,C,D (texto de matriz), num/den (potências decrescentes), canon, T, polos e sugestão de módulo.
"""
from __future__ import annotations

EXEMPLOS = {
    "A1 · Q1 e Q2 — Sistema 1 (3ª ordem)": dict(
        fonte="Atividade 1, questões 1 (Sistema 1) e 2 (itens a, b)",
        tipo="Contínuo", entrada="Matrizes",
        A="-4 3 0; -1 -4 -1; 0 0 -2", B="0; 0; 3", C="1 0 0", D="0",
        poles_K="-2+1j; -2-1j; -10", poles_L="-1+1j; -1-1j; -10", poles_F="",
        x0="1 0 0", xh0="0 0 0", modulo="Realimentação de estados",
        nota="Polinômio característico s³+10s²+35s+38; U com det = −81; V com det = −9; K = [8/3, −22/3, −4/3]; L = [2, −25/3, −16/3]ᵀ."),
    "A1 · Q1 e Q3 — Sistema 2 (3ª ordem)": dict(
        fonte="Atividade 1, questões 1 (Sistema 2) e 3 (itens c, d)",
        tipo="Contínuo", entrada="Matrizes",
        A="-5 1 0; 2 7 -1; 0 -3 6", B="0; 0; 1", C="1 -2 3", D="0",
        poles_K="-2+1j; -2-1j; -10", poles_L="-1+1j; -1-1j; -10", poles_F="",
        x0="1 0 0", xh0="0 0 0", modulo="Polinômio característico e estabilidade",
        nota="CONFERÊNCIA: pelos menores principais, p(s)=s³−8s²−28s+207 (o coeficiente de s no manuscrito parece diferente entre as páginas; o correto é −28). Sistema instável; controlável (det U = −1) e observável."),
    "A1 · Q4 — G(s)=(3s²+4s−2)/(s³+3s²+7s+5), canônica controlável": dict(
        fonte="Atividade 1, questão 4 (também Atividade 2, questão 1)",
        tipo="Contínuo", entrada="Função de transferência",
        num="3 4 -2", den="1 3 7 5", canon="controlável",
        poles_K="-4; -4; -5", poles_L="-1; -1; -2", poles_F="-1±0.1j; -5; -12",
        x0="1 0 0", xh0="0 0 0", modulo="Cadastro e conversão",
        nota="A=[0 1 0;0 0 1;−5 −7 −3], B=[0 0 1]ᵀ, C=[−2 4 3]. det V = −723; K = [−75, −49, −10]. O produto (s+4)²(s+5) é s³+13s²+56s+80."),
    "A2 · Q1 — seguidor para G(s) (polos −1±0,1j, −5, −12)": dict(
        fonte="Atividade 2, questão 1 b",
        tipo="Contínuo", entrada="Função de transferência",
        num="3 4 -2", den="1 3 7 5", canon="controlável",
        poles_K="-4; -4; -5", poles_L="-1; -1; -2", poles_F="-1±0.1j; -5; -12",
        x0="0 0 0", xh0="0 0 0", modulo="Seguidor de referência (degrau)",
        nota="CONFERÊNCIA: (s²+2s+1,01)(s²+17s+60) = s⁴+19s³+95,01s²+137,17s+60,6 (o coeficiente 137,17 = 120+17,17)."),
    "A2 · Q2 — 3ª ordem em matrizes (seguidor)": dict(
        fonte="Atividade 2, questão 2",
        tipo="Contínuo", entrada="Matrizes",
        A="-6 1 0; -11 0 1; -6 0 0", B="2; 6; 2", C="1 0 0", D="0",
        poles_K="-4; -4; -5", poles_L="-1; -1; -2", poles_F="-1±0.1j; -5; -12",
        x0="0 0 0", xh0="0 0 0", modulo="Seguidor de referência (degrau)",
        nota="CONFERÊNCIA: o enunciado traz a21 = −11. No manuscrito, as matrizes A e A_aumentada foram copiadas com +11, mas os produtos de U (por ex. −20 na 2ª coluna) correspondem a −11, então só a transcrição está errada. O app usa −11. det(U_a) = 16."),
    "A2 · Q3 — G(s)=(2s+1)/(s(s+3)), seguidor com polos em −1": dict(
        fonte="Atividade 2, questão 3",
        tipo="Contínuo", entrada="Função de transferência",
        num="2 1", den="1 3 0", canon="controlável",
        poles_K="-2; -3", poles_L="-2; -3", poles_F="-1; -1; -1",
        x0="0 0", xh0="0 0", modulo="Seguidor de referência (degrau)",
        nota="A=[0 1;0 −3], B=[0 1]ᵀ, C=[1 2]; K_a = [−1, −1, 0] (k₁ = −1, K₂ = [−1, 0])."),
    "A2 · Q4 — G(s)=3/(s²+2s−4), canônica observável, seguidor −2,−2,−3": dict(
        fonte="Atividade 2, questão 4",
        tipo="Contínuo", entrada="Função de transferência",
        num="3", den="1 2 -4", canon="observável",
        poles_K="-2; -3", poles_L="-2; -3", poles_F="-2; -2; -3",
        x0="0 0", xh0="0 0", modulo="Seguidor de referência (degrau)",
        nota="A=[0 4;1 −2], B=[3 0]ᵀ, C=[0 1]; K_a ≈ [−4, −5/3, −10/3]."),
    "A3 · Q1 — discreto de 2ª ordem": dict(
        fonte="Atividade 3, questão 1",
        tipo="Discreto", entrada="Matrizes",
        A="1.0 0.1; 0.5 0.1", B="1.0; 0.1", C="1 1", D="0", T="1",
        poles_K="0.5; 0.5", poles_L="0.2; 0.3", poles_F="",
        x0="1 0", xh0="0 0", modulo="Análise e projeto discretos",
        nota="Polos de G: z²−1,1z+0,05 → ≈1,0525 e ≈0,0475 (instável, |z|>1). Controlável e observável."),
    "A3 · Q2 — discreto de 3ª ordem": dict(
        fonte="Atividade 3, questão 2",
        tipo="Discreto", entrada="Matrizes",
        A="0.0877 0.1786 -0.0691; -0.0595 0.0877 -0.1056; 0 0 0.3679",
        B="-0.0561; -0.1438; 0.9482", C="1 0 0", D="0", T="1",
        poles_K="-0.5+0.5j; -0.5-0.5j; 0.5", poles_L="0.9+0.1j; 0.9-0.1j; 0.1", poles_F="",
        x0="1 0 0", xh0="0 0 0", modulo="Análise e projeto discretos",
        nota="Dados com 4 casas decimais, conforme o enunciado. O app mantém as frações exatas das entradas informadas."),
    "A3 · Q3 — discretização (T = 1 s) e projeto discreto": dict(
        fonte="Atividade 3, questão 3 (a, b, c)",
        tipo="Contínuo", entrada="Matrizes",
        A="0 1; 0 0", B="0; 1", C="1 1", D="0", T="1",
        poles_K="0.25+0.5j; 0.25-0.5j", poles_L="-0.5; 0.5", poles_F="",
        x0="1 0", xh0="0 0", modulo="Discretização",
        nota="A é singular (nilpotente): G=[1 1;0 1], H=[1/2; 1]. K = [−13/16, −35/32]; L = [5/4, 3/4]ᵀ."),
}
