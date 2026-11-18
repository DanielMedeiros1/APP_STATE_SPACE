# Espaço de Estados — Caderno de cálculos

Aplicativo Python / Streamlit de Daniel Medeiros Soares Carneiro para DCA3701, UFRN. Interface refeita a partir do projeto fornecido, seguindo o fluxo do aplicativo de LGR: dados na lateral, botão Resolver questão e resolução aberta. O núcleo matemático do ZIP original foi aproveitado e revisado; a interface e o caderno de resultados foram reconstruídos. 
https://app-p3-daniel-medeiros-soares-carneiro.streamlit.app/

## Executar no computador

Requer Python 3.10 ou superior (validado com Python 3.12).

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Publicar no Streamlit Community Cloud

Extraia o ZIP. Coloque `app.py`, `report.py`, `inputs.py`, `requirements.txt` e a pasta `statespace` na raiz do seu repositório. Preserve a pasta `statespace` com seus arquivos dentro dela. Não mova os arquivos internos para a raiz. A pasta `.streamlit` configura o tema claro e é opcional; `tests` contém os testes.

No Streamlit, selecione a branch `main` e o arquivo principal `app.py`. Ao atualizar um repositório existente, substitua os arquivos correspondentes. O aplicativo não foi publicado automaticamente.

## Resolver e capturar

1. Carregue um dos dez exemplos das listas ou digite o sistema.
2. Escolha os itens da questão e os polos desejados.
3. Clique em **Resolver questão**. Todos os cálculos selecionados aparecem abertos.
4. Informe o número da questão no campo de identificação.
5. Escolha frações exatas ou decimais e ajuste o tamanho das equações.
6. Selecione um bloco em **Bloco para capturar**, quando quiser um print por item.
7. Ative **Modo de captura** para esconder os dados de entrada. Use **Voltar à edição** para retornar.
8. Alternativamente, baixe o HTML, abra-o com conexão à internet e use o botão de impressão para salvar em PDF. As equações são renderizadas por MathJax via CDN. Esse recurso não envia suas matrizes a uma API de cálculo: elas já estão no documento exportado.

Colunas das matrizes: espaços. Linhas: Enter ou ponto e vírgula. Decimais: ponto ou vírgula. Exemplo de A: `0 1; -2 -3`. Um vetor B pode ser `0; 1` ou `0 1`. Coeficientes de função de transferência: potências decrescentes, separados por espaço. Polos: `-1±0,1j; -5; -12` ou `-2; -2; -3`. Também aceita `i` para a parte imaginária.

Se alterar dados, a resolução antiga permanece identificada até clicar novamente em Resolver questão. Não há atualização silenciosa do sistema que já foi calculado.

## Conteúdo

- Polinômio característico, polos e estabilidade interna.
- Controlabilidade: B, AB, A²B etc., matriz U/Wc, determinante e posto.
- Observabilidade: C, CA, CA² etc., matriz V/Wo, determinante e posto.
- Formas canônicas da disciplina, inclusive funções próprias com D diferente de zero.
- Ackermann: polinômio desejado, potências, q(A), inversa, produto matricial, ganho e conferência exata do polinômio.
- Observador e sistema combinado com realimentação de estados estimados.
- Seguidor contínuo de degrau, com ganho integral separado do ganho dos estados.
- Discretização exata por $G(T)=e^{AT}$ e $H(T)=\int_0^T e^{A\tau}B\,d\tau$; para matrizes nilpotentes, mostra as potências de A, a série finita e a integral desenvolvida. Projeto de K e L diretamente no sistema discretizado.
- Seguidor discreto de degrau, com transformação de ganhos consistente com u = Kx.
- Simulação opcional de saída e plano de polos; exportação CSV.
- Exportação de resolução em HTML e Markdown/LaTeX.

## Convenções e referências

Referência principal: capítulo 7 de `SisContr(2).pdf`, especialmente seções 7.4–7.18; numeração impressa aproximadamente 65–87 (o PDF possui cinco páginas iniciais). As duas atividades resolvidas e a terceira lista foram usadas para os exemplos e a sequência de apresentação.

- Realimentação: **u = Kx**, matrizes fechadas A + BK e G + HK.
- Ackermann: K = −eₙᵀ U⁻¹ q(A).
- Observador: L = q(A)V⁻¹ eₙ, dinâmica do erro A − LC / G − LC.
- Seguidor contínuo: e = y − r, ξ̇ = y − r, u = k₁ξ + K₂x.
- No seguidor discreto, a implementação usa u = −k₂x + k₁v, como na seção 7.18. O ganho K̂ do sistema aumentado é calculado na convenção positiva. Por substituição na matriz fechada, [k₂ k₁] = ([0 … 1] − K̂)Q⁻¹. A conferência polinomial valida a conversão de sinais.

As resoluções manuscritas orientam a apresentação. Na questão 3 da Lista 3, a discretização é desenvolvida separadamente por $G(T)=e^{AT}$ e pela integral de $H(T)$; os projetos seguem $W_c=[H\ GH\ \cdots]$, $W_o=[C;CG;\cdots]$, $K=-e_n^TW_c^{-1}q_c(G)$ e $L=q_o(G)W_o^{-1}e_n$. O exemplo carrega também a interligação com realimentação dos estados estimados. Resultados são recalculados, não copiados dos manuscritos. Por exemplo, o sistema 2 da primeira lista tem polinômio s³ − 8s² − 28s + 207; o seguidor da segunda lista tem coeficiente 137,17 no termo de primeiro grau.

## Limites

Sistemas reais SISO de ordem 1 a 6. Não realiza alocação parcial de polos de sistemas que não são completamente controláveis/observáveis; apresenta a condição e impede inversões singulares. O seguidor discreto requer D = 0 e transformação Q não singular.

A exponencial simbólica é oferecida para sistemas até ordem 3 com autovalores expressáveis sem solução cúbica/quártica geral; para A nilpotente, a série finita é exata em todas as ordens suportadas. Nos demais casos usa exponencial numérica em blocos com 12 algarismos significativos, explicitamente identificada. O projeto após discretização numérica é exato relativamente a esses valores aproximados.

A estabilidade reportada é interna, baseada nos autovalores da realização, incluindo análise de polos na fronteira; não se confunde com estabilidade BIBO de uma função de transferência após cancelamentos.

## Verificação

```bash
python -m pip install pytest
python -m pytest -q
```

Os testes verificam os exemplos, conversões, ganhos, polos repetidos, fronteira de estabilidade, exponenciais simbólicas e singulares, seguidores e interações da interface via Streamlit AppTest.
