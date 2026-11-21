# PROVA 3 - APP ESPAÇO DE ESTADOS

Aplicativo Python / Streamlit desenvolvido por Daniel Medeiros Soares Carneiro para DCA3701, UFRN.

Link: https://app-prova3-daniel-medeiros-soares-carneiro.streamlit.app/

## Executar no computador

Requer Python 3.10 ou superior (validado com Python 3.12).

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Resolver e capturar

1. Carregue um dos doze exemplos das listas ou digite o sistema.
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
- Seguidor discreto de degrau, com transformação de ganhos consistente com a fórmula da apostila (seção 7.18).
- Simulação opcional de saída e plano de polos; exportação CSV.
- Exportação de resolução em HTML e Markdown/LaTeX.

## Convenções e referências

Referência principal: capítulo 7 de `SisContr(2).pdf`, especialmente seções 7.4–7.18; numeração impressa aproximadamente 65–87 (o PDF possui cinco páginas iniciais). As duas atividades resolvidas a terceira lista e o 4º Exercício da 3ª Unidade foram usados para os exemplos e a sequência de apresentação.

- Realimentação: **u = Kx**, matrizes fechadas A + BK e G + HK.
- Ackermann: K = −eₙᵀ U⁻¹ q(A).
- Observador: L = q(A)V⁻¹ eₙ, dinâmica do erro A − LC / G − LC.
- Seguidor contínuo: e = y − r, ξ̇ = y − r, u = k₁ξ + K₂x.
- No seguidor discreto, o relatório apresenta a fórmula de Ackermann da folha para $K_A$, a conversão $\widehat K=-K_A$ e a transformação $[k_2\;k_1]=(\widehat K+[0\;\cdots\;0\;1])Q^{-1}$; a lei implementada é $u_k=-k_2x_k+k_1v_k$.


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
