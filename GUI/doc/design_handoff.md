# Handoff Spec: SPLACE (Flutter) — correções de design

Revisão crítica dos pedidos das três últimas rodadas comparada com o estado
atual (capturas em `build/shots_v3` → `build/shots_v5`, 1280×720) e com o app
HTML (`docs/`, renderizado com Edge em `build/html_ref/`).

## 1. Diagnóstico

### 1.1 O que piorou

| # | Onde | Antes (v3) | Agora (v5) | Por que piorou |
|---|------|------------|------------|----------------|
| P1 | Switches (Etapa 1, MAFFT, campos) | switch M3 ~32 px | encolhido para 22 px (`FittedBox`) | alvo de clique pequeno, ilegível ao lado do texto de 14 px |
| P2 | Select "Tipo de genoma" | valor 16/600 (DropdownButton) | valor 14/600 + ícone DNA + chevron baixo | valor em negrito ao lado de campos com texto normal; ícone parece ampulheta; chevron desalinhado |
| P3 | Linha de busca | rótulos alinhados em 2 linhas | 4 campos, mas rótulo "Tipo de genoma" 2 px abaixo dos demais | alinhamento por baixo (`CrossAxisAlignment.end`) com alturas diferentes |
| P4 | Modais | larguras 480–900, altura pelo conteúdo | todos 1000×720 fixos | modal da chave NCBI com metade vazia; progresso com 1000 px sobre janela de 1280 |
| P5 | Tooltips | alinhados à esquerda | justificados em 320 px | "rios" de espaço entre palavras em linhas curtas |
| P6 | Tela inicial | número em círculo suave (azul-claro, número azul) | círculo azul cheio | o pedido era "como antes"; mudou o estilo |

### 1.2 O que não foi feito como pedido

| # | Pedido | Situação | Ação |
|---|--------|----------|------|
| N1 | Logo preenchendo tudo, sem bordas brancas | tela inicial corrigida, mas o ícone do app (janela/barra de tarefas) ainda é o logo do Flutter | gerar `app_icon.ico` (Windows), ícones macOS e Linux a partir do logo SPLACE, preenchendo o quadrado |
| N2 | Taxonomia: as 3 fontes via API, inclusive SALVE | SALVE continua por tabela | **bloqueado**: não há API pública documentada e o acesso ao site do SALVE foi negado neste ambiente. Precisa do endpoint (DevTools → Network ao abrir uma espécie) ou de permissão para consultar `salve.icmbio.gov.br` |
| N3 | Janela nunca abaixo de 1280×720 | trava nativa feita e testada no Release; o usuário roda um `Debug\splace.exe` de 00:14, anterior à trava | parar o `flutter run` e rodar de novo (hot reload não recompila código nativo) |
| N4 | Boas práticas de tipografia no app todo | só os modais e as dicas foram tratados; restam 22 tamanhos de fonte (10 a 20 px, com 12.4, 12.8, 13.6, 14.4…) e pesos 500/600/700/800 misturados | aplicar a escala da seção 2 em todo `lib/src/ui` |

## 2. Tokens

### 2.1 Tipografia (Montserrat; mono JetBrains Mono)

| Token | Tamanho/peso | Uso |
|-------|--------------|-----|
| `type-display` | 44 / 800 | "SPLACE" da tela inicial (único) |
| `type-brand` | 20 / 700 | "SPLACE" no cabeçalho |
| `type-title` | 16 / 600 | título de etapa, card, modal |
| `type-subtitle` | 14 / 600 | subtítulos, nomes em listas, rótulos de campo, botões |
| `type-body` | 14 / 400 | texto, **valor de campo e de select** |
| `type-body-sm` | 13 / 400 | descrições, ajuda, células de tabela |
| `type-caption` | 12 / 400 | metadados, contagens, rodapés |
| `type-kicker` | 11 / 700, CAIXA ALTA, espaçamento 0,08 em | rótulos de seção ("OPÇÃO 1") |
| `type-mono` | 12 / 400 | acessos, comandos, cabeçalhos FASTA |

Regras:
- Um só peso por frase. Ênfase dentro de uma frase usa cor (link `accent`, sublinhado), não negrito.
- Exceções permitidas: siglas da marca ("**SPL**it **A**lign…") e o par "Título — descrição" das etapas da tela inicial, que são dois elementos e não uma frase.
- Valores dentro de campos (TextField, select, número) sempre `type-body` (14/400). Rótulo acima do campo `type-subtitle`.
- Sem 800 fora do `type-display`.

### 2.2 Cores (já existentes em `SplaceTokens`)

| Token | Claro | Uso |
|-------|-------|-----|
| `accent` | #323795 | botões primários, números, links |
| `accentSoft` | azul-50 | fundo do número das etapas, seleção |
| `accentSoftBorder` | azul-200 | borda do número, trilho da timeline |
| `accentText` | azul-700 | texto sobre `accentSoft` |
| `surface2` | #F3F4F6 | fundo de cards da tela inicial e rodapés de modal |
| `textMuted` | cinza-500 | descrições |

### 2.3 Espaçamento e raios

| Token | Valor | Uso |
|-------|-------|-----|
| `space-xs` | 4 | ícone ↔ texto |
| `space-sm` | 8 | entre rótulo e campo |
| `space-md` | 12 | entre campos de uma linha |
| `space-lg` | 16 | entre cards |
| `space-xl` | 24 | padding interno de modal |
| `control-h` | 44 | altura de TextField, select, caixa de switch, botão |
| `radius-sm` | 8 | campos, botões |
| `radius-md` | 10–12 | cards |
| `radius-lg` | 18 | modais |

## 3. Componentes

### 3.1 Tela inicial (`chrome/splash.dart`)

- Hero como está (logo 136 px de altura preenchendo o quadro, sem margem branca).
- Barra colorida: gradiente `splaceRainbow` com 160 % da largura, deslizando 0 % → 100 % em 4 s, `easeInOut`, vai-e-volta (= `splashGradientShift` do HTML).
- 8 linhas, cada uma um card (`surface2`, borda `border`, raio 10, altura = espaço disponível ÷ 8, entre 34 e 58 px, margem vertical 3).
- **Dentro do card**, à esquerda: círculo de 30 px **no estilo de antes / do HTML** (`.splash-step-number-primary`): fundo `accentSoft`, borda 1 px `accentSoftBorder`, número 13/600 `accentText`.
- Timeline: trilho de 2 px `accentSoftBorder` no centro dos círculos, atravessando os cards e os intervalos, do centro do 1º ao centro do 8º.
- Texto numa linha: título `type-subtitle` + " — " + descrição `type-body-sm` `textMuted`; reticências se não couber.
- Botão "Iniciar" 50 px, largura total.

### 3.2 Linha de busca da Opção 1 (`sections/load_section.dart`)

- Quatro colunas na mesma linha, `flex` 3/3/3/3, espaço `space-md`.
- Rótulos alinhados no topo (`CrossAxisAlignment.start`); todos os controles com `control-h` (44).
- Nome: TextField com ícone de lupa à esquerda.
- Select (3.3) para o tipo de genoma.
- Escopo e banco: caixa de switch (3.4).

### 3.3 Select (`SplaceSelect`, `widgets/kit.dart`)

- Campo: só o texto do valor (`type-body`, 14/400), sem ícone; chevron 11 px centrado verticalmente à direita.
- Menu: largura do campo, raio 12, padding 6, itens 40 px com ícone 14 px + rótulo `type-body` + detalhe `type-caption`; item atual com fundo `accentSoft` e check à direita; hover `surface2`.
- Teclado: Enter/Espaço abre; setas navegam; Esc fecha.

### 3.4 Caixa de switch (`SwitchBox`)

- Altura 44, raio 8, padding horizontal 10.
- Switch M3 com 26 px de altura (≈ 42×26), trilho `accent` com check quando ligado.
- Caixa inteira clicável; ligada = fundo `accentSoft`, borda `accentSoftBorder`, texto `accentText`; o peso do texto **não muda** (14/600 sempre).

### 3.5 Tooltip (`SplaceTooltip`)

- Texto justificado (pedido), largura máxima **360** (reduz rios), padding 12×8, `type-body-sm` 13/500 branco sobre #323795, seta 6 px.
- Última linha alinhada à esquerda (padrão do justify).

### 3.6 Modais (`SplaceDialog`)

| Tamanho | Largura | Altura | Uso |
|---------|---------|--------|-----|
| `standard` | 880 (ou janela − 64) | acompanha o conteúdo, até janela − 64; rola por dentro | chave NCBI, ferramentas, editar registro, construtor de cabeçalho, progresso |
| `standard` + `fill` | 880 | sempre janela − 64 (máx. 760) | listas/tabelas: revisão de genomas, genes não reconhecidos |
| `compact` | 560 | conteúdo | confirmações |

- Cabeçalho do HTML (`.modal-header`): gradiente 135° #1E2263 → #323795, ícone 40 px em caixa branca 15 %, título `type-title` 16/700 branco, subtítulo 13/400 branco 72 %.
- Corpo: padding `space-xl` (24).
- Rodapé: `surface2`, borda superior; à esquerda um texto `type-caption` opcional; à direita as ações (secundária em TextButton, primária em FilledButton).

### 3.7 Modal da chave NCBI (`chrome/api_key_dialog.dart`)

Uma coluna, de cima para baixo:
1. Situação numa linha (TonePanel info/sucesso): "Nenhuma chave salva — até 3 requisições/s" ou "Chave salva — até 10 requisições/s".
2. Rótulo "Sua API key" + linha: campo (ícone de chave, colar, mostrar/ocultar) + botão "Salvar chave" (44 px) ao lado, como no HTML.
3. Validação logo abaixo do campo (`type-caption`): "Formato válido" (verde) / "A API key parece inválida" (âmbar); só aparece com texto digitado.
4. Separador, "Como obter uma chave" (`type-subtitle`) e os 3 passos numerados (círculo 22 px `accentSoft`, texto `type-body-sm`), e o link "Abrir conta do NCBI".
- Rodapé: à esquerda "Guardada só neste computador…"; à direita "Apagar chave" (texto vermelho, só se houver chave) e "Fechar".

### 3.8 Modal de revisão de genomas

Mantém o layout enxuto da v5 (filtro, menu "Selecionar", aviso UNVERIFIED de uma linha, tabela Seleção/Acesso/Organismo/Descrição/Situação, rodapé "N de M selecionados" + "Importar N"), com `standard + fill`.

### 3.9 Modal de progresso

`standard` (880), altura pelo conteúdo. Item atual `type-body`, taxa e chave numa linha `type-caption`, barra 8 px, contador e percentual `type-caption`.

### 3.10 Card "Taxonomia dos registros"

- Três cartões de fonte em `Row` com `Expanded` (mesma largura) e `IntrinsicHeight` (mesma altura), espaço `space-md`.
- Pílula do dataFishing com um só peso (13/400, cor `accentText`), sem negrito parcial.
- SALVE: ver N2 — quando houver endpoint, o cartão passa a "online", e a busca segue o mesmo fluxo do GBIF (lote de 5, modal de progresso, cache).

### 3.11 Ícone do app

- Windows `windows/runner/resources/app_icon.ico` (16, 24, 32, 48, 64, 128, 256), macOS `AppIcon.appiconset`, Linux (ícone da janela via GTK).
- Arte: o cartão do logo (cinza #EEEEEE com o desenho) **ocupando o quadrado inteiro**: fundo #EEEEEE até as bordas, desenho centrado com margem de 4 %, sem transparência nem borda branca.

## 4. Estados e interações

| Elemento | Estado | Comportamento |
|----------|--------|---------------|
| Select | hover | borda `textFaint` |
| Select | aberto | chevron vira para cima, borda `accent` 2 px |
| SwitchBox | hover | cursor de clique |
| SwitchBox | desabilitado | opacidade 50 %, sem cursor |
| Botão "Salvar chave" | chave inválida | desabilitado (azul 40 %) |
| Modal | Esc | fecha (exceto progresso em andamento) |
| Linha da tabela de revisão | clique | alterna a seleção; selecionada = `accentSoft` 50 % |

## 5. Comportamento por largura

| Largura da janela | Mudança |
|-------------------|---------|
| 1280 (mínimo) | barra lateral 264; log flutuante; modais 880 |
| ≥ 1500 | log fixo à direita (340) |
| ≥ 1600 | conteúdo continua ocupando toda a largura útil (sem limite) |

Abaixo de 1280×720 a janela não reduz (WM_GETMINMAXINFO + WM_SIZING + WM_WINDOWPOSCHANGING; GTK geometry hints; `contentMinSize` no macOS).

## 6. Casos-limite

- Textos longos em PT/ES: descrição das etapas com reticências numa linha; rótulos de switch até 2 linhas.
- Lista de revisão vazia após filtro: "Nenhum genoma corresponde ao filtro." centralizado.
- Sem chave NCBI: aviso informativo, não erro.
- SALVE sem endpoint: cartão continua indicando "tabela" até a integração.

## 7. Movimento

| Elemento | Gatilho | Animação | Duração | Curva |
|----------|---------|----------|---------|-------|
| Barra colorida | contínuo | gradiente desliza 60 % | 4000 ms, vai-e-volta | easeInOut |
| SwitchBox | toggle | cor de fundo/borda | 150 ms | linear |
| Tooltip | hover 250 ms | fade | 140 ms | linear |

## 8. Acessibilidade

- Ordem de foco na linha de busca: nome → tipo de genoma → escopo → banco → Buscar.
- Switches e checkboxes com rótulo clicável (a caixa inteira).
- Contraste: texto `textMuted` sobre `surface2` ≥ 4,5:1; tooltip branco sobre #323795.
- Modais: foco inicial no primeiro campo; Esc fecha; botão de fechar com tooltip "Fechar".

## 9. Critérios de aceite

1. Captura 1280×720 da tela inicial cabe inteira, com números em círculo claro dentro dos cards e trilho ligando-os.
2. Na linha de busca, os quatro rótulos têm o mesmo topo e os quatro controles 44 px; o valor do select tem o mesmo peso do texto do campo de nome.
3. Nenhum modal de conteúdo com área vazia maior que 1/3 da altura; todos com 880 de largura.
4. Ícone da janela e da barra de tarefas é o logo SPLACE, sem borda branca.
5. `grep -rhoE "fontSize: [0-9.]+" lib/src/ui` só retorna valores da escala (11, 12, 13, 14, 16, 18, 20, 44 e os tamanhos de ícone/gráfico documentados).
6. Analisador sem avisos, testes passando, Release compilado, janela não reduz abaixo de 1280×720.

## 10. Situação (06/11/2026)

| Item | Estado | Verificação |
|------|--------|-------------|
| P1 switches | feito (26 px, caixa ≥ 44 px, rótulo longo quebra linha) | `build/shots_v6/02_load.png`, `07c_mafft.png` |
| P2 select | feito (valor 14/400, sem ícone no campo, chevron centrado) | `02e_select.png` |
| P3 rótulos da busca | feito (alinhados no topo) | `02_load.png` |
| P4 modais | feito (880 de largura, altura pelo conteúdo; revisão com `fill`) | `02a_api_key.png`, `02c_review.png`, `02d_progress_modal.png` |
| P5 tooltips | feito (justificado, 360 px) | `02b_tooltip.png` |
| P6 números da tela inicial | feito (círculo claro dentro do card, trilho ligando) | `01_splash.png` |
| N1 ícone do app | feito (Windows .ico 16–256, macOS 16–1024, Linux via GTK) | ícone extraído do `splace.exe` Release |
| N2 SALVE via API | **bloqueado** (endpoint desconhecido; acesso ao site negado neste ambiente) | — |
| N3 tamanho mínimo | feito e testado no Release (encolher, SetWindowPos, minimizar/restaurar) | `build/minsize.ps1` |
| N4 tipografia | feito (escala 11/12/13/14/16 + marca 18/20/44; sem 800 fora do display; sem negrito no meio de frase) | `grep -rhoE "fontSize: [0-9.]+" lib/src/ui` |
