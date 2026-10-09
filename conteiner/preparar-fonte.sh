#!/usr/bin/env bash
#==============================================================================
# preparar-fonte.sh -- junta o que a definicao Apptainer precisa em %files.
#
# POR QUE NAO `COPY .` / `%files .`:  a arvore de desenvolvimento tem 4,7 GB de
# dados de corrida em higflow/ve, 2,8 GB de bibliotecas compiladas e 259 MB de
# .git.  `git archive` de um commit leva apenas o que esta' VERSIONADO naquele
# commit -- sem lista de exclusao para manter, e com a garantia de que a imagem
# corresponde a um ponto identificavel da historia.
#
# Uso:  conteiner/preparar-fonte.sh [commit]     (padrao: HEAD)
#==============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${1:-HEAD}"
OUT="$ROOT/conteiner/_build"
mkdir -p "$OUT"

cd "$ROOT"
SHA=$(git rev-parse --short "$REF")
DESC=$(git describe --always --dirty "$REF" 2>/dev/null || echo "$SHA")

# DUAS recusas, e a segunda eu descobri falhando.
#
# (a) arquivo rastreado MODIFICADO: a imagem nao corresponderia a commit nenhum,
#     e a diferenca nao apareceria em lugar algum depois.
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo ">> ERRO: ha' arquivos rastreados modificados." >&2
  git status --short --untracked-files=no | sed 's/^/     /' >&2
  echo ">> Commite ou guarde antes: a imagem precisa corresponder a um commit." >&2
  exit 1
fi
# (b) arquivo NAO RASTREADO e nao ignorado: `git archive` simplesmente NAO o
#     inclui, em silencio.  Na primeira execucao disto os proprios arquivos da
#     definicao Apptainer estavam nao rastreados, e o tarball saiu sem eles --
#     a imagem teria sido construida sem a receita que a descreve.  Verificar
#     apenas "modificados" nao pega isso, porque nao ha' nada modificado.
#     A lista e' CLASSIFICADA, senao a guarda recusa sempre e vira ruido: a
#     arvore tem dezenas de .amr e saidas GERADAS que nao sao versionadas de
#     proposito.  Recusa por FONTE nao rastreada; apenas conta o resto.
NAORASTR=$(git ls-files --others --exclude-standard)
FONTE=$(echo "$NAORASTR" | grep -E '\.(c|h|cpp|hpp|def|sh|py|mk)$|(^|/)Makefile$|(^|/)varsrc$' || true)
RESTO=$(echo "$NAORASTR" | grep -vE '\.(c|h|cpp|hpp|def|sh|py|mk)$|(^|/)Makefile$|(^|/)varsrc$' | grep -c . || true)
if [ -n "$FONTE" ]; then
  echo ">> ERRO: FONTE nao rastreada seria omitida do tarball em silencio:" >&2
  echo "$FONTE" | sed 's/^/     /' >&2
  echo ">> Adicione ao commit, ou ignore explicitamente no .gitignore." >&2
  echo ">> (para seguir ciente da omissao: PERMITIR_NAO_RASTREADOS=1)" >&2
  [ "${PERMITIR_NAO_RASTREADOS:-0}" = "1" ] || exit 1
  echo ">> AVISO: seguindo com omissao deliberada de fonte." >&2
fi
[ "${RESTO:-0}" -gt 0 ] && \
  echo ">> nota: $RESTO arquivo(s) nao rastreado(s) que NAO sao fonte (malhas e saidas geradas); nao entram no tarball, como deve ser."

echo ">> fonte do HiGFlow: $DESC ($SHA)"
git archive --format=tar.gz --prefix='' -o "$OUT/higflow-src.tar.gz" "$REF"

# O tarball do PETSc ja' esta' em arvore mas NAO e' versionado (e' binario de
# 16 MB); copia-se do diretorio de bibliotecas.
PETSC_TGZ="$ROOT/bibliotecas/petsc-3.25.4.tar.gz"
ESPERADO=12c990fb39a5764ac8311211d09c01ed80fb983136c75bf7b558312b2509dbbd
if [ ! -f "$PETSC_TGZ" ]; then
  echo ">> tarball do PETSc ausente; baixando"
  curl -fsSL -o "$PETSC_TGZ" \
    "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.25.4.tar.gz"
fi
echo "$ESPERADO  $PETSC_TGZ" | sha256sum -c - || {
  echo ">> ERRO: soma de verificacao do PETSc nao bate." >&2; exit 1; }
cp -f "$PETSC_TGZ" "$OUT/petsc-3.25.4.tar.gz"

# GUARDA DE VACUIDADE do tarball: um archive vazio ou truncado gera uma imagem
# que falha tres horas adiante, no `make`, com erro que nao aponta para aqui.
# A listagem sai UMA vez para um arquivo, e as buscas sao nele.
# `tar tzf ... | grep -q` sob `set -o pipefail` da' FALSO NEGATIVO: o grep -q sai
# no primeiro casamento, o tar morre de SIGPIPE (141), e o pipefail transforma o
# casamento em falha.  Foi assim que esta guarda acusou a ausencia de
# higtree/Makefile num tarball que tinha 883 entradas de higtree/.
LISTA="$OUT/.conteudo.txt"
tar tzf "$OUT/higflow-src.tar.gz" > "$LISTA"
for esperado in higtree/Makefile higflow/Makefile conteiner/higflow-euler.def \
                higflow/tests/front-tracking varsrc; do
  if ! grep -q "^$esperado" "$LISTA"; then
    echo ">> ERRO: o tarball nao contem '$esperado'." >&2
    echo ">> Verifique se o caminho esta' versionado no commit $SHA." >&2
    exit 1
  fi
done
echo ">> tarball conferido: $(wc -l < "$LISTA") entradas"
rm -f "$LISTA"

# A procedencia vai DENTRO da imagem: sem isso, descobrir meses depois de que
# commit saiu um .sif exige adivinhacao.
cat > "$OUT/PROCEDENCIA" <<EOF
commit    $(git rev-parse "$REF")
describe  $DESC
ramo      $(git rev-parse --abbrev-ref HEAD)
data      $(git show -s --format=%cI "$REF")
preparado $(date -Iseconds)
maquina   $(uname -srm)
EOF

echo ">> pronto em $OUT:"
ls -la "$OUT" | awk 'NR>3 {printf "     %-28s %7.1f MB\n",$NF,$5/1048576}'
echo ">> construir:  apptainer build --fakeroot higflow-euler.sif conteiner/higflow-euler.def"
