#!/usr/bin/env bash
# Contribution status report.
#
# Answers, in one command:
#   - where we are and what is uncommitted
#   - which roadmap stages have a branch, commits, and are merged into the trunk
#   - what exists locally and was never pushed
#   - what origin has that upstream does not (pull request candidates)
#   - what upstream gained since the last fetch (conflict risk)
#   - open pull request state, when the GitHub CLI is available
#
# Usage:  bash tools/contrib/status.sh [--fetch]
#         --fetch   contact the remotes before reporting (needs network)

set -uo pipefail

TRUNK="juniormar/main"
STAGE_PREFIX="juniormar/"
UPSTREAM_MAIN="upstream/master"
ORIGIN_MAIN="origin/master"

DO_FETCH=0
[ "${1:-}" = "--fetch" ] && DO_FETCH=1

if [ -t 1 ] && [ -z "${NO_COLOR:-}" ]; then
    B=$'\033[1m'; DIM=$'\033[2m'; R=$'\033[0m'
    GRN=$'\033[32m'; YEL=$'\033[33m'; RED=$'\033[31m'; CYN=$'\033[36m'
else
    B=""; DIM=""; R=""; GRN=""; YEL=""; RED=""; CYN=""
fi

cd "$(git rev-parse --show-toplevel 2>/dev/null)" || {
    echo "not inside a git repository" >&2; exit 1
}

hr() { printf '%s\n' "────────────────────────────────────────────────────────────────────────"; }
head2() { printf '\n%s%s%s\n' "$B" "$1" "$R"; hr; }

has_ref()  { git rev-parse --verify --quiet "$1" >/dev/null 2>&1; }
count()    { git rev-list --count "$1" 2>/dev/null || echo "?"; }

printf '%s%sHigFlow - contribution status%s   %s\n' "$B" "$CYN" "$R" "$(date '+%Y-%m-%d %H:%M')"

# ── 1. working tree ────────────────────────────────────────────────────────────
head2 "1. Working tree"

CUR=$(git rev-parse --abbrev-ref HEAD)
printf 'branch        %s%s%s\n' "$B" "$CUR" "$R"
printf 'head          %s  %s\n' "$(git rev-parse --short HEAD)" "$(git log -1 --format=%s)"

DIRTY=$(git status --porcelain | wc -l | tr -d ' ')
if [ "$DIRTY" -eq 0 ]; then
    printf 'changes       %sclean%s\n' "$GRN" "$R"
else
    printf 'changes       %s%s uncommitted%s\n' "$YEL" "$DIRTY" "$R"
    git status --short | sed 's/^/              /' | head -20
    [ "$DIRTY" -gt 20 ] && printf '              %s... and %s more%s\n' "$DIM" "$((DIRTY - 20))" "$R"
fi

STASH=$(git stash list | wc -l | tr -d ' ')
[ "$STASH" -gt 0 ] && printf 'stashes       %s%s%s\n' "$YEL" "$STASH" "$R"

# ── 2. remotes ─────────────────────────────────────────────────────────────────
head2 "2. Remotes"

git remote -v | awk '$3=="(fetch)"{printf "%-10s %s\n", $1, $2}'

if [ "$DO_FETCH" -eq 1 ]; then
    printf '\nfetching ... '
    if git fetch --all --prune --quiet 2>/dev/null; then
        printf '%sok%s\n' "$GRN" "$R"
    else
        printf '%sfailed (offline?)%s\n' "$RED" "$R"
    fi
else
    printf '\n%s(run with --fetch to refresh remote state)%s\n' "$DIM" "$R"
fi

# ── 3. roadmap stages ──────────────────────────────────────────────────────────
head2 "3. Roadmap stages"

if ! has_ref "$TRUNK"; then
    printf '%strunk %s does not exist yet%s\n' "$YEL" "$TRUNK" "$R"
fi

BASE="master"
has_ref "$BASE" || BASE="$UPSTREAM_MAIN"

printf '%-28s %7s %7s %7s  %s\n' "BRANCH" "COMMITS" "PUSHED" "MERGED" "LAST COMMIT"
printf '%s\n' "$DIM────────────────────────────────────────────────────────────────────────$R"

FOUND=0
while IFS= read -r br; do
    [ -z "$br" ] && continue
    [ "$br" = "$TRUNK" ] && continue
    FOUND=1

    ahead=$(count "$BASE..$br")

    if has_ref "origin/$br"; then
        unpushed=$(count "origin/$br..$br")
        if [ "$unpushed" = "0" ]; then pushed="${GRN}yes${R}"; else pushed="${YEL}-$unpushed${R}"; fi
    else
        pushed="${DIM}no${R}"
    fi

    if has_ref "$TRUNK" && git merge-base --is-ancestor "$br" "$TRUNK" 2>/dev/null; then
        merged="${GRN}yes${R}"
    else
        merged="${DIM}no${R}"
    fi

    last=$(git log -1 --format='%ad %s' --date=short "$br" 2>/dev/null | cut -c1-42)
    printf '%-28s %7s %16s %16s  %s\n' "$br" "$ahead" "$pushed" "$merged" "$last"
done < <(git for-each-ref --format='%(refname:short)' refs/heads | grep "^${STAGE_PREFIX}" || true)

[ "$FOUND" -eq 0 ] && printf '%sno stage branches yet%s\n' "$DIM" "$R"

# trunk itself
if has_ref "$TRUNK"; then
    printf '\n%-28s %7s commits ahead of %s\n' "$TRUNK" "$(count "$BASE..$TRUNK")" "$BASE"
fi

# ── 4. not yet pushed ──────────────────────────────────────────────────────────
head2 "4. Local work not on origin"

ANY=0
while IFS= read -r br; do
    [ -z "$br" ] && continue
    if has_ref "origin/$br"; then
        n=$(count "origin/$br..$br")
        [ "$n" != "0" ] && { printf '%s%-28s%s %s commit(s) not pushed\n' "$YEL" "$br" "$R" "$n"; ANY=1; }
    else
        n=$(count "$BASE..$br")
        [ "$n" != "0" ] && { printf '%s%-28s%s branch never pushed (%s commit(s))\n' "$YEL" "$br" "$R" "$n"; ANY=1; }
    fi
done < <(git for-each-ref --format='%(refname:short)' refs/heads)

[ "$ANY" -eq 0 ] && printf '%severything local is on origin%s\n' "$GRN" "$R"

# ── 5. pull request candidates ─────────────────────────────────────────────────
head2 "5. Pull request candidates (on origin, not in upstream/master)"

ANY=0
while IFS= read -r br; do
    [ -z "$br" ] && continue
    short=${br#origin/}
    case "$short" in "$STAGE_PREFIX"*) ;; *) continue ;; esac
    [ "$short" = "$TRUNK" ] && continue
    if has_ref "$UPSTREAM_MAIN"; then
        n=$(count "$UPSTREAM_MAIN..$br")
        [ "$n" != "0" ] && { printf '%-30s %s commit(s) ahead of upstream/master\n' "$short" "$n"; ANY=1; }
    fi
done < <(git for-each-ref --format='%(refname:short)' refs/remotes/origin || true)

[ "$ANY" -eq 0 ] && printf '%snone%s\n' "$DIM" "$R"

# ── 6. upstream drift ──────────────────────────────────────────────────────────
head2 "6. Upstream drift (conflict risk)"

if has_ref "$UPSTREAM_MAIN"; then
    behind=$(count "master..$UPSTREAM_MAIN")
    if [ "$behind" = "0" ]; then
        printf 'master vs upstream/master     %sin sync%s\n' "$GRN" "$R"
    else
        printf 'master vs upstream/master     %s%s commit(s) behind%s\n' "$YEL" "$behind" "$R"
        git log --oneline "master..$UPSTREAM_MAIN" | head -8 | sed 's/^/    /'
    fi
fi

printf '\n%sactive upstream branches%s\n' "$B" "$R"
for b in Kaina PC_Daniel_Mesh Daniel PC_ImproveDocumentation Castelo; do
    if has_ref "upstream/$b"; then
        ahead=$(count "master..upstream/$b")
        info=$(git log -1 --format='%ad  %s' --date=short "upstream/$b" | cut -c1-52)
        printf '  %-26s %5s ahead   %s\n' "$b" "$ahead" "$info"
    fi
done

# ── 7. pull requests ───────────────────────────────────────────────────────────
head2 "7. Pull requests"

if command -v gh >/dev/null 2>&1; then
    if gh auth status >/dev/null 2>&1; then
        gh pr list --repo antoniocastelofilho/HigFlow --author '@me' \
            --state all --limit 20 \
            --json number,title,state,headRefName \
            --template '{{range .}}{{printf "  #%-5v %-10v %-30v %v\n" .number .state .headRefName .title}}{{end}}' \
            2>/dev/null || printf '%s  could not query (no access to upstream?)%s\n' "$DIM" "$R"
    else
        printf '%s  gh present but not authenticated - run: gh auth login%s\n' "$DIM" "$R"
    fi
else
    printf '%s  gh not installed - pull request state unavailable%s\n' "$DIM" "$R"
fi

# ── 8. repository weight ───────────────────────────────────────────────────────
head2 "8. Repository weight"

WT=$(du -sh --exclude=.git . 2>/dev/null | cut -f1)
GD=$(du -sh .git 2>/dev/null | cut -f1)
printf 'working tree  %s\n' "${WT:-?}"
printf '.git          %s\n' "${GD:-?}"
printf 'tracked files %s\n' "$(git ls-files | wc -l | tr -d ' ')"

printf '\n'
