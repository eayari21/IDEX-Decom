#!/usr/bin/env bash
set -euo pipefail

# Convert these nested git repos/submodules into normal directories.
# (Quotes matter because of the '->' in folder names.)
SUBS=(
  "L0->L1A/idex-decom"
  "L1A->L1B/idex-decom"
  "L1B->L2A/idex-decom"
  "L2A->L2B/idex-decom"
  "L2B->L2C/Mapping/imap_processing"
)

# Ensure we're inside a git work tree and cd to repo root
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: not inside a git repository." >&2
  exit 1
fi
repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

echo "Repo root: $repo_root"
echo "Vendoring nested repos as normal directories…"

for p in "${SUBS[@]}"; do
  if [[ ! -e "$p" ]]; then
    echo "  - Skipping missing: $p"
    continue
  fi
  echo "  - Processing: $p"

  # 1) Deinit submodule entry (if it exists)
  git submodule deinit -f "$p" 2>/dev/null || true

  # 2) Remove from the index but keep the files on disk
  git rm -r --cached "$p" 2>/dev/null || true

  # 3) Remove submodule metadata
  rm -rf ".git/modules/$p" 2>/dev/null || true

  # 4) Remove the inner repo metadata so it's just a normal folder
  rm -rf "$p/.git" 2>/dev/null || true
done

# Clean up .gitmodules if needed (deinit usually edits it; remove if empty)
if [[ -f .gitmodules ]]; then
  if ! git config -f .gitmodules -l >/dev/null 2>&1; then
    rm -f .gitmodules
  fi
  git add -f .gitmodules 2>/dev/null || true
fi

# (Optional) Ignore known large archives to avoid GitHub 100MB limit
# Add patterns only if they aren't present already.
touch .gitignore
grep -qE '(^|/)\*.zip($|\s)' .gitignore || printf '\n# Ignore archives\n*.zip\n' >> .gitignore
grep -Fq "L2A->L2B/idex-decom.zip" .gitignore || printf 'L2A->L2B/idex-decom.zip\n' >> .gitignore
git add .gitignore 2>/dev/null || true

# Stage everything and commit if there are changes
git add -A
if ! git diff --cached --quiet; then
  git commit -m "Vendor nested repos as normal directories; remove submodule metadata"
else
  echo "No changes to commit."
fi

# Push to the current branch
current_branch="$(git rev-parse --abbrev-ref HEAD)"
echo "Pushing to origin/${current_branch}…"
git push -u origin "${current_branch}"

echo "Done."

