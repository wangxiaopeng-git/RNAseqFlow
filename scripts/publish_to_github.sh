#!/usr/bin/env bash
set -euo pipefail

# Run from the extracted RNAseqFlow repository root after `gh auth login`.
OWNER="${GH_OWNER:-wangxiaopeng-git}"
REPO="${GH_REPO:-RNAseqFlow}"
VISIBILITY="${GH_VISIBILITY:-public}"
[[ "${VISIBILITY}" == "public" || "${VISIBILITY}" == "private" ]] || { echo 'GH_VISIBILITY must be public or private' >&2; exit 2; }
command -v git >/dev/null || { echo 'git is required' >&2; exit 2; }
command -v gh >/dev/null || { echo 'GitHub CLI (gh) is required: https://cli.github.com/' >&2; exit 2; }
[[ -f pyproject.toml && -f README.md && -d src/rnaseqflow ]] || { echo 'Run from the RNAseqFlow project root' >&2; exit 2; }
gh auth status >/dev/null || { echo 'Run gh auth login first' >&2; exit 2; }
if gh repo view "${OWNER}/${REPO}" >/dev/null 2>&1; then
  echo "Repository ${OWNER}/${REPO} already exists. Stopping to avoid overwriting data." >&2
  exit 2
fi
if [[ ! -d .git ]]; then
  git init -b main
fi
git add README.md LICENSE CHANGELOG.md Makefile .gitignore pyproject.toml src workflow scripts config demo docs tests envs .github
if ! git diff --cached --quiet; then
  git commit -m "Initial RNAseqFlow v0.1.0 development preview"
fi
# Requires GitHub user authorization. No repository is created by ChatGPT.
gh repo create "${OWNER}/${REPO}" "--${VISIBILITY}" --source=. --remote=origin --push
