#!/bin/bash
# ============================================================
# push_to_github.sh
# Run this script from YOUR terminal once you have a new token.
# It will create the GitHub repo and push all project files.
#
# Usage:
#   export GITHUB_TOKEN=<your_token_here>
#   bash push_to_github.sh
# ============================================================

set -e

# ── Configuration ──────────────────────────────────────────
USERNAME="Ensiuc"
REPO="CpG_Methylation_ML"

# Token must be set as an environment variable — never hardcode it
if [ -z "$GITHUB_TOKEN" ]; then
  echo "ERROR: GITHUB_TOKEN is not set."
  echo "Run: export GITHUB_TOKEN=ghp_yourtoken"
  echo "Then re-run this script."
  exit 1
fi

# ── Create GitHub repository ───────────────────────────────
echo ">>> Creating GitHub repository: $USERNAME/$REPO ..."
curl -s -X POST \
  -H "Authorization: token $GITHUB_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  https://api.github.com/user/repos \
  -d "{
    \"name\": \"$REPO\",
    \"description\": \"ML pipeline for CpG methylation analysis across a shared genomic region in multiple early-life exposure datasets\",
    \"private\": false,
    \"auto_init\": false,
    \"has_wiki\": true,
    \"has_issues\": true
  }"

echo ""
echo ">>> Initializing local git repo..."
cd "$(dirname "$0")"

git init
git config user.name "$USERNAME"
git config user.email "$USERNAME@users.noreply.github.com"

git add .
git commit -m "Initial commit: project structure, ML tools overview, pipeline, GitHub Pages"

echo ""
echo ">>> Pushing to GitHub..."
git branch -M main
git remote add origin "https://$GITHUB_TOKEN@github.com/$USERNAME/$REPO.git"
git push -u origin main

echo ""
echo ">>> Enabling GitHub Pages (docs/ folder)..."
curl -s -X POST \
  -H "Authorization: token $GITHUB_TOKEN" \
  -H "Accept: application/vnd.github.v3+json" \
  "https://api.github.com/repos/$USERNAME/$REPO/pages" \
  -d '{"source": {"branch": "main", "path": "/docs"}}'

echo ""
echo "============================================================"
echo "Done! Your project is live at:"
echo "  Repo:  https://github.com/$USERNAME/$REPO"
echo "  Page:  https://$USERNAME.github.io/$REPO"
echo "  (GitHub Pages may take 1-2 minutes to deploy)"
echo "============================================================"
