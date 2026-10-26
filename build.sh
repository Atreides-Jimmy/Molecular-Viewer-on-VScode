#!/bin/bash
set -e

echo "========================================"
echo "  Molecular Viewer - Build Script"
echo "========================================"
echo ""

cd "$(dirname "$0")"

if ! command -v node &> /dev/null; then
    echo "ERROR: Node.js is not installed."
    echo "Please install Node.js from https://nodejs.org/ (LTS version)"
    echo ""
    echo "The dist/ folder already contains pre-compiled JS files."
    echo "You can package directly with: npx vsce package --no-dependencies"
    exit 1
fi

if [ ! -d "node_modules" ]; then
    echo "[1/3] Installing dependencies..."
    npm install
else
    echo "[1/3] Dependencies already installed."
fi

echo ""
echo "[2/3] Compiling TypeScript..."
npx tsc -p ./ --pretty

echo ""
echo "[3/3] Packaging extension..."
rm -f molecular-viewer-3d*.vsix
npx vsce package --no-dependencies 2>/dev/null || {
    echo ""
    echo "WARNING: vsce packaging failed."
    echo "Try: npm install -g @vscode/vsce"
    echo "Then run: vsce package --no-dependencies"
    echo ""
    echo "The dist/ folder is ready for manual deployment."
    exit 0
}

echo ""
echo "========================================"
echo "  BUILD SUCCESS!"
echo "========================================"
echo "Output: molecular-viewer-3d-*.vsix"
echo ""
echo "Install in VS Code / Trae:"
echo "  1. Ctrl+Shift+P"
echo "  2. Extensions: Install from VSIX..."
echo "  3. Select molecular-viewer-3d-*.vsix"
echo "  4. For Remote-SSH: Install on Remote"
echo "========================================"
