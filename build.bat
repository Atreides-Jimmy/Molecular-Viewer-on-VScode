@echo off
chcp 65001 >nul 2>&1
echo ========================================
echo   Molecular Viewer - Build Script
echo ========================================
echo.

cd /d "%~dp0"

where node >nul 2>&1
if errorlevel 1 (
    echo ERROR: Node.js is not installed or not in PATH.
    echo Please install Node.js from https://nodejs.org/ (LTS version)
    echo.
    echo Alternatively, the dist/ folder already contains pre-compiled JS files.
    echo You can package directly with: npx vsce package --no-dependencies
    pause
    exit /b 1
)

if not exist node_modules (
    echo [1/3] Installing dependencies...
    call npm install
    if errorlevel 1 (
        echo ERROR: npm install failed.
        pause
        exit /b 1
    )
) else (
    echo [1/3] Dependencies already installed.
)

echo.
echo [2/3] Compiling TypeScript...
call npx tsc -p ./ --pretty
if errorlevel 1 (
    echo ERROR: TypeScript compilation failed.
    pause
    exit /b 1
)

echo.
echo [3/3] Packaging extension...
if exist molecular-viewer-3d*.vsix del /q molecular-viewer-3d*.vsix
call npx vsce package --no-dependencies
if errorlevel 1 (
    echo.
    echo WARNING: vsce packaging failed.
    echo Try: npm install -g @vscode/vsce
    echo Then run: vsce package --no-dependencies
    echo.
    echo The dist/ folder is ready for manual deployment.
) else (
    echo.
    echo ========================================
    echo   BUILD SUCCESS!
    echo ========================================
    echo Output: molecular-viewer-3d-*.vsix
    echo.
    echo Install in Trae/VS Code:
    echo   1. Ctrl+Shift+P
    echo   2. Extensions: Install from VSIX...
    echo   3. Select molecular-viewer-3d-*.vsix
    echo   4. For Remote-SSH: Install on Remote
    echo ========================================
)

pause
