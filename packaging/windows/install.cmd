@echo off
"%SystemRoot%\System32\WindowsPowerShell\v1.0\powershell.exe" -NoProfile -ExecutionPolicy Bypass -File "%~dp0install.ps1" -AddToPath
if errorlevel 1 (echo Installation failed.) else (echo Installation complete. Open a new terminal to use fast-gfn1-xtb.)
pause
