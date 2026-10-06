@echo off
setlocal
if not defined XTBPATH set "XTBPATH=%~dp0share\xtb"
if not defined OMP_NUM_THREADS set "OMP_NUM_THREADS=8"
if not defined MKL_NUM_THREADS set "MKL_NUM_THREADS=1"
if not defined OMP_STACKSIZE set "OMP_STACKSIZE=64M"
if not defined KMP_STACKSIZE set "KMP_STACKSIZE=64M"
"%~dp0bin\fast-gfn1-xtb.exe" %*
exit /b %errorlevel%
