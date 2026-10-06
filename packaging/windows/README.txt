fast-gfn1-xtb Windows x64
========================

Requirements: Windows 10/11 x64, SSE4.2-capable x64 CPU.
No oneAPI, Visual Studio, Python, Conda, WSL, administrator rights or network
access is needed to run this package. Runtime libraries are included in bin/.

Portable use:
1. Extract the whole ZIP into a folder. Keep bin/, share/ and licenses/ together.
2. In PowerShell, from this folder:
   .\fast-gfn1-xtb.cmd --version
   .\fast-gfn1-xtb.cmd .\examples\water.xyz --gfn 1 --grad --norestart
   .\fast-gfn1-xtb.cmd .\examples\water.xyz --gfn 1 --hess --alpb water --norestart
   Calculation output is written in the current working directory.

Install for your Windows user:
Double-click install.cmd. It checks package hashes, copies the files into
%LOCALAPPDATA%\Programs\fast-gfn1-xtb and adds that folder to your user PATH.
Open a NEW terminal, then use:
   fast-gfn1-xtb input.xyz --gfn 1 --grad --norestart
The installer refuses to replace an existing installation folder.

Custom destination without changing PATH:
   powershell -NoProfile -ExecutionPolicy Bypass -File .\install.ps1 -Destination C:\tools\fast-gfn1-xtb
Add -AddToPath if you also want to register that folder in your user PATH.

The launcher supplies the bundled parameter path and defaults to 8 OpenMP
threads and sequential MKL. Existing XTBPATH/OMP_NUM_THREADS settings are
respected. Example: $env:OMP_NUM_THREADS='4' in PowerShell.
Always specify --gfn 1. The inherited CLI defaults to GFN2.

Uninstall: remove the installation folder and its entry in your user PATH
(Windows "Edit environment variables for your account").

COPYING and licenses/ contain the applicable notices. source/ contains the
corresponding project and pinned mctc-lib source, plus build instructions.
This is a modified xTB 6.7.1 distribution (GFN1-fast 4.7.0).
