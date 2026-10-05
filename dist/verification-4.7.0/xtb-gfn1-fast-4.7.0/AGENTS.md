# AGENTS.md

このファイルは、Codex/自動エージェントが Fast_GFN1_xTB リポジトリで作業するときの
Python 実行環境、標準コマンド、サブエージェント運用を共有するためのメモです。

## 1. 実行環境（前提）

- OS / shell: Windows + PowerShell
- リポジトリルート: `C:\GitHub_repository\Fast_GFN1_xTB`
- Codex Agent environment: Windows native
- Python環境: Mambaforge / Conda `fastxtb` (`C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe`)
- Python標準バージョン: 3.13.15

## 2. Pythonコマンドの扱い

- Python は Windows版 Mambaforge / Conda の仮想環境 `fastxtb` を使用する。
- Python バージョンは `3.13.15` を標準とする。
- 実行ファイルは `C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe` を使う。
- Codex / 自動エージェントは `conda activate fastxtb` の成否や現在の shell activation 状態に依存してはならない。
- Python を実行するときは、可能な限り上記 `python.exe` を直接指定する。
- PowerShell では次の変数を標準として使う。
  - `$PY = "C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe"`
- pip は単独の `pip` コマンドではなく、必ず `& $PY -m pip ...` の形で実行する。
- pytest などの Python モジュールも、必ず `& $PY -m <module>` の形で実行する。
- 新しいスクリプトで子プロセス Python を起動する場合は、`sys.executable` を優先する。
- 使用中の Python を確認するときは `& $PY -c "import sys; print(sys.executable)"` を使い、上記パスと一致することを確認する。
- system Python、Microsoft Store Python、WSL Python、WinPython、他の Conda 環境を Fast_GFN1_xTB の開発・テストに使用しない。

## 2.5 環境確認

作業開始時、Python実行や依存関係の変更を伴う場合は、必要に応じて次を確認する。

```powershell
$PY = "C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe"
& $PY --version
& $PY -c "import sys; print(sys.executable)"
```

期待する実行ファイル:

```text
C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe
```

