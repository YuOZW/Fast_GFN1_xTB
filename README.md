# fast-gfn1-xtb

`fast-gfn1-xtb` は、[xTB 6.7.1](https://github.com/grimme-lab/xtb) を基にした、GFN1-xTBのCPU高速化版です。現在の実装は **GFN1-fast 4.7.0**。エネルギー・勾配・解析的Hessianを対象とし、数値精度を維持しながら不要なproperty計算の省略、行列演算の整理、OpenMPによる並列化を行っています。

生成される実行ファイルは、Windowsでは **`fast-gfn1-xtb.exe`**、Linuxでは **`fast-gfn1-xtb`** です。計算時には **`--gfn 1` を明示してください**。xTB由来のCLIの既定値はGFN2のままです。

## 現在の機能

- SCC: SYEVDを既定の固有値ソルバーとし、中間反復ではcompact densityを使用。最終確認反復ではfull densityを構築します。
- Gradient: GFN1専用経路、積分の再利用、不要なproperty計算の省略。
- Hessian: 気相・GBSA・ALPBの解析的Hessianを既定で使用。未対応条件や応答計算の失敗時には数値Hessianへfallbackします。
- OpenMP: Hessianの主要ループを並列化。行列演算のまとめ処理、疎なBorn/SASA導関数処理も実装しています。

通常版xTBと同じproperty出力を保証する用途ではありません。Wiberg/Mayer bond orderなど、対象計算に不要な処理は省略します。詳細な検証条件と過去の測定結果は [4.7.0の資料](tests/gfn1_fast_4_7_0/README.md) を参照してください。保存済みの旧配布物・測定ログには旧実行ファイル名 `xtb` が含まれます。

## ソースと依存ライブラリ

```text
git clone https://github.com/YuOZW/Fast_GFN1_xTB.git
cd Fast_GFN1_xTB
```

必要なライブラリはBLAS/LAPACKとmctc-libです。mctc-libは、検証済みのコミット **`77f65c6f2cf6330d05d0757ca173da097096780e`** を使用します。ソース配布ZIPにはこの依存ソースを同梱しています。以下の構成ではtblite、CPCMX、上流の単体テスト、JSONを無効にします。GBSA/ALPBは利用できます。

## Windowsでのコンパイル

必要なもの:

- Visual Studio 2022のC++ Build ToolsとWindows SDK（`cl`）
- Intel oneAPIのFortranコンパイラー `ifx` とMKL
- CMake 3.17以降、Ninja、Git、PowerShell

Intel oneAPIのIntel 64用コマンドプロンプトで `powershell -NoProfile` を起動し、リポジトリのルートへ移動してください。コンパイラーと実行時DLLのPATHを引き継ぎます。

mctc-libがまだない場合は取得します。

```powershell
if (!(Test-Path 'subprojects/mctc-lib/CMakeLists.txt')) {
    git clone https://github.com/grimme-lab/mctc-lib.git subprojects/mctc-lib
    git -C subprojects/mctc-lib checkout --detach 77f65c6f2cf6330d05d0757ca173da097096780e
}
```

Releaseビルド:

```powershell
cmake -S . -B build/windows-release -G Ninja `
    -DCMAKE_BUILD_TYPE=Release `
    -DCMAKE_C_COMPILER=cl -DCMAKE_Fortran_COMPILER=ifx `
    -DMCTCLIB_FIND_METHOD=subproject `
    -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_TESTS=OFF -DWITH_JSON=OFF `
    -DWITH_OBJECT=OFF -DWITH_OpenMP=ON `
    -DBLA_VENDOR=Intel10_64lp_seq `
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF `
    '-DCMAKE_Fortran_FLAGS_RELEASE=/O3 /DNDEBUG' `
    '-DCMAKE_EXE_LINKER_FLAGS=/Qoption,link,/STACK:67108864'
cmake --build build/windows-release --target fast-gfn1-xtb --parallel 4
& ./build/windows-release/fast-gfn1-xtb.exe --version
```

Debugビルドは、ビルド先を `build/windows-debug`、`CMAKE_BUILD_TYPE` を `Debug` に変え、`CMAKE_Fortran_FLAGS_RELEASE` の行を省略します。大きな局所配列に対応するため、Windowsでは64 MiBのスタックを指定しています。

通常のPowerShellから、oneAPI環境の取り込みとスモークテストをまとめて行う既存スクリプトも利用できます。こちらは開発用Python環境 `C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe` が必要です。

```powershell
./tests/gfn1_fast_2_2_6/run_windows_build.ps1 -BuildType Release
./tests/gfn1_fast_2_2_6/run_windows_build.ps1 -BuildType Debug
./tests/gfn1_fast_4_7_0/run_windows_parallel.ps1
./tests/gfn1_fast_4_6_0/run_windows_default.ps1
```

スクリプトのフォルダ名は履歴を引き継いでいますが、ビルドするのは現在のソースです。既定の出力先は `build-gfn1-fast-current-windows-ifx-release/fast-gfn1-xtb.exe` とDebug側の対応ディレクトリです。実行時もoneAPIのDLLがPATHに必要です。

## Linuxでのコンパイル

GNU FortranとシステムBLAS/LAPACKを使用する構成です。Ubuntu/Debianで必要なパッケージを用意する例:

```bash
sudo apt update
sudo apt install build-essential gfortran cmake ninja-build git libblas-dev liblapack-dev
```

リポジトリのルートで、mctc-libを準備してReleaseビルドします。

```bash
bash tests/gfn1_fast_2_2_6/prepare_debug_deps.sh
cmake -S . -B build/linux-release -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran \
    -DMCTCLIB_FIND_METHOD=subproject \
    -DWITH_TBLITE=OFF -DWITH_CPCMX=OFF -DWITH_TESTS=OFF -DWITH_JSON=OFF \
    -DWITH_OBJECT=OFF -DWITH_OpenMP=ON \
    -DBLA_VENDOR=Generic \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF
cmake --build build/linux-release --target fast-gfn1-xtb --parallel 4
./build/linux-release/fast-gfn1-xtb --version
```

Debugビルドはビルド先を `build/linux-debug`、`CMAKE_BUILD_TYPE` を `Debug` に変えます。上記はシステムBLASを使用するため、WindowsのMKL構成と性能は同じ条件になりません。

## 実行例とスレッド設定

パラメータファイルをソースから読み込む場合は、`XTBPATH` にリポジトリのルートを設定します。出力ファイルが作られるため、計算ごとに作業ディレクトリを分けると便利です。

Windows:

```powershell
$env:XTBPATH = (Get-Location).Path
$env:OMP_NUM_THREADS = '8'
$env:MKL_NUM_THREADS = '1'
& ./build/windows-release/fast-gfn1-xtb.exe input.xyz --gfn 1 --norestart
& ./build/windows-release/fast-gfn1-xtb.exe input.xyz --gfn 1 --grad --norestart
& ./build/windows-release/fast-gfn1-xtb.exe input.xyz --gfn 1 --hess --norestart
& ./build/windows-release/fast-gfn1-xtb.exe input.xyz --gfn 1 --hess --alpb water --norestart
```

Linux:

```bash
export XTBPATH="$PWD"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
./build/linux-release/fast-gfn1-xtb input.xyz --gfn 1 --grad --norestart
./build/linux-release/fast-gfn1-xtb input.xyz --gfn 1 --hess --gbsa water --norestart
```

`OMP_NUM_THREADS` は外側のOpenMP処理に使用します。BLASは1スレッドを基準にして二重並列を避けます。上記Windows構成は **sequential MKL** にリンクするため、`MKL_NUM_THREADS` を増やしてもMKL内部は並列化されません。スレッド版MKLの性能は別途検証が必要です。

Hessian用OpenMPは、既定では64 AO以上で使用し、最大8スレッド、作業配列の上限1 GiBで制限します。環境変数で調整できます。

| 環境変数 | 用途 |
| --- | --- |
| `XTB_GFN1_FAST_PROFILE=1` | 高速化経路のプロファイルを表示 |
| `XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN=1` | 数値Hessianを使用 |
| `XTB_GFN1_FAST_DISABLE_HESSIAN_OPENMP=1` | HessianのOpenMP経路を無効化 |
| `XTB_GFN1_FAST_HESSIAN_OMP_MIN_NAO` | OpenMP開始AO数 |
| `XTB_GFN1_FAST_HESSIAN_OMP_MAX_THREADS` | Hessianで使用する最大スレッド数 |
| `XTB_GFN1_FAST_ANALYTIC_HESSIAN_MAX_MIB` | 解析的Hessianの作業配列メモリ上限 |

## インストール

Linuxでユーザー領域にインストールする例:

```bash
cmake --build build/linux-release --parallel 4
cmake --install build/linux-release --prefix "$HOME/.local"
export PATH="$HOME/.local/bin:$PATH"
export XTBPATH="$HOME/.local/share/xtb"
fast-gfn1-xtb --version
```

インストール前には、上記のようにターゲットを限定せずビルドし、依存ライブラリの補助プログラムも生成してください。

Windowsでは次のようにインストールできます。

```powershell
cmake --build build/windows-release --parallel 4
cmake --install build/windows-release --prefix C:/fast-gfn1-xtb
```

実行ファイルは `bin/fast-gfn1-xtb.exe`、パラメータは `share/xtb` に入ります。oneAPIの実行時DLL環境は引き続き必要です。

CMakeの実行ファイル用ターゲット名とMesonの実行ファイル名も `fast-gfn1-xtb` です。ライブラリ名、API、パラメータのファイル名はxTB由来の名称を維持しています。

## ライセンスと引用

本リポジトリはxTBを改変した派生版です。元の著作権表示と [LGPL-3.0-or-later](COPYING) を引き継ぎます。GFN1-xTBを使った研究では、[GFN1-xTBの論文](https://doi.org/10.1021/acs.jctc.7b00118) と [xTBの引用資料](assets/references.bib) を参照してください。
