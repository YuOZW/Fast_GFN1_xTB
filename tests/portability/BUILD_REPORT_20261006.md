# fast-gfn1-xtb ビルド・実行ファイル名の検証（2026-10-06）

対象: 現在のGFN1-fast 4.7.0ソースのREADME更新と実行ファイル名変更。
配布済みZIP・保存済みベンチマークは再作成していない。

## ビルド

- Windows native: Intel ifx 2025.2、MSVC 19.44、MKL 2025.2 sequential、CMake 4.1.1、Ninja 1.12.1。
- Windows Release: READMEのCMakeコマンドを新規 `build/windows-release` に実行し成功。
- Windows Release/Debug: 更新した `run_windows_build.ps1` でビルド・source guards・policy検証・Energy/Gradientスモークテスト成功。
- Linux: WSL Ubuntu 22.04のGNU Fortran 11.4.0、CMake 3.22.1、システムBLAS/LAPACKでREADMEのRelease構成が成功。WSL内のPythonは使用していない。
- 新規Windows/Linuxビルドに旧名 `xtb(.exe)` は生成されない。
- Windows/Linuxとも `--help` に `Usage: fast-gfn1-xtb` が表示され、`--version` が正常終了。
- CMakeの全ターゲットをビルド後、ワークスペース内の別prefixへのインストールとインストール済み実行ファイルの起動が成功。
- Mesonの実行ファイル宣言も新名に変更したが、Mesonによるコンパイルは今回実施していない。

## 数値確認

Windowsの既存CLI検証で、OpenMPの実ワーカー・直列切り替え・Hessian一致および解析的Hessianの既定/fallback方針を確認。命名変更後のCMake新規Releaseビルドでも、同じ入力・制御ファイルをWindows/Linuxで実行して比較した。

| ケース | OpenMPスレッド数 | Energy差 (Eh) | Hessian最大差 (Eh/bohr²) |
| --- | ---: | ---: | ---: |
| gas_serial | 1 | 0.000e+00 | 0.000e+00 |
| gas_threads2 | 2 | 0.000e+00 | 0.000e+00 |
| alpb_threads2 | 2 | 0.000e+00 | 0.000e+00 |

差はファイルの出力精度内での比較。Linuxの気相1/2スレッド間でもHessianの出力値が一致した。

## Linuxで発見した既存不具合

水分子のHessian計算後の振動解析で、Cの対称性判定fallbackが `MaxRotAxis[2]` に `C2 ` を連結し、GNUのfortificationでSIGABRTになった。返却バッファをFortran側の6文字バッファに合わせ、最大の回転軸を境界付き書き込みで一つだけ返すように修正した。Energy/Gradient/Hessianの計算処理は変更していない。

専用回帰検証:

```bash
gcc -O2 -D_FORTIFY_SOURCE=2 tests/portability/test_symmetry_fallback.c -lm -o build/test_symmetry_fallback
build/test_symmetry_fallback
```

C4/C2が併存する場合の最高軸選択、C20、軸なしの各ケースが成功。修正後はLinuxの気相・ALPB Hessianが振動解析を含め正常終了した。WindowsのDebugスモークテストも再ビルド後に成功。

実行ログと比較値はローカルの `build/` およびWindows検証ビルドディレクトリに保存されている。
