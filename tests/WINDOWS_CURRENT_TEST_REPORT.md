# 最新作業ツリーのWindowsビルド・実行確認

## 4.7.0：解析Hessian並列化・最適化（2026-10-05）

主要な解析応答と積分／SASA微分をOpenMPで並列化。Windows nativeの
Release／Debugで実workerと全行列parity、非線形SCC Gradient-FD、fallback、
入れ子API、workerごとのscratchとメモリ制御を検証した。
気相／GBSA／ALPBのtaxol 5回中央値は8 threadsで4.27／4.51／4.55秒。
固定した4.6.0・1 threadの30.83／35.11／33.86秒に対し7.22／7.78／7.44倍。
直列の最適化と並列化の寄与、完全な試行値、source／binary SHAは
`gfn1_fast_4_7_0/RESULTS.md` と同フォルダのJSONに記録した。
配布ZIPの958 source／124依存ファイルを照合し、再展開から独立した
Release／Debugビルド、各19 parallel CLI／21 default／38 native FD条件もPASS。
最終記録は `../dist/REEXTRACTION_4.7.0.json`。配布物は作成時のsnapshotで、
再展開後の最終検証記録をZIP外に保存している。
以下の4.6.0以前は当時のビルド・検証履歴である。

## 最終4.6.0のWindows native検証（2026-10-05）

対応GFN1の解析Hessianを既定で有効化。SCF-screeningした非ゼロoverlapの
production Gradient整合性、D3／CN／repulsionのcutoff窓、halogenのcutoff／
近傍分岐を確認し、拒否時は数値微分へ戻す。ENABLE=0／DISABLE=1は維持。

- Release／Debugビルド・実行、既存Energy／Gradient smokeがPASS。
- 既定／fallbackの実CLI 21条件が両構成でPASS。解析11、fallback 9、thread parity 1。
- 非線形SCCの溶媒分子38条件を両構成で再検証し、全座標の未投影HessianがPASS。
- 実際のscreening点とD3 cutoff／有限差分窓のcomponent検証も両構成でPASS。
- 気相／溶媒CLI、電子／Coulomb応答、固定bare-charge溶媒応答、部分列の加算APIもPASS。
- 4.6.0のtaxol GBSA／ALPB CLI監査は1 GiB既定上限内で解析を使用。保存した
  4.5.1の全339×339行列との差は0、標準数値との差は7.2e-9／8.4e-9 Eh/bohr²。
- 4.6.0の7回反復・profileなし・1スレッド比較は、標準数値Hessianに対し
  気相ジシラン4.412倍、ベンゼン8.020倍、溶媒小分子3.788～4.331倍。
  標準側は既知のWindows ifx割当descriptor不具合を避けるためHessianの
  OpenMP指示だけを除去した比較用ビルド。複数スレッドの標準版比較ではない。

| 最終構成 | SHA-256 |
|---|---|
| release | `3f042588d183c64af8ab268ce8b6cbb27e40b91fd481142a237ed06a53916587` |
| debug | `e990092088820af0e369d80b32de85ec04e591a1a513d5e0cc7b5efc55c21645` |

測定済みバイナリのhashは各benchmark JSONに固定。上表の最終ビルドはpolicy
コメントだけを明確化した後の再ビルドで、最終CLI 21条件を再実行した。
ソース・実行ファイルを `tests/gfn1_fast_4_6_0/SOURCE_MANIFEST_20261005.json`
で固定した。配布ZIP再展開の結果はZIP外の `dist/REEXTRACTION_4.6.0.json`
に記録する。以下の4.5.1以前の記録は各当時のビルド・hashの履歴である。

## 今回のWindows native再確認（2026-10-05 13:02 JST）

現在の未コミット変更を含む4.5.1作業ツリーを、Release／Debugともに
`cmake --build ... --target xtb-exe --clean-first --parallel 4` でクリーンビルドした。
両構成で786ステップが完了し、コンパイル・リンク・実行に合格した。
今回の作業では計算ソースを変更していない。WSLは使用していない。

- 環境: Windows x64、Intel ifx 2025.2.0、MSVC 19.44.35217、oneMKL LP64 sequential、CMake 4.1.1、Ninja 1.12.1。
- Python: AGENTS.md指定のConda `fastxtb`、3.13.15を実行ファイルの絶対パスで使用。
- Release: 水／taxol 350 AOの気相・ALPBと従来density・汎用Gradient比較を含む7ケースPASS。
- Debug: `/check:all /fpe:0`付きで水の気相・ALPB、2ケースPASS。
- 両構成でsource guards、compact-density数学検証、コンパイルしたFortran dispatcher検証がPASS。
- 従来経路とのEnergy差は表示桁で0、Gradient最大差は `2.00534e-15 Eh/bohr`。
- 水のRelease／Debug間のEnergy差は0、Gradient最大差は `1.00614e-15 Eh/bohr`。
- GBSA／ALPBの実CLI Hessian検証も両構成で各15条件PASS。解析10条件の数値微分との差は最大 `1.94e-8 Eh/bohr²`。分岐近傍・メモリ制限・旧Gradient・塩のfallback 5条件は数値経路と一致。
- 実行ファイルはAMD64／PE32+、stack reserveは64 MiB。

| 構成 | 今回生成した実行ファイルのSHA-256 |
|---|---|
| Release | `bcb61b376e0612704901e51f44319db458b1656b4e077685ec29e4d40080048e` |
| Debug | `c767c06d8f65b729189420cbfe096d49e13eba5102200a0f499c9d95544e853f` |

クリーンビルドログは各buildフォルダの `native-recheck-20261005-clean-build.log`。
Energy／Gradientは `smoke/summary.json` と各実行の `run.log`、Hessianは
`solvent-cli-test/summary.json` と各実行のログに保存した。バイナリ形式・hash・
数値比較・現ソースhashをReleaseフォルダの `native-recheck-20261005.json` に記録した。
同じフォルダの以前の `benchmark_summary.json` は今回の測定結果ではない。

再確認はリポジトリルートのPowerShellから以下を実行する。
ビルドスクリプトがoneAPI環境と実行時DLLの探索パスを設定する。

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
.\tests\gfn1_fast_4_5_0\run_windows_regression.ps1 -BuildType Release
.\tests\gfn1_fast_4_5_0\run_windows_regression.ps1 -BuildType Debug
```

今回はWindowsでのビルド・動作確認で、速度ベンチマークではない。
前回の594 AO座標入力がないため、WSLでの同一計算の再現と速度比較は未確認。
解析Hessianの既定は引き続き無効。大分子溶媒Hessianや追加の開発目標全体を
今回の小分子CLI検証だけで完了したとは扱わない。
その後の同一ソースでの大分子溶媒全339列・CLIメモリ監査と小分子thread
crossoverは `gfn1_fast_4_5_1/RESULTS.md` に追記した。標準版との大分子溶媒
反復比較は実行中で、速度結果はまだ確定していない。
以下のhashと結果は各時点の履歴として保存する。

## 溶媒SCC Hessian接続後の検証（2026-10-05）

4.5.0の溶媒応答をSCC・calculator／CLIの解析Hessianへ接続し、Windows native
Release／Debugでビルド・実行を確認した。通常Energy／Gradientの値と経路間比較は
従来どおりPASS。GBSA／ALPB、Still／P16、零温度／有限温度／開殻を含む分子38条件、
全座標の未投影Hessianと部分加算API、CLIの解析10条件・fallback 5条件を両構成でPASS。
1e-4 bohrのGradient有限差分との差は最大1.172e-8 Eh/bohr²。

- 当時のRelease SHA-256: `193cae00930589ac20108c184d65d6377a2366f77ba8527ee945e31fac7dab78`
- 当時のDebug SHA-256: `aacf3e40afe5dd587923a6c1f10555abd803ec2560147048d797ded05fc6662b`

profileなし・1スレッド・warmup後7回の交互測定で、ジシランとBr複合体の溶媒
解析Hessianは標準数値の3.75～4.37倍。標準版のHessian OpenMP指示のみを外した
Windows ifx互換性ビルドとの比較で、SCC accuracy 1e-7／刻み1e-4 bohr。
全行列とEnergy／Gradientを全試行で比較し、別profile監査で解析dispatchを確認。
解析Hessianの既定は引き続き無効で、大分子溶媒・既定刻み・thread crossoverの
速度を示す結果ではない。表面分岐の保守的拒否と数値fallbackを維持する。

全記録は `gfn1_fast_4_5_0/RESULTS.md`、`BENCHMARK.md`、
`SOLVENT_RESULTS_20261005.json`、`SOURCE_SOLVENT_MANIFEST_20261005.json`。
以前の部品だけのバイナリ・hash・source snapshotは別に保存した。
以下の時刻付き記録とhashは各時点の履歴で、現在の生成物とは区別する。

## 溶媒微分部品追加後の再確認（2026-10-05 08:00 JST）

4.5.0のBorn半径／角度グリッド表面積／重み付きCM5の二次微分部品を含む最新
ソースをRelease／Debugで再ビルドした。通常Energy／GradientのRelease 7ケース、
Debug 2ケース、既存source／数学／dispatcher検証は全PASS。ゼロoverlapの
気相／ALPB／GBSA Energy微分と実CLI Hessian／fallbackも両構成で再確認した。
溶媒微分部品は独立した全座標一次微分の有限差分と比較して両構成でPASS。
検証範囲と詳細は `gfn1_fast_4_5_0/IMPLEMENTATION_STATUS.md` と
`COMPONENT_RESULTS_20261005.json` を参照。溶媒の解析Hessianへの接続は未完了。

- 当時のRelease SHA-256: `1271b860ded6e01ddbba176e5dbceddf41f724faf465dd446c1d1391a158fa20`
- 当時のDebug SHA-256: `9bff65b7bf9501c4dc264d356679852f70f71276a90bc6fa3f373c112ba2025d`

4.4.0のtaxol全Hessian測定は保存したバイナリで先に完了した。1スレッド・
profileなしの5回中央値は標準数値250.238602秒／解析30.070810秒、8.32倍。
全339×339行列差は全試行2.991e-7 Eh/bohr²。標準版のWindows互換性対応、
profile監査と固定hashを含む詳細は `gfn1_fast_4_4_0/BENCHMARK.md` を参照。
この速度測定は現在の4.5部品追加後バイナリの測定ではない。

## 今回の再確認（2026-10-05 07:15 JST）

現在の未コミット変更を含む全ソースをWindows nativeでRelease／Debugともに
クリーンビルドした。両構成で782ステップが完了し、コンパイル・リンク・実行に合格。
WSLは使用していない。今回の作業では計算ソースを変更していない。

- ifx 2025.2.0、MSVC 19.44.35217.0、oneMKL LP64 sequential、CMake 4.1.1／Ninja。
- Pythonは指定の `C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe`、3.13.15。
- Release: 水・taxol 350 AOの気相／ALPB、従来density／汎用Gradient比較を含む7ケースPASS。
- Debug: `/check:all /fpe:0`付きで水の気相／ALPB、2ケースPASS。
- 両構成でsource guards、compact-density数学検証、Fortran dispatcher検証もPASS。
- 経路間Energy差は表示桁で0、Gradient最大差は `2.005e-15 Eh/bohr`。
- 水のRelease／Debug比較はEnergy差0、Gradient最大差 `1.006e-15 Eh/bohr`。
- 生成物はAMD64／PE32+、stack reserveは64 MiB。

| 構成 | 今回生成した実行ファイルのSHA-256 |
|---|---|
| Release | `dd0f95b2e58758bad5be1e9a0987b66c039b27b7822af7086647febff21c74d5` |
| Debug | `df26018d4826156a8564813c0bd795935b67c255c7f5c39c57b00c331a4674da` |

各buildフォルダの `native-verification-clean-build.log` にクリーンビルドログ、
`smoke/<case>/run.log` に実行ログを保存した。数値比較・バイナリ形式・hashは
Releaseフォルダの `native-verification.json` に保存した。
以下のhashと結果は過去の記録であり、現在の生成物とは区別する。

今回確認した範囲はWindowsでのビルドと通常GFN1 Energy／Gradientの動作。
Hessian、GBSAの微分整合性、WSLとの速度比較は今回の対象外。
前回の594 AO座標入力がないため、その計算の再現は未確認。

再実行はリポジトリルートのPowerShellから:

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
```

## 過去の検証記録

4.3.0のhalogen Hessian接続後もRelease／Debugを再ビルドし、通常Energy／Gradient
検証に加えてhalogen補正、分子応答、CLIとfallbackを通した。今回の詳細・最新
実行ファイルhashは `gfn1_fast_4_3_0/RESULTS.md` と
`gfn1_fast_4_3_0/SOURCE_MANIFEST_20261005.json` を参照。
以下の06:16のhashは、その時点のクリーンビルドの記録である。

## 最新版のクリーン再確認（2026-10-05 06:16 JST）

Windows native / PowerShellで、現在の未コミット変更を含む全ソースをRelease／Debug
ともに `cmake --build ... --target xtb-exe --clean-first --parallel 4` で再コンパイルした。
両構成とも782ステップを完了し、リンクと実行テストに合格した。WSLは使用していない。
今回は実装の追加・変更は行っていない。

- 環境: Intel ifx 2025.2.0、MSVC 19.44.35217.0、oneMKL LP64 sequential、
  CMake 4.1.1、Ninja 1.12.1、指定Conda Python 3.13.15。
- Release: 水の気相／ALPB、taxol 350 AOの気相／ALPBと従来density／汎用Gradient比較、計7ケースPASS。
- Debug: allocation／bounds checksとfloating-point traps付きで、水の気相／ALPB、計2ケースPASS。
- 両構成のsource guards、compact-density数学検証、Fortran dispatcher検証もPASS。
- Energy差は表示桁で0。taxolの経路間Gradient成分最大差は `2.005e-15 Eh/bohr`、
  水のRelease／Debug間は `1.006e-15 Eh/bohr` 以下。
- PE形式は両方AMD64／PE32+、stack reserveは64 MiB。

クリーンビルドログは各buildフォルダの `clean-build.log`、実行ログは
`smoke/<case>/run.log`、数値比較と最新SHA-256はReleaseフォルダの
`release_debug_comparison.json` に保存した。

| 構成 | 今回生成した実行ファイルのSHA-256 |
|---|---|
| Release | `378F37C1CF9DA4B45502B8FA2EBD19B8B1153B1E671598FCFF9EEF13CE9AF879` |
| Debug | `147A4A1E9EC6DF3AE6695BC4FD1D3850DA1BAADAEB00C29501C817D3D1120086` |

今回の確認はビルドとEnergy／Gradientの動作検証。Hessianや性能測定については
以下の以前の検証記録を参照。前回の594 AO入力とWSL実行ファイルは今回比較していない。

## 以前の検証記録

確認日: 2026-10-05。Windows native / PowerShellで実施し、WSLは使用していない。
Git HEADは `c2ccf660ffe2ddc5f96971b77fdb325c82a8fbeb`（2.2.6）。
未コミットの2.3以降の開発変更と分子Hessian応答モジュールを含む現時点の全ソースを、
新しい `build-gfn1-fast-current-windows-ifx-release` / `-debug` フォルダでビルドした。

## 環境

- Intel Fortran ifx / IntelLLVM 2025.2.0、MSVC 19.44.35217.0
- oneMKL LP64 sequential BLAS/LAPACK、CMake 4.1.1、Ninja 1.12.1
- Python 3.13.15: `C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe`
- mctc-lib: `77f65c6f2cf6330d05d0757ca173da097096780e`
- Release `/O3`、Debug `/Od /debug:full /check:all /fpe:0`
- IntelのWindows用dialect flags、64 MiB main-thread stack、実行時1スレッド
- tblite / CPCM-X / JSON / upstream unit-test targets / LTOは無効

Windows用フラグ、HOME不在時のUSERPROFILE fallback、ライブラリ構成、
stack設定は既存のWindows対応を使用した。今回はrepulsion／D3の二次微分と
電子・電荷・座標応答を結合する分子Hessianモジュールも含めて再ビルドした。
追加のAPI経由Debug検証で、EEQ charge-only呼び出しが未確保の微分配列を
必須引数へ渡す既存の不具合を検出した。微分引数をoptionalにし、charge-onlyでは
省略するよう修正した。energy accumulatorもゼロ初期化し、微分を要求する呼び出しは
必要な引数の存在を検証する。修正後のRelease／Debugで以下の全通常テストを再実行した。
4.1.0の独立した厳密検証は `gfn1_fast_4_1_0/RESULTS.md` を参照。
互換性対応の詳細は `gfn1_fast_2_2_6/WINDOWS_TEST_REPORT.md` を参照。

## 結果

ReleaseとDebugの全ソースのコンパイル・リンク: **PASS**。
両方の実行ファイルがPE x64であり、stack reserveが67,108,864 bytesであることも確認した。
両構成で既存のsource guards、compact-density数学テスト、Fortran dispatcherテストに合格した。

ReleaseのEnergy/Gradient計算7ケースが合格:

| 入力 | NAO | Energy / Eh | Gradient norm / Eh per bohr |
|---|---:|---:|---:|
| 水、気相 | 8 | -5.764209595203 | 0.094423990378 |
| 水、ALPB water | 8 | -5.787247944321 | 0.075935279788 |
| taxol、気相 | 350 | -195.911950900123 | 0.149903893916 |
| taxol、ALPB water | 350 | -195.958435141829 | 0.136167479182 |

taxol気相についてproduction / legacy density / generic gradientの3経路、
ALPBについてproduction / legacy densityの2経路を実行した。
全ケースでSCC収束、有限なGradient成分、Gradientファイルと表示normの一致を確認した。
taxolは11 SCC反復で収束し、productionは10反復でcompact densityを使用した。

| 比較 | Energy絶対差 / Eh | Gradient成分最大絶対差 / Eh per bohr |
|---|---:|---:|
| taxol production / legacy density | 表示桁で0 | 2.005e-15 |
| taxol production / generic gradient | 表示桁で0 | 9.992e-16 |
| taxol ALPB production / legacy density | 表示桁で0 | 1.998e-15 |
| 水 気相 Release / Debug | 表示桁で0 | 9.992e-16 |
| 水 ALPB Release / Debug | 表示桁で0 | 1.006e-15 |

許容値はEnergy `1e-11 Eh`、Gradient成分 `1e-10 Eh/bohr`。
Debugは水の気相・ALPBの2ケースをallocation/bounds checksとfloating-point traps付きで実行した。
大分子のDebug実行、旧WSLバイナリとの比較は今回の検証範囲に含めない。
前回の594 AO入力はフォルダにないため、その計算の再現は未確認。
Hessian応答の独立した厳密実行検証もRelease／Debug双方で合格した。
水の0 K・30000 K・開殻3000 Kとdisilaneの300 Kについて、実際のSCC計算と
全座標のGradient有限差分を比較した。気相Hessianには電子応答、Coulomb、
repulsion、pairwise D3を含め、最大差は刻み2.5e-4 bohrで約4.17e-8 Eh/bohr²。
刻みを半分にすると差が約1/4になること、対称性と並進不変性も確認した。
4.2.0ではこのAPIを明示的に有効化できるCLI経路へ接続した。
水の温度・開殻4ケース、disilane、benzeneで数値Hessianとの比較を通した。
ALPB、HCl、拘束、無効化、メモリ制限で数値微分へfallbackし、4スレッドの数値
微分も1スレッドと一致した。Debugの小分子CLI検証と、両構成での部分列加算・
restart維持・polarizability要求fallbackのAPI検証も合格した。
通常CLIの既定Hessianは数値微分。溶媒・halogen等の解析応答は未完成。
詳しい検証範囲は `gfn1_fast_4_2_0/RESULTS.md`、同一Windows環境での標準版
との速度比較は `performance/BENCHMARK.md` と
`gfn1_fast_4_2_0/BENCHMARK.md` を参照。

コンパイラには既存のstructure alignment、無効化された機能のINTENT(OUT)、
`/machine:x64`を無視する警告がある。ビルドは正常終了し、生成物のx64形式を独立に確認した。

## 再実行

リポジトリルートのPowerShellから:

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
```

スクリプトはoneAPI環境初期化、ビルド、テストを実行し、呼出元の環境を復元する。
実行ファイルは各buildフォルダの `xtb.exe`、ビルドログは `build.log`、
実行ログは `smoke/<case>/run.log`、結果は `smoke/summary.json` に保存される。
Release/Debug比較はReleaseフォルダの `release_debug_comparison.json` に保存した。

追加の分子Hessian応答検証:

```powershell
.\tests\gfn1_fast_4_1_0\run_windows_molecular.ps1 -LibraryBuildType Release
.\tests\gfn1_fast_4_1_0\run_windows_molecular.ps1 -LibraryBuildType Debug
```

各buildフォルダの `molecular-response-test/molecular-response-test.log` に結果を保存した。

通常の入力を直接実行する場合は、oneAPIのDLL探索環境を初期化する:

```powershell
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_STACKSIZE = '64M'
cmd /d /c 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul && "C:\GitHub_repository\Fast_GFN1_xTB\build-gfn1-fast-current-windows-ifx-release\xtb.exe" input.xyz --gfn 1 --grad --norestart'
```

4.2.0 CLI接続とWindows数値Hessianのworker初期化修正後の実行ファイルSHA-256:

- Release: `0F36AB0A4DE856726A604353F9973CCB8B04C80BBCA3784D2F251FDA45B24BF2`
- Debug: `7A0DF1A1190CE68B64FD3A54AC359274077C531EB60357FD36DBE38B844966E8`

buildフォルダは再ビルドで更新される。固定された版の配布物として扱わない。

### 4.5.1 大分子溶媒の反復比較完了（2026-10-05）

保存した4.5.1実行ファイルで、taxol（113原子・350 AO）の全339×339
Hessianを標準数値微分と比較。1スレッド、profileなし、warmup 1回＋5回反復。
GBSA中央値は標準229.775900 s／解析35.167327 s（6.533789倍）、
ALPBは229.952206 s／33.692745 s（6.824977倍）。全行列差は最大
7.20e-9／1.00e-8 Eh/bohr²。標準側はWindows ifxの既知の割当descriptor
クラッシュを回避するためHessianのOpenMP指示のみ除去した比較用ビルド。
詳細は `gfn1_fast_4_5_1/LARGE_BENCHMARK_20261005.json`。

## 配布ソースの再展開確認完了（2026-10-05）

`dist/xtb-gfn1-fast-4.6.0.zip` のSHA-256は
`47094abd1d849ea8be4f91ea3afb848ee43343102263669d5fdb8bca8b6ebce9`。
新規再展開先で921ファイル、mctc-lib依存124ファイルと最終計算ソース594ファイル
を照合した。Git情報を同梱せず、親リポジトリ検索も遮断してRelease／Debugを
各786ステップでビルド。両構成でEnergy／Gradient smoke、実CLI既定21条件、
screening／D3 cutoff窓、溶媒全座標Hessian 38条件がすべてPASS。
ZIPは変更せず、新しいバイナリhash・ログhash・結果を外部の
`dist/REEXTRACTION_4.6.0.json` と `gfn1_fast_4_6_0/REEXTRACTION_20261005.json`
に固定した。ビルドパスとGit情報が異なるので実行ファイルのhashは元ビルドと
異なるが、計算ソースは594ファイルすべて一致している。
