# Fast_GFN1_xTB 作業状況（2026-10-05）

## 最新：4.7.0 解析HessianのOpenMP並列化

susceptibility、座標応答、AO積分微分、SASA微分を並列化した。
参照cacheは読み取り専用、scratchは各worker内で確保し、列の分離と
幾何行列のworker別集計で処理する。4列のBLAS収縮、overlap変換の再利用、
Born／SASAの6座標更新で直列経路も高速化した。
既定は64 AO以上、最大8 worker、解析workspace上限1 GiB。上限に応じて
workerを減らし、直列でも収まらない場合や従来のguard失敗時は数値微分へ戻す。

Release／Debugで実際の2／4／8 worker、runtimeの2-thread上限、並列fallbackを
19条件ずつ検証した。既定CLI 21条件、溶媒非線形SCCの全座標FD 38条件、
気相のdensity／shell電荷応答、部分列加算、入れ子API、OpenMPなしのresource
コンパイルもPASS。taxolの8-thread未投影Hessian全339列は、GBSAで最大
9.42e-9、ALPBで1.60e-8 Eh/bohr²の独立Gradient-FD差でPASS。

固定した4.6.0／4.7.0バイナリをwarmup後5回ずつランダムな順序で測定した。
taxol（113原子、350 AO）の全Hessian計算の中央値は次のとおり。

| モデル | 4.6.0・1 thread | 4.7.0・1 thread | 4.7.0・8 threads | 4.6.0比 |
|---|---:|---:|---:|---:|
| 気相 | 30.83 s | 15.60 s | 4.27 s | 7.22倍 |
| GBSA | 35.11 s | 16.21 s | 4.51 s | 7.78倍 |
| ALPB | 33.86 s | 16.00 s | 4.55 s | 7.44倍 |

全試行の出力全行列差は0、Energy／Gradient normは1e-10以内で一致。
同じ4.7.0の1 threadに対して8 threadsは3.51～3.65倍。MKLはsequential。
これは高速版の旧解析経路との比較であり、標準xTBの複数thread比較ではない。
記録は `tests/gfn1_fast_4_7_0/RESULTS.md` とhash付きJSONを参照。
`dist/xtb-gfn1-fast-4.7.0.zip` を再展開し、958 source／124依存ファイルの
SHA照合、Git metadataなしのRelease／Debugビルド、各19 parallel CLI／
21 default／38 native FD条件とresource単体検証を完了した。
最終再展開記録は `dist/REEXTRACTION_4.7.0.json` に保存した。
以下は各段階の履歴であり、4.6.0以前のsnapshotとZIPは維持している。

## 確認したベース

Git HEAD は `c2ccf660ffe2ddc5f96971b77fdb325c82a8fbeb`、コミット名は
`version 2.2.6`。前回の「高速化戦略とコード確認」の引き継ぎと実装を照合した。
CPU専用のGFN1 Energy／Gradientを対象とし、最終full diagonalization、
compact density、SYEVD既定、汎用Gradientへのfallbackを維持する。

594 AOの前回ログは確認したが、その座標入力はこのフォルダにない。
前回のWSL時間と今回のWindows時間を直接比較していない。

## 現在のWindows実行環境

Intel ifx 2025.2、Visual Studio 2022、oneMKL LP64 sequential、CMake／Ninjaで
Windows x64のRelease／Debugビルドと実行を確認した。WSLは使用していない。
PythonはAGENTS.md指定のConda `fastxtb`（3.13.15）を直接呼び出す。

Windows用コンパイラフラグ、HOME不在時のUSERPROFILE fallback、64 MiBの
stack設定を追加した。mctc-libは指定コミットに固定して配置済み。
Releaseで水／taxolのgas・ALPB、従来density／汎用Gradientとの比較を通した。
Debugで水／disilaneの小分子検証を通した。

詳しい環境設定と初回検証結果は
`tests/gfn1_fast_2_2_6/WINDOWS_TEST_REPORT.md`を参照。

## 実装・検証した開発変更

| 段階 | 現状 |
|---|---|
| 2.3.2 subspace reuse | 実装・数値検証済み。性能目標は未達のため明示的な有効化が必要 |
| 2.4.0 Gradient | packed H0利用とCartesian primitive収縮を実装・検証済み |
| 3.0.0 diagonalization-free SCC | 密度・energy-weighted density直接更新の検証用実装。既定では無効。高温の一部はfull solveへfallback |
| 4.0.0 Hessian展開 | 有限温度・開殻P/W応答、連立shell電荷応答、Coulomb応答kernelを実装。分子全体の解析Hessianは未完成 |
| 4.1.0 座標微分・分子応答 | H0／overlap／CN／Coulomb／repulsion／D3の二次微分を電子・電荷応答と結合。小分子気相で検証済み。CLI接続は4.2.0、溶媒・halogen対応は未完成 |
| 4.2.0 解析Hessian CLI | 対応済みの気相系で明示的に有効化可能。未対応モデル・失敗・メモリ超過は数値微分へfallback。小分子で標準版との速度比較を完了 |
| 4.3.0 halogen Hessian | Cl／Br／I／Atの気相系を解析経路へ接続。補正単独96ケース、分子7ケース、CLIとbranch fallbackをRelease／Debugで検証済み |
| 4.4.0／2.4.1 ゼロoverlapのGradient | H0/S再構築による実分子の微分欠落を直接prefactorで修正。気相／ALPB／GBSAのEnergy微分、Hessian、回転整合性をRelease／Debugで検証。taxol全Hessianの5回測定で標準版の8.32倍 |
| 4.5.0 溶媒Hessian | GBSA／ALPBのStill／P16、CM5、表面積、H-bond、ALPB形状微分をSCC応答・CLIへ接続。分子38条件とCLI／fallbackをRelease／Debugで検証。小分子溶媒Hessianは標準数値の3.75～4.37倍 |
| 4.5.1 溶媒分岐・メモリ配置 | 深く埋没した表面点の不要な分岐拒否を除去し、溶媒とAOの大きな一時配列を順番に確保。taxolのGBSA／ALPB全339列を検証し、CLIでも既定1 GiB上限内で解析経路を確認。小分子1／4／8スレッド比較と大分子溶媒の標準版との反復比較を完了。GBSA 6.53倍／ALPB 6.82倍 |
| 4.6.0 既定採用・分岐fallback | Gradient整合性とcutoff窓の拒否を追加。既定の解析Hessianとfallbackの実CLI 21条件をRelease／Debugで検証 |

4.4.0では、実際に収束した非対称水で従来Gradientの `0.0468 Eh/bohr` の誤差を
再現し、Energy有限差分で独立確認した。修正後は解析Gradientとの差が
`3e-16 Eh/bohr` 以下。小さい非ゼロoverlapのSCF screeningは区別し、通常の
水／taxol値を維持する。legacy Gradient指定では解析Hessianを数値微分へ戻す。
GBSAの既存SASA閾値をまたぐ有限差分の不連続も標準版で再現した。滑らかな
領域でゼロoverlapを維持するGBSA分子は厳密Energy微分検証を通している。
詳細は `tests/gfn1_fast_4_4_0/RESULTS.md` を参照。

7回の反復測定で、修正後のtaxol Energy／Gradientは標準版の1.49倍（気相）、
1.47倍（ALPB）。直前の高速版との差は約1%以内で、追加の速度向上は主張しない。
直接prefactorは微分の正しさのために採用する。4.3のtaxol長時間測定はこの
実分子の不具合を検出して中断し、保存した4.4バイナリで別フォルダへ再測定を完了。
taxolの5回中央値は標準版250.238602秒／解析版30.070810秒、8.32倍。
全339×339行列の最大差は全試行で2.991e-7 Eh/bohr²。profileなし・1スレッド、
別profile監査で解析dispatchを確認した。標準数値HessianのWindows互換性対応と
比較条件は `tests/gfn1_fast_4_4_0/BENCHMARK.md` に明記した。
4.5.0では溶媒中の全座標Hessianも独立Gradient有限差分で検証し、1e-4 bohrで
最大差1.172e-8 Eh/bohr²。塩・メモリ超過・旧Gradient経路・表面の分岐近傍は
数値微分へ戻す。ジシラン／Br複合体のGBSA／ALPBはprofileなし・1スレッド・
warmup後7回中央値で標準数値の3.75～4.37倍。全行列とEnergy／Gradientを比較した。
有限温度の既存電子数許容誤差によるEnergy微分差も標準版で再現し、処理は維持。
詳細は `tests/gfn1_fast_4_5_0/RESULTS.md` と `BENCHMARK.md` を参照。
4.5.1では最終メモリ配置でtaxolの未投影Hessian全339列を独立した非線形SCC
Gradient有限差分と比較し、GBSA最大2.854e-8／ALPB最大1.592e-8 Eh/bohr²でPASS。
CLI監査でも解析経路を使用し、推定workspace 880428680 bytes、記録したpeak
pagefile usageは両モデルとも約858 MBで既定上限内。1／4／8スレッドの小分子
7回反復比較では、同スレッド数の高速版数値微分に対して1.15～5.79倍。
通常の既定刻みの水・ジシランも改善。標準版の数値Hessianとの大分子溶媒
5回反復比較を保存した4.5.1バイナリで完了した。GBSAは、標準
229.775900秒／解析35.167327秒の中央値で6.53倍。全行列の最大差は全試行
7.20e-9 Eh/bohr²以下。ALPBは標準229.952206 s／解析33.692745 s、6.824977倍で、全行列差は最大1.00e-8 Eh/bohr²。
詳細は `tests/gfn1_fast_4_5_1/RESULTS.md` と各hash付きJSONを参照。
4.6.0では小さい非ゼロoverlapのscreening点で0.0468 Eh/bohrのGradient差を
実際に再現し、production Gradientとの整合性確認で数値微分へ戻すようにした。
D3の60 bohr cutoff上とその有限差分窓内も拒否。CN／repulsion／halogenの
距離・近傍分岐も窓を確認する。対応GFN1解析Hessianを既定で有効化し、
ENABLE=0とDISABLE=1を維持。実CLIの既定／fallback 21条件と溶媒分子38条件を
Release／Debugで検証した。詳細は `tests/gfn1_fast_4_6_0/README.md`。
配布ZIPを新しいフォルダへ再展開し、921ファイルと依存124ファイルのSHA照合、
Git情報なしのRelease／Debug各786ステップのビルド、既定CLI 21条件と溶媒
全座標Hessian 38条件、screening／D3境界検証を両構成で完了した。
`dist/xtb-gfn1-fast-4.6.0.zip` と外部検証記録 `dist/REEXTRACTION_4.6.0.json`
を保存した。2.3.2／2.4.x／3.x／4.xの実装・比較・採用判断を完了。
遅かったsubspace reuseとfermi-operator SCCは検証用のまま既定で無効とする。

2.3.2は同一Hamiltonianでfull解と固有値・密度・電荷・band free energyを
照合し、最終反復をfull経路に戻す。零温度、300 K、1000 K、開殻、ALPBの
Energy／Gradient比較を通した。taxolでfull solveを11回から9回に減らせるが、
強制再利用の全体時間は複数回の測定で10～22%遅い。目標のfull solve 4～6回は
達成していない。既定の通常実行はSYEVDを使用する。

2.4.0は従来・個別有効化・汎用Gradient・4スレッド・ALPB・開殻との比較を
通した。s/p/d全9組合せの積分微分と全Energyの有限差分も検証した。
taxolの5回中央値では全体約0.7～1.3%短縮、電子Gradient kernelは
0.036→0.035秒。小さい差であり、より広い入力での性能検証が必要。

これらは未コミットの変更で、検証済みソースZIPは上記のとおり保存している。
各段階の詳細は`tests/gfn1_fast_2_3_2/RESULTS.md`、
`tests/gfn1_fast_2_4_0/RESULTS.md`、
`tests/gfn1_fast_3_0_0/RESULTS.md`、
`tests/gfn1_fast_4_0_0/RESULTS.md`、
`tests/gfn1_fast_4_1_0/RESULTS.md`を参照。

## 再実行

リポジトリルートのPowerShellから:

```powershell
.\tests\gfn1_fast_2_4_0\run_windows_regression.ps1
.\tests\gfn1_fast_2_4_0\run_windows_regression.ps1 -BuildType Debug
```

最新版作業ツリー全体のWindows確認は以下で再実行する:

```powershell
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Release
.\tests\gfn1_fast_2_2_6\run_windows_build.ps1 -BuildName build-gfn1-fast-current-windows-ifx -BuildType Debug
```

最新の実行ファイルは`build-gfn1-fast-current-windows-ifx-release/xtb.exe`。
直接実行時もoneAPIの環境を初期化する。再現用スクリプトが環境初期化、
ビルド、検証をまとめて行う。buildフォルダは名前にかかわらず最後にビルドした
ソースの内容を持つため、バージョン別の保存済み配布物として扱わない。

## 次の作業

3.0.0は零温度のpurificationと有限温度の行列Fermi関数からP／Wを計算し、
電子数・HP=SW・定期full auditを検証する。disilaneの零温度～5000 K、
taxolの零温度・300 Kで固有値分解を省く反復を確認した。
taxolの1000 K・5000 Kはauditで棄却し、full solveで最終結果を維持する。
現在のdense演算は遅く、検証用の明示的な有効化が必要。
ALPB・開殻・強制fallbackを含む回帰と5回中央値の測定を完了した。
taxolのFOE強制経路は通常の約3.4倍の時間を要する。高温では従来の
occupation電子数誤差（約5e-10～1.2e-9）との差を確認したが、auditを緩和しない。

4.0.0ではfinite-temperature／開殻P/W応答、shell電荷の連立応答を追加し、
独立した一般化固有値計算と非線形SCCモデルの有限差分で検証した。
実際のCoulomb／三次電荷項の応答kernelも既存addShiftと照合した。
次は分子座標でのH0／overlap、CN／solvent／dispersion等の二次微分と結合し、
分子全体のCartesian Hessianを組み立てる。数値Hessian比較とfallbackを維持する。

4.1.0は実GFN1基底の水／disilaneでH0／overlapの全座標微分と固定密度Hessianを
既存ビルダー・Gradientの有限差分に照合した。CNと固定電荷Coulombの二次微分も
Release／Debugで検証した。厳密にゼロのoverlapで既存H0/S Gradientと任意の
試験密度が一致しない制限を検出し、新しい式は独立したEnergy二次差分で検証した。
詳細と未完了の項目は4.1.0のRESULTSに記録した。4.2.0では新APIをCLIへ接続し、
数値Hessian fallbackを維持する。`XTB_GFN1_FAST_ENABLE_ANALYTIC_HESSIAN=1`で
対応済みの気相モデルに対して有効にできる。既定では無効。

追加の分子Hessian APIは水の零温度・高温・開殻とdisilaneで全Gradient有限差分に
照合し、Release／Debugで合格した。CLIの既定は引き続き数値Hessianを使用する。
2026-10-05に最新版全ソースのWindowsビルドと通常実行を再確認し、API初期電荷の
EEQ charge-only経路の未確保引数と未初期化energyを修正した。再現手順と検証範囲は
`tests/WINDOWS_CURRENT_TEST_REPORT.md`に記録した。

標準xTB 6.7.1（公式commit `26b28010e805f7d1aeeef39813feb473e69cc4be`）も
同じWindows ifx／oneMKLでビルド済み。各段階のpolicyを同一バイナリで切り替え、
1スレッド・profileなし・warmup後7回の交互測定と標準版との数値比較を完了した。
最新版2.4経路のtaxol気相／ALPBは標準版の約1.45／1.44倍の速度だった。
reuse／dense FOEは追加すると遅いため、引き続き既定では無効とする。
詳細は`tests/performance/BENCHMARK.md`を参照。

4.2.0の解析Hessianはdisilane／benzeneで標準の数値微分の4.48／7.87倍、
現在の高速版数値微分の3.37／5.71倍を確認した。標準版はWindows ifxのOpenMP
descriptor停止を避けるためHessian内の並列指示のみ外した1スレッドの比較用ビルド。
元の数値微分式・SCC・Gradientを使用する。溶媒・halogen・大分子の解析Hessianの
速度を示す結果ではない。詳細は`tests/gfn1_fast_4_2_0/BENCHMARK.md`を参照。

4.3.0でhalogen補正の解析二次微分を分子全体へ接続した。実SCCの全座標Gradient
有限差分との差は刻み2e-4 bohrで約1.99e-8 Eh/bohr²以下。最近接原子の同距離点や
cutoff上では数値微分へfallbackする。CH3Br／CH3I + waterの8原子複合体では、
profileなし・warmup後7回の交互測定で標準数値Hessianの3.99／3.98倍を確認した。
同じ標準版の1スレッド互換性対応を使用した比較で、溶媒モデルは使用していない。
詳細は`tests/gfn1_fast_4_3_0/RESULTS.md`と`BENCHMARK.md`を参照。
溶媒応答、より大きい分子の性能・メモリ検証、screened/nodal積分の整合性と
正式配布物の作成は未完了。解析Hessianの既定は引き続き無効。

taxol（113原子／350 AO）の解析Hessianも完了し、標準版の339×339数値Hessianと
最大差2.991e-7 Eh/bohr²で一致した。SCC accuracy 1e-7、刻み5e-4 bohr。
保存した4.4バイナリのwarmup後5回・profileなし交互測定も完了し、標準250.238602秒／
解析30.070810秒、8.32倍を確認した。これは気相の結果で、4.5の大分子溶媒測定は未完了。
