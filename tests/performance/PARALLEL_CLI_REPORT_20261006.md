# --parallel の挙動調査（2026-10-06）

対象: GFN1-fast 4.7.0 Windows Releaseの配布ZIP内の実行ファイル。
CPU: Intel(R) Core(TM) i7-14700F。MKL sequential、SYEVD、oneAPIのPATHを使用しない実行。
実行ファイルSHA-256: `008fb7cfcd23e3658c1bc2ae85cd63f083e4bb13aac4358ce68891bd7b1a3085`。
実装・配布ZIPは変更していない。

## 指定の意味と優先順位

`--parallel N` または `-P N` はOpenMPのスレッド数を設定する。
`src/prog/main.F90` の引数処理で `omp_set_num_threads(N)` を呼ぶ。
`OMP_NUM_THREADS=1` + CLI 4で4、環境変数8 + CLI 1で1になった。
繰り返して指定した場合は最後の指定を採用する。
配布ランチャーは環境変数未指定時に8を設定するが、CLIで1/4に上書きできた。
ヘルプの「number of parallel processes」は実装のスレッドという意味に合っていない。

起動ログの `omp threads` は `src/printout.f90` の診断用parallel regionで
実際に作成されたチームの人数。Hessianのチーム人数を表すものではない。
`OMP_THREAD_LIMIT=2` とCLI 8では起動診断が2、Hessian requested workersが8、
actual response threadsが2だった。

## 計算経路ごとの制限

- MKLは直列。現在のビルドの前処理済みmainにも `mkl_set_num_threads` は存在せず、`--parallel` で固有値ソルバーやBLASを並列化しない。
- 解析的Hessianは既定で64 AO以上から並列化。最大8ワーカーに加えて、座標数・shell数・推定作業配列メモリ上限1 GiBを考慮する。小さい水分子ではCLI 8でもHessianの実ワーカーは1。
- 検証用にAO閾値を1にすると水でCLI 2/4の実ワーカーが2/4になった。
- Taxol（113原子、350 AO）ではCLI 1/4/8の実ワーカーが1/4/8。CLI 16でも起動診断は16、Hessianは8。
- ALPBのdisilaneでも、検証用のAO閾値1とCLI 4で実ワーカー4になった。
- 数値Hessian fallbackは別の方針。24 AO以上かつ座標ジョブ数が指定数の2倍以上の場合に外側の変位計算を並列化する。解析的Hessianの最大8制限はこの経路には適用されない。
- 水の数値HessianはCLI 8で外側並列なし。DisilaneではCLI 8で外側並列あり、CLI 16ではジョブ数24 < 32なので外側並列なし。外側並列なしでも、各SCC/Gradient内部のOpenMP設定は残る。

## 気相Taxolの速度

1回のwarm-up後、1/4/8の順序を各回ランダム化して5回測定。
表はプロセス起動・出力を含むwall timeの中央値。タイミング測定では追加profileを無効化し、実ワーカーの確認は別に実施した。
`--gfn 1 --norestart`、`$hess sccacc=1e-7, step=1e-5`。
初期SCCの `--acc 1e-7` 指定は既存CLIが `1e-4` に補正している。

| 指定 | Gradient (s) | 対1スレッド | Hessian (s) | 対1スレッド |
| --- | ---: | ---: | ---: | ---: |
| `--parallel 1` | 0.348 | 1.00倍 | 15.647 | 1.00倍 |
| `--parallel 4` | 0.244 | 1.43倍 | 5.993 | 2.61倍 |
| `--parallel 8` | 0.217 | 1.60倍 | 4.582 | 3.42倍 |

このPC・分子・設定での結果であり、全入力に同じ加速率を保証するものではない。

## 数値結果

Taxolの解析的Hessian全成分はCLI 1/4/8/16間、全測定反復とも出力値で一致。
Energyも出力桁で一致した。Gradient全339成分を全Gradient測定反復で比較した最大差は
`3.997e-15 Eh/bohr`。
Disilaneの数値HessianはCLI 1/8で出力値が一致し、CLI 16は最大差約
`1.0e-10 Eh/bohr²`（出力最終桁1単位）だった。この差の発生箇所の特定は今回の対象外。

## 引数検証の弱点

現在のIntel OpenMPランタイムを使う配布版では、以下がすべて終了コード0になった。

| 指定 | 環境変数OMP=8での起動チーム人数 | 挙動 |
| --- | ---: | --- |
| `--parallel 0` | 1 | 自動選択ではなく1スレッド扱い |
| `--parallel -1` | 1 | 1スレッド扱い |
| `--parallel abc` | 8 | 解析失敗の警告、指定を無視して継続 |
| `--parallel` のみ | 8 | このオプションの値省略を無視して継続 |

正の整数を指定すること。今後修正するなら、`N >= 1` と値の存在をCLIで検証し、
不正入力を非ゼロ終了コードで拒否するのが明確。今回は調査のみで挙動は変更していない。

## 再現と保存ログ

```powershell
$PY = 'C:\Users\f3r1i\mambaforge\envs\fastxtb\python.exe'
& $PY tests/performance/investigate_cli_parallel.py `
    --exe '<配布ZIP展開先>/bin/fast-gfn1-xtb.exe' `
    --output build/cli-parallel-new-investigation
```

主な実測・追加fallback/ランチャー監査のJSONと各run.log:
`build/cli-parallel-investigation-final/RESULTS.json`。
追加の数値Hessian監査は `XTB_GFN1_FAST_DISABLE_ANALYTIC_HESSIAN=1` と
CLI 1/8/16で実行し、`outer OpenMP` のprofileと行列全成分を比較した。
起動条件、raw timing samples、requested/actual worker数をJSONに保存している。
