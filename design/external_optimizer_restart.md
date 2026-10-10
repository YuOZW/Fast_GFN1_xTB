# 外部構造最適化向けrestart設計案

2026-10-07。対象はFast GFN1-xTB 4.7.0。実装前の設計案であり、以下のSession APIやCLIオプションはまだ存在しない。

外部ライブラリがCartesian座標を渡し、エネルギーとGradientを受け取る用途を想定する。採用・棄却の通知がないライブラリにも対応する。まず既存のWindows配布exeを利用できるCLI backendを作り、同じSessionを共有ライブラリ／常駐worker backendに拡張する。

## 現在のコードから分かること

- `src/restart.f90` のGFN1ファイルrestartはshell電荷 `qsh` を保存する。座標、元素の順序、溶媒条件は保存しない。読み込みの主判定は原子数・shell数であり、method／電子数／不対電子数の不一致は警告になる。
- `src/api/interface.f90` はResults内の `TRestart` を再利用する。`src/api/results.f90` の `xtb_copyResults` で試行用Resultsをコピーできる。従って採用済みResultsを保護しながら試行点を評価する土台はある。
- APIの `xtb_getCharges` は原子電荷を返す。shell電荷を取得／設定する公開APIは現在なく、予測を本体へ渡すための拡張が必要である。
- Broyden履歴はSCC呼び出しごとに初期化される。MOを保持しても、既定の直接SYEVDがそのMOから反復を開始するわけではない。
- 現在のportable配布はexeを使う構成。既存xtb-pythonをインストールするだけでは、このforkのDLLが使用されることを保証できない。

## 外部ライブラリとの契約

呼び出し側にSessionを置き、restart、結果cache、checkpointの管理をSessionに集約する。エンジンは各座標で通常のSCCと最終確認を実行する。

```text
session = Session(model, backend, restart_policy)
evaluation = session.evaluate(positions, parent_id=None)
# evaluation: id, energy, gradient, converged, restart_id, diagnostics

session.accept(evaluation.id)     # 通知可能なoptimizerでのみ使用
session.reject(evaluation.id)     # 通知可能なoptimizerでのみ使用
session.checkpoint(path)
```

`evaluate` のSCC成功と、optimizerによる構造の採用は別の事実として保存する。評価に成功しただけでは採用済み構造を更新しない。エネルギーが下がったことだけから採用を推測しない。拘束付き最適化や遷移状態探索にも同じ契約を使う。

採用通知がない場合は、採用軌道に基づく予測を無効にして、座標に基づく近傍restart選択を使う。呼び出し順やジョブの完了順を軌道と解釈しない。これは通常の `energy_gradient(R)` 型インターフェースで使用できる。

## 保存する状態

| 状態 | 内容と用途 |
|---|---|
| Model identity | 元素と順序、原子数、shell定義、電荷、スピン、GFN1パラメータ、電子温度、溶媒と設定、外部電荷／場、境界条件、エンジンと形式のversion |
| Evaluation snapshot | 一意の評価ID、座標と単位、収束済みqsh、収束条件、親ID、収束状態、SCC反復数・最終残差 |
| Result cache | 同一座標のenergy・Gradient、使用条件・単位、達成した精度。Hessianは必要な場合だけ追加 |
| Accepted history | optimizerが明示的に採用した評価ID。予測用の短い履歴を保持 |
| Optional memory state | API backendで使うResults／TRestart。濃密行列を含むため数と総メモリを制限 |

Model identityが変わった場合は別Sessionとして扱う。外部電荷の位置などが変化する用途は初版では履歴を無効にする。将来、周辺環境も予測変数として扱う拡張を検討する。

同一座標の結果cacheと、近傍座標の初期電荷は別々に扱う。結果cacheは丸めた座標や距離の近さで命中させない。座標を規定の単位・配列順に正規化して完全一致を判定し、要求精度が保存結果より厳しい場合は再計算する。並進・回転した構造のGradientをそのまま返すこともしない。

初版では、長期保存する電子状態はqshを基本とする。全MOや密度行列、Broyden履歴の永続化は測定後に判断する。

## 初期値選択と予測

1. 同一座標・条件・十分な精度の成功結果があればenergyとGradientをcacheから返す。
2. 指定された親、採用済み点、近傍の収束済み評価からrestart候補を選ぶ。座標RMS変位と最大原子変位を評価し、遠すぎる候補は使わない。
3. 基本は近傍点のqshを初期値にする。棄却された構造でもSCCが収束していれば、この候補cacheには残せる。ただし採用軌道の履歴には入れない。
4. 採用履歴があり、変位方向と歩幅に十分な連続性があるときだけ電荷予測を試す。
5. 予測が不適切な場合は通常warm startへ戻し、それでも収束しなければ標準初期guessを使う。

最初の予測候補は、2点の採用履歴を用いた方向付き線形予測とする。共通Cartesian座標系で `d = Rk-R(k-1)`、`t = Rnew-Rk` として、

```text
alpha = dot(t,d) / dot(d,d)
q_predict = qk + alpha * (qk-q(k-1))
```

ただし小さい `dot(d,d)`、大きい直交成分、方向反転、大きい変位、過大なalphaでは予測を無効にする。閾値は設定可能とし、実際の最適化軌道で決める。等間隔・同方向の場合だけ、前回実験の `2*qk-q(k-1)` に一致する。

孤立分子では、距離判定や予測用座標の剛体位置合わせを検討する。外部場、外部電荷、周期系などで向き・位置に物理的意味がある場合には、自動的に剛体運動を除かない。元素順序の自動並べ替えは初版では行わない。

予測電荷は有限値、次元、総電荷を検査し、変化量を制限する。電荷保存の微小な誤差は補正し、大きい不整合は棄却する。予測の初期残差に上限を設け、通常warm startの実際の残差と比較するauditを限定的に行う。比較のために毎回余分な固有値計算を行う方法は、その費用まで測定して判断する。

将来は短い履歴から正則化した局所応答を推定する方法を検討する。既に解析Hessian応答が得られている場合は `q_predict = q + (dqsh/dR)*dR` も候補となる。予測のためだけに毎回全座標の応答を計算することは初期方針にしない。

## CLI backend

各評価を独立した作業ディレクトリで実行する。親snapshotから従来形式の `xtbrestart` を作業ディレクトリにコピーし、そこでexeを実行する。コピー先を本体が上書きしても親snapshotは保持される。

```text
session/
  manifest.json
  evaluations/<id>/input.xyz, run.log, gradient, xtbrestart, metadata.json
  snapshots/<id>/restart.bin, metadata.json
  accepted-history.json
```

既存exeとの初版連携では、wrapperが座標とModel identityをsidecarに保存し、従来restartの対応エンジンも記録する。本体のunformatted binaryをPythonで直接加工する方法は、以前の実験に限定する。製品用の予測電荷入力は本体側の検証付きI/Oとして追加する。

将来のCLI拡張案は `--restart-in`、`--restart-out`、`--result-json`。入力snapshotは読み取り専用とし、出力snapshotは別のファイルへ保存する。形式version、明示された単位、サイズ・checksum検査、収束状態を機械可読にする。これらのオプションは現在未実装である。

正常終了、明示的なSCC収束、有限なenergy・全Gradient、期待する配列サイズ、restartの整合性を確認した場合だけsnapshotを公開する。未収束・途中終了の出力は親を置き換えない。保存は同一filesystem内の一時ファイルから原子的に公開し、snapshotは不変とする。評価と採用履歴の更新は別のtransactionにする。

失敗時はseedを通常warm start、標準guessの順に変えて限定回数再試行する。失敗が残ればoptimizerへ明示的に返し、未収束のenergy・Gradientを成功値として返さない。

## API／常駐worker backend

Calculator、basisなどの静的状態を再利用する。採用済み／選択元Resultsから試行Resultsをコピーして座標を更新し、試行用Resultsでsinglepointを実行する。失敗しても選択元Resultsを更新しない。xtb-pythonの `singlepoint(res, copy=True)` はこのコピー方式に対応するが、本forkの共有ライブラリに接続したbindingの検証が別途必要である。

shell電荷のget/set、Model identityの検証、収束・残差の機械可読取得をAPI拡張として設ける。汎用ラッパーからFortran内部配列を直接操作しない。

並列評価では各workerに独立したCalculator、Environment、Molecule、Resultsを持たせる。共有の可変Resultsを複数スレッドから使わない。初版は独立プロセスworkerを優先し、スレッド安全性を仮定しない。完了順でaccepted historyを更新しない。worker数と各workerのOpenMPスレッド数を合わせてCPU資源を割り当てる。

数値Hessianの正負変位は、共通の基準点からそれぞれseedする。変位点を順番につないでrestart履歴にしない。

## 精度と検証

SCC収束条件、有限温度occupation、最終full-density確認は維持する。履歴によって複数のSCF解の別の枝へ収束する可能性は、残差検査だけでは排除できない。低gap・状態変化が疑われるケースでは標準guessとの比較を行い、対象電子状態の方針を確認する。

最適化全体を、cold、単純warm、近傍選択、予測ありで比較する。成功した点だけでなく、予測準備、コピー、I/O、再試行を含む全SCC反復数と総wall timeを測定する。最終energy・Gradient・構造、optimizerの採用／棄却回数も比較する。

検証対象は同一座標の反復要求、歩幅・方向が変わる採用軌道、棄却後のrollback、採用通知のないcallback、有限差分変位、順不同の並列評価、再開checkpoint、未収束・破損ファイル、元素順序・電荷・スピン・温度・溶媒条件の変更を含める。気相とGBSA／ALPBを検証する。前回の19→15→13反復の実験は一つの同方向変位であり、実際の最適化全体の性能保証には使わない。

実装順は、Sessionと不変snapshot・同一座標cache・近傍warm start、次に採用通知と守られた電荷予測、その後API／常駐worker backendと高度な応答再利用を推奨する。

参考: [geomeTRICのEngineと最適化フロー](https://geometric.readthedocs.io/en/latest/how-it-works.html)、[xtb-python Resultsとsinglepointの仕様](https://xtb-python.readthedocs.io/en/latest/general-api.html)。
