# GFN1-fast 4.6.0 results

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
