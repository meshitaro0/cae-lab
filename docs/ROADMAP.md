# 学習ロードマップ

学習の階層は **Part → Phase** とします。各Partは独立した一連の学習として扱い、Phase番号はPartの中で付けます。パスの `<part>` には `mechanics` や `cad` の識別名を使います。

| Part | 状態 | 入口 |
| --- | --- | --- |
| 固体力学・CAE Part（mechanics） | Phase 1〜6修了 | [教材一覧](mechanics/README.md)・[詳細ロードマップ](mechanics/ROADMAP.md) |
| CAD Part（cad） | 次の学習テーマ。使用ソフト・到達目標・構成を検討 | 初稿作成時に専用ディレクトリとロードマップを作成 |

## 教材の配置

- `docs/<part>/texts/`：本文。
- `docs/<part>/feedback/`：採点・詳細フィードバック。
- `docs/<part>/learning-log.md`：学習者の感想・質問・反映判断。
- `docs/<part>/ROADMAP.md`：Part内の順序・到達条件・現在地。
- `docs/<part>/glossary.md`：Part内の用語集。必要になった時点で作成。
- `notebooks/<part>/`：答案欄付き演習Notebook。

CAD初稿の作成時に `docs/cad/` と `notebooks/cad/` を作ります。固体力学・CAE PartのPhase 7としては扱いません。
共通の設計原則は [PROFILE.md](PROFILE.md)、運用は [TEACHING_LOOP.md](TEACHING_LOOP.md)、ログの入口は [learning-log.md](learning-log.md) に置きます。
