# 学習ログの索引

Partごとの学習・採点・教材改訂は、そのPartのログへ記録します。Partをまたぐ構成変更はこのファイルに記録します。

| Part | 学習ログ | 採点・詳細フィードバック |
| --- | --- | --- |
| 固体力学・CAE Part（Phase 1〜6） | [mechanics/learning-log.md](mechanics/learning-log.md) | [mechanics/feedback/](mechanics/feedback/) |
| CAD Part | 初稿作成時に `docs/cad/learning-log.md` を新設 | 必要時に `docs/cad/feedback/` を新設 |

## 2026-09-27 — 学習単位の名称をPartへ統一

- **学習者の指定（表現）**：Phaseを束ねる学習単位の表現をPartに統一する。
- **反映**：AGENTS、README、全体とmechanicsのロードマップ、教材改善ループ、ライブラリ案内、ログの説明を統一。階層を **Part → Phase**、パスのプレースホルダーを `<part>` と明記した。過去の学習者発言を直接引用した箇所は原文を保持する。
- **変更範囲・次回**：PROFILE、教材本文、演習、実ディレクトリ名は変更なし。以後の教材作成・報告ではPartを使用する。名称統一に関する保留事項はない。

## 2026-09-27 — 教材作業のGit操作タイミングをルール化

- **学習者の指定（運用）**：Phase・Partの初稿作成前に新しいブランチを作る。初稿完成直後にはコミットせず、初稿レビュー完了または演習開始の連絡を受けてコミットする。解答と学習ログの記録完了、再提出の修正完了の連絡を受けたら、先にコミットしてから採点・確認・レビュー・教材修正へ進む。Phase完了時は完了作業後に差分を提示し、学習者の確認を得てからコミット・push・PR作成を行う。
- **反映**：正本を [AGENTS.md](../AGENTS.md#ブランチコミットpushprのタイミング) に追加し、[教材改善ループ](TEACHING_LOOP.md)へ順序を反映。既にコミット済みならそのIDを使い、採点前の空コミットを作らない。Phase完了時は未コミット差分とPR全体の変更概要を提示する。
- **変更範囲**：運用文書と本ログを更新。教材本文・演習・PROFILE・ROADMAPは変更なし。Partごとに複数Phaseを持つ構成を維持する。
- **次回確認**：次の初稿レビュー・提出・再提出・Phase完了時に、この順序で運用する。今回のルール編集自体は教材初稿やPhase完了に該当せず、自動コミット・pushの対象としない。

## 2026-09-27 — Part単位の教材構成・依存管理・CI

- **学習者の希望（範囲・運用）**：CADをPhase 1〜6と同じ粒度の独立したPartにしたい。既存教材・NotebookをPartごとにまとめ、リンクを付け替える。フィードバックとログもPartごとに分離し、`libs` を `src` へ変更する。extras・CIの目的と実装はチャットで説明し、CIを学びながら導入したい。
- **判断・変更**：`docs/mechanics/` の下に本文・フィードバック・ログ・ロードマップ・用語集を集約し、演習は `notebooks/mechanics/` へ移動。上位ロードマップと本ログをPartの索引に変更。過去の記録は省略せず移し、Notebookは参照パスのみ変更する。
- **実装**：`src/mechlite/` へ移動し、パッケージ探索設定を更新。SciPy・Matplotlibを `pyproject.toml` の `notebooks` extraにまとめ、`requirements.txt` はその参照に変更。GitHub Actionsでpush・PR・手動実行時にPython 3.12のLinux／Windowsでパッケージをインストールし、既存unittestを実行する設定を追加。
- **教材・運用への反映**：README、Partの入口、AGENTS、TEACHING_LOOP、ライブラリ案内を新配置へ更新。PROFILEの学習者特性・設計原則は変更なし。新しいCAD教材は初稿作成時に作る。
- **検証**：131個のローカルリンク先を確認。6冊のNotebookはMarkdownの参照パスのみ変更され、コード・出力・メタデータは移動前と一致。8個のライブラリソースも一致。既存 `.venv` へのeditable install成功後に17件のunittestが成功。GitHub Actions上の通常インストール・Linux／Windowsジョブは未実行で、push後に確認する。
- **保留・次回**：CADの使用ソフトと到達目標を決める。GitHub上でのCI初回実行は、変更をpushした後に確認する。自動公開・デプロイは今回の設定に含めない。
