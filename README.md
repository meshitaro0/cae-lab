# cae-lab

CAE の結果を力学から読み解き、最小のソルバ実装で確かめる学習リポジトリです。

教材は公式の復習ではなく、FEMで目にする応力・ひずみ・主応力・塑性量を起点に、連続体力学と離散化へ逆算して理解するよう設計しています。

- [学習者プロフィール](docs/PROFILE.md)
- [学習ロードマップ（Partの一覧）](docs/ROADMAP.md)
- [教材改善ループ](docs/TEACHING_LOOP.md)
- [学習ログの索引](docs/learning-log.md)
- [固体力学・CAE Part：Phase 1〜6の教材・演習・フィードバック](docs/mechanics/README.md)
- [固体力学・CAE Partの用語集](docs/mechanics/glossary.md)
- [計算ライブラリ mechlite：利用例・Notebook対応表・構成改善案](docs/mechlite.md)

CAD Partは独立したPartとして、初稿作成時に `docs/cad/` と `notebooks/cad/` を追加します。

## Setup

```bash
python -m pip install -e .
```

`src/mechlite/` にはテンソル変換、等方弾性、部材公式、1D棒FEM、J2塑性、結果レビューをまとめています。ライブラリはNumPyのみで動作します。演習用のSciPy・Matplotlibも入れる場合は `python -m pip install -e ".[notebooks]"` を使います。`python -m pip install -r requirements.txt` も同じ依存宣言を参照します。コマンドはリポジトリ直下で、利用するNotebookカーネルと同じPython環境で実行してください。

```python
from mechlite import IsotropicElastic
from mechlite.tensors import mises

material = IsotropicElastic(E=210000.0, nu=0.3)
stress = material.stress([[0.001, 0, 0], [0, 0, 0], [0, 0, 0]])
print(mises(stress))
```

回帰・物理検証テスト：上記のインストール後、リポジトリ直下で `python -m unittest discover -s tests -v`。
提出Notebookは学習記録として保持し、新しい計算ではライブラリを利用します。

[GitHub Actionsの設定](.github/workflows/tests.yml)はpush・PR・手動実行に対応し、Linux／WindowsのPython 3.12で通常インストールとテストを行います。
