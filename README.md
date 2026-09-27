# cae-lab

CAE の結果を力学から読み解き、最小のソルバ実装で確かめる学習リポジトリです。

教材は公式の復習ではなく、FEMで目にする応力・ひずみ・主応力・塑性量を起点に、連続体力学と離散化へ逆算して理解するよう設計しています。

- [学習者プロフィール](docs/PROFILE.md)
- [学習ロードマップ](docs/ROADMAP.md)
- [教材改善ループ](docs/TEACHING_LOOP.md)
- [学習ログ](docs/learning-log.md)
- [用語集](docs/glossary.md)
- [計算ライブラリ mechlite：利用例・Notebook対応表・構成改善案](docs/mechlite.md)
- [Phase 1: CAE結果を読むための材料力学](docs/texts/phase1-stress-strain.md)
- [Phase 5: J2塑性と材料点更新](docs/texts/phase5-j2-plasticity.md) / [演習](notebooks/Phase5.ipynb)
- [Phase 6: CAE結果レビューと検証計画](docs/texts/phase6-cae-review.md) / [演習](notebooks/Phase6.ipynb)

## Setup

```bash
pip install -e .
```

`libs/mechlite/` にはテンソル変換、等方弾性、部材公式、1D棒FEM、J2塑性、結果レビューをまとめています。ライブラリはNumPyのみで動作します。既存Notebookの描画などには `pip install -r requirements.txt` を使います。

```python
from mechlite import IsotropicElastic
from mechlite.tensors import mises

material = IsotropicElastic(E=210000.0, nu=0.3)
stress = material.stress([[0.001, 0, 0], [0, 0, 0], [0, 0, 0]])
print(mises(stress))
```

回帰・物理検証テスト：リポジトリ直下で `python -m unittest discover -s tests -v`。
提出Notebookは学習記録として保持し、新しい計算ではライブラリを利用します。
