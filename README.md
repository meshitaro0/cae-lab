# cae-lab

[![tests](https://github.com/meshitaro0/cae-lab/actions/workflows/tests.yml/badge.svg)](https://github.com/meshitaro0/cae-lab/actions/workflows/tests.yml)

CAE の結果を力学から読み解き、最小のソルバ実装で確かめる学習リポジトリです。

教材は公式の復習ではなく、FEMで目にする応力・ひずみ・主応力・塑性量を起点に、連続体力学と離散化へ逆算して理解するよう設計しています。

## 共通

- [トラック一覧](docs/ROADMAP.md)
- [学習者プロフィール](docs/PROFILE.md)
- [教材改善ループ](docs/TEACHING_LOOP.md)
- [用語集](docs/glossary.md)

## 固体力学トラック（Phase 1〜6 完了）

[ロードマップ](docs/solid-mechanics/ROADMAP.md) / [学習ログ](docs/solid-mechanics/learning-log.md)

| Phase | 本文 | 演習 |
| --- | --- | --- |
| 1 | [CAE結果を読むための材料力学：応力・ひずみは何を表すのか](docs/solid-mechanics/texts/phase1-stress-strain.md) | [Phase1.ipynb](notebooks/solid-mechanics/Phase1.ipynb) |
| 2 | [テンソル・主値・不変量：von Mises 応力を座標に依らず読む](docs/solid-mechanics/texts/phase2-tensors-invariants.md) | [Phase2.ipynb](notebooks/solid-mechanics/Phase2.ipynb) |
| 3 | [材料力学の再導出：3次元応力から部材の軸・曲げ・ねじり・座屈へ](docs/solid-mechanics/texts/phase3-member-mechanics.md) | [Phase3.ipynb](notebooks/solid-mechanics/Phase3.ipynb) |
| 4 | [最小FEM：釣合いを節点変位の連立方程式へ変える](docs/solid-mechanics/texts/phase4-minimal-fem.md) | [Phase4.ipynb](notebooks/solid-mechanics/Phase4.ipynb) |
| 5 | [相当塑性ひずみから材料点更新へ](docs/solid-mechanics/texts/phase5-j2-plasticity.md) | [Phase5.ipynb](notebooks/solid-mechanics/Phase5.ipynb) |
| 6 | [CAEの出力から、判断の根拠と検証計画へ](docs/solid-mechanics/texts/phase6-cae-review.md) | [Phase6.ipynb](notebooks/solid-mechanics/Phase6.ipynb) |

採点・フィードバックの記録は [docs/solid-mechanics/reviews/](docs/solid-mechanics/reviews/) にあります。

## CADトラック（Phase 1 完了）

[ロードマップ](docs/cad/ROADMAP.md) / [学習ログ](docs/cad/learning-log.md)

受け取ったCAD形状の表現（B-spline・NURBS・B-rep）を自作実装で確かめ、パラメトリック形状・メッシュ化・CAE向け形状処理・自動化パイプラインを経て、等幾何解析で固体力学トラックと合流する全8 Phaseの構成です。

| Phase | 本文 | 演習 |
| --- | --- | --- |
| 1 | [曲線：STEPの `B_SPLINE_CURVE_WITH_KNOTS` を自分で評価する](docs/cad/texts/phase1-curves.md) | [Phase1.ipynb](notebooks/cad/Phase1.ipynb) |

## Setup

```bash
pip install -e ".[notebooks]"          # ライブラリ本体 + Notebook用（matplotlib, ipykernel）
python -m unittest discover -s tests -v  # ライブラリのテストと文書リンクの検査
```

ライブラリだけを使うなら `pip install -e .` で足ります。テストは標準ライブラリの unittest なので追加のインストールは不要です。

## mechlite

演習Notebookで実装した関数を、再利用できる形に整理したライブラリです（[src/mechlite/](src/mechlite/)）。Notebookは提出物として当時のまま残し、ライブラリは独立に保守します。テストは提出Notebookの関数定義を読み込み、ライブラリと結果が一致するかも比較します。

| モジュール | 公開API | 由来 |
| --- | --- | --- |
| `tensors` | `stress_to_voigt`, `strain_to_voigt`, `voigt_to_stress`, `voigt_to_strain`, `invariants`, `mises` | Phase 1, 2, 5, 6 |
| `elasticity` | `IsotropicElastic`, `elastic_stress`, `elasticity_matrix_3d`, `stress_from_strain` | Phase 1, 5 |
| `members` | `axial_response`, `cantilever_tip_load`, `euler_buckling` | Phase 3 |
| `bar` | `bar_element`, `assemble_bar`, `solve_dirichlet`, `recover_bar`, `uniform_bar_mesh` | Phase 4 |
| `plasticity` | `J2Material`, `J2State`, `j2_update`, `accumulate_alpha`, `integrate_strain_path` | Phase 5 |
| `review` | `checked_history`, `audit_history`, `mesh_changes` | Phase 6 |

各関数の入出力の形状・単位・前提は docstring にあります。Notebookから移さなかったのは、特定の境界条件に依存する解析解・誤差評価（Phase 4）、未実装の任意課題（Phase 3 `recover_section_resultants`）、採点用チェック・描画・演習データ生成です。

```python
import numpy as np
from mechlite import J2Material, J2State

mat = J2Material(E=210e3, nu=0.3, sigma_y0=250.0, H=1000.0)
state, dalpha = mat.update(np.diag([0.003, 0.0, 0.0]), J2State())
```
