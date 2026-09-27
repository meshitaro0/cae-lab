# mechlite の使い方と移植範囲

`src/mechlite/` は、Phase 1〜6 の演習から抽出した NumPy ベースの計算ライブラリです。
提出済みNotebookは答案・当時の実装・実験の記録として保持します。今後の再利用コードはこのライブラリ側で管理し、新しいNotebookやスクリプトから import します。
学習用の小規模計算が対象で、実ソルバの代替や自動的な安全性判定は行いません。

## セットアップと検証

リポジトリ直下で実行します。ライブラリの実行依存は NumPy のみです。

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

テストは標準ライブラリの unittest を使い、インストールした `mechlite` を読み込みます。
Notebookとの比較では関数定義のみをメモリ上で読み込み、描画や長い演習セルは実行しません。
Notebookそのものの全セル再実行とは別の検証です。

## モジュールと移植元

| モジュール | 公開API | 移植元・整理内容 |
| --- | --- | --- |
| `tensors` | `strain_to_voigt`, `voigt_to_strain`, `stress_to_voigt`, `voigt_to_stress`, `invariants`, `mises` | Phase 1のVoigt変換、Phase 5・6の不変量計算。Misesを一本化 |
| `elasticity` | `IsotropicElastic`, `elastic_stress`, `elasticity_matrix_3d`, `stress_from_strain` | Phase 1・5の等方弾性。E・nuと派生係数を材料クラスに集約 |
| `members` | `axial_response`, `cantilever_tip_load`, `euler_buckling` | Phase 3の軸・曲げ・Euler座屈公式 |
| `bar` | `bar_element`, `assemble_bar`, `solve_dirichlet`, `recover_bar`, `uniform_bar_mesh` | Phase 4の要素積分・組立・拘束・回復。メッシュ生成を荷重指定から分離 |
| `plasticity` | `J2Material`, `J2State`, `j2_update`, `accumulate_alpha`, `integrate_strain_path` | Phase 5の材料点更新・累積量・履歴積分 |
| `review` | `checked_history`, `audit_history`, `mesh_changes` | Phase 6の台帳検証・エネルギー／力積残差・メッシュ間差 |

純粋な計算はモジュール関数にし、材料定数と材料点履歴だけをクラス化しました。
関数ごとの入出力形状、単位、前提は docstring でも確認できます。

## 使い方

```python
import numpy as np
from mechlite import IsotropicElastic, J2Material, J2State
from mechlite.tensors import mises
from mechlite.bar import uniform_bar_mesh, assemble_bar, solve_dirichlet, recover_bar
from mechlite.plasticity import integrate_strain_path

# 3D等方弾性：ひずみは無次元、Eと応力はMPa
elastic = IsotropicElastic(E=210000.0, nu=0.3)
sigma = elastic.stress(np.diag([1e-4, 0.0, 0.0]))
print(mises(sigma))

# 1D棒：座標mm、荷重N、断面積mm²
x, conn = uniform_bar_mesh(L=100.0, ne=4)
E, A, p = np.full(4, 210000.0), np.full(4, 100.0), np.zeros(4)
loads = np.zeros(5)
loads[-1] = 1000.0
K, f = assemble_bar(x, conn, E, A, p, loads)
d, residual = solve_dirichlet(K, f, {0: 0.0})
strain, stress, axial_force = recover_bar(x, conn, E, A, d)

# 試行更新を採用するかは呼び出し側で決める
material = J2Material(E=210000.0, nu=0.3, sigma_y0=250.0, H=1000.0)
committed = J2State()
trial, dalpha = material.update(np.diag([0.003, 0.0, 0.0]), committed)
committed = trial

# 経路の始点に対応する初期状態を渡す。既定値は応力・履歴ゼロ
path = np.array([0.0, 0.004, -0.004, 0.0])[:, None, None] * np.diag([1., -.5, -.5])
sigma_h, ep_h, alpha_h, da_h = integrate_strain_path(path, material, subdivisions=40)
```

J2Stateは入力配列をコピーし、通常の要素代入を禁止します。試行更新で確定状態を上書きしません。
履歴積分は各折れ線区間を指定数で分割し、初期値込みで応力・塑性ひずみ・累積量を返します。
初期状態と経路始点の物理的整合性は呼び出し側の責任です。

## 元実装からの修正・変更

- **せん断応力変換の修正**：Phase 1の `stress_from_strain` は応力Voigtベクトルを工学ひずみ用の変換へ渡しており、せん断応力が半分になる。ライブラリは応力用・ひずみ用の逆変換を区別し、弾性応力はテンソル式で計算する。旧名 `voigt_to_tensor` は意味が曖昧なため採用しない。答案自体は変更しない。
- **Voigt規約の固定**：順序は `[xx, yy, zz, yz, xz, xy]`。ひずみのみせん断成分を2倍にする。応力とひずみの内積が仕事に一致することをテストする。
- **入力検証**：実数・有限値・形状・正値条件・対称性を共通化。対称性は絶対差 `1e-12`、相対許容差0。単一テンソル `(3,3)` が基本で、バッチ入力は履歴関数だけが受ける。
- **FEMの堅牢化**：EAと分布荷重は各Gauss点で一度ずつスカラー評価。NaNを拒否し、整数荷重でも浮動小数点で組み立てる。非ゼロ拘束、全拘束、接続・自由度番号を検証する。拘束不足などで特異な自由剛性は NumPy の `LinAlgError` を返す。
- **メッシュと荷重の分離**：`generate_mesh` の座標近傍への荷重丸め込みは採用しない。`uniform_bar_mesh` は節点・接続だけを返し、集中荷重は節点番号で明示する。
- **塑性の依存関係**：Notebookの単位行列・材料定数のグローバル変数を除去。`history_driver` と `history_driver_general` は任意の3Dひずみ折れ線を扱う `integrate_strain_path` に統合。旧 `j2_update` の引数・戻り値形式は維持するが、弾性時も塑性ひずみを独立コピーで返す。
- **許容値**：J2更新の弾性判定は元の `1e-7 MPa` 相当を既定値とし、`J2Material.yield_tolerance` で指定可能。単位系を変える場合は同じ応力単位へ換算する。
- **部材公式**：片持ち梁は元docstringの契約どおり応力・たわみの絶対値を返す。負の荷重でも大きさは正になる。
- **収支計算**：Phase 6の台帳規約と正規化を保持。非ゼロ初期K/U/P、非等間隔時刻に対応。特定ソルバのエネルギー列対応や採否基準を自動推定しない。

## 今回抽出しない処理

- Phase 2のセル内計算は、既存のNumPy線形代数と今回の不変量APIで扱えるため専用クラスを増やさない。
- Phase 3の `recover_section_resultants` は任意発展の未実装関数。移植済みとして公開しない。
- Phase 4の `analytical_u` / `analytical_stress` / `evaluate_error` は特定境界条件の解析解に依存するため、汎用APIにしない。補間ヘルパーも右端の扱いが未整備。場の収束性は独立したテストで確認する。
- `checkpoint`、`check_ex*`、`check_when_ready`、描画、Phase 6の `force_signal`、演習データ生成は学習・検証シナリオ側へ残す。

有限変形、整合接線、全体非線形FEM、疎行列、梁・シェル要素、CAD入出力は今回の範囲外です。

## リポジトリ構成（2026-09-27更新）

| 項目 | 状態 | 配置・運用 |
| --- | --- | --- |
| Part単位の教材・演習 | 実施 | [固体力学・CAE Part](mechanics/README.md)を `docs/mechanics/` と `notebooks/mechanics/` に集約。CADは初稿作成時に `docs/cad/` と `notebooks/cad/` を新設 |
| フィードバック | 実施 | [mechanics/feedback/](mechanics/feedback/) に採点文書を移動。本文は `mechanics/texts/` |
| 学習ログ | 実施 | [Part別ログ](mechanics/learning-log.md)と[全体の索引](learning-log.md)に分離。過去の記録を保持し、参照パスを更新 |
| 依存宣言 | 実施 | `pyproject.toml` にNumPyと `notebooks` extra（SciPy・Matplotlib）を集約。`requirements.txt` は `-e .[notebooks]` を参照 |
| CI | 設定追加 | [GitHub Actions](../.github/workflows/tests.yml)でpush・PR・手動実行時にLinux／Windows・Python 3.12の通常インストールとunittestを実行。GitHub上の初回実行はpush後に確認 |
| パッケージ配置 | 実施 | `src/mechlite/` に移動し、setuptoolsの探索先を更新。テストはインストール済みパッケージを利用 |

CADエンジンを使う場合は、その依存を力学のNumPy計算へ混在させず、形状処理・メッシュへの受け渡しを担うモジュールを分けるのが適切です。具体的なCADソフトと学習範囲は次の教材設計時に決めます。
