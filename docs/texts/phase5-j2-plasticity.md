# Phase 5 — 相当塑性ひずみから材料点更新へ

Phase 4では、節点変位からひずみを求め、弾性行列を掛けて応力へ進んだ。塑性では、同じ現在ひずみでも過去の負荷によって応力が変わる。この履歴を積分点に保存し、次の増分へ引き継ぐことが今回の中心である。

- 前提：[Phase 2の不変量](phase2-tensors-invariants.md)、[Phase 4の最小FEM](phase4-minimal-fem.md)。
- 演習・答案：[Phase5.ipynb](../../notebooks/Phase5.ipynb)。本文には答案欄を置かない。
- 到達点：降伏条件・流れ則・硬化則の役割を分けて説明し、三次元J2材料点更新をNumPyで実装・検証できる。
- 必須範囲：微小ひずみ、等方線形弾性、関連流れ則、線形等方硬化、速度非依存、等温、損傷なし。全体非線形FEMの実装と有限変形塑性は発展範囲とする。

## 1. 除荷したのに塑性ひずみが残るのはなぜか

CAEで荷重を下げると、Mises応力は低下しても相当塑性ひずみが残ることがある。応力は現在の弾性変形に対応し、相当塑性ひずみは過去の塑性変形の蓄積を表すためである。塑性変形したことと、今この瞬間も塑性流動していることを区別しよう。

全ひずみ $\boldsymbol\varepsilon$、弾性ひずみ $\boldsymbol\varepsilon^e$、塑性ひずみ $\boldsymbol\varepsilon^p$ は対称な二階テンソル（無次元）で、微小ひずみでは

$$
\boldsymbol\varepsilon=\boldsymbol\varepsilon^e+\boldsymbol\varepsilon^p,
\qquad \boldsymbol\sigma=\mathbb C:(\boldsymbol\varepsilon-\boldsymbol\varepsilon^p).
$$

$\boldsymbol\sigma$ はCauchy応力（二階テンソル、MPa）、$\mathbb C$ は弾性剛性（四階テンソル、MPa）。コロンは二重縮約で、二階同士なら $\mathbf A:\mathbf B=\sum_{i,j=1}^3 A_{ij}B_{ij}$、四階と二階なら $(\mathbb C:\mathbf A)_{ij}=\sum_{k,l=1}^3 C_{ijkl}A_{kl}$ とする。添字は三つの直交座標成分を表す。

弾性だけなら現在の全ひずみで応力が決まる。塑性では $\boldsymbol\varepsilon^p$ が追加で必要になる。これが[内部変数](../glossary.md#内部変数)を保存する理由である。微視的な転位構造などをこの教材で直接計算するわけではなく、その履歴の効果を少数の変数でモデル化する。

### 等方弾性を体積変化と形状変化へ分ける

単位テンソル $\mathbf I$、トレース演算 $\operatorname{tr}\mathbf A=\sum_i A_{ii}$、偏差演算 $\operatorname{dev}\mathbf A=\mathbf A-\operatorname{tr}(\mathbf A)\mathbf I/3$ を使う。体積弾性率 $K$ とせん断弾性率 $G$ はスカラーである。

$$
\boldsymbol\sigma=K\operatorname{tr}(\boldsymbol\varepsilon^e)\mathbf I
+2G\operatorname{dev}(\boldsymbol\varepsilon^e),\qquad
K=\frac{E}{3(1-2\nu)},\quad G=\frac{E}{2(1+\nu)}.
$$

$E>0$ はYoung率（MPa）、$-1<\nu<1/2$ はPoisson比（無次元）。$K$ はPhase 3の有効長係数とは異なる。体積ひずみを $\operatorname{tr}\boldsymbol\varepsilon^e$ と定義すると、その係数が $K$ となる。等方性により体積成分と偏差成分を独立に応答させ、一軸応力でのYoung率と横収縮、および純せん断での応答に合わせると上の係数関係になる。

演習は常に3×3配列を使う。対称行列の両側のせん断成分を含めて縮約するので、せん断項は2回数える。工学せん断ひずみ $\gamma_{xy}=2\varepsilon_{xy}$ を、そのまま行列の $xy$ 成分へ入れない。

## 2. 降伏条件だけでは次の応力は決まらない

平均応力 $p=\operatorname{tr}\boldsymbol\sigma/3$（引張正、MPa）、偏差応力 $\mathbf s=\boldsymbol\sigma-p\mathbf I$（二階テンソル）を定義する。ここで $p$ は圧縮正の圧力ではない。相当応力 $q$ と偏差応力の第2不変量 $J_2$ はスカラーで、Phase 2のMises応力とつながる。

$$
J_2=\frac12\mathbf s:\mathbf s,\qquad
q=\sqrt{3J_2}=\sqrt{\frac32\mathbf s:\mathbf s}.
$$

一軸応力 $\operatorname{diag}(\sigma,0,0)$ では $q=|\sigma|$、純せん断応力 $\sigma_{xy}=\sigma_{yx}=\tau$ では $q=\sqrt3|\tau|$、静水圧応力 $p\mathbf I$ では $q=0$。$\operatorname{diag}$ は指定成分を対角に並べる演算、$\sigma,\tau$ はここだけのスカラー応力である。

[降伏関数](../glossary.md#降伏関数) $f$（MPa）を

$$
f(\boldsymbol\sigma,\alpha)=q-\sigma_y(\alpha),\qquad
\sigma_y(\alpha)=\sigma_{y0}+H\alpha
$$

とする。$\alpha$ は累積相当塑性ひずみ（非負の無次元スカラー）、$\sigma_{y0}>0$ は初期降伏応力、$H\ge0$ は塑性ひずみに対する硬化係数（ともにMPa）。$\sigma_y$ は現在の降伏応力。$f<0$ は降伏面の内側、$f=0$ は面上、$f>0$ はこの速度非依存モデルで許容しない応力である。

平均応力を足しても $q$ は変わらない。主応力空間では静水圧軸に沿う円柱状の降伏面となる。圧力で強度が変わる材料や引張・圧縮で降伏が大きく違う材料には、この仮定をそのまま使えない。

降伏条件は「許容範囲」を決めるだけである。塑性ひずみがどちらへ、どれだけ進むかを決めるため、次に流れ則と履歴の更新則が要る。

## 3. 流れ則から相当塑性ひずみを導く

時間 $t$ による微分を上付きドットで表す。塑性乗数 $\dot\lambda\ge0$ は塑性流動の速さを決めるスカラー（1/時間）。[関連流れ則](../glossary.md#関連流れ則)は降伏関数の応力勾配を塑性ひずみ速度の方向に取る。

$$
\dot{\boldsymbol\varepsilon}^p=\dot\lambda\frac{\partial f}{\partial\boldsymbol\sigma}.
$$

$q>0$ で $q^2=3\mathbf s:\mathbf s/2$ を微分すると、$2q\,dq=3\mathbf s:d\mathbf s$。$\operatorname{tr}\mathbf s=0$ より $\mathbf s:d\mathbf s=\mathbf s:d\boldsymbol\sigma$ だから

$$
\mathbf n:=\frac{\partial f}{\partial\boldsymbol\sigma}=\frac{3\mathbf s}{2q},\qquad
\dot{\boldsymbol\varepsilon}^p=\dot\lambda\mathbf n.
$$

$\mathbf n$ は無次元の二階テンソルであり、空間内の面法線ベクトルではない。$\mathbf n:\mathbf n=3/2$ なので単位ノルムでもない。トレースはゼロで、$\operatorname{tr}\dot{\boldsymbol\varepsilon}^p=0$、すなわち塑性体積変化がないことが導ける。弾性による体積変化は残る。

[累積相当塑性ひずみ](../glossary.md#累積相当塑性ひずみ)の定義を

$$
\dot\alpha=\sqrt{\frac23\dot{\boldsymbol\varepsilon}^p:\dot{\boldsymbol\varepsilon}^p},\qquad
\alpha(t)=\alpha(0)+\int_0^t\dot\alpha\,dt
$$

とすると、流れ則を代入して $\dot\alpha=\dot\lambda$ が得られる。この係数は、一軸塑性引張で $\dot\alpha=|\dot\varepsilon^p_{xx}|$ となるように選んでいる。単調引張では $\dot{\boldsymbol\varepsilon}^p=\dot\lambda\operatorname{diag}(1,-1/2,-1/2)$ である。

負荷反転すると塑性ひずみの成分は減少し得るが、$\alpha$ は減少しない。一般に $\alpha\ne\sqrt{2\boldsymbol\varepsilon^p:\boldsymbol\varepsilon^p/3}$。右辺は現在のテンソルの大きさ、左辺はそこへ至る経路の長さに相当する。向きが固定された単調な塑性経路をゼロ状態から進む場合には一致する。

### 降伏面上でも除荷は弾性になり得る

塑性の条件を相補条件としてまとめる。

$$
f\le0,\quad\dot\lambda\ge0,\quad\dot\lambda f=0.
$$

塑性流動中は $f=0$ を保ち、$\dot f=0$（整合条件）。面上でも内向きの負荷なら $\dot\lambda=0$ で弾性除荷する。したがって「$q$ が降伏応力に等しい」だけで現在の塑性増分を決められない。

塑性仕事率密度は $\boldsymbol\sigma:\dot{\boldsymbol\varepsilon}^p=q\dot\alpha\ge0$（MPa/時間）となる。ただし全量が直ちに熱へ変わるとは限らない。例えば硬化の保存エネルギー密度を $H\alpha^2/2$ とするモデルなら、塑性流動中の散逸率は $\sigma_{y0}\dot\alpha$ であり、硬化への蓄積分を差し引く。

## 4. 硬化係数は応力–全ひずみ曲線の傾きと同じか

等方硬化は降伏面を拡大する。$H=0$ は完全塑性、$H>0$ は線形等方硬化である。移動硬化は背応力（二階テンソル）により降伏面の中心を移動させ、逆方向降伏の変化を表す。繰返し変形のBauschinger効果を等方硬化だけで再現しようとしない。

一軸**応力**の単調引張に限定すると、軸方向スカラーについて

$$
\varepsilon=\frac\sigma E+\varepsilon^p,\qquad
\alpha=\varepsilon^p,\qquad d\sigma=H\,d\varepsilon^p.
$$

代入すると $d\varepsilon=(1/E+1/H)d\sigma$。塑性域の全ひずみに対する接線係数 $E_t=d\sigma/d\varepsilon$（MPa）は

$$
E_t=\frac{EH}{E+H},\qquad H=\frac{EE_t}{E-E_t}\quad(0\le E_t<E).
$$

$H=0$ は極限で $E_t=0$。測定曲線の傾き $E_t$ をそのまま $H$ として入力しない。弾性除荷の傾きはこのモデルでは $E$ へ戻る。

一軸**ひずみ** $\boldsymbol\varepsilon=\operatorname{diag}(e,0,0)$ では横応力が生じる。処女材の弾性域で $q=2G|e|$ となり、初回降伏は $|e|=\sigma_{y0}/(2G)$。一軸応力の初回降伏ひずみ $\sigma_{y0}/E$ や接線係数 $E_t$ と取り違えない。

## 5. 弾性予測から降伏面へ戻す

全ひずみ増分 $\Delta\boldsymbol\varepsilon=\boldsymbol\varepsilon_{n+1}-\boldsymbol\varepsilon_n$ を与える。添字 $n,n+1$ は確定済みと次の増分端の時点、$\Delta$ は差分である。入力履歴は $\boldsymbol\sigma_n,\boldsymbol\varepsilon^p_n,\alpha_n$。整合した初期状態から始める。

### 5.1 まず塑性ひずみを固定する

塑性増分がゼロならどうなるかを計算する。上付き $\mathrm{tr}$ はtrial（試行）で、トレース演算とは異なる。

$$
\boldsymbol\sigma^{\mathrm{tr}}=\boldsymbol\sigma_n+\mathbb C:\Delta\boldsymbol\varepsilon,
\quad \mathbf s^{\mathrm{tr}}=\operatorname{dev}\boldsymbol\sigma^{\mathrm{tr}},\quad
q^{\mathrm{tr}}=\sqrt{3\mathbf s^{\mathrm{tr}}:\mathbf s^{\mathrm{tr}}/2},
$$
$$
f^{\mathrm{tr}}=q^{\mathrm{tr}}-(\sigma_{y0}+H\alpha_n).
$$

$f^{\mathrm{tr}}\le0$ なら試行応力を採用し、内部変数を変えない。$q^{\mathrm{tr}}=0$ もここで処理できるため、流れ方向のゼロ除算を避けられる。

### 5.2 塑性の場合は増分端で流れ則を満たす

後退Euler法は流れ方向を新しい時点で評価する。$\Delta\alpha=\Delta\lambda$ として

$$
\Delta\boldsymbol\varepsilon^p=\Delta\alpha\frac{3\mathbf s_{n+1}}{2q_{n+1}},\quad
\mathbf s_{n+1}=\mathbf s^{\mathrm{tr}}-2G\Delta\boldsymbol\varepsilon^p.
$$

代入すると $(1+3G\Delta\alpha/q_{n+1})\mathbf s_{n+1}=\mathbf s^{\mathrm{tr}}$。両偏差応力は同じ方向なので、ノルムを取って

$$
q_{n+1}=q^{\mathrm{tr}}-3G\Delta\alpha.
$$

一方、降伏面上の条件は $q_{n+1}=\sigma_{y0}+H(\alpha_n+\Delta\alpha)$。二式を等置して

$$
\boxed{\Delta\alpha=\frac{f^{\mathrm{tr}}}{3G+H}},\qquad
\Delta\boldsymbol\varepsilon^p=\frac{3\Delta\alpha}{2q^{\mathrm{tr}}}\mathbf s^{\mathrm{tr}},
$$
$$
\boldsymbol\sigma_{n+1}=\boldsymbol\sigma^{\mathrm{tr}}-2G\Delta\boldsymbol\varepsilon^p,\quad
\boldsymbol\varepsilon^p_{n+1}=\boldsymbol\varepsilon^p_n+\Delta\boldsymbol\varepsilon^p,\quad
\alpha_{n+1}=\alpha_n+\Delta\alpha.
$$

これが[radial return](../glossary.md#radial-return)である。縮めるのは偏差応力で、平均応力は試行値のまま。応力テンソル全体を同じ比率で縮めると、塑性体積変化なしというモデルを壊す。非線形硬化なら $q^{\mathrm{tr}}-3G\Delta\alpha-\sigma_y(\alpha_n+\Delta\alpha)=0$ というスカラー方程式の局所反復が必要になる場合がある。

この閉形式の更新は、一定弾性係数・線形等方硬化に対する離散方程式の解である。任意の連続負荷経路を一増分で厳密に積分できるという意味ではない。途中の除荷や方向変化を増分端だけから復元することはできない。

### 5.3 実装の契約と検証

演習では `j2_update(deps, sigma_n, ep_n, alpha_n, E, nu, sigma_y0, H)` が `(sigma_new, ep_new, alpha_new, dalpha)` を返す。テンソルは有限・対称な3×3配列、その他はスカラー。入力を変更せず、確定済み履歴と試行状態を分ける。異常な形状・非有限値・範囲外の材料定数は `ValueError` とする。入力状態の構成則・履歴全体の整合性は呼出し側が保証する。

弾塑性分岐は理論上 $f^{\mathrm{tr}}\le0$ で行う。丸め誤差対策で許容値を導入するならMPaのスケールを明記し、降伏残差評価にも反映する。検証は一つの値の一致だけで終えず、次を別々に調べる。

- 静水圧負荷で $\Delta\alpha=0$、純せん断で $q=\sqrt3|\tau|$。
- 塑性ステップで降伏残差 $|q_{n+1}-\sigma_y(\alpha_{n+1})|$ が小さい。
- $\operatorname{tr}\Delta\boldsymbol\varepsilon^p=0$、$\alpha$ が非減少、平均応力が試行値と一致。
- 全ひずみを別に累積して、$\boldsymbol\sigma=\mathbb C:(\boldsymbol\varepsilon-\boldsymbol\varepsilon^p)$ と照合。
- 全テンソルを同じ直交行列で回転したとき、応力・塑性ひずみは同じ規則で回転し、$\alpha$ は不変。
- 同じ入力の再呼出しで同じ結果となり、入力配列は変更されない。

## 6. 積分点更新と全体釣合いをつなぐ

Phase 4の $\mathbf B$（ひずみ–節点変位の行列）、節点変位ベクトル $\mathbf d$、内力・外力ベクトル $\mathbf f_{\mathrm{int}},\mathbf f_{\mathrm{ext}}$ を再利用する。微小変形の固定形状で、材料更新後の応力を、仕事に整合するVoigt規約のベクトル $\boldsymbol\sigma_V$ に並べれば

$$
\mathbf r(\mathbf d)=\mathbf f_{\mathrm{int}}(\mathbf d)-\mathbf f_{\mathrm{ext}},\qquad
\mathbf f_{\mathrm{int}}=\int_\Omega\mathbf B^T\boldsymbol\sigma_V\,d\Omega.
$$

$\Omega$ は解析領域、$\mathbf r$ は残差ベクトル。未知変位を仮定するたびに積分点のひずみが変わり、応力も変わるため、弾性のように一定の剛性を一度解くだけでは済まない。

陰解法のNewton反復では、反復番号 $k$、自由自由度を示す添字 $f$、変位修正ベクトル $\delta\mathbf d_f$ を使い、$\mathbf K_{T,ff}\delta\mathbf d_f=-\mathbf r_f$ を解く。拘束自由度の残差は反力になる。接線剛性行列 $\mathbf K_T$ は、この材料更新をひずみで微分した[整合接線](../glossary.md#整合接線)の行列表現 $\mathbf D_{\mathrm{alg}}$ を使って

$$
\mathbf K_T=\int_\Omega\mathbf B^T\mathbf D_{\mathrm{alg}}\mathbf B\,d\Omega
$$

と組み立てる。$\mathbf D_{\mathrm{alg}}$ は単なる割線 $\sigma/\varepsilon$ でも、常に弾性行列でもない。本Phaseでは役割の理解を必須、導出・全体実装は任意とする。有限差分で確かめる場合も、確定済み状態を固定して入力ひずみ増分だけを摂動する。降伏の分岐点では一意の滑らかな微分を期待しない。

```text
増分 n の確定履歴を保持
  → 次の節点変位を仮定
  → 各積分点の「n からの全ひずみ増分」を求める
  → 毎回 n の履歴から材料更新して試行応力・試行内部変数を作る
  → 内力と自由残差を組み立て、必要なら変位を修正して反復
  → 全体収束後だけ試行履歴を n+1 の確定履歴にする
```

Newton反復ごとに $\alpha$ を確定履歴へ加算すると、同じ負荷増分を何度も塑性変形した扱いになる。陽解法でも材料点更新と履歴管理は必要だが、運動方程式を時間積分するため通常の全体Newton反復とは手順が異なる。材料点の降伏残差、静解析の自由残差、動解析の運動量・エネルギー収支は異なる検証である。

## 7. 数値実験で何を比較するか

共通定数を $E=210000$ MPa、$\nu=0.3$、$\sigma_{y0}=250$ MPa、$H=1000$ MPaとする。

1. 一軸ひずみ $\operatorname{diag}(0.003,0,0)$ をゼロ状態から与え、手計算と3D更新を照合する。$p=525$ MPa、$\Delta\alpha\approx0.0009642744$、$q\approx250.964274$ MPaが目印。
2. 全ひずみ $\operatorname{diag}(a,-a/2,-a/2)$ の経路を $a=0\to0.004\to-0.004\to0$ と進める。各区間を40増分に分割し、$\sigma_{xx}$ 対 $a$、$q$・$\alpha$ 対増分番号を描く。これは体積一定のひずみ指定であり、一軸応力試験ではない。全ひずみゼロの終点に応力ゼロを要求しない。
3. 非比例経路として $\operatorname{diag}(0.004,-0.002,-0.002)$ まで進み、その対角成分を保持して工学せん断 $\gamma_{xy}=0.008$ まで増す。二つの直線区間を必ず保持して各区間10、20、40、80分割し、640分割の結果に対する最終応力・$\alpha$ の差を表にする。640分割は厳密解でなく比較用の細分解。降伏面への復帰が正確でも経路の時間離散化誤差は残る。

これは材料点だけの実験であり、要素寸法を変えるメッシュ収束とは異なる。比例単調負荷では分割による差がほぼ出ないこともあるため、方向が変わる経路も検査する。

## 8. CAEの表示値へ戻る

| 観察 | このモデルでの説明・次に確認すること |
| --- | --- |
| 相当塑性ひずみが正、現在のMises応力は初期降伏応力以下 | 過去に塑性化して弾性除荷した可能性。現在の $\Delta\alpha$ と時刻歴を確認 |
| Mises応力が250 MPaを超えた | 硬化後は $250+H\alpha$ と比較。初期値との比較だけでは誤りといえない |
| 塑性ひずみ成分が負 | 圧縮方向ならあり得る。累積量 $\alpha$ の負値とは区別 |
| 平滑化コンターの $q$ と $\alpha$ が降伏則に合わない | 同じ積分点・時刻・層で比較。平均応力テンソルからの $q$ と、各点の $q$ の平均も別物 |
| 細かいメッシュほど塑性ひずみが集中 | 時間増分・要素形式・拘束・局所化・損傷則を分けて検証。材料点更新だけで原因を断定しない |

LS-DYNA等では、材料カード、硬化の種類、速度・温度依存、シェル厚さ方向の積分点、出力される塑性量の定義を確認する。本教材の $\alpha$ を、すべての材料モデルの「effective plastic strain」と無条件に同一視しない。材料更新と構造全体の安全性・破断判定は別段階である。

### 仮定を外すと何が変わるか

- **平面ひずみ**：全ひずみの面外成分がゼロでも面外応力・面外塑性ひずみは一般にゼロでない。3D内部状態を保持する。
- **平面応力**：面外応力ゼロを満たす面外ひずみを同時に求める必要がある。3D更新後に応力成分をゼロへ上書きする方法は不整合。
- **有限変形**：[物質座標・応力測度](../glossary.md#物質座標と応力測度)を区別する。基準位置ベクトル $\mathbf X$ と現在位置ベクトル $\mathbf x$ から変形勾配 $\mathbf F=\partial\mathbf x/\partial\mathbf X$（二階テンソル）を定義する。有限塑性では $\mathbf F=\mathbf F^e\mathbf F^p$ の乗算分解などを用いる。Cauchy応力は現在面積基準で、基準面積に対応する応力測度とは異なる。微小ひずみの加算分解と固定座標の増分更新を大回転へそのまま拡張しない。
- **異方性・圧力依存・速度依存**：降伏面・流れ則・時間発展式が変わるため、上の閉形式を共用できるとは限らない。
- **損傷・軟化**：本教材の $H\ge0$ とは別のモデル設計と検証を要する。破壊・局所化・メッシュ依存の扱いをPhase 6以降のレビュー課題へつなぐ。

## 9. 学習後の確認と参照先

演習は6問・100点。導出、実装、検証の根拠、CAE解釈を評価する。数値だけの一致を完了条件にしない。解答後は[学習ログ](../learning-log.md)のPhase 5記入枠へ、つまずいた導出・検証結果・節間のつながりを記録する。

式の照合・追加の実装例：

- [J. Bleyerほか：von Mises塑性のFEniCSx実装](https://bleyerj.github.io/comet-fenicsx/tours/nonlinear_problems/plasticity/plasticity.html) — 後退Eulerの材料更新と全体Newton反復。参照先の累積塑性ひずみ $p$ は、本教材では $\alpha$ に相当する。
- [LS-DYNA Support：Radial Return](https://www.dynasupport.com/tutorial/computational-plasticity/radial-return) — 偏差応力空間での復帰の説明。個別材料カードの仕様確認とは分けて読む。

上記は理論・アルゴリズムの照合先であり、演習の実行にFEniCSやLS-DYNAは不要。
