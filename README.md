このプログラムは、**分数組成で書かれた化学式**から、**総原子数 (N)** を固定して、それに近い**整数比の化学式候補**を列挙するためのものです。
たとえば `Ag0.33Cl0.67` に対して、`N=3` なら `Ag1Cl2` のような整数組成候補を返します。

---

## 何をするプログラムか

中心の関数は `enumerate_ratio_formulas(...)` です。

これは、

* 入力化学式: `Ag0.33Cl0.67`
* 総原子数: `N=3`
* 許容誤差: `eps`

を与えると、

* `AgCl2`
* `Ag2Cl4`
* などのような整数組成候補

のうち、**元の分率に十分近いもの**を返します。

さらに各候補について、

* 還元化学式 (`reduced_formula`)
* 非還元化学式 (`unreduced_formula`)
* 組成分率 (`frac`)
* 誤差 (`err_l1`, `err_linf`)
* 酸化数推定 (`ion_guesses`)

も持たせています。

---

## 各部の役割

### 1. `IntChemCandidate`

候補1件を表すデータクラスです。

主な属性は以下です。

* `N`
  総原子数
* `reduced_formula`
  約分後の化学式
* `unreduced_formula`
  約分前の整数化学式
* `integer_reduced_composition`
  約分後の `Composition`
* `integer_unreduced_composition`
  約分前の `Composition`
* `frac`
  候補の元素分率
* `err_l1`
  L1誤差
* `err_linf`
  L∞誤差
* `ion_guesses`
  `pymatgen` の酸化数推定結果

---

### 2. `_target_fraction(comp)`

`Composition` を元素分率辞書に変換します。

例:

```python
Composition("Ag1Cl2")
```

から

```python
{"Ag": 1/3, "Cl": 2/3}
```

のような辞書を作ります。

---

### 3. `_errors(target, got)`

目標分率と候補分率の差を計算します。

返り値は

* `err_l1 = Σ |差|`
* `err_linf = max |差|`

です。

---

### 4. `unreduced_formula_str(comp, order)`

元素順を保ったまま、約分しない化学式文字列を作ります。

例:

* `Composition({"Ag":2, "Cl":4})`
* 順序 `["Ag", "Cl"]`

なら

```python
"Ag2Cl4"
```

を返します。

`pymatgen` の `reduced_formula` だと約分されて `AgCl2` になるので、それとは別に非還元表記を残すための関数です。

---

### 5. `_all_integer_tuples_sum_N(m, N)`

長さ `m` の非負整数タプルで、和が `N` になるものを全列挙します。

例:

```python
m=2, N=3
```

なら

```python
(0,3), (1,2), (2,1), (3,0)
```

を生成します。

ただし通常は全探索はせず、まずはもっと軽い方法で候補を作っています。

---

## 中心関数の使い方

## `enumerate_ratio_formulas(...)`

### 引数

```python
enumerate_ratio_formulas(
    fraction_formula: str,
    N: int,
    eps: float,
    metric: str = "linf",
    max_candidate: int = 10,
    compute_ion_guesses: bool = True,
    all_oxi_states: bool = False,
    max_sites: Optional[int] = None,
    oxi_states_override: Optional[dict] = None,
) -> List[IntChemCandidate]
```

### 主要引数の意味

#### `fraction_formula`

分数組成の化学式文字列です。

例:

```python
"Ag0.33Cl0.67"
"Li0.2Mn0.4O0.4"
```

#### `N`

整数化したときの総原子数です。

例:

* `Ag0.33Cl0.67` に対して `N=3` なら `Ag1Cl2`
* `N=6` なら `Ag2Cl4`

#### `eps`

許容誤差です。
誤差がこれ以下の候補だけ返します。

#### `metric`

誤差判定に使う指標です。

* `"linf"`: 最大差
* `"l1"`: 差の総和

#### `max_candidate`

丸めで処理できない場合に全探索へ落ちますが、そのとき
`m >= max_candidate` なら計算量を避けるため空配列を返します。
ここでの `m` は元素数です。

#### `compute_ion_guesses`

`True` なら酸化数推定も行います。

#### `all_oxi_states`, `max_sites`, `oxi_states_override`

`pymatgen.Composition.oxi_state_guesses()` に渡す設定です。

---

## 候補生成の仕組み

この関数はまず、目標分率に `N` を掛けます。

たとえば

```python
fraction_formula = "Ag0.33Cl0.67"
N = 3
```

なら理想個数はだいたい

* Ag: 0.99
* Cl: 2.01

です。

これを床関数で丸めると

* Ag: 0
* Cl: 2

合計 2 個なので、あと 1 個必要です。
その 1 個をどの元素に足すかを全組合せで試します。

この例では

* Ag に 1 個足す → `(1,2)`
* Cl に 1 個足す → `(0,3)`

のような候補を作って誤差判定します。

---

## 使用例

### 例1: 基本

```python
cands = enumerate_ratio_formulas(
    fraction_formula="Ag0.33Cl0.67",
    N=3,
    eps=0.05,
)
```

### 表示例

```python
for c in cands:
    print("unreduced:", c.unreduced_formula)
    print("reduced:", c.reduced_formula)
    print("frac:", c.frac)
    print("err_linf:", c.err_linf)
    print("err_l1:", c.err_l1)
    print("ion_guesses:", c.ion_guesses)
    print()
```

期待される候補の一例は

```python
unreduced: AgCl2
reduced: AgCl2
frac: {'Ag': 0.3333333333333333, 'Cl': 0.6666666666666666}
err_linf: 0.0033333333333332993
err_l1: 0.006666666666666599
ion_guesses: ...
```

---

### 例2: 酸化数推定を切る

```python
cands = enumerate_ratio_formulas(
    fraction_formula="Li0.2Mn0.4O0.4",
    N=5,
    eps=0.1,
    compute_ion_guesses=False,
)
```

酸化数推定はやや重い場合があるので、不要なら `False` が便利です。

---

### 例3: L1誤差で判定

```python
cands = enumerate_ratio_formulas(
    fraction_formula="Ag0.33Cl0.67",
    N=3,
    eps=0.02,
    metric="l1",
)
```

---

## 返り値

返り値は `List[IntChemCandidate]` です。
誤差の小さい順に並んで返ります。

ソート順は

1. `err_linf`
2. `err_l1`
3. `reduced_formula`
4. `unreduced_formula`

です。

---

## 重複除去

候補は最後に

```python
key = c.unreduced_formula
```

で重複除去されています。

つまり、同じ非還元化学式なら、より誤差の小さいものだけ残します。

---

## 注意点

### 1. `fraction_formula` は `pymatgen.Composition` が読める形式である必要がある

たとえば

```python
"Ag0.33Cl0.67"
```

のような形式は読めます。

---

### 2. `N <= 0` なら空配列

```python
if N <= 0:
    return []
```

です。

---

### 3. 元素数が多いと全探索を避ける

丸めベースの簡便法で処理できないとき、全探索に入ります。
ただし元素数 `m` が `max_candidate` 以上だと

```python
return []
```

になります。

ここは変数名 `max_candidate` より、実際には
**「全探索を許す最大元素数」**
の意味に近いです。

---

### 4. 誤差は分率空間で計算している

個数差ではなく、**元素分率の差**で判定しています。

---

### 5. `if False:` の部分は今は無効

ここ:

```python
if False:
    err_l1 = err_l1/N
    err_linf = err_linf/N
```

は実行されません。
つまり現在は、誤差はそのままの分率差で評価しています。

---

## 最小実行例

```python
from pymatgen.core.composition import Composition

cands = enumerate_ratio_formulas(
    fraction_formula="Ag0.33Cl0.67",
    N=3,
    eps=0.05,
    metric="linf",
)

for c in cands:
    print(c)
```

---

## どういう場面で使うか

この関数は、たとえば

* 実験組成が小数で与えられている
* しかし結晶構造生成や列挙には整数組成が必要
* 総サイト数 `N` を仮定して候補式を出したい

という場面で便利です。

特に、

* CSP の入力式候補生成
* 分率表記から近い整数比組成の探索
* 酸化数推定込みの事前フィルタ

に向いています。

---

## 使い方まとめ

最も基本的にはこれです。

```python
cands = enumerate_ratio_formulas(
    fraction_formula="Ag0.33Cl0.67",
    N=3,
    eps=0.05,
)
```

そして各候補を確認します。

```python
for c in cands:
    print(c.unreduced_formula, c.reduced_formula, c.err_linf)
```
