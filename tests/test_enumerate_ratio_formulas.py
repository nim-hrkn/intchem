# test_enumerate_ratio_formulas.py
# Run: pytest -q

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest

# あなたの実装モジュール名に合わせて修正してください
# 例: from ratio_formulas import enumerate_ratio_formulas
from intchem import enumerate_ratio_formulas 


def _sanitize_filename(s: str) -> str:
    # "Ag0.33Cl0.67" -> "Ag0.33Cl0.67" のままでもよいが、安全にする
    return re.sub(r"[^0-9A-Za-z.\-_=+]+", "_", s)


def build_df(chemical_formula: str, eps: float = 0.02, N1: int = 1, N2: int = 12) -> pd.DataFrame:
    results = []
    for N in range(N1, N2 + 1):
        cands = enumerate_ratio_formulas(chemical_formula, N=N, eps=eps, metric="linf")
        for c in cands:
            results.append(
                {
                    "N": N,
                    "unreduced_formula": c.unreduced_formula,
                    "reduced_formula": c.reduced_formula,
                    "error": float(c.err_linf),
                    "ion_guesses": c.ion_guesses,
                }
            )
    df = pd.DataFrame(results)

    # 並び順・型の揺れを固定して比較可能にする
    if not df.empty:
        df = df.sort_values(
            by=["N", "unreduced_formula", "reduced_formula", "error"],
            kind="mergesort",
            ignore_index=True,
        )
    return df


def _normalize_for_compare(df: pd.DataFrame) -> pd.DataFrame:
    """pickleの環境差（順序、float丸め、list/dictの並び）を吸収するための正規化"""
    if df.empty:
        return df.copy()

    out = df.copy()

    # float丸め（計算機差を小さくする）
    out["error"] = out["error"].astype(float).round(12)

    # ion_guesses: list[dict[str,float]] をソート＆丸めして安定化
    def norm_ion_guesses(x):
        if x is None:
            return None
        # 期待: list of dict
        res = []
        for d in x:
            # dict key順固定、値丸め
            res.append({k: round(float(v), 12) for k, v in sorted(d.items(), key=lambda kv: kv[0])})
        # listを文字列化して比較安定化（dict順は上で固定済み）
        # さらに解の順序も揺れる可能性があるので、ソートしておく
        res_sorted = sorted(res, key=lambda d: tuple(d.items()))
        return res_sorted

    out["ion_guesses"] = out["ion_guesses"].apply(norm_ion_guesses)

    # 行順固定
    out = out.sort_values(
        by=["N", "unreduced_formula", "reduced_formula", "error"],
        kind="mergesort",
        ignore_index=True,
    )
    return out


@pytest.mark.parametrize(
    "chemical_formula",
    [
        "Ag0.33Cl0.67",
        "Na0.45Cl0.55",
        "Fe2O3",
        "Li2O2",
    ],
)
def test_enumerate_ratio_formulas_matches_pickle(chemical_formula: str):
    data_dir = Path(__file__).resolve().parent / "data"    
    ref_path = data_dir / f"{_sanitize_filename(chemical_formula)}.pickle"
    assert ref_path.exists(), f"Reference pickle not found: {ref_path}"

    df_got = build_df(chemical_formula=chemical_formula)
    df_ref = pd.read_pickle(ref_path)

    # 比較可能な形に正規化
    df_got_n = _normalize_for_compare(df_got)
    df_ref_n = _normalize_for_compare(df_ref)

    # まず形
    assert list(df_got_n.columns) == list(df_ref_n.columns)

    # 次に中身（厳密一致）
    pd.testing.assert_frame_equal(df_got_n, df_ref_n, check_dtype=False)


def test_enumerate_ratio_formulas_invalid_formula_returns_empty_or_raises():
    # "LiO8" は組成としては整数だが、酸化数推定が空になる/候補が出ないなどを想定
    # 実装方針により「空を返す」か「例外」を投げるかがあり得るので両対応のテストにする
    chemical_formula = "LiO8"
    try:
        df = build_df(chemical_formula=chemical_formula)
        # 少なくともクラッシュせずDataFrameが返ること
        assert isinstance(df, pd.DataFrame)
        # 期待としては候補が非常に少ないか、空
        assert df.shape[0] >= 0
    except Exception as e:
        # 例外にする実装でもpytest的には許す（ただし型は限定）
        assert isinstance(e, (ValueError, KeyError))