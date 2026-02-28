from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import itertools

from pymatgen.core.composition import Composition


@dataclass(frozen=True)
class IntChemCandidate:
    N: int
    reduced_formula: str
    unreduced_formula: str
    integer_reduced_composition: Composition
    integer_unreduced_composition: Composition
    frac: Dict[str, float]
    err_l1: float
    err_linf: float
    ion_guesses: List[Dict[str, float]]  


def _target_fraction(comp: Composition) -> Dict[str, float]:
    fc = comp.fractional_composition
    return {str(el): float(v) for el, v in fc.items()}


def _errors(target: Dict[str, float], got: Dict[str, float]) -> Tuple[float, float]:
    keys = sorted(set(target) | set(got))
    diffs = [abs(got.get(k, 0.0) - target.get(k, 0.0)) for k in keys]
    return (sum(diffs), max(diffs) if diffs else 0.0)


def unreduced_formula_str(comp: Composition, order: List[str]) -> str:
    d = comp.get_el_amt_dict()
    s = ""
    for el in order:
        if el not in d:
            continue
        n = int(round(d[el]))
        if n == 0:
            continue
        s += el if n == 1 else f"{el}{n}"
    return s


def _all_integer_tuples_sum_N(m: int, N: int):
    if m == 1:
        yield (N,)
        return
    for x in range(N + 1):
        for rest in _all_integer_tuples_sum_N(m - 1, N - x):
            yield (x,) + rest


def enumerate_ratio_formulas(
    fraction_formula: str,
    N: int,
    eps: float,
    metric: str = "linf",          # "linf" or "l1"
    max_candidate: int = 10,
    compute_ion_guesses: bool = True,      
    all_oxi_states: bool = False,          
    max_sites: Optional[int] = None,       
    oxi_states_override: Optional[dict] = None,  
) -> List[IntChemCandidate]:
    """
    Input:
      fraction_formula: e.g. "Ag0.33Cl0.67"
      N: total atom count (single value)
      eps: tolerance parameter (accept if err <= eps)
    Output:
      List[intChemCandidate] that satisfy acceptance criterion.

    Acceptance:
      err = max_i |x_int_i - x_target_i| (metric="linf")
      err = sum_i |x_int_i - x_target_i|  (metric="l1")
      accept if err <= eps
    """
    if N <= 0:
        return []

    target_comp = Composition(fraction_formula)
    target_frac = _target_fraction(target_comp)

    elems = list(target_frac.keys())
    m = len(elems)

    raw = [N * target_frac[e] for e in elems]
    floors = [int(x // 1) for x in raw]
    base_sum = sum(floors)
    need = N - base_sum

    candidates: List[IntChemCandidate] = []

    def add_candidate(counts: List[int]):
        if sum(counts) != N:
            return
        if all(c == 0 for c in counts):
            return

        int_comp = Composition({elems[i]: counts[i] for i in range(m) if counts[i] != 0})
        frac = _target_fraction(int_comp)
        err_l1, err_linf = _errors(target_frac, frac)
        if False:
            err_l1 = err_l1/N
            err_linf = err_linf/N
        err = err_linf if metric.lower() == "linf" else err_l1

        if err <= eps:
            red_comp = Composition(int_comp.reduced_formula)

            if compute_ion_guesses:
                ion_guesses = red_comp.oxi_state_guesses(
                    oxi_states_override=oxi_states_override,
                    all_oxi_states=all_oxi_states,
                    max_sites=max_sites,
                )
            else:
                ion_guesses = []

            candidates.append(
                IntChemCandidate(
                    N=N,
                    reduced_formula=int_comp.reduced_formula,
                    unreduced_formula=unreduced_formula_str(int_comp, elems),
                    integer_reduced_composition=red_comp,
                    integer_unreduced_composition=int_comp,
                    frac=frac,
                    err_l1=err_l1,
                    err_linf=err_linf,
                    ion_guesses=ion_guesses,
                )
            )

    if 0 <= need <= m:
        for idxs in itertools.combinations(range(m), need):
            counts = floors.copy()
            for i in idxs:
                counts[i] += 1
            add_candidate(counts)
    else:
        if m >= max_candidate:
            return []
        for counts in _all_integer_tuples_sum_N(m, N):
            add_candidate(list(counts))

    uniq = {}
    for c in candidates:
        key = c.unreduced_formula
        if key not in uniq or (c.err_linf, c.err_l1) < (uniq[key].err_linf, uniq[key].err_l1):
            uniq[key] = c

    out = sorted(uniq.values(), key=lambda x: (x.err_linf, x.err_l1, x.reduced_formula, x.unreduced_formula))
    return out

