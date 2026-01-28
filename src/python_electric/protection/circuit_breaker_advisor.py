from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Sequence
import warnings

from .. import Quantity, Q_
from .circuit_breaker import CircuitBreaker

if TYPE_CHECKING:
    from ..network.components.cable import Cable


__all__ = [
    "CircuitBreakerSuggestion",
    "CircuitBreakerAdvisor",
]


@dataclass(frozen=True, slots=True)
class CircuitBreakerSuggestion:
    """
    Result of an advisor evaluation for a single circuit-breaker candidate.
    """
    standard: CircuitBreaker.Standard
    category: CircuitBreaker.Category
    I_cu: Quantity
    k_m: float | None
    passed: bool
    reasons: list[str]
    warnings: list[str]
    sort_key: tuple

    def to_kwargs(self) -> dict:
        """
        Convert the suggestion to keyword arguments accepted by
        Cable.connect_circuit_breaker(...).

        Notes
        -----
        - `I_sc_max` and `I_sc_min` are intentionally not included here, since
          Cable.connect_circuit_breaker(...) accepts them separately, and you
          typically pass the cable's stored values.
        """
        return {
            "standard": self.standard,
            "category": self.category,
            "I_cu": self.I_cu,
            "k_m": self.k_m,
        }

    def __str__(self) -> str:
        s = [
            f"standard: {self.standard}",
            f"category: {self.category}",
            f"I_cu: {self.I_cu.to('kA'):~P.0f}",
            f"k_m: {self.k_m}"
        ]
        if self.warnings:
            s.append("\nWarnings:")
            s.append("---------")
            s.extend([str(w) for w in self.warnings])
        return "\n".join(s)


class CircuitBreakerAdvisor:
    """
    Suggest suitable values for circuit-breaker parameters (category, I_cu, k_m)
    for a given Cable, based on currents already computed/stored in the Cable.

    Philosophy
    ----------
    This advisor does NOT "auto-install" a circuit breaker. It enumerates
    candidates and validates them by re-using the same CircuitBreaker checks
    that Cable.connect_circuit_breaker(...) uses, returning transparent
    recommendations (with reasons and warnings).

    Required Cable data (minimal)
    -----------------------------
    - cable.I_sc_max and cable.I_sc_min must be available (short-circuit calc done)
    - cable.I_b, cable.I_n, cable.I_z and cable.I2t should be available (sizing done)
    """

    def __init__(
        self,
        cable: Cable,
        standard: CircuitBreaker.Standard,
        *,
        safety_factor_Icu: float = 1.0,
        icu_series: Sequence[Quantity | float] | None = None,
        km_grid: Sequence[float] | None = None,
        prefer_adjustable: bool = False
    ) -> None:
        self.cable = cable
        self.standard = standard
        self.safety_factor_Icu = float(safety_factor_Icu)

        # Icu series can be provided as floats (interpreted as kA) or Quantities.
        self.icu_series = self._normalize_icu_series(
            icu_series
            if icu_series is not None
            else [6, 10, 15, 20, 25, 36, 50, 70]  # default: kA values
        )

        self.km_grid = list(km_grid) if km_grid is not None else [4, 5, 6, 8, 10, 12, 14]
        self.prefer_adjustable = bool(prefer_adjustable)

        self._validate_inputs()

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def suggest(self, top: int = 1) -> list[CircuitBreakerSuggestion]:
        """
        Return the best suggestion(s), ranked from best to worst.

        Parameters
        ----------
        top:
            Number of best suggestions to return. Use top=1 for a single best
            suggestion. Use top>1 to inspect alternatives.

        Returns
        -------
        list[CircuitBreakerSuggestion]
            Sorted list of best candidates. If no candidate passes, the returned
            list is empty (and you should inspect `suggest_all()` for diagnostics).
        """
        top = int(top)
        if top < 1:
            raise ValueError("top must be >= 1.")

        all_suggestions = self.suggest_all()
        passed = [s for s in all_suggestions if s.passed]
        # passed.sort(key=lambda s: s.sort_key)
        return passed[:top]

    def suggest_all(self) -> list[CircuitBreakerSuggestion]:
        """
        Generate and evaluate ALL candidates, return them sorted.

        Useful for diagnostics (why did everything fail?).
        """
        candidates = self._generate_candidates()
        results: list[CircuitBreakerSuggestion] = []
        for category, I_cu, k_m in candidates:
            results.append(self._evaluate_candidate(category, I_cu, k_m))
        results.sort(key=lambda s: s.sort_key)
        return results

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------

    def _validate_inputs(self) -> None:
        # Ensure short-circuit data exists
        if (getattr(self.cable, "I_sc_max", None) is None
                or getattr(self.cable, "I_sc_min", None) is None):
            raise ValueError(
                "Cable short-circuit currents are not available "
                "(I_sc_max / I_sc_min). Run the short-circuit "
                "calculation first and propagate results to the "
                "cable."
            )

        # Ensure sizing/current data exists (these should normally exist in Cable)
        for attr in ("I_b_ph", "I_n_ph", "I_z_ph"):
            if getattr(self.cable, attr, None) is None:
                raise ValueError(
                    f"Cable attribute '{attr}' is required but not set."
                )

        # I2t may be None for some cases; we allow it because CircuitBreaker
        # accepts I2t=None. But the check might be less strict in that case.
        _ = getattr(self.cable, "I2t_ph", None)

        # Ensure series not empty
        if not self.icu_series:
            raise ValueError("icu_series is empty.")
        if not self.km_grid:
            raise ValueError("km_grid is empty.")

    # -------------------------------------------------------------------------
    # Candidate generation
    # -------------------------------------------------------------------------

    def _generate_candidates(
        self
    ) -> list[tuple[CircuitBreaker.Category, Quantity, float | None]]:
        """
        Returns a list of (category, I_cu, k_m) candidates.
        """
        I_sc_max = self.cable.I_sc_max.to("A")
        I_sc_min = self.cable.I_sc_min.to("A")

        # Determine minimal Icu requirement and choose a small set
        # of Icu candidates
        I_required = (self.safety_factor_Icu * I_sc_max.m) * Q_(1.0, "A")
        icu_candidates = self._select_icu_candidates(I_required)

        # Category ordering (can be refined later with "motor feeder" hints, etc.)
        if self.standard == CircuitBreaker.Standard.RESIDENTIAL:
            cat_order = [
                CircuitBreaker.Category.C,
                CircuitBreaker.Category.B,
                CircuitBreaker.Category.D,
            ]
        else:
            # Industrial: include adjustable as fallback (unless user prefers it)
            if self.prefer_adjustable:
                cat_order = [
                    CircuitBreaker.Category.ADJUSTABLE,
                    CircuitBreaker.Category.C,
                    CircuitBreaker.Category.B,
                    CircuitBreaker.Category.D,
                ]
            else:
                cat_order = [
                    CircuitBreaker.Category.C,
                    CircuitBreaker.Category.B,
                    CircuitBreaker.Category.D,
                    CircuitBreaker.Category.ADJUSTABLE,
                ]

        candidates: list[tuple[CircuitBreaker.Category, Quantity, float | None]] = []

        for I_cu in icu_candidates:
            for cat in cat_order:
                if cat != CircuitBreaker.Category.ADJUSTABLE:
                    candidates.append((cat, I_cu, None))
                    continue

                # Adjustable only meaningful for INDUSTRIAL
                if self.standard != CircuitBreaker.Standard.INDUSTRIAL:
                    continue

                # Limit k_m grid using the fundamental feasibility constraint:
                # For industrial adjustable, CircuitBreaker uses I_r = I_b and
                # I_m_min = 0.8 * k_m * I_r; we need I_sc_min >= I_m_min.
                I_b = self.cable.I_b_ph.to("A")
                km_max_allowed = (I_sc_min.m / (0.8 * I_b.m)) if I_b.m > 0 else 0.0

                for km in self.km_grid:
                    if km <= km_max_allowed + 1e-12:
                        candidates.append((cat, I_cu, float(km)))

        return candidates

    def _select_icu_candidates(self, I_required: Quantity) -> list[Quantity]:
        """
        Pick a small set of Icu candidates from the series:
        - the smallest value >= I_required
        - optionally 1-2 higher values (gives alternatives)
        """
        series = sorted([v.to("A") for v in self.icu_series], key=lambda x: x.m)

        # find first >= required
        idx = None
        for i, v in enumerate(series):
            if v.m >= I_required.to("A").m - 1e-12:
                idx = i
                break
        if idx is None:
            # required is above series max -> return just the max (will fail breaking-capacity check anyway)
            return [series[-1]]

        # return the first + next two (if exist)
        out = [series[idx]]
        if idx + 1 < len(series):
            out.append(series[idx + 1])
        if idx + 2 < len(series):
            out.append(series[idx + 2])
        return out

    # -------------------------------------------------------------------------
    # Candidate evaluation
    # -------------------------------------------------------------------------

    def _evaluate_candidate(
        self,
        category: CircuitBreaker.Category,
        I_cu: Quantity,
        k_m: float | None,
    ) -> CircuitBreakerSuggestion:
        reasons: list[str] = []
        warn_msgs: list[str] = []

        I_sc_max = self.cable.I_sc_max.to("A")
        I_sc_min = self.cable.I_sc_min.to("A")

        # Quick breaking capacity pre-check
        if I_cu.to("A").m + 1e-12 < I_sc_max.m:
            reasons.append(
                f"Fail: I_cu ({I_cu.to('kA').m:.1f} kA) < I_sc_max ({I_sc_max.to('kA').m:.1f} kA)."
            )
            sort_key = self._make_sort_key(
                passed=False, I_cu=I_cu, category=category, k_m=k_m, margin_sc=0.0
            )
            return CircuitBreakerSuggestion(
                standard=self.standard,
                category=category,
                I_cu=I_cu,
                k_m=k_m,
                passed=False,
                reasons=reasons,
                warnings=warn_msgs,
                sort_key=sort_key,
            )
        reasons.append(
            f"Pass: I_cu ({I_cu.to('kA').m:.1f} kA) >= I_sc_max ({I_sc_max.to('kA').m:.1f} kA)."
        )

        # Evaluate by re-using CircuitBreaker checks (same input set as Cable.connect_circuit_breaker)
        passed = False
        margin_sc = 0.0

        with warnings.catch_warnings(record=True) as wrec:
            warnings.simplefilter("always")

            cb = CircuitBreaker(
                standard=self.standard,
                category=category,
                I_b=self.cable.I_b_ph,
                I_n=self.cable.I_n_ph,
                I_z=self.cable.I_z_ph,
                I2t=getattr(self.cable, "I2t_tot", None),
                I_cu=I_cu,
                E_t=None,
                k_m=k_m,
                t_m=None,
            )

            c_over = cb.check_overload_protection()
            c_sc = cb.check_shortcircuit_protection(I_sc_max, I_sc_min)

            # Capture warnings emitted by the CB checks
            for w in wrec:
                warn_msgs.append(str(w.message))

            if c_over:
                reasons.append("Pass: overload protection check succeeded.")
            else:
                reasons.append("Fail: overload protection check failed.")

            if c_sc:
                reasons.append("Pass: short-circuit protection check succeeded.")
            else:
                reasons.append("Fail: short-circuit protection check failed.")

            passed = bool(c_over and c_sc)

            # Add a simple margin indicator for ranking/diagnostics:
            # For non-adjustable categories, estimate margin using cb's magnetic threshold,
            # otherwise use I_sc_min - I_m_min based on k_m and I_r.
            try:
                # noinspection PyProtectedMember
                I_m_min = cb._get_min_magn_trip_threshold().to("A").m  # internal helper exists in your class
                margin_sc = I_sc_min.m - I_m_min
                reasons.append(
                    f"Info: magnetic trip margin at I_sc_min = {margin_sc:.1f} A "
                    f"(I_sc_min={I_sc_min.m:.1f} A, I_m_min={I_m_min:.1f} A)."
                )
            except Exception:
                # Keep margin at 0 if we cannot compute it for any reason
                margin_sc = 0.0

        sort_key = self._make_sort_key(
            passed=passed, I_cu=I_cu, category=category, k_m=k_m, margin_sc=margin_sc
        )

        return CircuitBreakerSuggestion(
            standard=self.standard,
            category=category,
            I_cu=I_cu,
            k_m=k_m,
            passed=passed,
            reasons=reasons,
            warnings=warn_msgs,
            sort_key=sort_key,
        )

    # -------------------------------------------------------------------------
    # Ranking / sorting
    # -------------------------------------------------------------------------

    def _make_sort_key(
        self,
        *,
        passed: bool,
        I_cu: Quantity,
        category: CircuitBreaker.Category,
        k_m: float | None,
        margin_sc: float,
    ) -> tuple:
        """
        Lower sort_key is better.

        Ranking policy (simple, explainable):
        1) Passed candidates first
        2) Lowest I_cu (avoid overkill)
        3) Prefer fixed curves over adjustable (unless prefer_adjustable=True)
        4) Category preference (general feeder heuristic)
        5) For adjustable: prefer higher k_m (more immunity) among passing candidates
        6) Prefer larger short-circuit margin (I_sc_min - I_m_min)
        """
        passed_rank = 0 if passed else 1
        icu_rank = I_cu.to("A").m

        is_adjustable = (category == CircuitBreaker.Category.ADJUSTABLE)
        if self.prefer_adjustable:
            adjustable_rank = 0 if is_adjustable else 1
        else:
            adjustable_rank = 1 if is_adjustable else 0

        # General curve preference (you can later inject "motor feeder" logic here):
        # C > B > D > ADJUSTABLE (adjustable already handled separately).
        cat_pref = {
            CircuitBreaker.Category.C: 0,
            CircuitBreaker.Category.B: 1,
            CircuitBreaker.Category.D: 2,
            CircuitBreaker.Category.ADJUSTABLE: 3,
        }.get(category, 99)

        # For adjustable: higher k_m is better (within passing set)
        km_rank = 0.0
        if is_adjustable and k_m is not None:
            km_rank = -float(k_m)  # negative so higher k_m becomes "smaller" => better

        # Higher margin_sc is better => negative for sorting
        margin_rank = -float(margin_sc)

        return passed_rank, icu_rank, adjustable_rank, cat_pref, km_rank, margin_rank

    # -------------------------------------------------------------------------
    # Helpers
    # -------------------------------------------------------------------------

    @staticmethod
    def _normalize_icu_series(values: Iterable[Quantity | float]) -> list[Quantity]:
        out: list[Quantity] = []
        for v in values:
            if isinstance(v, (int, float)):
                out.append(Q_(float(v), "kA"))
            else:
                out.append(v.to("A"))
        # Ensure strictly positive and unique-ish
        out = [x for x in out if x.to("A").m > 0.0]
        return out
