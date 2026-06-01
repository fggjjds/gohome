"""Reported experiment tables and comparison helpers for the reproduced papers.

The numbers in this module are transcribed from the paper screenshots supplied in
this task.  They are used to check whether the implementation follows the same
algorithmic trend as the papers: EvalComp improves bootstrapping time/precision,
and NewCompH3&5 reduces comparison complexity/runtime versus NewCompG/H.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ComparisonComplexityRow:
    """NewComp comparison-complexity row from the optimized sign paper."""

    alpha: int
    newcomp_g: int
    newcomp_h: int
    newcomp_h35: int

    @property
    def h35_vs_g_reduction(self) -> float:
        return 1 - self.newcomp_h35 / self.newcomp_g

    @property
    def h35_vs_h_reduction(self) -> float:
        return 1 - self.newcomp_h35 / self.newcomp_h


@dataclass(frozen=True)
class ComparisonRuntimeRow:
    """HEAAN runtime row for comparing NewCompG/H/H3&5."""

    alpha: int
    newcomp_g_ms: float
    newcomp_h_ms: float
    newcomp_h35_ms: float

    @property
    def h35_vs_g_speedup(self) -> float:
        return self.newcomp_g_ms / self.newcomp_h35_ms

    @property
    def h35_vs_h_speedup(self) -> float:
        return self.newcomp_h_ms / self.newcomp_h35_ms


@dataclass(frozen=True)
class EvalCompAccuracyRow:
    """EvalComp paper Table IV row for time and bootstrapping precision."""

    scheme: str
    parameter: str
    remaining_modulus_bits: int
    total_time_ms: float
    amortized_time_ms: float
    precision_bits: float


@dataclass(frozen=True)
class NewEvalCompH35Estimate:
    """Estimated combined EvalComp + NewCompH3&5 row.

    The estimate keeps EvalComp's bootstrapping precision/modulus-retention
    metrics and applies the NewCompH3&5 HCF-component runtime factor to a
    configurable share of EvalComp's total time.
    """

    family: str
    alpha: int
    hcf_share: float
    baseline_comparison: str
    evalcomp_time_ms: float
    estimated_time_ms: float
    evalcomp_precision_bits: float
    remaining_modulus_bits: int
    hcf_component_speedup: float

    @property
    def time_reduction_vs_evalcomp(self) -> float:
        return 1 - self.estimated_time_ms / self.evalcomp_time_ms


NEWCOMP_COMPLEXITY_ROWS = (
    ComparisonComplexityRow(4, 16, 14, 13),
    ComparisonComplexityRow(8, 24, 20, 19),
    ComparisonComplexityRow(12, 28, 26, 25),
    ComparisonComplexityRow(16, 36, 32, 30),
    ComparisonComplexityRow(20, 40, 38, 36),
)

NEWCOMP_RUNTIME_ROWS = (
    ComparisonRuntimeRow(4, 20.8, 17.9, 11.7),
    ComparisonRuntimeRow(8, 39.8, 26.9, 17.0),
    ComparisonRuntimeRow(12, 48.1, 38.1, 25.7),
    ComparisonRuntimeRow(16, 68.8, 54.7, 33.8),
    ComparisonRuntimeRow(20, 83.7, 66.0, 45.4),
)

EVALCOMP_ACCURACY_ROWS = (
    EvalCompAccuracyRow("EvalMod", "I-1", 510, 1646, 3.214, 14.4),
    EvalCompAccuracyRow("EvalMod", "II-1", 510, 13902, 3.394, 16.7),
    EvalCompAccuracyRow("EvalMod", "III-1", 515, 115081, 3.512, 20.1),
    EvalCompAccuracyRow("EvalRound", "I-2", 513, 1799, 3.514, 23.5),
    EvalCompAccuracyRow("EvalRound", "II-2", 582, 14843, 3.624, 25.5),
    EvalCompAccuracyRow("EvalRound", "III-2", 590, 121208, 3.699, 29.1),
    EvalCompAccuracyRow("EvalComp", "I-3", 633, 1369, 2.673, 25.7),
    EvalCompAccuracyRow("EvalComp", "II-3", 632, 11698, 2.856, 28.1),
    EvalCompAccuracyRow("EvalComp", "III-3", 615, 98435, 3.004, 32.2),
)


def evalcomp_time_precision_summary() -> list[dict[str, float | str]]:
    """Compare EvalComp rows against matching EvalMod/EvalRound rows."""

    rows: list[dict[str, float | str]] = []
    for family in ("I", "II", "III"):
        evalmod = _find_eval_row("EvalMod", f"{family}-1")
        evalround = _find_eval_row("EvalRound", f"{family}-2")
        evalcomp = _find_eval_row("EvalComp", f"{family}-3")
        rows.append(
            {
                "family": family,
                "evalcomp_vs_evalmod_time_reduction": 1
                - evalcomp.total_time_ms / evalmod.total_time_ms,
                "evalcomp_vs_evalround_time_reduction": 1
                - evalcomp.total_time_ms / evalround.total_time_ms,
                "evalcomp_vs_evalmod_precision_gain_bits": evalcomp.precision_bits
                - evalmod.precision_bits,
                "evalcomp_vs_evalround_precision_gain_bits": evalcomp.precision_bits
                - evalround.precision_bits,
            }
        )
    return rows


def newcomp_h35_summary() -> list[dict[str, float | int]]:
    """Return complexity/runtime improvements for NewCompH3&5."""

    rows: list[dict[str, float | int]] = []
    runtime_by_alpha = {row.alpha: row for row in NEWCOMP_RUNTIME_ROWS}
    for complexity in NEWCOMP_COMPLEXITY_ROWS:
        runtime = runtime_by_alpha[complexity.alpha]
        rows.append(
            {
                "alpha": complexity.alpha,
                "complexity_vs_g_reduction": complexity.h35_vs_g_reduction,
                "complexity_vs_h_reduction": complexity.h35_vs_h_reduction,
                "runtime_vs_g_speedup": runtime.h35_vs_g_speedup,
                "runtime_vs_h_speedup": runtime.h35_vs_h_speedup,
            }
        )
    return rows


def newevalcomp_h35_estimates(
    alpha: int = 8, hcf_share: float = 1.0, baseline_comparison: str = "H"
) -> list[NewEvalCompH35Estimate]:
    """Estimate the combined NewEvalCompH3&5 scheme.

    ``baseline_comparison`` is the comparison component assumed inside the
    original EvalComp HCF path: ``"H"`` uses the NewCompH row and ``"G"`` uses
    the NewCompG row.  ``hcf_share`` is the fraction of EvalComp total time
    attributed to HCF comparisons; use 1.0 for a component-level upper-bound
    estimate, or a smaller value to model non-HCF transforms.
    """

    if not 0 <= hcf_share <= 1:
        raise ValueError("hcf_share must lie in [0, 1]")
    runtime = _find_runtime_row(alpha)
    normalized = baseline_comparison.strip().upper()
    if normalized == "H":
        hcf_factor = runtime.newcomp_h35_ms / runtime.newcomp_h_ms
    elif normalized == "G":
        hcf_factor = runtime.newcomp_h35_ms / runtime.newcomp_g_ms
    else:
        raise ValueError("baseline_comparison must be 'H' or 'G'")
    hcf_speedup = 1 / hcf_factor

    estimates: list[NewEvalCompH35Estimate] = []
    for family in ("I", "II", "III"):
        evalcomp = _find_eval_row("EvalComp", f"{family}-3")
        estimated_time = evalcomp.total_time_ms * (
            1 - hcf_share * (1 - hcf_factor)
        )
        estimates.append(
            NewEvalCompH35Estimate(
                family=family,
                alpha=alpha,
                hcf_share=hcf_share,
                baseline_comparison=normalized,
                evalcomp_time_ms=evalcomp.total_time_ms,
                estimated_time_ms=estimated_time,
                evalcomp_precision_bits=evalcomp.precision_bits,
                remaining_modulus_bits=evalcomp.remaining_modulus_bits,
                hcf_component_speedup=hcf_speedup,
            )
        )
    return estimates


def format_paper_report(
    combined_alpha: int = 8, hcf_share: float = 1.0, baseline_comparison: str = "H"
) -> str:
    """Return a human-readable report comparing against the paper tables."""

    lines = ["EvalComp paper comparison (reported Table IV):"]
    for row in evalcomp_time_precision_summary():
        lines.append(
            "  family {family}: time -{mod:.1%} vs EvalMod, -{rnd:.1%} vs "
            "EvalRound; precision +{pmod:.1f}b vs EvalMod, +{prnd:.1f}b "
            "vs EvalRound".format(
                family=row["family"],
                mod=row["evalcomp_vs_evalmod_time_reduction"],
                rnd=row["evalcomp_vs_evalround_time_reduction"],
                pmod=row["evalcomp_vs_evalmod_precision_gain_bits"],
                prnd=row["evalcomp_vs_evalround_precision_gain_bits"],
            )
        )

    lines.append("NewCompH3&5 paper comparison (reported Tables II/III):")
    for row in newcomp_h35_summary():
        lines.append(
            "  alpha={alpha}: complexity -{cg:.1%} vs G, -{ch:.1%} vs H; "
            "runtime {sg:.2f}x vs G, {sh:.2f}x vs H".format(
                alpha=row["alpha"],
                cg=row["complexity_vs_g_reduction"],
                ch=row["complexity_vs_h_reduction"],
                sg=row["runtime_vs_g_speedup"],
                sh=row["runtime_vs_h_speedup"],
            )
        )

    lines.append(
        "NewEvalCompH3&5 combined estimate "
        f"(alpha={combined_alpha}, HCF share={hcf_share:.0%}, "
        f"baseline={baseline_comparison.upper()}):"
    )
    for row in newevalcomp_h35_estimates(
        alpha=combined_alpha,
        hcf_share=hcf_share,
        baseline_comparison=baseline_comparison,
    ):
        lines.append(
            "  family {family}: est. time {time:.1f}ms "
            "(-{reduction:.1%} vs EvalComp), precision {precision:.1f}b, "
            "remaining modulus 2^{modulus}".format(
                family=row.family,
                time=row.estimated_time_ms,
                reduction=row.time_reduction_vs_evalcomp,
                precision=row.evalcomp_precision_bits,
                modulus=row.remaining_modulus_bits,
            )
        )
    return "\n".join(lines)


def _find_eval_row(scheme: str, parameter: str) -> EvalCompAccuracyRow:
    for row in EVALCOMP_ACCURACY_ROWS:
        if row.scheme == scheme and row.parameter == parameter:
            return row
    raise ValueError(f"missing reported row: {scheme} {parameter}")


def _find_runtime_row(alpha: int) -> ComparisonRuntimeRow:
    for row in NEWCOMP_RUNTIME_ROWS:
        if row.alpha == alpha:
            return row
    raise ValueError(f"missing NewComp runtime row for alpha={alpha}")


__all__ = [
    "EVALCOMP_ACCURACY_ROWS",
    "NEWCOMP_COMPLEXITY_ROWS",
    "NEWCOMP_RUNTIME_ROWS",
    "ComparisonComplexityRow",
    "ComparisonRuntimeRow",
    "EvalCompAccuracyRow",
    "NewEvalCompH35Estimate",
    "evalcomp_time_precision_summary",
    "format_paper_report",
    "newcomp_h35_summary",
    "newevalcomp_h35_estimates",
]
