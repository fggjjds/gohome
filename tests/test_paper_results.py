import pytest

from homcomp.paper_results import (
    evalcomp_time_precision_summary,
    format_paper_report,
    newcomp_h35_summary,
    newevalcomp_h35_estimates,
)


def test_evalcomp_reported_table_improves_time_and_precision():
    for row in evalcomp_time_precision_summary():
        assert row["evalcomp_vs_evalmod_time_reduction"] > 0
        assert row["evalcomp_vs_evalround_time_reduction"] > 0
        assert row["evalcomp_vs_evalmod_precision_gain_bits"] > 0
        assert row["evalcomp_vs_evalround_precision_gain_bits"] > 0


def test_newcomp_h35_reported_table_improves_complexity_and_runtime():
    for row in newcomp_h35_summary():
        assert row["complexity_vs_g_reduction"] > 0
        assert row["complexity_vs_h_reduction"] >= 0
        assert row["runtime_vs_g_speedup"] > 1
        assert row["runtime_vs_h_speedup"] > 1


def test_format_paper_report_mentions_both_papers():
    report = format_paper_report()

    assert "EvalComp" in report
    assert "NewCompH3&5" in report


def test_newevalcomp_h35_combined_estimate_improves_evalcomp_time():
    estimates = newevalcomp_h35_estimates(alpha=8, hcf_share=1.0)

    assert [row.family for row in estimates] == ["I", "II", "III"]
    for row in estimates:
        assert row.estimated_time_ms < row.evalcomp_time_ms
        assert row.time_reduction_vs_evalcomp == pytest.approx(
            1 - 17.0 / 26.9
        )
        assert row.evalcomp_precision_bits > 0


def test_newevalcomp_h35_share_can_model_partial_hcf_time():
    full = newevalcomp_h35_estimates(alpha=8, hcf_share=1.0)[0]
    half = newevalcomp_h35_estimates(alpha=8, hcf_share=0.5)[0]

    assert full.estimated_time_ms < half.estimated_time_ms < half.evalcomp_time_ms
