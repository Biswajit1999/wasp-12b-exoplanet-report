from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import analyze_orbital_decay as decay


def test_weighted_quadratic_fit_recovers_curvature():
    x = np.arange(-20, 21, dtype=float)
    sigma = np.full_like(x, 2e-6)
    truth = 1e-4 + 2e-7 * x - 4e-10 * x**2
    fit = decay.weighted_fit(x, truth, sigma, degree=2)
    assert np.allclose(fit["coefficients"], [1e-4, 2e-7, -4e-10], rtol=1e-8)


def test_tidal_quality_factor_is_positive():
    assert 1e4 < decay.tidal_quality_factor(-1e-9) < 1e7


def test_published_timing_table_reproduces_decay_scale():
    data = decay.load_published_timings()
    linear = decay.ephemeris_fit(data, degree=1)
    quadratic = decay.ephemeris_fit(data, degree=2)
    coefficient = quadratic["coefficients"][2]
    period = quadratic["coefficients"][1]
    period_dot_ms_per_year = 2 * coefficient / period * 86400 * 1000 * 365.25
    assert len(data["time"]) == 158
    assert -33 < period_dot_ms_per_year < -26
    assert linear["bic"] - quadratic["bic"] > 100


def test_2024_vizier_compilation_reproduces_decay_and_source_robustness():
    data = decay.load_leonardi_transits()
    assert len(data["time"]) == 379
    assert data["epoch"].min() == -2833
    assert data["epoch"].max() == 1869
    linear = decay.ephemeris_fit(data, degree=1)
    quadratic = decay.ephemeris_fit(data, degree=2)
    estimate, formal_error = decay.period_derivative(quadratic)
    assert -32 < estimate < -30
    assert 0 < formal_error < 2
    assert linear["bic"] - quadratic["bic"] > 1_000
    sensitivity = decay.source_jackknife(data)
    assert len(sensitivity) >= 10
    assert min(row["period_dot_ms_per_year"] for row in sensitivity) > -33
    assert max(row["period_dot_ms_per_year"] for row in sensitivity) < -29


def test_yee_timing_set_prefers_decay_over_apsidal_proxy():
    data = decay.load_published_timings()
    quadratic = decay.ephemeris_fit(data, degree=2)
    apsidal = decay.apsidal_precession_fit(data)
    assert apsidal["bic"] - quadratic["bic"] > 10
    assert 5 <= apsidal["precession_period_years"] <= 100


def test_full_analysis_writes_new_machine_readable_outputs():
    result = decay.main()
    assert len(result["supported"]) >= 110
    assert result["scatter_scaled_period_dot_error_ms_per_year"] > 0
    assert decay.SENSITIVITY_FILE.stat().st_size > 500
    assert decay.SENSITIVITY_FIGURE.stat().st_size > 10_000
