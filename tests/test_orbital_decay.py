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
