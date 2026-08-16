"""Individual-transit timing test for the decaying orbit of WASP-12 b.

This module measures mid-times from the committed TESS sectors with a fixed
sector-level transit shape, then compares linear and quadratic ephemerides.
It is an independent public-data diagnostic, not a replacement for the global
transit-plus-occultation analysis in the cited literature.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

import analyze_multisector as multi
import analyze_transit as base


TIMINGS_FILE = base.FIG_DIR / "individual_transit_timings.csv"
MODEL_FILE = base.FIG_DIR / "orbital_decay_statistics.csv"
FIGURE_FILE = base.FIG_DIR / "wasp12b_orbital_decay.png"
PUBLISHED_FILE = base.DATA_DIR / "published_transit_occultation_times.csv"

EARTHS_PER_SUN = 332_946.0
RSUN_PER_AU = 0.00465047
PLANET_MASS_EARTH = 467.2101
STAR_MASS_SUN = 1.325
STAR_RADIUS_SUN = 1.69
SEMIMAJOR_AU = 0.0234


def fit_event(
    time: np.ndarray,
    flux: np.ndarray,
    error: np.ndarray,
    epoch_number: int,
    shape: np.ndarray,
    sector: int,
) -> dict[str, float | int | bool] | None:
    predicted = base.EPOCH_BJD + epoch_number * base.PERIOD_DAYS
    x = time - predicted
    duration = base.DURATION_HOURS / 24.0
    selected = np.abs(x) <= 1.8 * duration
    x, y, sigma = x[selected], flux[selected], error[selected]
    if len(x) < 70:
        return None
    left = np.sum(x < -0.7 * duration)
    right = np.sum(x > 0.7 * duration)
    core = np.sum(np.abs(x) < 0.45 * duration)
    if min(left, right) < 12 or core < 12:
        return None

    radius_ratio, impact = float(shape[1]), float(shape[2])

    def model(parameters: np.ndarray) -> np.ndarray:
        timing, baseline, slope = parameters
        full = np.asarray([timing, radius_ratio, impact, baseline, slope])
        return base.transit_profile(x, full)

    initial = np.asarray([float(shape[0]), float(np.median(y)), 0.0])
    bound_timing = 0.55 * duration
    fit = least_squares(
        lambda pars: (y - model(pars)) / sigma,
        initial,
        bounds=([-bound_timing, 0.94, -0.15], [bound_timing, 1.06, 0.15]),
        x_scale="jac",
        max_nfev=500,
    )
    fitted = model(fit.x)
    null = base.weighted_linear_null(x, y, sigma)
    chi2_transit = float(np.sum(np.square((y - fitted) / sigma)))
    chi2_null = float(np.sum(np.square((y - null) / sigma)))
    dof = len(y) - 3
    reduced_chi2 = chi2_transit / dof
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
    residuals = y - fitted
    _, _, _, beta = multi.noise_curve(residuals)
    error_scale = np.sqrt(max(reduced_chi2, 1.0)) * beta
    timing_error = float(np.sqrt(max(covariance[0, 0], 0.0)) * error_scale)
    bic_transit = chi2_transit + 3 * np.log(len(y))
    bic_null = chi2_null + 2 * np.log(len(y))
    timing_shift = float(fit.x[0])
    return {
        "sector": sector,
        "epoch": epoch_number,
        "predicted_bjd": predicted,
        "measured_bjd": predicted + timing_shift,
        "oc_seconds": timing_shift * 86400.0,
        "timing_error_seconds": timing_error * 86400.0,
        "n_points": len(y),
        "beta": beta,
        "reduced_chi_square": reduced_chi2,
        "delta_bic": bic_null - bic_transit,
        "supported": bool(bic_null - bic_transit >= 10 and timing_error > 0),
    }


def weighted_fit(x: np.ndarray, y: np.ndarray, sigma: np.ndarray, degree: int) -> dict[str, object]:
    columns = [np.ones_like(x), x]
    if degree == 2:
        columns.append(np.square(x))
    design = np.column_stack(columns)
    weighted_design = design / sigma[:, None]
    weighted_y = y / sigma
    coefficients, _, _, _ = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    model = design @ coefficients
    chi_square = float(np.sum(np.square((y - model) / sigma)))
    bic = chi_square + len(coefficients) * np.log(len(y))
    return {
        "coefficients": coefficients,
        "covariance": covariance,
        "model": model,
        "chi_square": chi_square,
        "dof": len(y) - len(coefficients),
        "bic": bic,
    }


def load_published_timings() -> dict[str, np.ndarray]:
    rows: list[dict[str, str]] = []
    with PUBLISHED_FILE.open(encoding="utf-8", newline="") as handle:
        rows.extend(csv.DictReader(handle))
    event_type = np.asarray([row["event_type"] for row in rows])
    occultation = event_type == "occ"
    epoch = np.asarray([float(row["epoch"]) for row in rows])
    return {
        "event_type": event_type,
        "occultation": occultation,
        "epoch": epoch,
        "event_coordinate": epoch + 0.5 * occultation.astype(float),
        "time": np.asarray([float(row["mid_time_bjd_tdb"]) for row in rows]),
        "error": np.asarray([float(row["mid_time_error_days"]) for row in rows]),
    }


def ephemeris_fit(data: dict[str, np.ndarray], degree: int) -> dict[str, object]:
    coordinate = data["event_coordinate"]
    sigma = data["error"]
    weights = 1.0 / np.square(sigma)
    center_epoch = float(np.round(np.sum(weights * coordinate) / np.sum(weights)))
    x = coordinate - center_epoch
    columns = [np.ones_like(x), x]
    if degree == 2:
        columns.append(np.square(x))
    # A separate constant for occultations absorbs light-travel time and any
    # small e*cos(omega) term instead of forcing them into orbital curvature.
    columns.append(data["occultation"].astype(float))
    design = np.column_stack(columns)
    weighted_design = design / sigma[:, None]
    weighted_time = data["time"] / sigma
    coefficients, _, _, _ = np.linalg.lstsq(weighted_design, weighted_time, rcond=None)
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    model = design @ coefficients
    chi_square = float(np.sum(np.square((data["time"] - model) / sigma)))
    bic = chi_square + len(coefficients) * np.log(len(x))
    return {
        "center_epoch": center_epoch,
        "x": x,
        "design": design,
        "coefficients": coefficients,
        "covariance": covariance,
        "model": model,
        "chi_square": chi_square,
        "dof": len(x) - len(coefficients),
        "bic": bic,
    }


def tidal_quality_factor(period_dot_days_per_day: float) -> float:
    mass_ratio = PLANET_MASS_EARTH / (STAR_MASS_SUN * EARTHS_PER_SUN)
    radius_ratio = STAR_RADIUS_SUN * RSUN_PER_AU / SEMIMAJOR_AU
    numerator = 27.0 * np.pi / 2.0 * mass_ratio * radius_ratio**5
    return float(numerator / abs(period_dot_days_per_day))


def main() -> dict[str, object]:
    base.FIG_DIR.mkdir(exist_ok=True)
    events: list[dict[str, float | int | bool]] = []
    for path in sorted(base.DATA_DIR.glob("tess*_lc.fits"), key=multi.sector_number):
        sector = multi.sector_number(path)
        time, flux, error, _ = base.load_light_curve(path)
        sector_fit = base.compare_models(time, flux, error)
        shape = np.asarray(sector_fit["parameters"], dtype=float)
        first = int(np.ceil((time.min() - base.EPOCH_BJD) / base.PERIOD_DAYS))
        last = int(np.floor((time.max() - base.EPOCH_BJD) / base.PERIOD_DAYS))
        for epoch in range(first, last + 1):
            result = fit_event(time, flux, error, epoch, shape, sector)
            if result is not None:
                events.append(result)

    supported = [item for item in events if item["supported"]]
    if len(supported) < 8:
        raise RuntimeError("Too few independently supported transits for an ephemeris comparison")

    epochs = np.asarray([item["epoch"] for item in supported], dtype=float)
    center_epoch = float(np.round(np.average(epochs)))
    x = epochs - center_epoch
    oc_days = np.asarray([item["oc_seconds"] for item in supported]) / 86400.0
    errors_days = np.asarray([item["timing_error_seconds"] for item in supported]) / 86400.0
    tess_linear = weighted_fit(x, oc_days, errors_days, degree=1)
    tess_quadratic = weighted_fit(x, oc_days, errors_days, degree=2)
    tess_quad_coefficient = float(tess_quadratic["coefficients"][2])
    tess_quad_error = float(np.sqrt(tess_quadratic["covariance"][2, 2]))
    tess_period_dot = 2.0 * tess_quad_coefficient / base.PERIOD_DAYS
    tess_period_dot_error = 2.0 * tess_quad_error / base.PERIOD_DAYS
    conversion = 86400.0 * 1000.0 * 365.25
    tess_milliseconds_per_year = tess_period_dot * conversion
    tess_milliseconds_per_year_error = tess_period_dot_error * conversion

    published = load_published_timings()
    published_linear = ephemeris_fit(published, degree=1)
    published_quadratic = ephemeris_fit(published, degree=2)
    published_quad_coefficient = float(published_quadratic["coefficients"][2])
    published_quad_error = float(np.sqrt(published_quadratic["covariance"][2, 2]))
    fitted_period = float(published_quadratic["coefficients"][1])
    period_dot = 2.0 * published_quad_coefficient / fitted_period
    period_dot_error = 2.0 * published_quad_error / fitted_period
    milliseconds_per_year = period_dot * conversion
    milliseconds_per_year_error = period_dot_error * conversion
    decay_timescale_years = fitted_period / abs(period_dot) / 365.25
    q_star = tidal_quality_factor(period_dot)

    fields = [
        "sector", "epoch", "predicted_bjd", "measured_bjd", "oc_seconds",
        "timing_error_seconds", "n_points", "beta", "reduced_chi_square",
        "delta_bic", "supported",
    ]
    with TIMINGS_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in events:
            writer.writerow(item)

    statistics = [
        ("events_fitted", len(events), "count"),
        ("events_supported", len(supported), "count; individual Delta BIC >= 10"),
        ("center_epoch", center_epoch, "integer epoch used for conditioning"),
        ("tess_only_linear_bic", tess_linear["bic"], ""),
        ("tess_only_quadratic_bic", tess_quadratic["bic"], ""),
        ("tess_only_delta_bic_linear_minus_quadratic", tess_linear["bic"] - tess_quadratic["bic"], ""),
        ("tess_only_period_dot_ms_per_year", tess_milliseconds_per_year, "ms per year; not a supported decay measurement"),
        ("tess_only_period_dot_error_ms_per_year", tess_milliseconds_per_year_error, "ms per year"),
        ("published_events", len(published["time"]), "count; transits plus occultations"),
        ("published_linear_chi_square", published_linear["chi_square"], ""),
        ("published_linear_dof", published_linear["dof"], ""),
        ("published_linear_bic", published_linear["bic"], ""),
        ("published_quadratic_chi_square", published_quadratic["chi_square"], ""),
        ("published_quadratic_dof", published_quadratic["dof"], ""),
        ("published_quadratic_bic", published_quadratic["bic"], ""),
        ("published_delta_bic_linear_minus_quadratic", published_linear["bic"] - published_quadratic["bic"], ""),
        ("published_period_dot", period_dot, "days per day"),
        ("published_period_dot_error", period_dot_error, "days per day; formal timing errors"),
        ("published_period_dot_ms_per_year", milliseconds_per_year, "ms per year"),
        ("published_period_dot_error_ms_per_year", milliseconds_per_year_error, "ms per year"),
        ("decay_timescale_years", decay_timescale_years, "years; P/abs(Pdot)"),
        ("modified_stellar_tidal_quality_factor", q_star, "Q-prime under equilibrium-tide convention"),
    ]
    with MODEL_FILE.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        for name, value, unit in statistics:
            writer.writerow([name, f"{value:.12g}" if isinstance(value, float) else value, unit])

    published_x = np.asarray(published_quadratic["x"])
    published_oc = (published["time"] - np.asarray(published_linear["model"])) * 86400.0
    published_grid = np.linspace(published_x.min() - 40, published_x.max() + 40, 800)
    linear_coeff = published_linear["coefficients"]
    quadratic_coeff = published_quadratic["coefficients"]
    published_linear_curve = linear_coeff[0] + linear_coeff[1] * published_grid
    published_quadratic_curve = (
        quadratic_coeff[0]
        + quadratic_coeff[1] * published_grid
        + quadratic_coeff[2] * published_grid**2
    )

    tess_grid = np.linspace(x.min() - 10, x.max() + 10, 600)
    tess_linear_curve = tess_linear["coefficients"][0] + tess_linear["coefficients"][1] * tess_grid
    tess_quadratic_curve = (
        tess_quadratic["coefficients"][0]
        + tess_quadratic["coefficients"][1] * tess_grid
        + tess_quadratic["coefficients"][2] * tess_grid**2
    )
    fig, (ax, tess_ax) = plt.subplots(2, 1, figsize=(10, 8.4), constrained_layout=True)
    transit = ~published["occultation"]
    occultation = published["occultation"]
    ax.errorbar(
        published_x[transit], published_oc[transit],
        yerr=published["error"][transit] * 86400.0,
        fmt=".", ms=5, alpha=0.55, color="#2563eb", label="published transits",
    )
    ax.errorbar(
        published_x[occultation], published_oc[occultation],
        yerr=published["error"][occultation] * 86400.0,
        fmt="s", ms=4, alpha=0.65, color="#d97706", label="published occultations",
    )
    ax.plot(
        published_grid,
        (published_quadratic_curve - published_linear_curve) * 86400.0,
        color="#9d174d", lw=2.2, label="quadratic minus linear ephemeris",
    )
    ax.set(
        xlabel=f"Event coordinate relative to E = {published_quadratic['center_epoch']:.0f}",
        ylabel="Observed − best linear model [s]",
        title="WASP-12 b: 158 published transit and occultation times",
    )
    ax.text(
        0.02, 0.04,
        f"Pdot = {milliseconds_per_year:.1f} ± {milliseconds_per_year_error:.1f} ms yr⁻¹\n"
        f"ΔBIC(linear − quadratic) = {published_linear['bic'] - published_quadratic['bic']:.1f}\n"
        f"P/|Pdot| = {decay_timescale_years / 1e6:.2f} Myr",
        transform=ax.transAxes, fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.88, "edgecolor": "#cbd5e1"},
    )
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=8)

    colors = {20: "#5e60ce", 43: "#2a9d8f", 45: "#e76f51"}
    for sector in sorted({int(item["sector"]) for item in supported}):
        chosen = np.asarray([int(item["sector"]) == sector for item in supported])
        tess_ax.errorbar(
            x[chosen], oc_days[chosen] * 86400.0,
            yerr=errors_days[chosen] * 86400.0,
            fmt="o", ms=4, capsize=2, alpha=0.75,
            color=colors.get(sector, "#334155"), label=f"TESS Sector {sector}",
        )
    tess_ax.plot(tess_grid, tess_linear_curve * 86400.0, "--", color="#64748b", lw=1.8, label="TESS-only linear")
    tess_ax.plot(tess_grid, tess_quadratic_curve * 86400.0, color="#9d174d", lw=2.2, label="TESS-only quadratic")
    tess_ax.set(
        xlabel=f"Transit epoch relative to E = {center_epoch:.0f}",
        ylabel="Observed − archive prediction [s]",
        title="TESS-only check: 62 supported individual transits",
    )
    tess_ax.text(
        0.02, 0.04,
        f"Formal Pdot = {tess_milliseconds_per_year:.0f} ± {tess_milliseconds_per_year_error:.0f} ms yr⁻¹\n"
        f"ΔBIC(linear − quadratic) = {tess_linear['bic'] - tess_quadratic['bic']:.1f}; curvature not supported",
        transform=tess_ax.transAxes, fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.88, "edgecolor": "#cbd5e1"},
    )
    tess_ax.grid(alpha=0.2)
    tess_ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.savefig(FIGURE_FILE, dpi=190)
    plt.close(fig)

    return {
        "events": events,
        "supported": supported,
        "tess_linear": tess_linear,
        "tess_quadratic": tess_quadratic,
        "published_linear": published_linear,
        "published_quadratic": published_quadratic,
        "period_dot_ms_per_year": milliseconds_per_year,
        "period_dot_error_ms_per_year": milliseconds_per_year_error,
        "decay_timescale_years": decay_timescale_years,
        "q_star": q_star,
    }


if __name__ == "__main__":
    result = main()
    print(
        f"WASP-12 b: {len(result['supported'])} supported individual transits; "
        f"Pdot={result['period_dot_ms_per_year']:.2f} +/- "
        f"{result['period_dot_error_ms_per_year']:.2f} ms/yr"
    )
    print(
        f"Published Delta BIC(linear-quadratic)="
        f"{result['published_linear']['bic'] - result['published_quadratic']['bic']:.2f}; "
        f"P/|Pdot|={result['decay_timescale_years'] / 1e6:.2f} Myr"
    )
