#!/usr/bin/env python3
"""Search E24 resistor combinations for the Candidate 2 current-sense circuit.

Circuit naming
--------------
R1 : reference-divider resistor from VREF node to GND
R2 : reference-divider resistor from +19 V to VREF node
R3 : resistor from VSENSE to the weighted-summing node
R4 : resistor from buffered VREF to the weighted-summing node
Rg : feedback resistor from the inverting input to GND
Rf : feedback resistor from VOUT to the inverting input
Rsense : 0.1-ohm current-sense resistor

Ideal circuit model
-------------------
VREF = Vsupply * R1 / (R1 + R2)
G    = 1 + Rf / Rg

VOUT(0 A) = G * [R3 / (R3 + R4)] * VREF
VOUT(full) = VOUT(0 A)
             + G * [R4 / (R3 + R4)] * Imax * Rsense

The script ranks resistor sets in two ways:
1. smallest nominal maximum endpoint error;
2. smallest worst-case endpoint error under independent resistor tolerances.

For tolerance checking, each of the six ordinary resistors and Rsense is
placed at its upper or lower tolerance limit, giving 2^7 = 128 corners.
No probability distribution is assumed.

Example
-------
python3 resistor_optimizer.py --top 10 --output-dir results
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np


# Standard E24 base values. Decade scaling is applied later.
E24_BASE_VALUES = (
    10, 11, 12, 13, 15, 16, 18, 20, 22, 24, 27, 30,
    33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91,
)

RESISTOR_NAMES = ("R1", "R2", "R3", "R4", "Rg", "Rf", "Rsense")


@dataclass(frozen=True)
class SearchConfig:
    """Electrical targets and search settings."""

    supply_voltage: float = 19.0
    max_load_current: float = 3.5
    target_output_zero: float = 0.25
    target_output_full: float = 3.3
    sense_resistance: float = 0.1
    resistor_tolerance: float = 0.05
    sense_resistor_tolerance: float = 0.05
    top_count: int = 10


def generate_e24_resistor_values() -> list[int]:
    """Return all E24 resistor values from 1 kOhm to 100 kOhm."""
    values = {
        base * decade
        for base in E24_BASE_VALUES
        for decade in (100, 1000, 10000)
        if 1_000 <= base * decade <= 100_000
    }
    return sorted(values)


def generate_unique_ratio_pairs(
    resistor_values: list[int],
    ratio_kind: str,
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Return unique exact resistor ratios and one representative pair each.

    For every pair (lower, upper), the stored ratio is upper/lower.

    ratio_kind = "divider"  -> R2/R1
    ratio_kind = "mixing"   -> R4/R3, limited to 2 ... 50
    ratio_kind = "feedback" -> Rf/Rg

    Electrically equivalent resistor scales are collapsed so that the search
    does not repeatedly evaluate identical ratios.
    """
    representatives: dict[Fraction, tuple[int, int]] = {}

    for lower_resistor in resistor_values:
        for upper_resistor in resistor_values:
            ratio = Fraction(upper_resistor, lower_resistor)

            if ratio_kind == "mixing" and not 2 <= ratio <= 50:
                continue

            representatives.setdefault(
                ratio,
                (lower_resistor, upper_resistor),
            )

    sorted_items = sorted(representatives.items())
    ratio_values = np.array([float(ratio) for ratio, _ in sorted_items])
    resistor_pairs = [pair for _, pair in sorted_items]
    return ratio_values, resistor_pairs


def calculate_pre_gain_bounds(
    mixing_ratio: float,
    divider_ratios: np.ndarray,
    config: SearchConfig,
) -> tuple[np.ndarray, ...]:
    """Calculate nominal and tolerance-bounded node voltages before gain.

    The second op-amp multiplies these node voltages by (1 + Rf/Rg).
    Vectorization over all divider ratios makes the search much faster.
    """
    tol = config.resistor_tolerance

    # Divider ratio = R2/R1, with R1 connected to ground.
    vref_nominal = config.supply_voltage / (1.0 + divider_ratios)
    vref_low = config.supply_voltage / (
        1.0 + divider_ratios * (1.0 + tol) / (1.0 - tol)
    )
    vref_high = config.supply_voltage / (
        1.0 + divider_ratios * (1.0 - tol) / (1.0 + tol)
    )

    # Mixing ratio = R4/R3.
    sense_weight_nominal = mixing_ratio / (1.0 + mixing_ratio)
    mixing_ratio_low = mixing_ratio * (1.0 - tol) / (1.0 + tol)
    mixing_ratio_high = mixing_ratio * (1.0 + tol) / (1.0 - tol)
    sense_weight_low = mixing_ratio_low / (1.0 + mixing_ratio_low)
    sense_weight_high = mixing_ratio_high / (1.0 + mixing_ratio_high)

    vsense_full_nominal = config.max_load_current * config.sense_resistance
    vsense_full_low = vsense_full_nominal * (1.0 - config.sense_resistor_tolerance)
    vsense_full_high = vsense_full_nominal * (1.0 + config.sense_resistor_tolerance)

    # Before final non-inverting gain.
    node_zero_nominal = vref_nominal * (1.0 - sense_weight_nominal)
    node_full_nominal = (
        sense_weight_nominal * vsense_full_nominal
        + (1.0 - sense_weight_nominal) * vref_nominal
    )

    node_zero_low = (1.0 - sense_weight_high) * vref_low
    node_zero_high = (1.0 - sense_weight_low) * vref_high

    node_full_low = np.minimum(
        sense_weight_low * vsense_full_low + (1.0 - sense_weight_low) * vref_low,
        sense_weight_high * vsense_full_low + (1.0 - sense_weight_high) * vref_low,
    )
    node_full_high = np.maximum(
        sense_weight_low * vsense_full_high + (1.0 - sense_weight_low) * vref_high,
        sense_weight_high * vsense_full_high + (1.0 - sense_weight_high) * vref_high,
    )

    return (
        node_zero_nominal,
        node_full_nominal,
        node_zero_low,
        node_zero_high,
        node_full_low,
        node_full_high,
    )


def endpoint_error_for_feedback_ratio(
    feedback_ratios: np.ndarray,
    pre_gain_bounds: tuple[np.ndarray, ...],
    config: SearchConfig,
    use_worst_case_tolerance: bool,
) -> np.ndarray:
    """Return maximum endpoint error for each Rf/Rg ratio."""
    (
        node_zero_nominal,
        node_full_nominal,
        node_zero_low,
        node_zero_high,
        node_full_low,
        node_full_high,
    ) = pre_gain_bounds

    if not use_worst_case_tolerance:
        gain = 1.0 + feedback_ratios
        zero_error = np.abs(gain * node_zero_nominal - config.target_output_zero)
        full_error = np.abs(gain * node_full_nominal - config.target_output_full)
        return np.maximum(zero_error, full_error)

    tol = config.resistor_tolerance
    gain_low = 1.0 + feedback_ratios * (1.0 - tol) / (1.0 + tol)
    gain_high = 1.0 + feedback_ratios * (1.0 + tol) / (1.0 - tol)

    return np.maximum.reduce(
        [
            config.target_output_zero - gain_low * node_zero_low,
            gain_high * node_zero_high - config.target_output_zero,
            config.target_output_full - gain_low * node_full_low,
            gain_high * node_full_high - config.target_output_full,
        ]
    )


def find_best_feedback_candidates(
    mixing_ratio: float,
    divider_ratios: np.ndarray,
    feedback_ratios: np.ndarray,
    config: SearchConfig,
    use_worst_case_tolerance: bool,
) -> list[tuple[float, int, int]]:
    """Find the best feedback ratio for every divider ratio efficiently.

    For fixed divider and mixing ratios, endpoint error is unimodal with
    feedback ratio. A discrete binary search locates the minimum; nearby
    values are then checked so the global top-N candidates are retained.
    """
    pre_gain_bounds = calculate_pre_gain_bounds(mixing_ratio, divider_ratios, config)

    left_index = np.zeros(len(divider_ratios), dtype=int)
    right_index = np.full(len(divider_ratios), len(feedback_ratios) - 1, dtype=int)

    while np.any(left_index < right_index):
        middle_index = (left_index + right_index) // 2
        next_index = np.minimum(middle_index + 1, len(feedback_ratios) - 1)

        middle_error = endpoint_error_for_feedback_ratio(
            feedback_ratios[middle_index],
            pre_gain_bounds,
            config,
            use_worst_case_tolerance,
        )
        next_error = endpoint_error_for_feedback_ratio(
            feedback_ratios[next_index],
            pre_gain_bounds,
            config,
            use_worst_case_tolerance,
        )

        move_right = middle_error > next_error
        active = left_index < right_index

        left_index = np.where(active & move_right, middle_index + 1, left_index)
        right_index = np.where(active & ~move_right, middle_index, right_index)

    # Check a small neighborhood around each discrete minimum.
    nearby_indices = np.clip(
        left_index[:, None]
        + np.arange(-config.top_count, config.top_count + 1),
        0,
        len(feedback_ratios) - 1,
    )

    expanded_bounds = tuple(value[:, None] for value in pre_gain_bounds)
    nearby_errors = endpoint_error_for_feedback_ratio(
        feedback_ratios[nearby_indices],
        expanded_bounds,
        config,
        use_worst_case_tolerance,
    )

    unique_candidate_ids = (
        np.arange(len(divider_ratios))[:, None] * len(feedback_ratios)
        + nearby_indices
    )

    _, first_occurrence = np.unique(
        unique_candidate_ids.ravel(),
        return_index=True,
    )

    flat_errors = nearby_errors.ravel()
    best_positions = first_occurrence[
        np.argsort(flat_errors[first_occurrence], kind="stable")[: config.top_count]
    ]

    results: list[tuple[float, int, int]] = []
    candidates_per_divider = nearby_indices.shape[1]

    for flat_position in best_positions:
        divider_index = int(flat_position // candidates_per_divider)
        feedback_index = int(nearby_indices.ravel()[flat_position])
        error = float(flat_errors[flat_position])
        results.append((error, divider_index, feedback_index))

    return results


def calculate_nominal_response(
    resistor_values: tuple[float, float, float, float, float, float, float],
    config: SearchConfig,
) -> tuple[float, float, float, float]:
    """Return VREF, closed-loop gain, VOUT at 0 A, and VOUT at full scale."""
    r1, r2, r3, r4, rg, rf, rsense = resistor_values

    vref = config.supply_voltage * r1 / (r1 + r2)
    closed_loop_gain = 1.0 + rf / rg

    reference_weight = r3 / (r3 + r4)
    sense_weight = r4 / (r3 + r4)

    vout_zero = closed_loop_gain * reference_weight * vref
    vout_full = (
        vout_zero
        + closed_loop_gain
        * sense_weight
        * config.max_load_current
        * rsense
    )

    return vref, closed_loop_gain, vout_zero, vout_full


def enumerate_tolerance_corners(
    resistor_values: tuple[float, float, float, float, float, float, float],
    config: SearchConfig,
) -> list[dict]:
    """Evaluate all 128 upper/lower resistor-tolerance corners."""
    corner_results: list[dict] = []

    for directions in product((-1, 1), repeat=7):
        actual_values = []

        for index, (nominal_value, direction) in enumerate(
            zip(resistor_values, directions)
        ):
            tolerance = (
                config.sense_resistor_tolerance
                if index == 6
                else config.resistor_tolerance
            )
            actual_values.append(nominal_value * (1.0 + direction * tolerance))

        vref, gain, vout_zero, vout_full = calculate_nominal_response(
            tuple(actual_values),
            config,
        )

        corner_results.append(
            {
                "directions": dict(zip(RESISTOR_NAMES, directions)),
                "actual_ohm": dict(zip(RESISTOR_NAMES, actual_values)),
                "Vref_V": vref,
                "gain": gain,
                "Vout_0A_V": vout_zero,
                "Vout_full_V": vout_full,
            }
        )

    return corner_results


def build_result_record(
    candidate: tuple[float, float, float, float, float, float, float],
    config: SearchConfig,
) -> dict:
    """Convert one resistor candidate into a complete result row."""
    _, r1, r2, r3, r4, rg, rf = candidate
    resistor_values = (
        r1,
        r2,
        r3,
        r4,
        rg,
        rf,
        config.sense_resistance,
    )

    vref, gain, vout_zero, vout_full = calculate_nominal_response(
        resistor_values,
        config,
    )

    tolerance_corners = enumerate_tolerance_corners(resistor_values, config)

    zero_outputs = [row["Vout_0A_V"] for row in tolerance_corners]
    full_outputs = [row["Vout_full_V"] for row in tolerance_corners]

    vout_zero_min = min(zero_outputs)
    vout_zero_max = max(zero_outputs)
    vout_full_min = min(full_outputs)
    vout_full_max = max(full_outputs)

    nominal_max_error = max(
        abs(vout_zero - config.target_output_zero),
        abs(vout_full - config.target_output_full),
    )

    worst_case_error = max(
        abs(vout_zero_min - config.target_output_zero),
        abs(vout_zero_max - config.target_output_zero),
        abs(vout_full_min - config.target_output_full),
        abs(vout_full_max - config.target_output_full),
    )

    result = {
        f"{name}_ohm": value
        for name, value in zip(RESISTOR_NAMES, resistor_values)
    }

    result.update(
        Vref_V=vref,
        gain=gain,
        Vout_0A_V=vout_zero,
        Vout_full_V=vout_full,
        nominal_max_error_mV=nominal_max_error * 1000.0,
        Vout_0A_min_V=vout_zero_min,
        Vout_0A_max_V=vout_zero_max,
        Vout_full_min_V=vout_full_min,
        Vout_full_max_V=vout_full_max,
        worst_case_error_mV=worst_case_error * 1000.0,
        resistor_tolerance=config.resistor_tolerance,
        sense_tolerance=config.sense_resistor_tolerance,
    )

    return result


def search_resistor_configurations(
    config: SearchConfig,
    resistor_values: list[int] | None = None,
) -> dict:
    """Search Candidate 2 resistor ratios and return nominal/robust rankings."""
    resistor_values = (
        generate_e24_resistor_values()
        if resistor_values is None
        else resistor_values
    )

    divider_ratios, divider_pairs = generate_unique_ratio_pairs(
        resistor_values,
        "divider",
    )
    feedback_ratios, feedback_pairs = generate_unique_ratio_pairs(
        resistor_values,
        "feedback",
    )
    mixing_ratios, mixing_pairs = generate_unique_ratio_pairs(
        resistor_values,
        "mixing",
    )

    nominal_pool: list[tuple] = []
    robust_pool: list[tuple] = []
    start_time = time.perf_counter()

    for index, (mixing_ratio, (r3, r4)) in enumerate(
        zip(mixing_ratios, mixing_pairs),
        start=1,
    ):
        for use_worst_case, candidate_pool in (
            (False, nominal_pool),
            (True, robust_pool),
        ):
            feedback_candidates = find_best_feedback_candidates(
                mixing_ratio,
                divider_ratios,
                feedback_ratios,
                config,
                use_worst_case,
            )

            for error, divider_index, feedback_index in feedback_candidates:
                r1, r2 = divider_pairs[divider_index]
                rg, rf = feedback_pairs[feedback_index]
                candidate_pool.append((error, r1, r2, r3, r4, rg, rf))

            candidate_pool[:] = sorted(candidate_pool)[: config.top_count]

        if index % 100 == 0:
            print(
                f"Searched {index}/{len(mixing_ratios)} mixing ratios",
                flush=True,
            )

    rankings = []

    for ranking_name, candidate_pool in (
        ("theoretical", nominal_pool),
        ("robust", robust_pool),
    ):
        result_rows = [
            build_result_record(candidate, config)
            for candidate in candidate_pool
        ]

        # Confirm that the fast search metric agrees with the explicit
        # endpoint/corner calculation used in the exported results.
        error_key = (
            "worst_case_error_mV"
            if ranking_name == "robust"
            else "nominal_max_error_mV"
        )

        for candidate, row in zip(candidate_pool, result_rows):
            assert abs(candidate[0] * 1000.0 - row[error_key]) < 1e-8

        for rank, row in enumerate(result_rows, start=1):
            row["rank"] = rank

        rankings.append(result_rows)

    return {
        "config": asdict(config),
        "search": {
            "E24_count": len(resistor_values),
            "divider_ratios": len(divider_ratios),
            "feedback_ratios": len(feedback_ratios),
            "mixing_ratios": len(mixing_ratios),
            "ratio_designs": (
                len(divider_ratios)
                * len(feedback_ratios)
                * len(mixing_ratios)
            ),
            "seconds": time.perf_counter() - start_time,
            "corners_per_design": 128,
        },
        "theoretical": rankings[0],
        "robust": rankings[1],
    }


def write_results(search_results: dict, output_directory: Path) -> None:
    """Write ranking tables to CSV and complete metadata to JSON."""
    output_directory.mkdir(parents=True, exist_ok=True)

    for ranking_name in ("theoretical", "robust"):
        rows = search_results[ranking_name]
        field_names = ["rank"] + [key for key in rows[0] if key != "rank"]

        csv_path = output_directory / f"{ranking_name}_best.csv"
        with csv_path.open("w", newline="", encoding="utf-8-sig") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=field_names)
            writer.writeheader()
            writer.writerows(rows)

    json_path = output_directory / "analysis.json"
    json_path.write_text(
        json.dumps(search_results, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.05,
        help="Tolerance of ordinary resistors (default: 0.05 = 5%%).",
    )
    parser.add_argument(
        "--sense-tolerance",
        type=float,
        default=None,
        help="Tolerance of Rsense. Defaults to --tolerance.",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of results to retain in each ranking.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV and JSON output files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    sense_tolerance = (
        args.tolerance
        if args.sense_tolerance is None
        else args.sense_tolerance
    )

    if not 0 <= args.tolerance < 1:
        raise SystemExit("--tolerance must be in [0, 1).")
    if not 0 <= sense_tolerance < 1:
        raise SystemExit("--sense-tolerance must be in [0, 1).")
    if args.top < 1:
        raise SystemExit("--top must be at least 1.")

    config = SearchConfig(
        resistor_tolerance=args.tolerance,
        sense_resistor_tolerance=sense_tolerance,
        top_count=args.top,
    )

    search_results = search_resistor_configurations(config)
    write_results(search_results, args.output_dir)

    for ranking_name in ("theoretical", "robust"):
        best_result = search_results[ranking_name][0]
        print(
            ranking_name,
            json.dumps(best_result, ensure_ascii=False),
        )


if __name__ == "__main__":
    main()
