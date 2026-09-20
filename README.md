# MIE366 Current-Sense Resistor Optimizer

This repository contains the numerical resistor-selection tool and
supporting files developed for MIE366 Design Assignment 1 at the
University of Toronto.

## Design Target

The current-sense circuit is designed to approximate

\[
V_{OUT} = 8.714 V_{SENSE} + 0.25
\]

with:

- \(R_{SENSE} = 0.1 \Omega\)
- \(I_{LOAD} = 0\) to \(3.5\text{ A}\)
- \(V_{SENSE} = 0\) to \(0.35\text{ V}\)
- \(V_{OUT} = 0.25\) to \(3.3\text{ V}\)
- Standard E24 5% resistor values
- Ordinary resistors restricted to 1 kΩ to 100 kΩ

## Optimization Method

The Python script:

1. Generates permitted E24 resistor values.
2. Removes electrically duplicate resistor-ratio combinations.
3. Evaluates the nominal circuit transfer function.
4. Calculates maximum endpoint error.
5. Evaluates all \(2^7 = 128\) tolerance corners for the six ordinary
   resistors and the current-sense resistor.
6. Ranks candidate resistor configurations.

## Selected Configuration

| Component | Value |
|---|---:|
| Rsense | 0.1 Ω |
| R1 | 2.4 kΩ |
| R2 | 91 kΩ |
| R3 | 3 kΩ |
| R4 | 51 kΩ |
| Rg | 6.2 kΩ |
| Rf | 51 kΩ |

Nominal maximum endpoint error: approximately **0.236 mV**.

## Repository Contents

- `resistor_optimizer_v2.py` — optimization script
- `results/` — generated ranking tables
- `schematics/` — candidate and final circuit schematics

## Interactive Falstad Models

- [Candidate 1](https://www.falstad.com/s.php?s=Vyss9r)
- [Candidate 2](https://www.falstad.com/s.php?s=RiAqLA)

## Notes

This repository is provided as supporting material for the design report.
The report contains the engineering analysis; this repository provides
the computational implementation and reproducibility files.
