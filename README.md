# Configuration Space Temporality (CST)

**A framework proposing that time emerges from energy-driven transitions through configuration space.**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

Configuration Space Temporality (CST) is a theoretical physics framework in which temporal evolution is not a fundamental parameter but an emergent phenomenon arising from transitions between configurations in a high-dimensional space. The theory introduces a **generator operator G** — a Witten Laplacian on configuration space — whose spectral properties encode the dynamical structure of physical systems.

This repository contains the computational code and results from the first empirical test of CST predictions on real quantum hardware data from IBM Quantum processors.

## Key Results

### Witten Laplacian for Coupled Transmons

We construct the CST operator **G** explicitly for a system of two coupled transmon qubits on the configuration torus T² = S¹ × S¹:

$$H_W = -\Delta + |\nabla W|^2 - \Delta W$$

where W = Φ/ℏ_eff and Φ is the transmon potential energy landscape including Josephson coupling.

### Empirical Findings

| Observable | Correlation | p-value |
|---|---|---|
| Spectral gap vs asymmetry ratio (parametric) | r = −0.95 | < 10⁻¹² |
| Fidelity susceptibility vs asymmetry (parametric) | r = +0.99 | < 10⁻²² |
| gap₀₃/gap₀₁ vs IBM gate length (772 real qubit pairs) | r = +0.14 | < 0.01 |
| Raw asymmetry ratio vs IBM gate length (baseline) | r = 0.004 | 0.90 |

The Witten spectrum extracts a statistically significant signal (r = 0.14) from IBM calibration data where raw parameters show zero correlation — an 18× improvement over the naive predictor. Level spacing statistics confirm integrable (Poisson) behavior of the operator.

### Honest Assessment

The r = 0.14 correlation is weak (explains ~2% of variance). The dominant factor in gate performance — the coupling constant g — is not available in public IBM calibration data. A definitive test requires access to coupling strengths, available in IBM Heron processors with tunable couplers.

## Repository Contents

### Code

| File | Description |
|---|---|
| `p8b_test.py` | Tests CST prediction P8b (E_gate ∝ g·√(E_J,A/E_J,B)) using IBM FakeProvider backends. Extracts E_J from qubit frequency/anharmonicity, computes correlations with gate lengths across 772 qubit pairs from 6 IBM processors. |
| `generate_paper.py` | Generates the academic paper PDF using ReportLab. Contains all analysis results, tables, and references. |

### Results

| File | Description |
|---|---|
| `p8b_real_results.json` | Raw numerical results: per-backend correlations, all qubit pair data, statistical summaries. |

### Figures

| File | Description |
|---|---|
| `G_final_test.png` | Witten Laplacian spectral correlations with IBM data (772 pairs, 6 backends). Main result figure. |
| `G_witten_analysis.png` | Parametric analysis: spectral gaps vs asymmetry ratio (r = −0.95), fidelity susceptibility (r = +0.99), Poisson level spacing statistics. |
| `G_operator_spectrum.png` | Initial Fokker-Planck operator spectrum (non-self-adjoint, unstable — included for completeness). |
| `p8b_real_results.png` | P8b prediction test: R_ij vs gate length showing null result (r = 0.004). |
| `p8b_scales.png` | Alternative scale analysis (log, power law, interactions) — all null. |
| `p8b_validation.png` | Split-half validation protocol results. |
| `p8b_results.png` | Summary of P8b test across all backends. |

## Requirements

```
python >= 3.9
numpy
scipy
qiskit >= 1.0
qiskit-ibm-runtime
matplotlib
reportlab  # for paper generation only
```

## Quick Start

```bash
# Run the P8b test (uses Qiskit FakeProvider, no IBM account needed)
python p8b_test.py

# Generate the paper PDF
pip install reportlab
python generate_paper.py
```

To test with live IBM backends (optional), set your IBM Quantum token:
```bash
export IBM_TOKEN="your_token_here"
python p8b_test.py
```

## Related Publications

- **F.D. Blum, Claude (Anthropic), Catalyst AI.** "Witten Laplacian on the Configuration Space of Coupled Transmons: Spectral Structure and Correlations with IBM Quantum Calibration Data" (2025). [Zenodo](https://doi.org/10.5281/zenodo.XXXXXXX)
- **F.D. Blum.** "Configuration Space Temporality v14: Emergent Time from Configuration Transitions" (2025). [Zenodo](https://doi.org/10.5281/zenodo.XXXXXXX)
- **F.D. Blum.** "Configuration Space Temporality v13" (2025). [Zenodo](https://doi.org/10.5281/zenodo.XXXXXXX)

## Theory Summary

CST proposes that for any physical system with configuration space **C**:

1. **States** are probability distributions ρ on C (points on a statistical manifold)
2. **Time** emerges from the generator **G** of transitions between configurations
3. **G** takes the form of a **Witten Laplacian**: H_W = −Δ + |∇W|² − ΔW, where W encodes the energy landscape
4. The **spectrum of G** determines transition rates, gate speeds, and dynamical timescales
5. **Spectral gaps** of G predict which transitions are fast (small gap → slow) or slow (large gap → fast)

For transmon qubits, the configuration space is the torus T² with coordinates (φ_A, φ_B) representing the superconducting phase differences. The potential W incorporates Josephson energies E_J and the coupling E_J^(c).

## Future Directions

- **Coupling constant access**: Test with IBM Heron tunable couplers where g is a known control parameter
- **Direct spectral observables**: Correlate G spectrum with Rabi frequencies, ZZ interaction rates
- **Multi-qubit systems**: Extend to T^n for n > 2 transmons
- **Analytic predictions**: Derive closed-form expressions for spectral gaps in limiting regimes

## Authors

- **Frédéric David Blum** — Theory, physics, direction ([Catalyst AI](https://catalyst-ai.co))
- **Claude** (Anthropic) — Computation, analysis, validation, writing
- **Catalyst AI** — Conceptual analysis, framework development

## License

MIT License

## Citation

```bibtex
@article{blum2025witten,
  title={Witten Laplacian on the Configuration Space of Coupled Transmons: 
         Spectral Structure and Correlations with IBM Quantum Calibration Data},
  author={Blum, Fr{\'e}d{\'e}ric David and Claude (Anthropic) and Catalyst AI},
  year={2025},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXXXX}
}
```
