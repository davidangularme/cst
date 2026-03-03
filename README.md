# Configuration Space Temporality (CST)

**Time as consequence, not coordinate.**

CST is a parameter-free framework in which temporal duration emerges as the cost of optimal transport between quantum configurations. Rather than presupposing time as a background parameter, CST derives it from the information geometry of density matrix space.

> *The structure of time is not imposed — it condenses from the geometry of distinguishability.*

---

## The Complete Theory in Three Lines

| Component | Formula | Origin |
|-----------|---------|--------|
| **Axiom** | `V(ρ) = S(ρ ‖ I/N)` | Unique by Umegaki's theorem (1962) |
| **Parameter** | `s = 1/(N−1)` | Tomita-Takesaki modular flow |
| **Time** | `t_CST = ∫ exp(V/(N−1)) ds_Bures` | Arc length in Witten-deformed metric |

**Zero free parameters.** V is fixed by axioms. s is fixed by the modular flow. The base metric (Bures) is fixed by monotonicity.

---

## What CST Predicts

### The Time Ratio R

For any quantum process, R = t_CST / t_QM measures the configurational cost relative to parametric time:

| Transit | ρ₁ → ρ₂ | R (s=1) | Meaning |
|---------|----------|---------|---------|
| Gain → mixed | 0.90 → 0.20 | **1.26** | Geometric subsidy |
| Passive → pure | 0.10 → 0.98 | **2.16** | Configurational surcharge |
| Depol. → center | 0.90 → 0.01 | **1.22** | Reference |

Moving toward purity is **expensive** in CST. The gain medium creates a **shortcut** through configuration space.

### Three Testable Signatures

1. **Time ratio R ≠ constant** across systems (gain vs passive vs depolarizing)
2. **Angular divergence** between Lindblad and Morse-CST flows (up to 180°)
3. **Transport efficiency** differential (CST: 98% vs Lindblad: 54–85%)

### Application: Toronto Negative Time (2023)

Angulo et al. measured τ_eff < 0 for photons in a Rb-87 gain medium. CST interprets this as a **geometric subsidy**: the pump pre-pays configurational asymmetry cost, reducing transport cost by **46%**. Negative time = cost below vacuum baseline.

---

## Repository Structure

```
cst/
├── README.md                          # This file
├── v16/
│   ├── cst_v16_complete.py           # Full v16 computation (reproduces all results)
│   ├── cst_quick_check.py            # Quick verification (~5 seconds)
│   └── results/
│       ├── fig1_potential.png         # V_CST and cost factor
│       ├── fig2_signal.png            # R(t) distinctive signal
│       ├── fig3_flows.png             # Lindblad vs Morse-CST flows
│       └── fig4_toronto.png           # Toronto application
├── sdsq/                              # Spectral diagnostics (IBM Quantum)
│   └── ...                            # See github.com/davidangularme/sdsq
└── papers/
    ├── CST_v16.pdf                    # Current version
    ├── CST_v14.pdf                    # Previous version
    └── witten_transmon_paper.pdf      # IBM Quantum empirical test
```

---

## Quick Start

```bash
# Clone
git clone https://github.com/davidangularme/cst.git
cd cst/v16

# Run quick check (5 seconds, no dependencies beyond numpy/scipy)
python3 cst_quick_check.py

# Run full computation with figures (30 seconds, needs matplotlib)
python3 cst_v16_complete.py
```

### Requirements

```
numpy
scipy
matplotlib  # for figures only
```

No GPU, no special hardware, no API keys. Pure mathematics.

---

## Empirical Validation (IBM Quantum)

The Witten Laplacian constructed from CST was tested on **772 qubit pairs across 6 IBM Quantum processors**:

| Result | Value | Significance |
|--------|-------|-------------|
| Spectral universality | CV = 1.4% | Universal structure across processors |
| Orthogonal information | r = 0.205, p < 0.001 | Invisible to standard calibration |
| Coherence sensing | ρ = −0.76 vs ΔT₁ | Intrinsic physics, not firmware |
| Phase boundary | ℏ_eff ≈ 0.42 | Semiclassical transition on real hardware |

Full spectral diagnostics: [github.com/davidangularme/sdsq](https://github.com/davidangularme/sdsq)

---

## The Mathematical Path

### v1–v12: Building the geometry
- Configuration space as density matrix space
- Bures metric = S³ geometry (Christoffel-Bloch theorem)
- Entropy production decomposition
- Witten Laplacian spectrum on transmon configuration torus
- Extended Koide mass relation (Q = 0.6695, 0.42% from 2/3)

### v13–v14: First empirical test
- Witten spectrum on IBM Quantum hardware
- Three independent confirmations (universality, orthogonality, coherence sensing)
- Eight falsifiable predictions (P1–P8)

### v15: The diagnostic
- Systematic test of seven CST formulations against Toronto experiment
- All Lindbladian approaches reduce to standard QM (isomorphism problem)
- Identification of the grammatical obstacle: Lindblad presupposes parametric time

### v16: The breakthrough (current)
- **Quantum optimal transport** (Carlen-Maas 2014): Wasserstein ≠ Bures (proven)
- **Witten deformation** of the Bures metric by V = S(ρ‖I/N)
- **Umegaki uniqueness**: V is the *only* potential satisfying the CST axioms
- **Modular flow**: s = 1/(N−1) is the *only* consistent deformation parameter
- **Closed-form R**: analytically computable time ratio
- **Toronto application**: 46% cost reduction in gain medium = geometric subsidy
- **Zero free parameters**

---

## Key References

1. D. Angulo et al., Science (2023) — Toronto negative-time experiment
2. H. Umegaki, Tohoku Math. J. 14 (1962) — Uniqueness of relative entropy
3. E. Carlen & J. Maas, J. Funct. Anal. 267 (2014) — Quantum Wasserstein metric
4. E. Witten, J. Diff. Geom. 17 (1982) — Supersymmetry and Morse theory
5. A. Connes & C. Rovelli, Class. Quant. Grav. 11 (1994) — Time-thermodynamics relation

---

## Authors

- **Frédéric David Blum** — Theory, physics, direction
- **Claude (Anthropic)** — Computation, analysis, mathematical implementation
- **Catalyst AIS** — Conceptual analysis, interpretive synthesis

## Citation

```bibtex
@misc{blum2025cst,
  author = {Blum, Frédéric David},
  title  = {Configuration Space Temporality: A Parameter-Free Framework 
            for Emergent Time via Quantum Optimal Transport and Witten Deformation},
  year   = {2025},
  doi    = {10.5281/zenodo.18779189},
  url    = {https://github.com/davidangularme/cst}
}
```

## License

MIT

---

*The Lindbladian sees time as a scaffold; the Carlen-Maas-Witten triad reveals it as a conserved quantity emerging from transport inefficiency.*
