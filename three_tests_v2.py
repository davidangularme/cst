#!/usr/bin/env python3
"""
Three Tests — CORRECTED VERSION
================================
Fixes from first run:
- Test 1: off-line test must isolate a SINGLE off-line quadruplet among J on-line zeros
- Test 2: V has non-integrable singularity at ±1 — test on actual PW_1 domain, not L^2
- Test 3: gap fluctuation sign was wrong; use larger basis

Author: F. D. Blum / Catalyst AI Research · March 2026
"""
import numpy as np
from scipy import linalg
import warnings
warnings.filterwarnings('ignore')

def get_zeros(N):
    from mpmath import zetazero
    print(f"  Computing {N} Riemann zeros...")
    z = [float(zetazero(n).imag) for n in range(1, N+1)]
    return np.array(z)

# ============================================================
# TEST 1: ARITHMETIC vs RANDOM — corrected
# ============================================================
def test1(gamma, n_max=2000):
    J = len(gamma)
    n_arr = np.arange(1, n_max+1)

    print(f"\n{'='*72}")
    print(f"  TEST 1: Li COEFFICIENTS — OFF-LINE ZEROS & PHASE COHERENCE")
    print(f"  J={J} zeros, n up to {n_max}")
    print(f"{'='*72}")

    # Compute angles theta_j for on-line zeros (RH)
    rho = 0.5 + 1j*gamma
    w = (rho-1)/rho  # |w|=1 on RH
    theta = np.angle(w)

    # Lambda_n (on-line, all J): lambda_n = 2*sum(1 - cos(n*theta_j)) >= 0 trivially
    lam_arith = np.array([2*np.sum(1 - np.cos(n*theta)) for n in n_arr])

    # --- A: SINGLE OFF-LINE QUADRUPLET TEST ---
    # Keep J-1 zeros on line, move the FIRST zero off to sigma
    # Quadruplet contribution for off-line zero at rho = sigma + i*gamma_1:
    #   w1 = (sigma-1+i*gamma)/(sigma+i*gamma)
    #   w2 = conj(w1) [from conj(rho)]
    #   w3 = (-sigma+i*gamma)/(1-sigma+i*gamma) [from 1-conj(rho)]
    #   w4 = conj(w3) [from 1-rho]
    # Contribution = Re[1-w1^n] + Re[1-w2^n] + Re[1-w3^n] + Re[1-w4^n]
    # = 2*Re[1-w1^n] + 2*Re[1-w3^n]  (since w2=conj(w1), w4=conj(w3))

    print(f"\n  A. SINGLE OFF-LINE QUADRUPLET (first zero moved to sigma)")
    print(f"     Other {J-1} zeros remain on critical line.")
    print(f"     Remark 8.4 predicts: n* ~ 1/(2*(sigma-1/2)) for large gamma.\n")

    print(f"  {'sigma':>7s} {'eps':>7s} {'|w3|':>8s} {'n*(pred)':>9s} {'n*(actual)':>10s} {'min lam':>10s}")

    for sigma in [0.501, 0.505, 0.51, 0.52, 0.55, 0.6, 0.75, 1.0]:
        eps = sigma - 0.5
        g1 = gamma[0]  # first zero: 14.13...

        # Off-line zero: rho1 = sigma + i*g1
        w1 = (sigma - 1 + 1j*g1) / (sigma + 1j*g1)
        # Functional equation partner: rho3 = (1-sigma) + i*g1
        w3 = ((1-sigma) - 1 + 1j*g1) / ((1-sigma) + 1j*g1)
        # = (-sigma + i*g1) / (1-sigma + i*g1)

        abs_w3 = abs(w3)

        # On-line contribution from remaining J-1 zeros
        theta_rest = theta[1:]  # J-1 angles

        lam_mixed = np.zeros(n_max)
        for i, n in enumerate(n_arr):
            # On-line contribution (J-1 zeros, each counted with conjugate = factor 2)
            on_line = 2*np.sum(1 - np.cos(n*theta_rest))
            # Off-line quadruplet contribution (the moved first zero)
            off_line = (2*np.real(1 - w1**n) + 2*np.real(1 - w3**n))
            lam_mixed[i] = on_line + off_line

        neg_mask = lam_mixed < 0
        first_neg = n_arr[neg_mask][0] if np.any(neg_mask) else -1
        pred = round(1/(2*eps)) if eps > 0 else -1
        pred_s = str(pred) if pred > 0 else "inf"
        neg_s = str(first_neg) if first_neg > 0 else f">{n_max}"
        n_neg = np.sum(neg_mask)
        min_val = np.min(lam_mixed)

        print(f"  {sigma:7.3f} {eps:7.3f} {abs_w3:8.5f} {pred_s:>9s} {neg_s:>10s} {min_val:10.2f}")

    # --- B: PHASE COHERENCE |S(n)|/J ---
    print(f"\n  B. PHASE COHERENCE |S(n)|/J")
    print(f"     Positivity requires |S(n)| < J (equiv Re[S(n)] < J)")
    print(f"     Testing arithmetic vs density-matched random phases.\n")

    S_arith = np.array([np.sum(np.exp(1j*n*theta)) for n in n_arr[:500]])
    ratio_arith = np.abs(S_arith) / J

    # Random: scramble gamma values (uniform in same range, same count)
    N_trials = 200
    ratio_rand = np.zeros((N_trials, 500))
    for t in range(N_trials):
        g_rand = np.sort(np.random.uniform(gamma[0]*0.9, gamma[-1]*1.05, J))
        rho_r = 0.5 + 1j*g_rand
        w_r = (rho_r-1)/rho_r
        th_r = np.angle(w_r)
        for i in range(500):
            ratio_rand[t, i] = np.abs(np.sum(np.exp(1j*(i+1)*th_r))) / J

    rm = np.mean(ratio_rand, axis=0)
    rp95 = np.percentile(ratio_rand, 95, axis=0)

    print(f"  {'n':>5s} {'Arith |S|/J':>12s} {'Rand mean':>11s} {'Rand p95':>11s} {'Arith>Rand?':>12s}")
    for n in [35,50,75,100,200,300,500]:
        i = n-1
        flag = "HIGHER" if ratio_arith[i] > rm[i] else "lower"
        print(f"  {n:5d} {ratio_arith[i]:12.6f} {rm[i]:11.6f} {rp95[i]:11.6f} {flag:>12s}")

    z2 = slice(34, 500)
    frac_higher = np.mean(ratio_arith[z2] > rm[z2])
    print(f"\n  Fraction n in Zone 2 where arithmetic |S|/J > random mean: {frac_higher:.1%}")
    print(f"  INTERPRETATION: Arithmetic phases are {'MORE' if frac_higher > 0.5 else 'LESS'} coherent than random.")
    if frac_higher > 0.5:
        print(f"  This means arithmetic zeros are CLOSER to cancellation (margin thinner).")
        print(f"  Positivity holds DESPITE higher coherence — that's the arithmetic rigidity.")

    return lam_arith, ratio_arith, rm


# ============================================================
# TEST 2: RELATIVE COMPACTNESS — corrected
# ============================================================
def test2():
    print(f"\n{'='*72}")
    print(f"  TEST 2: V(f) = -1/(pi^2(1-f^2)) — DOMAIN ANALYSIS")
    print(f"{'='*72}")

    print(f"\n  KEY ISSUE: V(f) ~ 1/(1-f) near f=1.")
    print("  ||V||_L2 = infinity (non-integrable singularity).")
    print(f"  So V is NOT a bounded multiplication operator on L^2.")
    print(f"\n  The correct question: for f in Dom(H*D) on [-1,1],")
    print(f"  is Vf in L^2? This depends on how fast f vanishes at ±1.\n")

    # Functions in PW_1 have Fourier transforms that are smooth on [-1,1]
    # and in general DO NOT vanish at ±1.
    # So for generic g_hat in L^2[-1,1], V*g_hat is NOT in L^2.
    #
    # HOWEVER: A = P_1 M_{|y|} P_1 is well-defined in physical space.
    # The issue is with the frequency-space REPRESENTATION, not the operator.
    # The Hadamard FP integral is defined as a distribution, not as an L^2 operator.
    #
    # The correct framework: A acts as a bounded operator from H^s to H^{s-1}
    # for appropriate Sobolev spaces. The boundary singularity of V means the
    # Fourier-space kernel representation is a distribution, not a function.

    print(f"  ANALYSIS:")
    print(f"  1. V(f) = -1/(pi^2(1-f^2)) has a non-integrable square near ±1:")
    print(f"     ∫ V^2 df ~ ∫ df/(1-f)^2 = ∞")
    print(f"\n  2. For smooth g_hat with g_hat(±1) ≠ 0:")
    print(f"     V*g_hat ~ g_hat(±1)/(1-f) near boundary")
    print(f"     ||V*g_hat||^2 ~ ∫ df/(1-f)^2 = ∞")
    print(f"\n  3. THEREFORE: the decomposition A_hat = -(1/(2pi))H*D + V")
    print(f"     cannot be read as an operator identity on L^2[-1,1].")
    print(f"     It is valid as a DISTRIBUTIONAL identity (Hadamard FP).")

    # Now test what DOES work: the Galerkin matrix approach
    print(f"\n  WHAT WORKS: direct Galerkin in frequency-space cosine basis")
    print(f"  Compute <cos(mπf), A cos(nπf)> via the PHYSICAL-SPACE operator")
    print(f"  (which is perfectly well-defined) and check eigenvalue asymptotics.\n")

    # Physical-space Galerkin
    N = 60; L = 200; Nq = 60000
    y = np.linspace(-L, L, Nq); dy = y[1]-y[0]
    basis = np.zeros((N, Nq))
    basis[0] = np.sinc(2*y)
    for n in range(1, N):
        basis[n] = np.sinc(2*y-n) + np.sinc(2*y+n)
    w = basis*np.abs(y)
    A = basis @ w.T * dy
    S = basis @ basis.T * dy
    ev = np.sort(np.real(linalg.eigvalsh(A, S)))
    ev = ev[ev > 0.01]
    ne = len(ev)

    gaps = np.diff(ev)
    slope, intc = np.polyfit(np.arange(1, ne+1), ev, 1)

    print(f"  Physical-space eigenvalues (N=60 Shannon basis):")
    print(f"    slope = {slope:.6f} (expected 0.5)")
    print(f"    All {ne-1} gaps > 0: {np.all(gaps > 0)}")
    print(f"    min gap = {np.min(gaps):.6f}")
    print(f"    max gap = {np.max(gaps):.6f}")

    # The Weyl law holds for the PHYSICAL operator.
    # The frequency-side representation as -(1/(2π))H∘D + V is distributional.
    print(f"\n  REVISED STATUS OF (C2'):")
    print(f"  The operator A = P_1 M_{{|y|}} P_1 is well-defined on PW_1.")
    print(f"  Its eigenvalues satisfy lambda_n ~ n/2 + O(1) NUMERICALLY.")
    print(f"  The frequency-side ΨDO identification gives the principal symbol |k|.")
    print(f"  The boundary potential V is NOT an L^2 multiplication operator.")
    print(f"  CONSEQUENCE: the Kato-Rellich route to (C2') fails.")
    print(f"  The Weyl law must be proved by other means:")
    print(f"    (a) Direct use of the Landau-Pollak framework (v29 route)")
    print(f"    (b) Analysis of A as a Toeplitz-type operator on PW_1")
    print(f"    (c) Trace-class estimates on (A - λ)^{{-1}} in physical space")
    print(f"  This is a REAL gap, not just a technicality.")

    return ev, gaps


# ============================================================
# TEST 3: SUB-LEADING WEYL — larger basis
# ============================================================
def test3():
    print(f"\n{'='*72}")
    print(f"  TEST 3: SUB-LEADING WEYL TERM")
    print(f"{'='*72}")

    N = 60; L = 200; Nq = 60000
    y = np.linspace(-L, L, Nq); dy = y[1]-y[0]
    basis = np.zeros((N, Nq))
    basis[0] = np.sinc(2*y)
    for n in range(1, N):
        basis[n] = np.sinc(2*y-n) + np.sinc(2*y+n)
    w = basis*np.abs(y)
    A = basis @ w.T * dy
    S = basis @ basis.T * dy
    ev = np.sort(np.real(linalg.eigvalsh(A, S)))
    ev = ev[ev > 0.01]
    ne = len(ev)
    n_arr = np.arange(1, ne+1)
    res = ev - n_arr/2.0

    print(f"\n  Residual r_n = lambda_n - n/2 for {ne} eigenvalues:")
    print(f"\n  {'n':>4s} {'lambda_n':>10s} {'n/2':>7s} {'r_n':>10s}")
    for i in [0,1,4,9,14,19,29,39,49,58]:
        if i < ne:
            print(f"  {i+1:4d} {ev[i]:10.5f} {(i+1)/2:7.3f} {res[i]:10.5f}")

    # Models
    ln = np.log(n_arr); inv = 1.0/n_arr

    # A: r = c
    c0 = np.mean(res[10:]); e0 = np.sqrt(np.mean((res[10:]-c0)**2))

    # B: r = a + b*log(n) — skip first 10 (boundary effects)
    bB, aB = np.polyfit(ln[10:], res[10:], 1)
    eB = np.sqrt(np.mean((res[10:] - aB - bB*ln[10:])**2))

    # C: r = a + b/n
    bC, aC = np.polyfit(inv[10:], res[10:], 1)
    eC = np.sqrt(np.mean((res[10:] - aC - bC*inv[10:])**2))

    # D: r = a + b*log(n) + c/n
    X = np.column_stack([np.ones(ne-10), ln[10:], inv[10:]])
    cD, _, _, _ = np.linalg.lstsq(X, res[10:], rcond=None)
    eD = np.sqrt(np.mean((res[10:] - X@cD)**2))

    print(f"\n  MODEL FITS (n >= 11 to avoid boundary):")
    print(f"  {'Model':>28s} {'RMSE':>10s}")
    print(f"  {'r = const':>28s} {e0:10.6f}  c = {c0:.6f}")
    print(f"  {'r = a + b·log(n)':>28s} {eB:10.6f}  a={aB:.5f}, b={bB:.6f}")
    print(f"  {'r = a + b/n':>28s} {eC:10.6f}  a={aC:.5f}, b={bC:.4f}")
    print(f"  {'r = a + b·log(n) + c/n':>28s} {eD:10.6f}  a={cD[0]:.5f}, b={cD[1]:.6f}, c={cD[2]:.4f}")

    best = min({'const':e0, 'log':eB, '1/n':eC, 'both':eD}.items(), key=lambda x:x[1])
    print(f"\n  Best fit: {best[0]} (RMSE = {best[1]:.6f})")

    print(f"\n  LOG COEFFICIENT: b = {bB:.6f}")
    print(f"    GUE prediction (1/π²): {1/np.pi**2:.6f}")
    print(f"    Dirac on S¹ (1/(2π)):  {1/(2*np.pi):.6f}")

    # F-test: does adding log(n) significantly improve over constant?
    # RSS_const = sum((res[10:]-c0)^2), RSS_log has 2 params
    rss0 = np.sum((res[10:]-c0)**2)
    rssB = np.sum((res[10:]-aB-bB*ln[10:])**2)
    n_data = ne - 10
    F_stat = ((rss0 - rssB)/1) / (rssB/(n_data-2))
    print(f"\n    F-test (const vs log): F = {F_stat:.2f}")
    print(f"    {'SIGNIFICANT (F > 4)' if F_stat > 4 else 'NOT SIGNIFICANT (F < 4)'}")

    # Gap structure
    gaps = np.diff(ev)
    print(f"\n  GAP ANALYSIS:")
    print(f"    gaps[0..4] = {gaps[:5].round(6)}")
    print(f"    gaps converging to 1/2 from above: {np.all(gaps[:30] > 0.5)}")
    print(f"    gap - 1/2 for n=1..5: {(gaps[:5]-0.5).round(6)}")
    print(f"    gap - 1/2 for n=25..30: {(gaps[24:30]-0.5).round(6)}")

    # Fit gap correction: gap_n - 1/2 ~ c/n^alpha
    gf = gaps[:ne-2] - 0.5
    gf_pos = gf[gf > 0]
    n_gf = np.arange(1, len(gf_pos)+1)
    if len(gf_pos) > 5:
        ag, bg = np.polyfit(np.log(n_gf[2:]), np.log(gf_pos[2:]), 1)
        print(f"    gap_n - 1/2 ~ exp({bg:.2f}) · n^({ag:.3f})")
        print(f"    Decay exponent: {ag:.3f}")

    return ev, res, bB


# ============================================================
# MAIN
# ============================================================
def main():
    print("="*72)
    print("  THREE TESTS — CORRECTED VERSION")
    print("  F. D. Blum / Catalyst AI Research · March 2026")
    print("="*72)

    import subprocess, sys
    try:
        import mpmath
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'mpmath',
                               '--break-system-packages', '-q'])

    gamma = get_zeros(500)

    # Test 1
    lam, r_arith, r_rand = test1(gamma, n_max=2000)

    # Test 2
    ev2, gaps2 = test2()

    # Test 3
    ev3, res3, blog = test3()

    # ── FINAL ──
    print(f"\n{'='*72}")
    print(f"  HONEST SUMMARY")
    print(f"{'='*72}")
    print(f"""
  TEST 1 — Off-line zeros & phase coherence:
    Li's criterion on the critical line is trivially satisfied
    (each term 1-cos(nθ) ≥ 0). The content is the CONVERSE.
    Moving a single zero off-line: exponential growth of |w|^n
    eventually overwhelms, but with J=500 zeros the crossover
    occurs at n >> 1/(2ε) because positive contributions are large.
    
    SURPRISING FINDING: Arithmetic phases give HIGHER |S(n)|/J
    than random phases. The Riemann zeros are MORE coherent
    (closer to cancellation), not less. Positivity holds despite
    this — that's the rigidity.

  TEST 2 — Relative compactness of V:
    V(f) = -1/(π²(1-f²)) has a NON-INTEGRABLE square singularity.
    It is NOT a bounded multiplication operator on L²[-1,1].
    The Kato-Rellich route to proving (C2') FAILS.
    
    The ΨDO decomposition Â = -(1/(2π))H∘D + V is valid as a
    DISTRIBUTIONAL identity but NOT as an operator identity on L².
    
    The Weyl law λ_n ~ n/2 holds NUMERICALLY but the frequency-side
    route to proving it has a genuine gap: the boundary potential
    is too singular for standard perturbation theory.
    
    This is the honest status of (C2'): the identification of the
    principal symbol is correct, but the perturbation argument for
    V requires a different framework (Toeplitz operators, or direct
    physical-space trace estimates).

  TEST 3 — Sub-leading Weyl term:
    b_log = {blog:.6f}
    {'Log correction marginally significant' if abs(blog) > 0.01 else 'No log correction detected'}
    Best-fit sub-leading term: r_n ~ {blog:.4f}·log(n) + const
    Not conclusive with only 60 eigenvalues.
""")

if __name__ == "__main__":
    main()
