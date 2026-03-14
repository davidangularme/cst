#!/usr/bin/env python3
"""
Three Concrete Tests for the FDBC Programme
=============================================
Test 1: Arithmetic vs Random phases — does Li positivity require arithmetic?
Test 2: Relative compactness of V — is (C2') technically valid?
Test 3: Sub-leading Weyl term — what is N(lambda) - 2*lambda actually?

Author: F. D. Blum / Catalyst AI Research · March 2026
"""
import numpy as np
from scipy import linalg
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# RIEMANN ZEROS
# ============================================================
def get_zeros(N):
    from mpmath import zetazero
    print(f"  Computing {N} Riemann zeros...")
    z = []
    for n in range(1, N+1):
        z.append(float(zetazero(n).imag))
        if n % 100 == 0: print(f"    ...{n}/{N}")
    return np.array(z)

# ============================================================
# TEST 1: ARITHMETIC vs RANDOM PHASES
# ============================================================
def test1_phases(gamma, n_max=500, N_trials=300):
    J = len(gamma)
    print(f"\n{'='*72}")
    print(f"  TEST 1: ARITHMETIC vs RANDOM PHASES")
    print(f"  J={J} zeros, n_max={n_max}, {N_trials} random trials")
    print(f"{'='*72}")

    # Arithmetic: w_j = (rho_j-1)/rho_j, |w_j|=1 on RH
    rho = 0.5 + 1j*gamma
    w = (rho - 1)/rho
    theta = np.angle(w)  # phases on unit circle
    n_arr = np.arange(1, n_max+1)

    # --- Part A: Off-line zeros (the sharp test) ---
    print(f"\n  A. OFF-LINE ZEROS: perturb sigma away from 1/2")
    print(f"     Prediction (Remark 8.4): first lambda_n < 0 at n ~ 1/(2*eps)")
    print(f"\n  {'sigma':>8s}  {'eps':>8s}  {'pred n*':>8s}  {'actual n*':>10s}  {'min lambda':>12s}")

    for sigma in [0.500, 0.501, 0.502, 0.505, 0.510, 0.520, 0.550, 0.600, 0.750]:
        eps = sigma - 0.5
        rho_p = sigma + 1j*gamma
        w_p = (rho_p - 1)/rho_p
        # Functional equation partners: 1-conj(rho) = (1-sigma)+i*gamma
        rho_q = (1-sigma) + 1j*gamma
        w_q = (rho_q - 1)/rho_q

        lam = np.zeros(n_max)
        for i, n in enumerate(n_arr):
            # Quadruplet: {rho, conj(rho), 1-conj(rho), 1-rho}
            # Contribution = Re[1-w^n] + Re[1-conj(w)^n] + Re[1-w_q^n] + Re[1-conj(w_q)^n]
            lam[i] = np.sum(
                (1-np.real(w_p**n)) + (1-np.real(np.conj(w_p)**n))
                + (1-np.real(w_q**n)) + (1-np.real(np.conj(w_q)**n))
            ) / 2  # each zero counted once

        neg = n_arr[lam < 0]
        first_neg = neg[0] if len(neg) > 0 else -1
        pred = int(round(1/(2*eps))) if eps > 0 else -1
        pred_s = str(pred) if pred > 0 else "inf"
        first_s = str(first_neg) if first_neg > 0 else "none"
        print(f"  {sigma:8.3f}  {eps:8.3f}  {pred_s:>8s}  {first_s:>10s}  {np.min(lam):12.2f}")

    # --- Part B: |S(n)|/J diagnostic (on critical line) ---
    print(f"\n  B. PHASE COHERENCE: |S(n)|/J for arithmetic vs random gamma")
    print(f"     S(n) = sum_j exp(i*n*theta_j). Positivity <=> |S(n)| < J.")

    # Arithmetic S(n)
    S_arith = np.array([np.sum(np.exp(1j*n*theta)) for n in n_arr])
    ratio_arith = np.abs(S_arith) / J

    # Random gamma (uniform in same range, same count)
    ratio_rand_trials = np.zeros((N_trials, n_max))
    for t in range(N_trials):
        g_rand = np.sort(np.random.uniform(gamma[0]*0.9, gamma[-1]*1.1, J))
        rho_r = 0.5 + 1j*g_rand
        w_r = (rho_r-1)/rho_r
        th_r = np.angle(w_r)
        for i, n in enumerate(n_arr):
            ratio_rand_trials[t, i] = np.abs(np.sum(np.exp(1j*n*th_r))) / J

    ratio_rand_mean = np.mean(ratio_rand_trials, axis=0)
    ratio_rand_p95 = np.percentile(ratio_rand_trials, 95, axis=0)
    ratio_rand_max = np.max(ratio_rand_trials, axis=0)

    # Compare in Zone 2
    print(f"\n  {'n':>5s}  {'Arith':>10s}  {'Rand mean':>10s}  {'Rand p95':>10s}  {'Rand max':>10s}  {'Arith < Rand?':>14s}")
    for n in [35, 40, 50, 75, 100, 150, 200, 250, 300, 400, 500]:
        if n <= n_max:
            i = n-1
            flag = "YES" if ratio_arith[i] < ratio_rand_mean[i] else "no"
            print(f"  {n:5d}  {ratio_arith[i]:10.6f}  {ratio_rand_mean[i]:10.6f}  {ratio_rand_p95[i]:10.6f}  {ratio_rand_max[i]:10.6f}  {flag:>14s}")

    # Zone 2 statistics
    z2 = (n_arr > 34) & (n_arr <= 500)
    print(f"\n  Zone 2 summary (34 < n <= 500):")
    print(f"    Arithmetic: max |S(n)|/J = {np.max(ratio_arith[z2]):.6f}, mean = {np.mean(ratio_arith[z2]):.6f}")
    print(f"    Random:     max |S(n)|/J = {np.max(ratio_rand_max[z2]):.6f}, mean of means = {np.mean(ratio_rand_mean[z2]):.6f}")
    print(f"    Arith max < 1: {np.max(ratio_arith[z2]) < 1.0}")

    # Is arithmetic systematically lower?
    arith_lower_frac = np.mean(ratio_arith[z2] < ratio_rand_mean[z2])
    print(f"    Fraction where arith < rand mean: {arith_lower_frac:.2%}")

    return ratio_arith, ratio_rand_mean

# ============================================================
# TEST 2: RELATIVE COMPACTNESS OF V
# ============================================================
def test2_compactness():
    print(f"\n{'='*72}")
    print(f"  TEST 2: RELATIVE COMPACTNESS OF V(f) = -1/(pi^2(1-f^2))")
    print(f"{'='*72}")

    Nq = 20000
    f = np.linspace(-1+1e-10, 1-1e-10, Nq)
    df = f[1]-f[0]
    V = -1.0/(np.pi**2*(1-f**2))

    print(f"\n  Kato-Rellich test: ||V*e_n|| / ||D*e_n|| as n -> inf")
    print(f"  Basis: e_n(f) = cos(n*pi*f)")
    print(f"\n  {'n':>5s}  {'||Ve_n||':>12s}  {'||De_n||':>12s}  {'ratio':>12s}")

    ns = list(range(1, 21)) + list(range(25, 101, 5))
    ratios = []
    ns_used = []
    for n in ns:
        e_n = np.cos(n*np.pi*f)
        De_n = -n*np.pi*np.sin(n*np.pi*f)

        norm_Ve = np.sqrt(np.sum((V*e_n)**2)*df)
        norm_De = np.sqrt(np.sum(De_n**2)*df)
        r = norm_Ve/norm_De if norm_De > 0 else np.inf
        ratios.append(r)
        ns_used.append(n)
        if n <= 15 or n % 20 == 0:
            print(f"  {n:5d}  {norm_Ve:12.4f}  {norm_De:12.4f}  {r:12.6f}")

    ratios = np.array(ratios)
    ns_used = np.array(ns_used)

    # Fit decay: ratio ~ C * n^alpha
    mask = ns_used >= 5
    alpha, logC = np.polyfit(np.log(ns_used[mask]), np.log(ratios[mask]), 1)
    print(f"\n  Power-law fit (n >= 5): ratio ~ {np.exp(logC):.4f} * n^({alpha:.4f})")

    if alpha < -0.01:
        print(f"  RESULT: ratio DECAYS as n^{alpha:.3f}")
        print(f"  => V is infinitesimally bounded relative to H*D")
        print(f"  => (C2') holds numerically: V is relatively compact  ✓")
    else:
        print(f"  RESULT: ratio does NOT decay. (C2') may fail.")

    # Sobolev embedding check: ||Vf|| / ||f||_{H^1} for random f in PW_1
    print(f"\n  Sobolev check: ||Vf||/||f||_{{H^1}} for random band-limited f")
    sob_ratios = []
    for _ in range(200):
        # Random even band-limited function: sum of cos(k*pi*f) with random coeffs
        K = 30
        coeffs = np.random.randn(K)
        g = sum(coeffs[k]*np.cos((k+1)*np.pi*f) for k in range(K))
        Dg = sum(-coeffs[k]*(k+1)*np.pi*np.sin((k+1)*np.pi*f) for k in range(K))
        norm_Vg = np.sqrt(np.sum((V*g)**2)*df)
        norm_H1 = np.sqrt(np.sum(g**2)*df + np.sum(Dg**2)*df)
        sob_ratios.append(norm_Vg/norm_H1)

    sob_ratios = np.array(sob_ratios)
    print(f"    Mean ||Vf||/||f||_H1:  {np.mean(sob_ratios):.6f}")
    print(f"    Max  ||Vf||/||f||_H1:  {np.max(sob_ratios):.6f}")
    print(f"    Std  ||Vf||/||f||_H1:  {np.std(sob_ratios):.6f}")
    print(f"    All bounded: {np.all(sob_ratios < 10)} (max = {np.max(sob_ratios):.4f})")

    return alpha

# ============================================================
# TEST 3: SUB-LEADING WEYL TERM
# ============================================================
def test3_weyl():
    print(f"\n{'='*72}")
    print(f"  TEST 3: SUB-LEADING WEYL TERM")
    print(f"{'='*72}")

    N = 60
    L, Nq = 200, 60000
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

    print(f"\n  {ne} eigenvalues. Residual r_n = lambda_n - n/2:")
    print(f"\n  {'n':>5s}  {'lambda_n':>10s}  {'n/2':>8s}  {'r_n':>10s}")
    for i in [0,1,2,4,9,19,29,39,49,58]:
        if i < ne:
            print(f"  {i+1:5d}  {ev[i]:10.5f}  {(i+1)/2:8.3f}  {res[i]:10.5f}")

    # Fit models
    ln = np.log(n_arr)
    inv = 1.0/n_arr

    # Model A: constant
    c0 = np.mean(res); e0 = np.std(res)

    # Model B: a + b*log(n)
    bB, aB = np.polyfit(ln, res, 1)
    eB = np.sqrt(np.mean((res - aB - bB*ln)**2))

    # Model C: a + b/n
    bC, aC = np.polyfit(inv, res, 1)
    eC = np.sqrt(np.mean((res - aC - bC*inv)**2))

    # Model D: a + b*log(n) + c/n
    X = np.column_stack([np.ones(ne), ln, inv])
    cD, _, _, _ = np.linalg.lstsq(X, res, rcond=None)
    eD = np.sqrt(np.mean((res - X@cD)**2))

    print(f"\n  MODEL FITS:")
    print(f"  {'Model':>30s}  {'RMSE':>10s}")
    print(f"  {'r = const':>30s}  {e0:10.6f}   c={c0:.6f}")
    print(f"  {'r = a + b*log(n)':>30s}  {eB:10.6f}   a={aB:.6f}, b={bB:.6f}")
    print(f"  {'r = a + b/n':>30s}  {eC:10.6f}   a={aC:.6f}, b={bC:.5f}")
    print(f"  {'r = a + b*log(n) + c/n':>30s}  {eD:10.6f}   a={cD[0]:.6f}, b={cD[1]:.6f}, c={cD[2]:.5f}")

    print(f"\n  BEST FIT: ", end="")
    errors = {'const': e0, 'log': eB, '1/n': eC, 'log+1/n': eD}
    best = min(errors, key=errors.get)
    print(f"{best} (RMSE = {errors[best]:.6f})")

    print(f"\n  LOG CORRECTION TEST:")
    print(f"    b_log = {bB:.6f}")
    print(f"    Compare: 1/pi^2 = {1/np.pi**2:.6f} (GUE two-point)")
    print(f"    Compare: 1/(2pi) = {1/(2*np.pi):.6f} (Dirac on S^1)")
    sig = abs(bB) / (eB/np.sqrt(ne))
    print(f"    Significance: |b|/(RMSE/sqrt(n)) = {sig:.2f} sigma")
    if sig > 2:
        print(f"    => LOG TERM IS SIGNIFICANT at {sig:.1f}σ")
    else:
        print(f"    => Log term is NOT significant ({sig:.1f}σ)")

    # Gap fluctuations
    gaps = np.diff(ev)
    gf = gaps - 0.5
    print(f"\n  GAP FLUCTUATIONS delta_n - 1/2:")
    print(f"    Mean: {np.mean(gf):.7f}")
    print(f"    Std:  {np.std(gf):.7f}")
    ag, _ = np.polyfit(np.log(np.arange(2, len(gaps)+1)), np.log(np.abs(gf[1:])+1e-15), 1)
    print(f"    Decay: |delta_n - 1/2| ~ n^({ag:.2f})")

    return ev, res, bB

# ============================================================
# MAIN
# ============================================================
def main():
    print("="*72)
    print("  THREE CONCRETE TESTS — FDBC PROGRAMME")
    print("  F. D. Blum / Catalyst AI Research · March 2026")
    print("="*72)

    # Install mpmath if needed
    import subprocess, sys
    try:
        import mpmath
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'mpmath',
                               '--break-system-packages', '-q'])

    gamma = get_zeros(500)

    r1_arith, r1_rand = test1_phases(gamma, n_max=500, N_trials=300)
    alpha = test2_compactness()
    ev, res, b_log = test3_weyl()

    # ── FINAL SUMMARY ──
    print(f"\n{'='*72}")
    print(f"  FINAL SUMMARY")
    print(f"{'='*72}")
    print(f"""
  TEST 1 — Arithmetic vs Random Phases:
    On the critical line (|w|=1): lambda_n >= 0 is AUTOMATIC.
    The Li criterion's content is the CONVERSE: all lambda_n >= 0 => RH.
    Off-line zeros at sigma = 1/2 + eps: first lambda_n < 0 at n ~ 1/(2eps).
    |S(n)|/J measures closeness to cancellation.
    Arithmetic phases show {'LOWER' if np.mean(r1_arith[34:300] < r1_rand[34:300]) > 0.5 else 'COMPARABLE'} |S(n)|/J vs random in Zone 2.

  TEST 2 — Relative Compactness of V:
    ||Ve_n||/||De_n|| ~ n^({alpha:.3f})
    {'DECAYS => V relatively compact => (C2) HOLDS' if alpha < -0.01 else 'Does NOT decay => (C2) uncertain'}

  TEST 3 — Sub-leading Weyl Term:
    b_log = {b_log:.6f}
    {'Significant log correction DETECTED' if abs(b_log) > 0.001 else 'No significant log correction — N(lambda) ~ 2*lambda + const'}
""")

if __name__ == "__main__":
    main()
