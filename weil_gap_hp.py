#!/usr/bin/env python3
"""
HIGH-PRECISION SPECTRAL GAP — WEIL QUADRATIC FORM
Uses mpmath (50 decimal digits) to resolve the smallest eigenvalues.
"""
import mpmath
import time

mpmath.mp.dps = 50  # 50 decimal digits

def compute_zeros(num_zeros):
    """Compute imaginary parts of first num_zeros Riemann zeros."""
    print(f"  Computing {num_zeros} zeros at {mpmath.mp.dps} digits...")
    t0 = time.time()
    gammas = []
    for n in range(1, num_zeros + 1):
        g = mpmath.im(mpmath.zetazero(n))
        gammas.append(g)
        if n % 50 == 0:
            elapsed = time.time() - t0
            print(f"    {n}/{num_zeros} ({elapsed:.1f}s)")
    print(f"    Done in {time.time()-t0:.1f}s")
    return gammas


def even_weil_matrix_mp(lam, N, gammas):
    """
    Build even Weil quadratic form matrix at mpmath precision.
    
    Basis: cos(k·t), k = 0,...,N on [-L, L], L = log(λ)
    M_{jk} = 2 Σ_{γ>0} Φ_j(γ) Φ_k(γ)
    """
    L = mpmath.log(lam)
    size = N + 1
    num_g = len(gammas)
    
    # Build Phi matrix
    Phi = mpmath.matrix(size, num_g)
    for gi in range(num_g):
        gamma = gammas[gi]
        # k = 0
        Phi[0, gi] = 2 * mpmath.sin(gamma * L) / gamma
        # k >= 1
        for k in range(1, size):
            km = mpmath.mpf(k) - gamma
            kp = mpmath.mpf(k) + gamma
            if mpmath.fabs(km) < mpmath.mpf('1e-40'):
                p1 = L
            else:
                p1 = mpmath.sin(km * L) / km
            if mpmath.fabs(kp) < mpmath.mpf('1e-40'):
                p2 = L
            else:
                p2 = mpmath.sin(kp * L) / kp
            Phi[k, gi] = p1 + p2
    
    # M = 2 * Phi * Phi^T
    PhiT = Phi.T
    M = 2 * Phi * PhiT
    
    return M


def main():
    print("=" * 78)
    print("  HIGH-PRECISION SPECTRAL GAP (50 decimal digits)")
    print("  Weil Quadratic Form — Even Sector")
    print("=" * 78)
    print()
    
    # Compute zeros
    num_zeros = 200
    gammas = compute_zeros(num_zeros)
    
    # ─── Main scan ────────────────────────────────────
    print()
    print("=" * 78)
    print("  SPECTRAL GAP SCAN")
    print("=" * 78)
    
    configs = [
        # (λ, N)
        (3, 10),
        (3, 15),
        (3, 25),
        (5, 10),
        (5, 15),
        (5, 25),
        (10, 10),
        (10, 15),
        (10, 25),
        (20, 10),
        (20, 15),
        (50, 10),
        (50, 15),
        (100, 10),
    ]
    
    results = []
    
    for lam, N in configs:
        t0 = time.time()
        print(f"\n  λ={lam}, N={N} (matrix {N+1}×{N+1})...")
        
        M = even_weil_matrix_mp(lam, N, gammas)
        
        # Eigenvalues via mpmath.eigsy (symmetric solver)
        t1 = time.time()
        print(f"    Matrix built in {t1-t0:.1f}s. Computing eigenvalues...")
        
        eigvals_mp, eigvecs_mp = mpmath.eigsy(M)
        
        # Convert to sorted list
        eigvals = sorted([eigvals_mp[i] for i in range(len(eigvals_mp))])
        
        elapsed = time.time() - t0
        
        mu0 = eigvals[0]
        mu1 = eigvals[1]
        mu2 = eigvals[2] if len(eigvals) > 2 else mpmath.mpf(0)
        gap = mu1 - mu0
        tr = sum(M[i,i] for i in range(N+1))
        mu_max = eigvals[-1]
        
        # Ratio Δ/μ₁
        if mpmath.fabs(mu1) > mpmath.mpf('1e-100'):
            ratio = gap / mu1
        else:
            ratio = mpmath.mpf('inf')
        
        results.append({
            'lam': lam, 'N': N, 'mu0': mu0, 'mu1': mu1, 'mu2': mu2,
            'gap': gap, 'ratio': ratio, 'tr': tr, 'mu_max': mu_max,
            'eigvals': eigvals, 'time': elapsed
        })
        
        print(f"    Done in {elapsed:.1f}s")
        print(f"    trace = {mpmath.nstr(tr, 8)}")
        print(f"    μ₀    = {mpmath.nstr(mu0, 15)}")
        print(f"    μ₁    = {mpmath.nstr(mu1, 15)}")
        print(f"    μ₂    = {mpmath.nstr(mu2, 15)}")
        print(f"    Δ     = {mpmath.nstr(gap, 15)}")
        print(f"    Δ/μ₁  = {mpmath.nstr(ratio, 8)}")
        print(f"    μ_max = {mpmath.nstr(mu_max, 8)}")
    
    # ─── Summary table ────────────────────────────────
    print()
    print("=" * 78)
    print("  SUMMARY TABLE")
    print("=" * 78)
    print()
    print(f"  {'λ':>4s}  {'N':>3s}  {'μ₀':>20s}  {'μ₁':>20s}  "
          f"{'Δ=μ₁−μ₀':>20s}  {'Δ/μ₁':>10s}")
    print(f"  {'─'*4}  {'─'*3}  {'─'*20}  {'─'*20}  {'─'*20}  {'─'*10}")
    
    for r in results:
        print(f"  {r['lam']:4d}  {r['N']:3d}  {mpmath.nstr(r['mu0'],12):>20s}  "
              f"{mpmath.nstr(r['mu1'],12):>20s}  "
              f"{mpmath.nstr(r['gap'],12):>20s}  "
              f"{mpmath.nstr(r['ratio'],6):>10s}")
    
    # ─── Scaling analysis ─────────────────────────────
    print()
    print("=" * 78)
    print("  SCALING ANALYSIS: How does Δ depend on λ?")
    print("=" * 78)
    
    for N in [10, 15]:
        print(f"\n  N = {N}:")
        subset = [(r['lam'], float(r['gap']), float(r['mu0']), float(r['mu1'])) 
                  for r in results if r['N'] == N and float(r['gap']) > 0]
        if len(subset) >= 2:
            import numpy as np
            lams = np.array([s[0] for s in subset])
            gaps = np.array([s[1] for s in subset])
            mu0s = np.array([abs(s[2]) for s in subset])
            mu1s = np.array([s[3] for s in subset])
            
            # Gap scaling
            if np.all(gaps > 0):
                alpha = np.polyfit(np.log(lams), np.log(gaps), 1)[0]
                print(f"    Δ(λ) ~ λ^{alpha:.3f}")
            
            # μ₀ scaling
            if np.all(mu0s > 0):
                beta = np.polyfit(np.log(lams), np.log(mu0s), 1)[0]
                print(f"    |μ₀(λ)| ~ λ^{beta:.3f}")
            
            # μ₁ scaling
            if np.all(mu1s > 0):
                gamma_sc = np.polyfit(np.log(lams), np.log(mu1s), 1)[0]
                print(f"    μ₁(λ) ~ λ^{gamma_sc:.3f}")
    
    # ─── Eigenvalue spectrum at λ=3, N=25 ─────────────
    print()
    print("=" * 78)
    print("  EIGENVALUE SPECTRUM: λ=3, N=25")
    print("=" * 78)
    
    r325 = [r for r in results if r['lam'] == 3 and r['N'] == 25]
    if r325:
        eigvals = r325[0]['eigvals']
        print(f"\n  {'k':>3s}  {'μ_k':>25s}  {'log₁₀|μ_k|':>12s}")
        print(f"  {'─'*3}  {'─'*25}  {'─'*12}")
        for k in range(min(15, len(eigvals))):
            v = eigvals[k]
            if v != 0:
                log10 = float(mpmath.log(mpmath.fabs(v), 10))
            else:
                log10 = float('-inf')
            print(f"  {k:3d}  {mpmath.nstr(v,18):>25s}  {log10:12.2f}")
    
    # ─── Verdict ──────────────────────────────────────
    print()
    print("=" * 78)
    print("  VERDICT")
    print("=" * 78)
    
    print("""
  The spectral gap Δ(λ) = μ₁ − μ₀ measures the isolation of the
  Weil minimiser ξ_λ from the next eigenfunction.
  
  For Path C (variational contraction):
    - Need Δ(λ) > 0 for all λ (gap stays open)
    - Ideally Δ(λ) ≥ c · λ^{-α} with α < ∞ (polynomial decay OK)
    - If Δ(λ) → 0 faster than μ₀(λ) → 0, the PL inequality fails
  
  For Path B (Fisher information):
    - The gap is related to I(s) = −ζ''/ζ + (ζ'/ζ)² on Re(s) = 1/2
    - Δ > 0 ↔ I(s) finite ↔ Bures curvature bounded
    """)
    
    print("  Computation complete.")


if __name__ == '__main__':
    main()
