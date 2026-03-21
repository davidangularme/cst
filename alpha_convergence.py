#!/usr/bin/env python3
"""
ALPHA(N) CONVERGENCE TEST
=========================
Compute α(N) = log₁₀(μ₁/μ₀) / N for increasing N.
If α(N) → α∞ > 0, the exponential dominance lemma holds.

200 digits, 500 zeros, λ=3 and λ=5
"""
import mpmath
import time

DPS = 200
mpmath.mp.dps = DPS

print("=" * 78)
print(f"  ALPHA(N) CONVERGENCE TEST — {DPS} digits")
print("=" * 78)
print()

# ═══════════════════════════════════════════════════════════════
# ZEROS
# ═══════════════════════════════════════════════════════════════
NUM_ZEROS = 500
print(f"  Computing {NUM_ZEROS} zeros at {DPS} digits...")
t0 = time.time()
gammas = []
for n in range(1, NUM_ZEROS + 1):
    g = mpmath.im(mpmath.zetazero(n))
    gammas.append(g)
    if n % 100 == 0:
        print(f"    {n}/{NUM_ZEROS} ({time.time()-t0:.0f}s)")
print(f"    Done in {time.time()-t0:.0f}s")

# ═══════════════════════════════════════════════════════════════
# WEIL MATRIX
# ═══════════════════════════════════════════════════════════════
def build_weil(lam, N, gammas):
    L = mpmath.log(lam)
    size = N + 1
    Phi = mpmath.matrix(size, len(gammas))
    for gi in range(len(gammas)):
        gamma = gammas[gi]
        Phi[0, gi] = 2 * mpmath.sin(gamma * L) / gamma
        for k in range(1, size):
            km = mpmath.mpf(k) - gamma
            kp = mpmath.mpf(k) + gamma
            p1 = L if mpmath.fabs(km) < mpmath.power(10, -DPS+10) else mpmath.sin(km*L)/km
            p2 = L if mpmath.fabs(kp) < mpmath.power(10, -DPS+10) else mpmath.sin(kp*L)/kp
            Phi[k, gi] = p1 + p2
    return 2 * Phi * Phi.T

# ═══════════════════════════════════════════════════════════════
# MAIN SCAN
# ═══════════════════════════════════════════════════════════════
for LAM in [3, 5]:
    print()
    print("=" * 78)
    print(f"  λ = {LAM}")
    print("=" * 78)
    
    N_values = [5, 10, 15, 20, 25, 30, 35, 40, 45]
    
    print()
    print(f"  {'N':>4s}  {'log₁₀|μ₀|':>12s}  {'log₁₀|μ₁|':>12s}  "
          f"{'log₁₀(μ₁/μ₀)':>14s}  {'α(N)=log/N':>12s}  {'time':>6s}")
    print(f"  {'─'*4}  {'─'*12}  {'─'*12}  {'─'*14}  {'─'*12}  {'─'*6}")
    
    alphas = []
    Ns_done = []
    
    for N in N_values:
        size = N + 1
        # Check if eigenvalues will be resolvable
        # Rough estimate: μ₀ ~ 0.007^N → log₁₀|μ₀| ~ -2.1·N
        est_log = -2.1 * N
        if abs(est_log) > DPS - 20:
            print(f"  {N:4d}  {'SKIP — need more digits':>50s}")
            continue
        
        t1 = time.time()
        M = build_weil(LAM, N, gammas)
        eigvals_mp, _ = mpmath.eigsy(M)
        elapsed = time.time() - t1
        
        eigvals = sorted([eigvals_mp[i] for i in range(len(eigvals_mp))])
        mu0 = eigvals[0]
        mu1 = eigvals[1]
        
        # Handle potential negative eigenvalues (numerical noise)
        abs_mu0 = mpmath.fabs(mu0)
        abs_mu1 = mpmath.fabs(mu1)
        
        if abs_mu0 > 0 and abs_mu1 > 0:
            log_mu0 = float(mpmath.log(abs_mu0, 10))
            log_mu1 = float(mpmath.log(abs_mu1, 10))
            log_ratio = log_mu1 - log_mu0
            alpha = log_ratio / N
            
            alphas.append(alpha)
            Ns_done.append(N)
            
            print(f"  {N:4d}  {log_mu0:12.2f}  {log_mu1:12.2f}  "
                  f"{log_ratio:14.2f}  {alpha:12.4f}  {elapsed:6.1f}s")
        else:
            print(f"  {N:4d}  {'ZERO eigenvalue':>50s}  {elapsed:6.1f}s")
    
    # ─── Alpha convergence analysis ───────────────────
    print()
    print(f"  α(N) sequence for λ={LAM}:")
    for N_val, a_val in zip(Ns_done, alphas):
        bar = '█' * int(max(0, a_val * 100))
        print(f"    N={N_val:3d}:  α = {a_val:.4f}  {bar}")
    
    if len(alphas) >= 3:
        # Check convergence: is the last third stable?
        n3 = len(alphas) // 3
        last_third = alphas[-n3:] if n3 > 0 else alphas[-2:]
        mean_last = sum(last_third) / len(last_third)
        spread = max(last_third) - min(last_third)
        
        print()
        print(f"  Last {len(last_third)} values: mean α = {mean_last:.4f}, "
              f"spread = {spread:.4f}")
        
        if mean_last > 0.05 and spread < 0.1:
            print(f"  → α(N) CONVERGES to α∞ ≈ {mean_last:.3f} > 0  ✓")
            print(f"  → Exponential dominance lemma is SUPPORTED")
            print(f"  → μ₁/μ₀ ~ 10^({mean_last:.3f}·N) → ∞")
        elif mean_last > 0:
            print(f"  → α(N) is POSITIVE but not yet converged (spread = {spread:.3f})")
            print(f"  → Need larger N to confirm")
        else:
            print(f"  → α(N) ≤ 0 — exponential dominance FAILS ✗")

# ═══════════════════════════════════════════════════════════════
# FINAL
# ═══════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("  CONCLUSION")
print("=" * 78)
print(f"""
  The exponential dominance lemma states:
  
    If  lim inf_{{N→∞}} (1/N) · log₁₀(μ₁/μ₀) > 0,
    then ξ_λ → c_λ · k_λ exponentially in N.
    
  This is the first concrete, provable step toward Assumption A.
  
  Total time: {time.time()-t0:.0f}s
""")
