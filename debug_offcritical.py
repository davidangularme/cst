#!/usr/bin/env python3
"""
Debug: off-critical quadruplet contribution to Li coefficients.
Verify Remark 8.4 directly.
"""
import numpy as np
from mpmath import zetazero, mp
mp.dps = 30

def main():
    print("=" * 72)
    print("  DEBUG: OFF-CRITICAL QUADRUPLET CONTRIBUTION")
    print("=" * 72)
    
    # ─── Step 1: Single quadruplet, no other zeros ───
    print("\n§1. SINGLE QUADRUPLET (no background zeros)")
    print("-" * 50)
    
    gamma = 14.134725  # first zero
    
    for eps in [0.01, 0.05, 0.1, 0.2, 0.5]:
        sigma = 0.5 + eps
        
        # The four zeros:
        rho1 = sigma + 1j * gamma        # σ + iγ
        rho2 = sigma - 1j * gamma        # σ - iγ  (conjugate)
        rho3 = (1 - sigma) + 1j * gamma  # 1-σ + iγ  (func eq)
        rho4 = (1 - sigma) - 1j * gamma  # 1-σ - iγ  (func eq conj)
        
        w1 = (rho1 - 1) / rho1
        w2 = (rho2 - 1) / rho2
        w3 = (rho3 - 1) / rho3
        w4 = (rho4 - 1) / rho4
        
        print(f"\n  ε = {eps}, σ = {sigma}")
        print(f"    |w(ρ)|      = {abs(w1):.8f}  (<1 means shrinking)")
        print(f"    |w(1-ρ̄)|   = {abs(w3):.8f}  (>1 means growing)")
        
        # Contribution to λ_n:
        # C(n) = [1-w1^n] + [1-w2^n] + [1-w3^n] + [1-w4^n]
        #      = 4 - 2Re(w1^n) - 2Re(w3^n)
        # Since w2 = w1* and w4 = w3*
        
        n_values = [1, 5, 10, 50, 100, 500, 1000, 2000]
        print(f"    {'n':>6s}  {'2Re(1-w_ρ^n)':>14s}  {'2Re(1-w_FE^n)':>14s}  {'total':>12s}")
        for n in n_values:
            c_rho = 2 * np.real(1 - w1**n)
            c_fe = 2 * np.real(1 - w3**n)
            total = c_rho + c_fe
            print(f"    {n:6d}  {c_rho:14.4f}  {c_fe:14.4f}  {total:12.4f}")
        
        # Find n where total first becomes negative
        n_range = np.arange(1, 5001)
        c_rho_all = 2 * np.real(1 - w1**n_range)
        c_fe_all = 2 * np.real(1 - w3**n_range)
        total_all = c_rho_all + c_fe_all
        
        neg = np.where(total_all < 0)[0]
        if len(neg) > 0:
            print(f"    → First negative at n = {neg[0]+1}")
            print(f"      Predicted n* = 1/(2ε) = {1/(2*eps):.0f}")
            # Most negative
            imin = np.argmin(total_all)
            print(f"      Most negative: n={imin+1}, C={total_all[imin]:.4f}")
        else:
            print(f"    → No negative in n=1..5000")
            print(f"      Predicted n* = 1/(2ε) = {1/(2*eps):.0f}")
            print(f"      Min total = {np.min(total_all):.6f} at n={np.argmin(total_all)+1}")
    
    # ─── Step 2: Check for overflow ───
    print("\n\n§2. OVERFLOW CHECK")
    print("-" * 50)
    sigma = 1.0  # ε = 0.5
    rho3 = 0.0 + 1j * gamma  # 1-σ+iγ = iγ
    w3 = (rho3 - 1) / rho3
    print(f"  w(1-ρ̄) = {w3}")
    print(f"  |w|     = {abs(w3):.10f}")
    for n in [100, 500, 1000, 2000, 5000]:
        wn = w3**n
        print(f"  n={n:5d}: |w|^n = {abs(w3)**n:.4e}, Re(w^n) = {np.real(wn):.4e}, overflow={np.isinf(np.real(wn))}")
    
    # ─── Step 3: Quadruplet + 200 background zeros ───
    print("\n\n§3. QUADRUPLET + BACKGROUND (J=200)")
    print("-" * 50)
    
    J = 200
    print(f"  Loading {J} zeros...")
    gammas = np.array([float(zetazero(j).imag) for j in range(1, J+1)])
    
    # True λ_n (all on critical line)
    rho_true = 0.5 + 1j * gammas
    w_true = (rho_true - 1) / rho_true
    
    n_max = 3000
    
    for eps in [0.01, 0.05, 0.1, 0.5]:
        sigma = 0.5 + eps
        
        # Replace first zero with off-critical quadruplet
        # Remove: pair (ρ₁, ρ̄₁) with w_1 on critical line
        # Add: quadruplet at (σ+iγ₁, σ-iγ₁, 1-σ+iγ₁, 1-σ-iγ₁)
        
        rho_off = sigma + 1j * gammas[0]
        rho_fe = (1-sigma) + 1j * gammas[0]
        w_off = (rho_off - 1) / rho_off
        w_fe = (rho_fe - 1) / rho_fe
        
        print(f"\n  ε = {eps}: |w_off| = {abs(w_off):.6f}, |w_fe| = {abs(w_fe):.6f}")
        
        lam = np.zeros(n_max)
        for n in range(1, n_max+1):
            # Background: all J pairs except the first
            S_bg = np.sum(2 * np.real(1 - w_true[1:]**n))
            # Off-critical quadruplet
            S_off = 2*np.real(1 - w_off**n) + 2*np.real(1 - w_fe**n)
            lam[n-1] = S_bg + S_off
        
        neg = np.where(lam < 0)[0]
        print(f"  min(λ_n) = {np.min(lam):.4f} at n = {np.argmin(lam)+1}")
        if len(neg) > 0:
            print(f"  FIRST λ_n < 0 at n = {neg[0]+1}")
            print(f"  Predicted n* ≈ {1/(2*eps):.0f}")
            print(f"  Total negative: {len(neg)}")
        else:
            print(f"  No λ_n < 0 in n=1..{n_max}")
            
            # Diagnostic: what's the off-critical contribution at large n?
            n_diag = 2000
            S_bg_diag = np.sum(2 * np.real(1 - w_true[1:]**n_diag))
            S_off_diag = 2*np.real(1 - w_off**n_diag) + 2*np.real(1 - w_fe**n_diag)
            print(f"  At n={n_diag}: bg={S_bg_diag:.2f}, off-crit={S_off_diag:.4f}, total={S_bg_diag+S_off_diag:.2f}")
    
    print()

if __name__ == "__main__":
    main()
