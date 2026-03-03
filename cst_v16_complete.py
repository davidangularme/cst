#!/usr/bin/env python3
"""
CST v16 — Complete Computation
Reproduces all results and generates all figures from the paper.

Usage:
    python3 cst_v16_complete.py

Outputs:
    results/fig1_potential.png  — V_CST, cost factor, closed-form R
    results/fig2_signal.png     — R(t), ΔR signature, R vs s
    results/fig3_flows.png      — Bloch ball trajectories, purity, V(t)
    results/fig4_toronto.png    — Rb87 costs, subsidy, predictions

Requirements:
    numpy, scipy, matplotlib

Author: Frédéric David Blum
With: Claude (Anthropic) & Catalyst AIS
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad, solve_ivp
import os, sys

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(OUTDIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════
# CORE FRAMEWORK
# ═══════════════════════════════════════════════════════════

def V_CST(rho):
    """V(ρ) = S(ρ || I/2) = -S_vN + log2. Unique by Umegaki (1962)."""
    rho = np.clip(rho, 1e-12, 1 - 1e-12)
    lp, lm = (1 + rho) / 2, (1 - rho) / 2
    return lp * np.log(2 * lp) + lm * np.log(2 * lm)

def V_CST_N(pops):
    """V = S(ρ || I/N) for N-level system with given populations."""
    N = len(pops)
    return sum(p * np.log(N * p) for p in pops if p > 1e-15)

def R_ratio(rho1, rho2, s=1.0):
    """R = t_CST / t_QM for radial transit. Pure number, no dynamics."""
    a, b = min(rho1, rho2), max(rho1, rho2)
    if abs(b - a) < 1e-10: return 1.0
    num, _ = quad(lambda r: np.exp(s * V_CST(r)) / (2 * np.sqrt(1 - r**2)), a, b)
    den, _ = quad(lambda r: 1.0 / (2 * np.sqrt(1 - r**2)), a, b)
    return num / max(den, 1e-15)

# ═══════════════════════════════════════════════════════════
# LINDBLAD SIMULATION
# ═══════════════════════════════════════════════════════════

def bloch_to_rho(r):
    I = np.eye(2, dtype=complex)
    sx = np.array([[0,1],[1,0]], dtype=complex)
    sy = np.array([[0,-1j],[1j,0]], dtype=complex)
    sz = np.array([[1,0],[0,-1]], dtype=complex)
    return 0.5 * (I + r[0]*sx + r[1]*sy + r[2]*sz)

def bures_metric(r):
    rho2 = min(np.dot(r, r), 0.9999)
    return (np.eye(3) + np.outer(r, r) / (1 - rho2)) / 4.0

def lindblad_traj(r0, H, L_ops, t_max=10., n_steps=500):
    sx = np.array([[0,1],[1,0]], dtype=complex)
    sy = np.array([[0,-1j],[1j,0]], dtype=complex)
    sz = np.array([[1,0],[0,-1]], dtype=complex)
    paulis = [sx, sy, sz]
    def rhs(t, r):
        rho = bloch_to_rho(r)
        drho = -1j * (H @ rho - rho @ H)
        for Lk in L_ops:
            LdL = Lk.conj().T @ Lk
            drho += Lk @ rho @ Lk.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
        return np.array([np.real(np.trace(drho @ p)) for p in paulis])
    sol = solve_ivp(rhs, (0, t_max), r0,
                    t_eval=np.linspace(0, t_max, n_steps), rtol=1e-10, atol=1e-12)
    return sol.t, sol.y.T

def compute_R_traj(t_arr, traj, s=1.0):
    """Compute cumulative R(t) along a trajectory."""
    d_b, d_w = np.zeros(len(t_arr)), np.zeros(len(t_arr))
    for i in range(1, len(t_arr)):
        rm = 0.5 * (traj[i-1] + traj[i])
        dr = traj[i] - traj[i-1]
        rho = np.clip(np.sqrt(np.dot(rm, rm)), 1e-10, 0.9999)
        gb = (np.eye(3) + np.outer(rm, rm) / (1 - rho**2 + 1e-15)) / 4
        V = V_CST(rho)
        ds_b = np.sqrt(max(np.dot(dr, gb @ dr), 0))
        d_b[i] = d_b[i-1] + ds_b
        d_w[i] = d_w[i-1] + np.exp(s * V) * ds_b
    R = d_w / np.maximum(d_b, 1e-15)
    R[0] = R[1] if len(R) > 1 else 1
    return d_b, d_w, R

def cst_gradient_flow(r0, n_steps=2000, dt=0.003):
    """Morse-CST flow: dr/ds = -g_Bures^{-1} · ∇V_CST."""
    traj = [r0.copy()]; r = r0.copy()
    s_arr = [0.0]; s = 0.0
    for _ in range(n_steps):
        rho = np.clip(np.sqrt(np.dot(r, r)), 1e-10, 0.9999)
        grad_V = np.arctanh(rho) * r / rho
        gb = bures_metric(r)
        try: gi = np.linalg.inv(gb)
        except: break
        flow = -gi @ grad_V
        norm = np.sqrt(max(flow @ gb @ flow, 1e-20))
        if norm < 1e-12: break
        dr = flow * dt / norm
        r_new = r + dr
        if np.linalg.norm(r_new) > 0.99:
            r_new *= 0.99 / np.linalg.norm(r_new)
        r = r_new; s += dt
        traj.append(r.copy()); s_arr.append(s)
    return np.array(s_arr), np.array(traj)

# ═══════════════════════════════════════════════════════════
# MAIN COMPUTATION
# ═══════════════════════════════════════════════════════════

def main():
    print("CST v16 — Full computation")
    print("=" * 50)
    
    # ── Define systems ──
    H0 = np.zeros((2,2), dtype=complex)
    L_decay = np.zeros((2,2), dtype=complex); L_decay[0,1] = np.sqrt(1.0)
    L_pump  = np.zeros((2,2), dtype=complex); L_pump[1,0]  = np.sqrt(1.5)
    L_dep = [np.sqrt(0.5/3) * np.array([[0,1],[0,0]], dtype=complex),
             np.sqrt(0.5/3) * np.array([[0,0],[1,0]], dtype=complex),
             np.sqrt(0.5/3) * np.array([[1,0],[0,-1]], dtype=complex)]
    
    # ── Lindblad trajectories ──
    print("Computing Lindblad trajectories...")
    t_g, tr_g = lindblad_traj(np.array([0,0,0.9]), H0, [L_decay, L_pump], 8, 600)
    t_p, tr_p = lindblad_traj(np.array([0,0,0.1]), H0, [L_decay], 8, 600)
    t_d, tr_d = lindblad_traj(np.array([0,0,0.9]), H0, L_dep, 8, 600)
    
    db_g, dw_g, R_g = compute_R_traj(t_g, tr_g)
    db_p, dw_p, R_p = compute_R_traj(t_p, tr_p)
    db_d, dw_d, R_d = compute_R_traj(t_d, tr_d)
    
    # ── CST gradient flows ──
    print("Computing Morse-CST flows...")
    s_cg, tr_cg = cst_gradient_flow(np.array([0,0,0.9]))
    s_cp, tr_cp = cst_gradient_flow(np.array([0,0,0.1]))
    
    # ── Print results ──
    print(f"\nR(t→∞):  Gain={R_g[-1]:.4f}  Passive={R_p[-1]:.4f}  Depol={R_d[-1]:.4f}")
    print(f"Gain cost reduction vs passive: {(1-R_g[-1]/R_p[-1])*100:.1f}%")
    
    # ══════════════════════════════════════════
    # FIGURE 1: Potential, cost factor, R
    # ══════════════════════════════════════════
    print("\nGenerating Figure 1...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    rho_line = np.linspace(0.01, 0.99, 300)
    
    ax = axes[0]
    ax.plot(rho_line, [V_CST(r) for r in rho_line], 'b-', lw=2.5)
    ax.fill_between(rho_line, [V_CST(r) for r in rho_line], alpha=0.08, color='blue')
    for n, rs, c in [('Gain',0.20,'red'),('Passive',0.98,'green'),('Depol.',0.01,'orange')]:
        ax.plot(rs, V_CST(rs), 'o', color=c, ms=10, zorder=5, label=f'{n} (ρ={rs})')
    ax.set_xlabel('Purity ρ = |r|'); ax.set_ylabel('V(ρ) = S(ρ || I/2)')
    ax.set_title('(a) CST Potential — Unique (Umegaki)', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[1]
    ax.plot(rho_line, [np.exp(V_CST(r)) for r in rho_line], 'b-', lw=2.5)
    ax.fill_between(rho_line, 1, [np.exp(V_CST(r)) for r in rho_line], alpha=0.08, color='blue')
    for n, rs, c in [('Gain',0.20,'red'),('Passive',0.98,'green'),('Depol.',0.01,'orange')]:
        ax.plot(rs, np.exp(V_CST(rs)), 'o', color=c, ms=10, zorder=5, label=n)
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Purity ρ'); ax.set_ylabel('e^{V(ρ)} (cost factor)')
    ax.set_title('(b) Witten Deformation Factor (s=1)', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[2]
    rho2_scan = np.linspace(0.02, 0.98, 150)
    ax.plot(rho2_scan, [R_ratio(0.90, r2) for r2 in rho2_scan], 'r-', lw=2, label='From ρ₁=0.90')
    ax.plot(rho2_scan, [R_ratio(0.50, r2) for r2 in rho2_scan], 'b-', lw=2, label='From ρ₁=0.50')
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    ax.plot(0.20, R_ratio(0.90, 0.20), 'r*', ms=14, zorder=5, label='Toronto (gain→ρ=0.20)')
    ax.set_xlabel('Final purity ρ₂'); ax.set_ylabel('R = t_CST / t_QM')
    ax.set_title('(c) Time Ratio R (closed form)', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig1_potential.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # ══════════════════════════════════════════
    # FIGURE 2: R(t), ΔR, R vs s
    # ══════════════════════════════════════════
    print("Generating Figure 2...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    ax = axes[0]
    ax.plot(t_g, R_g, 'r-', lw=2.5, label='Gain (Toronto-like)')
    ax.plot(t_p, R_p, 'g-', lw=2.5, label='Passive decay')
    ax.plot(t_d, R_d, color='orange', lw=2.5, label='Depolarizing')
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Parametric time t'); ax.set_ylabel('R(t) = d_Witten / d_Bures')
    ax.set_title('(a) CST Time Ratio R(t)', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[1]
    tc = np.linspace(0.1, 7.5, 300)
    Rgi = np.interp(tc, t_g, R_g); Rpi = np.interp(tc, t_p, R_p); Rdi = np.interp(tc, t_d, R_d)
    ax.plot(tc, Rgi-Rdi, 'r-', lw=2.5, label='ΔR = R(gain) − R(depol)')
    ax.plot(tc, Rpi-Rdi, 'g-', lw=2.5, label='ΔR = R(passive) − R(depol)')
    ax.fill_between(tc, Rgi-Rdi, alpha=0.12, color='red')
    ax.fill_between(tc, Rpi-Rdi, alpha=0.12, color='green')
    ax.axhline(0, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Parametric time t'); ax.set_ylabel('ΔR')
    ax.set_title('(b) CST Signature: ΔR ≠ 0', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[2]
    s_scan = np.linspace(0.01, 3.0, 80)
    ax.plot(s_scan, [R_ratio(0.90,0.20,s) for s in s_scan], 'r-', lw=2, label='Gain (0.90→0.20)')
    ax.plot(s_scan, [R_ratio(0.10,0.98,s) for s in s_scan], 'g-', lw=2, label='Passive (0.10→0.98)')
    ax.plot(s_scan, [R_ratio(0.90,0.01,s) for s in s_scan], color='orange', lw=2, label='Depol. (0.90→0.01)')
    ax.axvline(1.0, color='blue', ls=':', lw=2, alpha=0.5, label='s=1 (N=2)')
    ax.axvline(0.5, color='purple', ls=':', lw=2, alpha=0.5, label='s=½ (N=3, Rb87)')
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    ax.set_xlabel('Deformation parameter s'); ax.set_ylabel('R = t_CST / t_QM')
    ax.set_title('(c) R vs s — The Single Dial', fontweight='bold')
    ax.legend(fontsize=8); ax.grid(alpha=0.15); ax.set_yscale('log'); ax.set_ylim(0.5, 30)
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig2_signal.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # ══════════════════════════════════════════
    # FIGURE 3: Flows
    # ══════════════════════════════════════════
    print("Generating Figure 3...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    ax = axes[0]
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(0.95*np.cos(theta), 0.95*np.sin(theta), 'k-', alpha=0.1)
    ax.plot(tr_g[:,0], tr_g[:,2], 'r-', lw=2.5, label='Gain Lindblad')
    ax.plot(tr_p[:,0], tr_p[:,2], 'g-', lw=2.5, label='Passive Lindblad')
    ax.plot(tr_d[:,0], tr_d[:,2], color='orange', lw=2.5, label='Depol. Lindblad')
    ax.plot(tr_cg[:,0], tr_cg[:,2], 'r--', lw=1.5, alpha=0.6, label='CST gradient')
    ax.plot(tr_cp[:,0], tr_cp[:,2], 'g--', lw=1.5, alpha=0.6)
    ax.plot(0, 0, 'k+', ms=12, mew=2, label='I/2 (center)')
    ax.set_xlabel('r_x'); ax.set_ylabel('r_z')
    ax.set_title('(a) Bloch ball: Lindblad (—) vs CST (--)', fontweight='bold')
    ax.set_aspect('equal'); ax.legend(fontsize=7, loc='lower left'); ax.grid(alpha=0.15)
    
    ax = axes[1]
    rho_g = np.sqrt(np.sum(tr_g**2, axis=1))
    rho_p = np.sqrt(np.sum(tr_p**2, axis=1))
    rho_d = np.sqrt(np.sum(tr_d**2, axis=1))
    ax.plot(t_g, rho_g, 'r-', lw=2.5, label='Gain')
    ax.plot(t_p, rho_p, 'g-', lw=2.5, label='Passive')
    ax.plot(t_d, rho_d, color='orange', lw=2.5, label='Depolarizing')
    ax.plot(s_cg, np.sqrt(np.sum(tr_cg**2,axis=1)), 'r--', lw=1.5, alpha=0.6)
    ax.plot(s_cp, np.sqrt(np.sum(tr_cp**2,axis=1)), 'g--', lw=1.5, alpha=0.6)
    ax.set_xlabel('Parameter (t or s)'); ax.set_ylabel('Purity ρ(t)')
    ax.set_title('(b) Purity evolution', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[2]
    ax.plot(t_g, [V_CST(rho_g[i]) for i in range(len(t_g))], 'r-', lw=2.5, label='Gain → descends V')
    ax.plot(t_p, [V_CST(rho_p[i]) for i in range(len(t_p))], 'g-', lw=2.5, label='Passive → ascends V')
    ax.plot(t_d, [V_CST(rho_d[i]) for i in range(len(t_d))], color='orange', lw=2.5, label='Depol. → descends V')
    ax.set_xlabel('Parametric time t'); ax.set_ylabel('V_CST(ρ(t))')
    ax.set_title('(c) Configurational potential along Lindblad', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig3_flows.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # ══════════════════════════════════════════
    # FIGURE 4: Toronto
    # ══════════════════════════════════════════
    print("Generating Figure 4...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    
    ax = axes[0]
    states = ['I/3', 'Gain\n(Toronto)', 'Moderate\ngain', 'Passive\nthermal', 'Pure |2⟩']
    pops_list = [[1/3,1/3,1/3],[.2,.6,.2],[.25,.5,.25],[.6,.3,.1],[.01,.98,.01]]
    esV3 = [np.exp(0.5 * V_CST_N(p)) for p in pops_list]
    colors_bar = ['orange','red','salmon','green','darkgreen']
    bars = ax.bar(range(len(states)), esV3, color=colors_bar, alpha=0.7, edgecolor='black')
    ax.set_xticks(range(len(states))); ax.set_xticklabels(states, fontsize=8)
    ax.set_ylabel('e^{sV} (configurational cost)')
    ax.set_title('(a) Rb87 (N=3, s=½): Cost per state', fontweight='bold')
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    for i, (bar, v) in enumerate(zip(bars, esV3)):
        ax.text(i, v + 0.015, f'{v:.3f}', ha='center', fontsize=9)
    ax.grid(alpha=0.15, axis='y')
    
    ax = axes[1]
    g_strengths = np.linspace(0.01, 0.99, 100)
    ax.plot(g_strengths, [np.exp(V_CST(g)-V_CST(0.98)) for g in g_strengths],
            'b-', lw=2.5, label='N=2 (qubit, s=1)')
    ax.plot(g_strengths,
            [np.exp(0.5*(V_CST_N([g/2,1-g,g/2])-V_CST_N([0.6,0.3,0.1]))) for g in g_strengths],
            'r-', lw=2.5, label='N=3 (Rb87, s=½)')
    ax.axhline(1, color='k', ls='--', alpha=0.3)
    ax.axvline(0.20, color='red', ls=':', alpha=0.4, label='Toronto ρ_ss≈0.20')
    ax.set_xlabel('Gain steady-state purity ρ_ss')
    ax.set_ylabel('Cost ratio: gain / passive')
    ax.set_title('(b) Geometric subsidy vs gain strength', fontweight='bold')
    ax.legend(fontsize=9); ax.grid(alpha=0.15)
    
    ax = axes[2]
    summary = (
        "CST PREDICTIONS FOR TORONTO\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        f"System         R(t→∞)   ΔR vs depol.\n"
        f"────────────── ──────── ────────────\n"
        f"Gain (Toronto)  {R_g[-1]:.4f}   {R_g[-1]-R_d[-1]:+.4f}\n"
        f"Passive decay   {R_p[-1]:.4f}   {R_p[-1]-R_d[-1]:+.4f}\n"
        f"Depolarizing    {R_d[-1]:.4f}   {0.0:+.4f}\n\n"
        f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"
        f"Gain cost reduction: {(1-R_g[-1]/R_p[-1])*100:.1f}%\n\n"
        f"The gain medium maintains atoms\n"
        f"in a LOW-V region of C-space.\n"
        f"Probe transit costs LESS in CST\n"
        f"time → geometric subsidy.\n\n"
        f"τ_eff < 0 is not an illusion:\n"
        f"it reflects the pre-paid\n"
        f"configurational asymmetry cost."
    )
    ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
    ax.axis('off')
    ax.set_title('(c) Quantitative predictions', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/fig4_toronto.png', dpi=200, bbox_inches='tight')
    plt.close()
    
    # ══════════════════════════════════════════
    # DONE
    # ══════════════════════════════════════════
    print(f"\nAll figures saved to {OUTDIR}/")
    for f in sorted(os.listdir(OUTDIR)):
        print(f"  {f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
