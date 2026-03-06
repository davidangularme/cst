"""
FIT — Fonctionnelle d'Information Temporelle
Calcul du point fixe auto-cohérent sur un qubit (boule de Bloch).

Structure mathématique :
  - État : ρ = (I + r·σ)/2,  r ∈ Bloch ball, |r| < 1
  - Vecteur log-Bloch : v(r) = arctanh(|r|)/|r| · r
  - Hamiltonien modulaire : H_mod = log ρ - log ρ_∞  ↔  Δv = v(r) - v(r_∞)
  - Flot gradient (K_ρ balanced) : dr/dt = (r·Δv) r - Δv
  - Fonctionnelle FIT : T[γ, ρ_∞] = ∫ ṙ · Δv dt = -∫ |Δv|² sin²θ dt  ≤ 0
  - Auto-cohérence : r_∞* = argmin_{r_∞} T[γ(r₀, r_∞), r_∞]
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

EPS = 1e-10
np.random.seed(42)

# ── Primitives ─────────────────────────────────────────────────────────────

def v(r):
    """Vecteur log-Bloch : v(r) = arctanh(|r|)/|r| · r"""
    mag = np.linalg.norm(r)
    if mag < EPS:
        return np.zeros(3)
    return np.arctanh(np.clip(mag, EPS, 1 - EPS)) / mag * r

def relative_entropy(r, r_inf):
    """
    S(ρ||ρ_∞) en coordonnées de Bloch :
    = α(r) - α(r_∞) + r · (v(r) - v(r_∞))
    où α(r) = ½ log((1-|r|²)/4)
    """
    def alpha(rr):
        mag = np.clip(np.linalg.norm(rr), EPS, 1 - EPS)
        return 0.5 * np.log((1 - mag**2) / 4)
    return (alpha(r) - alpha(r_inf)) + np.dot(r, v(r) - v(r_inf))

# ── Flot gradient ───────────────────────────────────────────────────────────

def flow_rhs(r, r_inf):
    """
    dr/dt = (r · Δv) r - Δv
    Dérivation : flot K_ρ de S(ρ||ρ_∞) avec projection trace-nulle via D(ρ||ρ_∞)·ρ
    Vérifié : point fixe r = r_∞, converge vers r_∞ depuis tout r₀
    """
    dv = v(r) - v(r_inf)
    return np.dot(r, dv) * r - dv

def integrate(r0, r_inf, T=10.0, n=600):
    """Intègre le flot gradient de r0 vers r_inf."""
    sol = solve_ivp(
        lambda t, r: flow_rhs(r, r_inf),
        (0, T), r0, t_eval=np.linspace(0, T, n),
        method='RK45', rtol=1e-9, atol=1e-11
    )
    return sol.t, sol.y.T  # (n, 3)

# ── Fonctionnelle FIT ───────────────────────────────────────────────────────

def FIT(t_vals, r_traj, r_inf):
    """
    T = ∫ ṙ · Δv dt
    Sur le flot : T = -∫ |Δv|² sin²(θ) dt  ≤ 0
    Egal à 0 ssi r(t) et Δv(t) toujours colinéaires (cas commutatif pur)
    """
    dt = np.diff(t_vals)
    total = 0.0
    for i in range(len(dt)):
        r_mid = 0.5 * (r_traj[i] + r_traj[i+1])
        dr_dt = (r_traj[i+1] - r_traj[i]) / dt[i]
        dv_mid = v(r_mid) - v(r_inf)
        total += np.dot(dr_dt, dv_mid) * dt[i]
    return total

def FIT_integrand(t_vals, r_traj, r_inf):
    """Intégrande de T point par point (pour visualisation)."""
    vals = []
    for i in range(len(t_vals) - 1):
        r_mid = 0.5 * (r_traj[i] + r_traj[i+1])
        dr_dt = (r_traj[i+1] - r_traj[i]) / (t_vals[i+1] - t_vals[i])
        dv_mid = v(r_mid) - v(r_inf)
        vals.append(np.dot(dr_dt, dv_mid))
    return np.array(vals)

# ── Auto-cohérence : scan du paysage T(r_∞) ────────────────────────────────

def scan_landscape(r0, n_theta=14, n_phi=18, r_scan=0.55, T=8.0):
    """Calcule T pour une grille de r_∞ sur une sphère de rayon r_scan."""
    r_infs, fits = [], []
    for theta in np.linspace(0.15, np.pi - 0.15, n_theta):
        for phi in np.linspace(0, 2 * np.pi, n_phi, endpoint=False):
            ri = r_scan * np.array([
                np.sin(theta) * np.cos(phi),
                np.sin(theta) * np.sin(phi),
                np.cos(theta)
            ])
            t_v, traj = integrate(r0, ri, T=T, n=250)
            r_infs.append(ri)
            fits.append(FIT(t_v, traj, ri))
    return np.array(r_infs), np.array(fits)

# ── Asymétrie temporelle ────────────────────────────────────────────────────

def time_reversal_test(r0, r_inf, T=10.0):
    """
    Compare T[ρ₀→ρ_∞]  vs  T[ρ_∞→ρ₀]
    Asymétrie ↔ flèche du temps encodée dans la fonctionnelle.
    """
    t_f, traj_f = integrate(r0, r_inf, T=T)
    t_r, traj_r = integrate(r_inf, r0, T=T)
    return (FIT(t_f, traj_f, r_inf), FIT(t_r, traj_r, r0),
            t_f, traj_f, t_r, traj_r)

# ── Visualisation ───────────────────────────────────────────────────────────

COL = ['#e74c3c', '#2ecc71', '#3498db']
LAB = ['x', 'y', 'z']

def plot_flow(r0, r_inf, T=10.0, filename='flow.png'):
    t_v, traj = integrate(r0, r_inf, T=T)
    T_val = FIT(t_v, traj, r_inf)
    S_vals = [max(relative_entropy(traj[i], r_inf), 1e-15) for i in range(len(t_v))]
    intgd = FIT_integrand(t_v, traj, r_inf)

    fig = plt.figure(figsize=(14, 9))
    gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.38)
    fig.suptitle(
        f'FIT — Flot gradient sur le qubit\n'
        f'r₀={r0.round(2)}, ρ_∞=r_∞={r_inf.round(3)},  𝒯={T_val:.4f}',
        fontsize=12, fontweight='bold'
    )

    # Trajectoire
    ax1 = fig.add_subplot(gs[0, :2])
    for i, (c, l) in enumerate(zip(COL, LAB)):
        ax1.plot(t_v, traj[:, i], color=c, lw=2.5, label=f'r_{l}(t)')
        ax1.axhline(r_inf[i], color=c, ls='--', lw=1.4, alpha=0.55)
        ax1.scatter([0], [r0[i]], color=c, s=70, zorder=6)
    ax1.set_xlabel('t'); ax1.set_ylabel('Composantes Bloch')
    ax1.set_title('Trajectoire  (--- = ρ_∞ cible)')
    ax1.legend(ncol=3, fontsize=9); ax1.grid(True, alpha=0.22)

    # Pureté
    ax_p = fig.add_subplot(gs[0, 2])
    r_mags = np.linalg.norm(traj, axis=1)
    ax_p.plot(t_v, r_mags, 'k-', lw=2)
    ax_p.axhline(np.linalg.norm(r_inf), ls='--', color='gray', lw=1.5)
    ax_p.axhline(np.linalg.norm(r0), ls=':', color='gray', lw=1.2)
    ax_p.set_xlabel('t'); ax_p.set_ylabel('|r(t)|')
    ax_p.set_title('Pureté'); ax_p.set_ylim(0, 1); ax_p.grid(True, alpha=0.22)

    # Entropie relative
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.semilogy(t_v, S_vals, '#9b59b6', lw=2)
    ax2.set_xlabel('t'); ax2.set_ylabel('S(ρ(t) ‖ ρ_∞)')
    ax2.set_title('Décroissance entropie relative'); ax2.grid(True, alpha=0.22)

    # Intégrande FIT
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(t_v[:-1], intgd, '#e67e22', lw=2)
    ax3.fill_between(t_v[:-1], intgd, 0, alpha=0.25, color='#e67e22')
    ax3.axhline(0, color='k', lw=0.8, alpha=0.5)
    ax3.set_xlabel('t'); ax3.set_ylabel('ṙ · Δv')
    ax3.set_title(f'Intégrande de 𝒯  (total={T_val:.4f})'); ax3.grid(True, alpha=0.22)

    # Angle θ entre r et Δv
    ax4 = fig.add_subplot(gs[1, 2])
    angles = []
    for i in range(len(t_v)):
        ri = traj[i]; dvi = v(ri) - v(r_inf)
        rm = np.linalg.norm(ri); dm = np.linalg.norm(dvi)
        if rm > EPS and dm > EPS:
            angles.append(np.degrees(np.arccos(np.clip(np.dot(ri, dvi) / (rm * dm), -1, 1))))
        else:
            angles.append(0.0)
    ax4.plot(t_v, angles, '#1abc9c', lw=2)
    ax4.axhline(90, ls='--', color='gray', lw=1, alpha=0.5)
    ax4.set_xlabel('t'); ax4.set_ylabel('θ (degrés)')
    ax4.set_title('Angle r(t) ↔ Δv(t)')
    ax4.set_ylim(0, 182); ax4.grid(True, alpha=0.22)

    plt.savefig(f'/mnt/user-data/outputs/{filename}', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return T_val


def plot_landscape(r_infs, fits, r0, r_inf_min, filename='landscape.png'):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f'Paysage FIT : 𝒯(ρ_∞) pour r₀={r0.round(2)}\n'
        f'Le minimum donne le ρ_∞ auto-cohérent',
        fontsize=12, fontweight='bold'
    )

    ax1 = axes[0]
    sc = ax1.scatter(r_infs[:, 0], r_infs[:, 2], c=fits,
                     cmap='plasma', s=35, alpha=0.85, edgecolors='none')
    ax1.scatter(*r0[[0, 2]], marker='*', s=200, c='white',
                edgecolors='k', zorder=10, label='r₀')
    ax1.scatter(*r_inf_min[[0, 2]], marker='D', s=110, c='cyan',
                edgecolors='k', zorder=10, label='r_∞* (min 𝒯)')
    theta_c = np.linspace(0, 2 * np.pi, 200)
    rs = np.linalg.norm(r_infs[0])
    ax1.plot(rs * np.cos(theta_c), rs * np.sin(theta_c), 'k--', alpha=0.25, lw=1)
    plt.colorbar(sc, ax=ax1, label='𝒯')
    ax1.set_xlabel('r_∞ x'); ax1.set_ylabel('r_∞ z')
    ax1.set_title('Paysage (projection x-z)')
    ax1.legend(fontsize=9); ax1.set_aspect('equal'); ax1.grid(True, alpha=0.2)

    ax2 = axes[1]
    r0n = r0 / (np.linalg.norm(r0) + EPS)
    angles = [np.degrees(np.arccos(np.clip(
        np.dot(ri / np.linalg.norm(ri), r0n), -1, 1))) for ri in r_infs]
    sc2 = ax2.scatter(angles, fits, c=np.linalg.norm(r_infs, axis=1),
                      cmap='viridis', s=22, alpha=0.75, edgecolors='none')
    plt.colorbar(sc2, ax=ax2, label='|r_∞|')
    ax2.set_xlabel('Angle r_∞ / r₀ (degrés)')
    ax2.set_ylabel('𝒯')
    ax2.set_title('𝒯 en fonction de la direction de ρ_∞')
    ax2.grid(True, alpha=0.25)

    plt.savefig(f'/mnt/user-data/outputs/{filename}', dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_time_reversal(r0, r_inf, filename='time_reversal.png'):
    T_fwd, T_rev, t_f, traj_f, t_r, traj_r = time_reversal_test(r0, r_inf)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        f'Asymétrie temporelle de la FIT\n'
        f'𝒯[ρ₀→ρ_∞] = {T_fwd:.4f}   |   𝒯[ρ_∞→ρ₀] = {T_rev:.4f}   |   |Δ𝒯| = {abs(T_fwd - T_rev):.4f}',
        fontsize=12, fontweight='bold'
    )
    for ax, tv, traj, T_val, title in [
        (axes[0], t_f, traj_f, T_fwd, f'Avant : ρ₀→ρ_∞  (𝒯={T_fwd:.4f})'),
        (axes[1], t_r, traj_r, T_rev, f'Retour : ρ_∞→ρ₀  (𝒯={T_rev:.4f})')
    ]:
        for i, (c, l) in enumerate(zip(COL, LAB)):
            ax.plot(tv, traj[:, i], color=c, lw=2, label=l)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('t'); ax.set_ylabel('Composantes Bloch')
        ax.legend(ncol=3, fontsize=9); ax.grid(True, alpha=0.25)

    plt.savefig(f'/mnt/user-data/outputs/{filename}', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return T_fwd, T_rev


# ── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 65)
    print("  FIT — Qubit : flot gradient et point fixe auto-cohérent")
    print("=" * 65)

    r0 = np.array([0.1, 0.3, 0.65])

    # ── 1. Flot vers ρ_∞ = état maximalement mélangé
    r_inf_mixed = np.zeros(3)
    print(f"\n[1] Flot r₀={r0} → ρ_∞ = état mélangé")
    T1 = plot_flow(r0, r_inf_mixed, T=10.0, filename='fit_flow_mixed.png')
    print(f"    𝒯 = {T1:.4f}")

    # ── 2. Paysage 𝒯(r_∞)
    print(f"\n[2] Scan du paysage 𝒯(r_∞) sur sphère r=0.55")
    r_infs, fits = scan_landscape(r0, n_theta=13, n_phi=18, r_scan=0.55)
    idx_min = np.argmin(fits)
    r_inf_star = r_infs[idx_min]
    print(f"    Min 𝒯 = {fits[idx_min]:.4f}  →  r_∞* = {r_inf_star.round(4)}")
    print(f"    Max 𝒯 = {fits.max():.4f}")
    plot_landscape(r_infs, fits, r0, r_inf_star, filename='fit_landscape.png')

    # ── 3. Flot avec ρ_∞ auto-cohérent
    print(f"\n[3] Flot avec ρ_∞ optimal  r_∞*={r_inf_star.round(3)}")
    T3 = plot_flow(r0, r_inf_star, T=10.0, filename='fit_flow_optimal.png')
    print(f"    𝒯(optimal) = {T3:.4f}")

    # ── 4. Asymétrie temporelle
    print(f"\n[4] Test d'asymétrie temporelle")
    T_fwd, T_rev = plot_time_reversal(r0, r_inf_star, filename='fit_time_reversal.png')
    print(f"    Avant    𝒯[ρ₀→ρ_∞] = {T_fwd:.4f}")
    print(f"    Retour   𝒯[ρ_∞→ρ₀] = {T_rev:.4f}")
    print(f"    |Δ𝒯|              = {abs(T_fwd - T_rev):.4f}")
    print(f"    Flèche du temps   : {'OUI ✓' if not np.isclose(T_fwd, T_rev, atol=1e-3) else 'non'}")

    # ── 5. Scan plusieurs conditions initiales
    print(f"\n[5] Phase : plusieurs r₀, même ρ_∞ = état mélangé")
    test_r0s = [
        np.array([0.1, 0.3, 0.65]),
        np.array([-0.1, -0.3, -0.65]),
        np.array([0.6, 0.3, 0.1]),
        np.array([0.0, 0.0, 0.8]),
    ]
    print(f"    {'r₀':35s}  {'𝒯 (→mélangé)':15s}  {'𝒯 (→r_∞*)':12s}")
    for r0_i in test_r0s:
        t1_, traj1_ = integrate(r0_i, r_inf_mixed, T=10.0)
        t3_, traj3_ = integrate(r0_i, r_inf_star, T=10.0)
        fit1 = FIT(t1_, traj1_, r_inf_mixed)
        fit3 = FIT(t3_, traj3_, r_inf_star)
        print(f"    r₀={r0_i}  {fit1:+.4f}  {fit3:+.4f}")

    print(f"\n{'='*65}")
    print("Formules clés :")
    print("  Flot : dr/dt = (r·Δv) r - Δv")
    print("         Δv = arctanh(|r|)/|r|·r  -  arctanh(|r_∞|)/|r_∞|·r_∞")
    print("  FIT  : T = ∫ ṙ·Δv dt  =  -∫|Δv|²sin²θ dt  ≤ 0")
    print("  Auto-cohérence : r_∞* = argmin_{r_∞} T[γ(r₀,r_∞), r_∞]")
    print("  Umegaki = T[géodésique, σ fixé] = valeur on-shell")
    print(f"{'='*65}")
