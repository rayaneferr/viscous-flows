"""
Écoulements visqueux — de la diffusion de quantité de mouvement aux pertes de charge.

Quatre démonstrations, toutes réglables en ligne de commande
(voir `uv run main.py --help`). Les figures sont écrites dans `figures/`.

    uv run main.py --mode stokes        # mise en mouvement d'une plaque (MF3-1/2/4)
    uv run main.py --mode poiseuille    # artère et citerne fêlée (MF3-7, MF3-8)
    uv run main.py --mode moody         # château d'eau et oléoduc (MF3-6, MF3-9)
    uv run main.py --mode check         # vérifications numériques

Convention : code en anglais, commentaires en français.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from src import (
    RE_LAMINAR,
    RE_TURBULENT,
    bend_loss_coefficient,
    boundary_layer_thickness,
    colebrook_friction_factor,
    couette_cylindrical_profile,
    couette_steady,
    couette_torque,
    diffusion_time,
    friction_factor,
    haaland_friction_factor,
    laminar_friction_factor,
    mean_velocity,
    mercury_pressure,
    poiseuille_flow_rate,
    poiseuille_profile,
    pressure_drop_for_flow,
    regular_head_loss,
    reynolds,
    roughness_from_velocity,
    similarity_variable,
    singular_head_loss,
    slit_flow_rate,
    slit_profile,
    solve_startup,
    stenosis_pressure_ratio,
    stokes_first_problem,
    torricelli_flow_rate,
    velocity_from_head,
    viscosity_from_sliding_block,
    viscosity_from_torque,
    wall_shear_stress,
)

FIGURES = Path(__file__).parent / "figures"

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREY = "0.45"


def style_axes(ax, *, grid=True):
    """Habillage commun à toutes les figures (sobre, lisible en clair comme en sombre)."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(alpha=0.25, linewidth=0.6)


# --------------------------------------------------------------------------- #
# Mode 1 — diffusion de quantité de mouvement (MF3-1, MF3-2, MF3-4)
# --------------------------------------------------------------------------- #
def run_stokes(args):
    nu = args.eta / args.rho
    tau = diffusion_time(args.gap, nu=nu)
    print(f"  ν = η/ρ = {nu:.2e} m²/s")
    print(f"  temps de diffusion τ = L²/ν = {tau:.3g} s  ({tau / 3600:.2f} h)")

    z, times, V = solve_startup(
        gap=args.gap,
        n_points=args.n_points,
        nu=nu,
        plate_speed=args.plate_speed,
        t_max=args.t_ratio * tau,
        n_steps=args.n_steps,
    )

    # ---- 1a. animation du transitoire ------------------------------------- #
    # Échantillonnage en √t : la diffusion avance en √(νt), un échantillonnage
    # linéaire écraserait tout le début intéressant du transitoire.
    frac = np.linspace(0.0, 1.0, args.frames) ** 2
    idx = np.unique((frac * (len(times) - 1)).astype(int))

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    style_axes(ax)
    (line,) = ax.plot([], [], color=BLUE, lw=2.2, label="Crank–Nicolson")
    (exact,) = ax.plot([], [], color=ACCENT, lw=1.4, ls="--", label=r"erfc (Stokes, milieu semi-infini)")
    (delta,) = ax.plot([], [], color=GREY, lw=1.0, ls=":",
                       label=r"couche limite $\delta = \sqrt{\nu t}$")
    steady = couette_steady(z, gap=args.gap, plate_speed=args.plate_speed)
    ax.plot(steady, z, color="0.7", lw=1.2, label="Couette stationnaire")
    ax.set_xlim(-0.03 * args.plate_speed, 1.05 * args.plate_speed)
    ax.set_ylim(0, args.gap)
    ax.set_xlabel("vitesse $v$ (m/s)")
    ax.set_ylabel("$z$ (m)")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    title = ax.set_title("")

    def update(k):
        t = times[k]
        line.set_data(V[k], z)
        exact.set_data(stokes_first_problem(z, t, plate_speed=args.plate_speed, nu=nu), z)
        d = boundary_layer_thickness(t, nu=nu)
        delta.set_data([0, args.plate_speed], [d, d])
        title.set_text(
            f"$t = {t / tau:.3f}\\,\\tau$   —   $\\delta=\\sqrt{{\\nu t}} = {d * 1e3:.1f}$ mm"
        )
        return line, exact, delta, title

    anim = FuncAnimation(fig, update, frames=idx, blit=False, interval=1000 / args.fps)
    out = FIGURES / "stokes_startup.gif"
    anim.save(out, writer="pillow", fps=args.fps, dpi=110)
    plt.close(fig)
    print(f"  → {out.name}")

    # ---- 1b. similitude et convergence ------------------------------------ #
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    # Panneau gauche : effondrement des profils sur la courbe de similitude.
    # Tant que δ(t) << L, le fluide « ignore » la plaque du haut et les profils
    # numériques se superposent tous à erfc(ζ) une fois tracés en ζ = z/2√(νt).
    for ratio, color in zip((0.002, 0.005, 0.01, 0.02), plt.cm.viridis(np.linspace(0.15, 0.8, 4))):
        k = int(ratio * (len(times) - 1) / args.t_ratio)
        t = times[k]
        zeta = similarity_variable(z[1:], t, nu=nu)
        ax1.plot(zeta, V[k][1:] / args.plate_speed, color=color, lw=1.8,
                 label=f"$t = {ratio:g}\\,\\tau$")
    zeta_ref = np.linspace(0, 2.5, 200)
    ax1.plot(zeta_ref, np.asarray([stokes_first_problem(np.array([2 * zz]), 1.0,
                                                        plate_speed=1.0, nu=1.0)[0]
                                   for zz in zeta_ref]),
             color=ACCENT, ls="--", lw=1.6, label=r"$\mathrm{erfc}(\zeta)$")
    ax1.set_xlim(0, 2.5)
    ax1.set_xlabel(r"variable de similitude  $\zeta = z / 2\sqrt{\nu t}$")
    ax1.set_ylabel("$v / U$")
    ax1.set_title("Aux temps courts : une diffusion auto-semblable")
    ax1.legend(fontsize=8)
    style_axes(ax1)

    # Panneau droit : convergence vers le profil de Couette linéaire.
    for ratio, color in zip((0.05, 0.15, 0.4, 1.0), plt.cm.viridis(np.linspace(0.15, 0.8, 4))):
        k = min(int(ratio * (len(times) - 1) / args.t_ratio), len(times) - 1)
        ax2.plot(V[k], z, color=color, lw=1.8, label=f"$t = {ratio:g}\\,\\tau$")
    ax2.plot(steady, z, color=ACCENT, ls="--", lw=1.6, label="Couette : $U(1-z/L)$")
    ax2.set_xlabel("vitesse $v$ (m/s)")
    ax2.set_ylabel("$z$ (m)")
    ax2.set_ylim(0, args.gap)
    ax2.set_title("Aux temps longs : le profil de Couette")
    ax2.legend(fontsize=8)
    style_axes(ax2)

    fig.tight_layout()
    out = FIGURES / "stokes_similarity.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")

    tau_w = wall_shear_stress(eta=args.eta, plate_speed=args.plate_speed, gap=args.gap)
    print(f"  contrainte pariétale en régime établi : τ = ηU/L = {tau_w:.3g} Pa")

    # ---- 1c. viscosimètre de Couette cylindrique (MF3-2) ------------------ #
    r = np.linspace(args.r_inner, args.r_outer, 300)
    v = couette_cylindrical_profile(r, r_inner=args.r_inner, r_outer=args.r_outer,
                                    omega=args.omega)
    gamma = couette_torque(eta=args.eta_couette, height=args.height,
                           r_inner=args.r_inner, r_outer=args.r_outer, omega=args.omega)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.plot(r * 1e2, v, color=BLUE, lw=2.2, label=r"$v(r) = Ar + B/r$")
    # rotation en bloc (B = 0) : le fluide tournerait comme un solide, sans cisaillement
    ax1.plot(r * 1e2, args.omega * r, color=GREY, ls=":", lw=1.6,
             label=r"rotation en bloc ($B=0$) : aucun couple")
    ax1.axvline(args.r_inner * 1e2, color="0.3", lw=1.0)
    ax1.axvline(args.r_outer * 1e2, color="0.3", lw=1.0)
    ax1.set_xlabel("$r$ (cm)")
    ax1.set_ylabel(r"vitesse orthoradiale $v_\theta$ (m/s)")
    ax1.set_title("Viscosimètre de Couette : profil entre les cylindres")
    ax1.legend(fontsize=8)
    style_axes(ax1)

    # Inversion : le couple est proportionnel à η, avec un facteur géométrique connu.
    etas = np.linspace(0.05, 2.0, 200)
    torques = [couette_torque(eta=e, height=args.height, r_inner=args.r_inner,
                              r_outer=args.r_outer, omega=args.omega) for e in etas]
    ax2.plot(etas, torques, color=BLUE, lw=2.2)
    ax2.plot([args.eta_couette], [gamma], "o", color=ACCENT, ms=7, zorder=5)
    ax2.annotate(f"mesure : Γ = {gamma * 1e3:.2f} mN·m\n→ η = {args.eta_couette:.2f} Pa·s",
                 xy=(args.eta_couette, gamma), xytext=(0.35, 0.75),
                 textcoords="axes fraction", fontsize=9, color=ACCENT,
                 arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.2))
    ax2.set_xlabel(r"viscosité dynamique $\eta$ (Pa·s)")
    ax2.set_ylabel(r"couple $\Gamma$ (N·m)")
    ax2.set_title(r"Étalonnage : $\Gamma = 4\pi\eta H\,\omega R_1^2R_2^2/(R_2^2-R_1^2)$")
    style_axes(ax2)

    fig.tight_layout()
    out = FIGURES / "couette_viscometer.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")

    eta_back = viscosity_from_torque(torque=gamma, height=args.height,
                                     r_inner=args.r_inner, r_outer=args.r_outer,
                                     omega=args.omega)
    print(f"  couple mesuré Γ = {gamma * 1e3:.3f} mN·m  →  η retrouvé = {eta_back:.4f} Pa·s")

    # ---- 1d. application numérique MF3-4 ---------------------------------- #
    eta_oil = viscosity_from_sliding_block(mass=10.0, area=0.1, angle=np.radians(25.0),
                                           speed=0.80, film=1e-3)
    re_film = reynolds(rho=900.0, velocity=0.80, diameter=1e-3, eta=eta_oil)
    t_estab = (1e-3) ** 2 / (eta_oil / 900.0)
    print("\n  MF3-4 — solide de 10 kg glissant sur un film d'huile de 1 mm :")
    print(f"    η = Mg·sinα·e/(S·v₀) = {eta_oil:.3f} Pa·s  (ordre de grandeur d'une huile moteur)")
    print(f"    Re = ρv₀e/η = {re_film:.2f} ≪ 2300  →  écoulement bien laminaire")
    print(f"    établissement e²/ν = {t_estab * 1e3:.2f} ms  →  régime stationnaire instantané")


# --------------------------------------------------------------------------- #
# Mode 2 — écoulements internes établis (MF3-7, MF3-8)
# --------------------------------------------------------------------------- #
def run_poiseuille(args):
    # ---- 2a. l'artère ------------------------------------------------------ #
    dp = pressure_drop_for_flow(flow_rate=args.flow_rate, radius=args.radius,
                                length=args.length, eta=args.eta_blood)
    dp_heart = mercury_pressure(args.heart_pressure)
    u_mean = mean_velocity(flow_rate=args.flow_rate, radius=args.radius)
    re = reynolds(rho=args.rho_blood, velocity=u_mean, diameter=2 * args.radius,
                  eta=args.eta_blood)

    print(f"  artère : R = {args.radius * 1e3:.1f} mm, L = {args.length:.2f} m, "
          f"Dv = {args.flow_rate * 1e6:.0f} cm³/s")
    print(f"    ΔP_artère = 8ηL·Dv/(πR⁴) = {dp:.0f} Pa  ({dp / dp_heart * 100:.0f} % "
          f"de ΔP_cœur = {dp_heart:.0f} Pa)")
    # Re ≈ 2700 : on est pile sur la transition laminaire/turbulent. Le profil de
    # Poiseuille reste une approximation raisonnable, mais ce n'est plus un résultat
    # exact — et le débit de l'énoncé (80 cm³/s dans UNE artère) est volontairement
    # généreux : c'est l'ordre du débit cardiaque total.
    regime = ("laminaire" if re < RE_LAMINAR
              else "à la limite de la transition" if re < RE_TURBULENT else "turbulent")
    print(f"    vitesse débitante U = {u_mean * 1e2:.1f} cm/s,  Re = {re:.0f} → {regime}")

    r = np.linspace(-args.radius, args.radius, 400)
    v = poiseuille_profile(np.abs(r), radius=args.radius, pressure_drop=dp,
                           length=args.length, eta=args.eta_blood)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.fill_betweenx(r * 1e3, 0, v * 1e2, color=BLUE, alpha=0.18)
    ax1.plot(v * 1e2, r * 1e3, color=BLUE, lw=2.2)
    # quelques flèches pour matérialiser le champ de vitesse
    for rr in np.linspace(-args.radius, args.radius, 11)[1:-1]:
        vv = poiseuille_profile(np.array([abs(rr)]), radius=args.radius, pressure_drop=dp,
                                length=args.length, eta=args.eta_blood)[0]
        ax1.annotate("", xy=(vv * 1e2, rr * 1e3), xytext=(0, rr * 1e3),
                     arrowprops=dict(arrowstyle="->", color=BLUE, lw=1.1, alpha=0.75))
    ax1.axhline(args.radius * 1e3, color="0.3", lw=2)
    ax1.axhline(-args.radius * 1e3, color="0.3", lw=2)
    ax1.axvline(u_mean * 1e2, color=ACCENT, ls="--", lw=1.4,
                label=f"vitesse débitante $U$ = {u_mean * 1e2:.0f} cm/s")
    ax1.set_xlabel("vitesse $v$ (cm/s)")
    ax1.set_ylabel("$r$ (mm)")
    ax1.set_title(r"Profil de Poiseuille  $v(r)=\Delta P(R^2-r^2)/4\eta L$")
    ax1.legend(fontsize=8, loc="lower right")
    style_axes(ax1)

    # Budget de pression : ce que coûte l'artère par rapport à ce que fournit le cœur.
    ax2.barh(["ΔP fourni\npar le cœur\n(« 12–8 »)", "ΔP dissipé\ndans l'artère"],
             [dp_heart, dp], color=["0.75", BLUE], height=0.55)
    for y, val in enumerate([dp_heart, dp]):
        ax2.text(val * 1.02, y, f"{val:.0f} Pa\n({val / 133.3:.0f} mmHg)",
                 va="center", fontsize=9)
    ax2.set_xlim(0, dp_heart * 1.35)
    ax2.set_xlabel("pression (Pa)")
    ax2.set_title("Le cœur a de la marge… tant que le rayon ne bouge pas")
    style_axes(ax2, grid=False)

    fig.tight_layout()
    out = FIGURES / "artery.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")

    # ---- 2b. sténose : la loi en R⁻⁴ -------------------------------------- #
    x = np.linspace(0.0, 0.5, 400)
    ratio = stenosis_pressure_ratio(x)
    r5 = float(stenosis_pressure_ratio(args.stenosis))

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(x * 100, ratio, color=BLUE, lw=2.4)
    ax.axhline(1.0, color=GREY, lw=1.0, ls=":")
    ax.plot([args.stenosis * 100], [r5], "o", color=ACCENT, ms=8, zorder=5)
    ax.annotate(f"−{args.stenosis * 100:.0f} % de rayon\n→ +{(r5 - 1) * 100:.1f} % de pression",
                xy=(args.stenosis * 100, r5), xytext=(0.22, 0.55),
                textcoords="axes fraction", fontsize=10, color=ACCENT,
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.3))
    # repère : à −50 % de rayon il faut 16 fois plus de pression
    ax.plot([50], [16], "o", color="0.4", ms=6)
    ax.annotate("−50 % → ×16", xy=(50, 16), xytext=(41, 19), fontsize=9, color="0.35")
    ax.set_xlim(0, 50)
    ax.set_ylim(0.9, 20)
    ax.set_yscale("log")
    ax.set_xlabel("réduction relative du rayon (%)")
    ax.set_ylabel(r"$\Delta P / \Delta P_0$  (à débit constant)")
    ax.set_title(r"Sténose artérielle : le prix de la loi en $R^{-4}$")
    style_axes(ax)
    fig.tight_layout()
    out = FIGURES / "stenosis.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")
    print(f"    sténose de {args.stenosis * 100:.0f} % : ΔP/ΔP₀ = (1−x)⁻⁴ = {r5:.4f} "
          f"→ +{(r5 - 1) * 100:.1f} %")

    # ---- 2c. la citerne fêlée (MF3-8) ------------------------------------- #
    q_real = slit_flow_rate(width=args.slit_width, span=args.slit_span,
                            thickness=args.wall, head=args.head, eta=args.eta_water)
    q_ideal = torricelli_flow_rate(width=args.slit_width, span=args.slit_span,
                                   head=args.head)
    u_slit = q_real / (args.slit_width * args.slit_span)
    re_slit = reynolds(rho=1000.0, velocity=u_slit, diameter=args.slit_width,
                       eta=args.eta_water)

    print(f"\n  citerne fêlée : fissure {args.slit_span * 1e2:.0f} cm × "
          f"{args.slit_width * 1e6:.0f} µm, paroi {args.wall * 1e2:.0f} cm")
    print(f"    Dv (visqueux)       = {q_real * 1e9:.2f} mm³/s, soit {q_real * 86400e3:.2f} L/jour")
    print(f"    Dv (fluide parfait) = {q_ideal * 86400e3:.0f} L/jour  "
          f"(×{q_ideal / q_real:.0f} !)")
    print(f"    Re = {re_slit:.3f} → écoulement rampant, la viscosité domine tout")

    x_slit = np.linspace(-args.slit_width / 2, args.slit_width / 2, 300)
    v_slit = slit_profile(x_slit, width=args.slit_width, thickness=args.wall,
                          head=args.head, eta=args.eta_water)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ax1.fill_between(x_slit * 1e6, 0, v_slit * 1e3, color=BLUE, alpha=0.18)
    ax1.plot(x_slit * 1e6, v_slit * 1e3, color=BLUE, lw=2.2)
    ax1.axvline(-args.slit_width / 2 * 1e6, color="0.3", lw=2)
    ax1.axvline(args.slit_width / 2 * 1e6, color="0.3", lw=2)
    ax1.set_xlabel(r"$x$ (µm)   —   largeur de la fissure $b$")
    ax1.set_ylabel("vitesse $v$ (mm/s)")
    ax1.set_title(r"Fissure : $v(x)=\frac{\rho g h}{2\eta e}\left(\frac{b^2}{4}-x^2\right)$")
    style_axes(ax1)

    ax2.bar(["fluide réel\n(visqueux)", "fluide parfait\n(Torricelli)"],
            [q_real * 86400e3, q_ideal * 86400e3], color=[BLUE, "0.75"], width=0.5)
    ax2.set_yscale("log")
    for i, val in enumerate([q_real * 86400e3, q_ideal * 86400e3]):
        ax2.text(i, val * 1.4, f"{val:.3g} L/jour", ha="center", fontsize=9)
    ax2.set_ylabel("fuite (L/jour, échelle log)")
    ax2.set_title(f"La viscosité freine la fuite d'un facteur {q_ideal / q_real:.0f}")
    ax2.set_ylim(top=q_ideal * 86400e3 * 8)
    style_axes(ax2, grid=False)

    fig.tight_layout()
    out = FIGURES / "slit_leak.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")


# --------------------------------------------------------------------------- #
# Mode 3 — pertes de charge : le diagramme de Moody (MF3-6, MF3-9)
# --------------------------------------------------------------------------- #
def run_moody(args):
    # ---- 3a. le château d'eau : Poiseuille s'effondre ---------------------- #
    diameter = 2 * args.pipe_radius
    q_lam = poiseuille_flow_rate(radius=args.pipe_radius,
                                 pressure_drop=1000.0 * 9.81 * args.head_tower,
                                 length=args.pipe_length, eta=args.eta_water)
    u_lam = mean_velocity(flow_rate=q_lam, radius=args.pipe_radius)
    re_lam = reynolds(rho=1000.0, velocity=u_lam, diameter=diameter, eta=args.eta_water)

    print(f"  château d'eau : H = {args.head_tower:.0f} m, L = {args.pipe_length:.0f} m, "
          f"D = {diameter * 1e2:.0f} cm")
    print(f"    1) si l'écoulement était laminaire : Dv = {q_lam * 1e3:.0f} L/s "
          f"(U = {u_lam:.0f} m/s)")
    print(f"    2) Re = {re_lam:.2e} ≫ 2300  →  l'hypothèse laminaire est absurde")

    eps, lam_meas, re_meas = roughness_from_velocity(
        velocity=args.measured_velocity, head=args.head_tower, length=args.pipe_length,
        diameter=diameter, rho=1000.0, eta=args.eta_water,
    )
    print(f"    3) vitesse mesurée U = {args.measured_velocity:.1f} m/s → Re = {re_meas:.2e}, "
          f"λ = {lam_meas:.4f}")
    print(f"       rugosité extraite de Colebrook : ε = {eps * 1e3:.2f} mm "
          f"(énoncé, lu sur l'abaque : 0,24 mm)")
    q_real = args.measured_velocity * np.pi * args.pipe_radius**2
    print(f"       débit réel = {q_real * 1e3:.1f} L/s, soit {q_lam / q_real:.0f} fois "
          f"moins que la prédiction laminaire")

    # ---- 3b. l'oléoduc d'Alaska ------------------------------------------- #
    d_oil = args.oil_diameter
    u_oil = args.oil_flow / (np.pi * d_oil**2 / 4)
    re_oil = reynolds(rho=args.rho_oil, velocity=u_oil, diameter=d_oil, eta=args.eta_oil)
    dp_poiseuille = pressure_drop_for_flow(flow_rate=args.oil_flow, radius=d_oil / 2,
                                           length=1.0, eta=args.eta_oil)
    rr_oil = args.oil_roughness / d_oil
    lam_oil = float(friction_factor(re_oil, rr_oil))
    dp_moody = regular_head_loss(lam=lam_oil, length=1.0, diameter=d_oil,
                                 rho=args.rho_oil, velocity=u_oil)

    print(f"\n  oléoduc d'Alaska : D = {d_oil:.1f} m, Dv = {args.oil_flow * 1e3:.0f} L/s, "
          f"U = {u_oil:.2f} m/s")
    print(f"    1) Poiseuille donnerait ΔP/L = {dp_poiseuille:.1f} Pa/m, "
          f"mais Re = {re_oil:.0f} > 4000 → hypothèse non valable")
    print(f"    2) Moody (ε = {args.oil_roughness * 1e3:.1f} mm) : λ = {lam_oil:.4f}, "
          f"ΔP/L = {dp_moody:.1f} Pa/m  (×{dp_moody / dp_poiseuille:.1f})")

    k_bend = bend_loss_coefficient(diameter=d_oil, curvature_radius=args.bend_radius,
                                   angle_deg=args.bend_angle)
    dp_bend = singular_head_loss(k=k_bend, rho=args.rho_oil, velocity=u_oil)
    print(f"    3) coude de {args.bend_angle:.0f}° (R_c = {args.bend_radius:.0f} m) : "
          f"Λ = {k_bend:.4f} → ΔP = {dp_bend:.1f} Pa")
    print(f"       c'est une perte de charge SINGULIÈRE, équivalente à "
          f"{dp_bend / dp_moody:.2f} m de conduite droite.")

    # ---- 3c. le diagramme de Moody, recalculé ----------------------------- #
    fig, ax = plt.subplots(figsize=(9.5, 6))

    re_lam_range = np.logspace(2.6, np.log10(RE_LAMINAR), 100)
    ax.plot(re_lam_range, laminar_friction_factor(re_lam_range), color=ACCENT, lw=2.4)
    ax.annotate(r"laminaire : $\lambda = 64/Re$", xy=(700, 64 / 700), xytext=(450, 0.115),
                fontsize=10, color=ACCENT)

    re_turb = np.logspace(np.log10(RE_TURBULENT), 8.3, 400)
    roughnesses = [0.0, 1e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
    for rr in roughnesses:
        lam = colebrook_friction_factor(re_turb, max(rr, 1e-12))
        ax.plot(re_turb, lam, color=BLUE if rr else "0.25", lw=1.9 if rr == 0 else 1.1,
                ls="--" if rr == 0 else "-", alpha=0.9)
        if rr:
            ax.text(1.05 * re_turb[-1], lam[-1], f"{rr:g}", fontsize=7.5,
                    va="center", color="0.25")
    ax.text(1.05 * re_turb[-1], 0.135, r"$\varepsilon/D$", fontsize=10, va="center",
            color="0.25", fontweight="bold")
    # le tuyau parfaitement lisse est la borne inférieure : aucune rugosité ne peut
    # faire descendre λ en dessous de cette courbe
    ax.annotate("tuyau lisse ($\\varepsilon = 0$)",
                xy=(2e6, float(colebrook_friction_factor(2e6, 1e-12))),
                xytext=(6e5, 0.0072), fontsize=9, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.25", lw=1.0))

    # zone critique : ni laminaire ni pleinement turbulent, λ non prédictible
    ax.axvspan(RE_LAMINAR, RE_TURBULENT, color="0.85", alpha=0.6)
    ax.text(3000, 0.075, "zone\ncritique", fontsize=8, ha="center", color="0.4")

    # nos deux cas d'étude, placés sur l'abaque
    ax.plot([re_meas], [lam_meas], "o", color=ACCENT, ms=9, zorder=6,
            markeredgecolor="white", markeredgewidth=1.2)
    ax.annotate(f"château d'eau\nRe = {re_meas:.1e}, λ = {lam_meas:.3f}\n→ ε = {eps * 1e3:.2f} mm",
                xy=(re_meas, lam_meas), xytext=(2.2e5, 0.062), fontsize=9, color=ACCENT,
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.2))
    ax.plot([re_oil], [lam_oil], "s", color="#1e8449", ms=8, zorder=6,
            markeredgecolor="white", markeredgewidth=1.2)
    ax.annotate(f"oléoduc d'Alaska\nRe = {re_oil:.0f}, λ = {lam_oil:.3f}",
                xy=(re_oil, lam_oil), xytext=(1.1e4, 0.058), fontsize=9, color="#1e8449",
                arrowprops=dict(arrowstyle="->", color="#1e8449", lw=1.2))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(400, 1e8)
    ax.set_ylim(0.0062, 0.17)
    ax.set_xlabel(r"nombre de Reynolds  $Re = \rho U D / \eta$")
    ax.set_ylabel(r"coefficient de frottement de Darcy  $\lambda$")
    ax.set_title("Diagramme de Moody — recalculé en résolvant Colebrook–White par Newton")
    ax.grid(which="both", alpha=0.22, linewidth=0.5)
    ax.set_yticks([0.008, 0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08, 0.1, 0.15])
    ax.set_yticklabels(["0.008", "0.01", "0.015", "0.02", "0.03", "0.04", "0.06",
                        "0.08", "0.10", "0.15"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    out = FIGURES / "moody_diagram.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\n  → {out.name}")

    # ---- 3d. le coût réel de la turbulence -------------------------------- #
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.bar(["Poiseuille\n(laminaire)", "réalité\n(turbulent)"],
            [q_lam * 1e3, q_real * 1e3], color=["0.75", BLUE], width=0.5)
    for i, val in enumerate([q_lam * 1e3, q_real * 1e3]):
        ax1.text(i, val + 3, f"{val:.1f} L/s", ha="center", fontsize=10)
    ax1.set_ylabel("débit au robinet (L/s)")
    ax1.set_title(f"Château d'eau : la turbulence divise\nle débit par {q_lam / q_real:.0f}")
    style_axes(ax1, grid=False)

    # Balayage : ce que coûte la rugosité sur l'oléoduc, à débit fixé
    eps_range = np.logspace(-6, -2.3, 200)
    lam_range = np.array([float(friction_factor(re_oil, e / d_oil)) for e in eps_range])
    dp_range = np.array([regular_head_loss(lam=l, length=1.0, diameter=d_oil,
                                           rho=args.rho_oil, velocity=u_oil)
                         for l in lam_range])
    ax2.plot(eps_range * 1e3, dp_range, color=BLUE, lw=2.2)
    ax2.axhline(dp_poiseuille, color=GREY, ls=":", lw=1.4)
    ax2.text(2e-3, dp_poiseuille * 1.06, "ce que prédirait Poiseuille", fontsize=8.5,
             color=GREY)
    ax2.plot([args.oil_roughness * 1e3], [dp_moody], "o", color=ACCENT, ms=8, zorder=5)
    ax2.annotate(f"aspérités de {args.oil_roughness * 1e3:.1f} mm\n→ {dp_moody:.0f} Pa/m",
                 xy=(args.oil_roughness * 1e3, dp_moody), xytext=(0.12, 0.3),
                 textcoords="axes fraction", fontsize=9, color=ACCENT,
                 arrowprops=dict(arrowstyle="->", color=ACCENT, lw=1.2))
    ax2.set_xscale("log")
    ax2.set_xlabel("rugosité $\\varepsilon$ (mm)")
    ax2.set_ylabel(r"$\Delta P / L$ (Pa/m)")
    ax2.set_title("Oléoduc : coût de la rugosité à débit imposé")
    style_axes(ax2)

    fig.tight_layout()
    out = FIGURES / "pipe_cases.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  → {out.name}")


# --------------------------------------------------------------------------- #
# Mode 4 — vérifications numériques
# --------------------------------------------------------------------------- #
def run_check(args):
    """Trois validations indépendantes des solveurs du dépôt."""
    nu = args.eta / args.rho
    tau = diffusion_time(args.gap, nu=nu)

    z, times, V = solve_startup(gap=args.gap, n_points=args.n_points, nu=nu,
                                plate_speed=args.plate_speed, t_max=args.t_ratio * tau,
                                n_steps=args.n_steps)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))

    # ---- (1) Crank–Nicolson vs solution analytique en erfc ---------------- #
    # Valable tant que la couche limite n'a pas atteint la plaque du haut.
    errors = []
    for ratio, color in zip((0.002, 0.005, 0.01, 0.02),
                            plt.cm.viridis(np.linspace(0.15, 0.8, 4))):
        k = int(ratio * (len(times) - 1) / args.t_ratio)
        t = times[k]
        exact = stokes_first_problem(z, t, plate_speed=args.plate_speed, nu=nu)
        err = np.sqrt(np.mean((V[k] - exact) ** 2)) / args.plate_speed
        errors.append((ratio, err, boundary_layer_thickness(t, nu=nu)))
        axes[0].plot(V[k], z * 1e2, color=color, lw=2.0, label=f"CN, $t={ratio:g}\\tau$")
        axes[0].plot(exact, z * 1e2, color="0.25", ls="--", lw=1.0)
    axes[0].set_xlabel("vitesse $v$ (m/s)")
    axes[0].set_ylabel("$z$ (cm)")
    axes[0].set_ylim(0, args.gap * 1e2 * 0.5)
    axes[0].set_title("(1) Crank–Nicolson vs erfc analytique\n(traits pointillés : exact)")
    axes[0].legend(fontsize=8)
    style_axes(axes[0])

    print("  (1) transitoire de Stokes — CN contre la solution exacte en erfc :")
    for ratio, err, d in errors:
        print(f"      t = {ratio:>6g} τ  (δ = {d * 1e3:5.1f} mm)  →  RMS = {err * 100:.4f} % de U")

    # ---- (2) convergence vers le profil de Couette ------------------------ #
    steady = couette_steady(z, gap=args.gap, plate_speed=args.plate_speed)
    residuals = np.sqrt(np.mean((V - steady) ** 2, axis=1)) / args.plate_speed
    axes[1].semilogy(times / tau, residuals, color=BLUE, lw=2.0)
    # la théorie prédit une décroissance en exp(-π²νt/L²) = exp(-π² t/τ) : le mode
    # fondamental de la boîte survit seul aux temps longs
    mask = times / tau > 0.05
    theory = residuals[mask][0] * np.exp(-np.pi**2 * (times[mask] - times[mask][0]) / tau)
    axes[1].semilogy(times[mask] / tau, theory, color=ACCENT, ls="--", lw=1.6,
                     label=r"$\propto e^{-\pi^2 t/\tau}$ (mode fondamental)")
    axes[1].set_xlabel(r"$t / \tau$")
    axes[1].set_ylabel("écart RMS au profil de Couette")
    axes[1].set_title("(2) Convergence vers le régime établi")
    axes[1].legend(fontsize=8)
    style_axes(axes[1])

    tail = residuals[times / tau > 0.3]
    rate = np.polyfit(times[times / tau > 0.3], np.log(tail), 1)[0] * tau
    print(f"\n  (2) taux de décroissance mesuré : {rate:.3f} / τ  "
          f"(théorie : -π² = {-np.pi**2:.3f})")
    print(f"      écart final au profil de Couette : {residuals[-1] * 100:.3f} % de U")

    # ---- (3) Colebrook : résidu de Newton et comparaison à Haaland -------- #
    re_test = np.logspace(3.7, 8, 300)
    rr_test = 1e-4
    lam_newton = colebrook_friction_factor(re_test, rr_test)
    lam_haaland = haaland_friction_factor(re_test, rr_test)
    # résidu de l'équation implicite : doit être à la précision machine
    x = 1.0 / np.sqrt(lam_newton)
    residual = np.abs(x + 2.0 * np.log10(rr_test / 3.7 + 2.51 * x / re_test))
    axes[2].loglog(re_test, np.maximum(residual, 1e-17), color=BLUE, lw=1.8,
                   label="résidu de Colebrook (Newton)")
    axes[2].loglog(re_test, np.abs(lam_haaland - lam_newton) / lam_newton, color=ACCENT,
                   lw=1.8, ls="--", label="écart relatif de Haaland")
    axes[2].set_xlabel("$Re$")
    axes[2].set_ylabel("résidu / écart relatif")
    axes[2].set_ylim(1e-17, 1e-1)
    axes[2].set_title(r"(3) Colebrook–White résolu par Newton ($\varepsilon/D=10^{-4}$)")
    axes[2].legend(fontsize=8, loc="center right")
    style_axes(axes[2])

    print(f"\n  (3) Colebrook résolu par Newton : résidu max = {residual.max():.2e} "
          "(précision machine)")
    print(f"      l'approximation explicite de Haaland s'en écarte de "
          f"{np.max(np.abs(lam_haaland - lam_newton) / lam_newton) * 100:.2f} % au pire")

    fig.tight_layout()
    out = FIGURES / "check.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\n  → {out.name}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Écoulements visqueux — diffusion de quantité de mouvement, "
                    "Poiseuille, pertes de charge.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--mode", choices=("stokes", "poiseuille", "moody", "check"), default="stokes",
        help="démonstration à lancer",
    )

    st = p.add_argument_group("plaque en mouvement / Couette (MF3-1, MF3-2, MF3-4)")
    st.add_argument("--gap", type=float, default=0.1, help="entrefer entre les plaques L (m)")
    st.add_argument("--rho", type=float, default=1000.0, help="masse volumique ρ (kg/m³)")
    st.add_argument("--eta", type=float, default=1e-3, help="viscosité dynamique η (Pa·s)")
    st.add_argument("--plate-speed", type=float, default=1.0, help="vitesse de la plaque U (m/s)")
    st.add_argument("--t-ratio", type=float, default=1.0, help="durée simulée (en τ = L²/ν)")
    st.add_argument("--n-points", type=int, default=240, help="points de discrétisation en z")
    st.add_argument("--n-steps", type=int, default=3000, help="pas de temps du schéma")
    st.add_argument("--frames", type=int, default=140, help="images de l'animation")
    st.add_argument("--fps", type=int, default=25, help="images/s du GIF")
    st.add_argument("--r-inner", type=float, default=0.04, help="rayon du cylindre intérieur R₁ (m)")
    st.add_argument("--r-outer", type=float, default=0.05, help="rayon du cylindre extérieur R₂ (m)")
    st.add_argument("--height", type=float, default=0.1, help="hauteur des cylindres H (m)")
    st.add_argument("--omega", type=float, default=10.0, help="vitesse angulaire ω (rad/s)")
    st.add_argument("--eta-couette", type=float, default=0.85,
                    help="viscosité du fluide testé au viscosimètre (Pa·s)")

    po = p.add_argument_group("écoulements internes (MF3-7, MF3-8)")
    po.add_argument("--radius", type=float, default=0.005, help="rayon de l'artère R (m)")
    po.add_argument("--length", type=float, default=1.0, help="longueur de l'artère L (m)")
    po.add_argument("--flow-rate", type=float, default=80e-6, help="débit sanguin Dv (m³/s)")
    po.add_argument("--eta-blood", type=float, default=4e-3, help="viscosité du sang (Pa·s)")
    po.add_argument("--rho-blood", type=float, default=1060.0, help="masse volumique du sang (kg/m³)")
    po.add_argument("--heart-pressure", type=float, default=0.04,
                    help="ΔP du cœur, en mètres de mercure (« 12–8 » = 4 cmHg)")
    po.add_argument("--stenosis", type=float, default=0.05, help="réduction relative du rayon")
    po.add_argument("--slit-width", type=float, default=10e-6, help="largeur de la fissure b (m)")
    po.add_argument("--slit-span", type=float, default=0.05, help="longueur de la fissure a (m)")
    po.add_argument("--wall", type=float, default=0.02, help="épaisseur de la paroi e (m)")
    po.add_argument("--head", type=float, default=1.0, help="hauteur d'eau dans la citerne h (m)")
    po.add_argument("--eta-water", type=float, default=1e-3, help="viscosité de l'eau (Pa·s)")

    mo = p.add_argument_group("pertes de charge (MF3-6, MF3-9)")
    mo.add_argument("--head-tower", type=float, default=60.0, help="hauteur du château d'eau (m)")
    mo.add_argument("--pipe-length", type=float, default=100.0, help="longueur de la conduite (m)")
    mo.add_argument("--pipe-radius", type=float, default=0.015, help="rayon de la conduite (m)")
    mo.add_argument("--measured-velocity", type=float, default=3.2,
                    help="vitesse débitante mesurée au robinet (m/s)")
    mo.add_argument("--oil-diameter", type=float, default=1.2, help="diamètre de l'oléoduc (m)")
    mo.add_argument("--oil-flow", type=float, default=3.4, help="débit de pétrole (m³/s)")
    mo.add_argument("--rho-oil", type=float, default=800.0, help="masse volumique du pétrole (kg/m³)")
    mo.add_argument("--eta-oil", type=float, default=0.3, help="viscosité du pétrole (Pa·s)")
    mo.add_argument("--oil-roughness", type=float, default=2e-4, help="rugosité des parois ε (m)")
    mo.add_argument("--bend-radius", type=float, default=10.0, help="rayon de courbure d'un coude (m)")
    mo.add_argument("--bend-angle", type=float, default=10.0, help="angle du coude (degrés)")

    return p


def main():
    args = build_parser().parse_args()
    FIGURES.mkdir(exist_ok=True)

    print(f"Mode : {args.mode}")
    {"stokes": run_stokes, "poiseuille": run_poiseuille,
     "moody": run_moody, "check": run_check}[args.mode](args)
    print("Terminé.")


if __name__ == "__main__":
    main()
