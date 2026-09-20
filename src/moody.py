"""
Pertes de charge : nombre de Reynolds, Colebrook–White et diagramme de Moody.

Poiseuille suppose l'écoulement laminaire. Dès que le nombre de Reynolds

    Re = ρ U D / η

dépasse ~2300, l'écoulement transite vers la turbulence et la loi en R⁴ s'effondre :
la dissipation devient bien plus forte que prévu. On abandonne alors le calcul exact
au profit d'une formulation empirique, la **perte de charge régulière** :

    ΔP = λ · (L/D) · ½ρU²

où le coefficient de frottement λ (Darcy–Weisbach) vaut 64/Re en laminaire, et est
donné en turbulent par l'équation **implicite** de Colebrook–White :

    1/√λ = −2 log₁₀( ε/(3,7 D) + 2,51/(Re √λ) )

Le diagramme de Moody n'est rien d'autre que l'abaque de cette équation. Plutôt que
de le lire, on le **recalcule** ici en résolvant Colebrook par la méthode de Newton.

Exercices MF3-6 (château d'eau) et MF3-9 (oléoduc d'Alaska).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]

RE_LAMINAR = 2300.0  # fin du régime laminaire (valeur usuelle)
RE_TURBULENT = 4000.0  # début du régime pleinement turbulent


def reynolds(*, rho: float, velocity: float, diameter: float, eta: float) -> float:
    """
    Nombre de Reynolds Re = ρUD/η : rapport des effets inertiels aux effets visqueux.

    C'est le seul paramètre sans dimension du problème (lois de similitude, MF3-12) :
    deux écoulements de même Re dans des conduites géométriquement semblables ont le
    même comportement, quelles que soient leurs tailles réelles.
    """
    return rho * velocity * diameter / eta


def laminar_friction_factor(re: float | Array) -> float | Array:
    """
    Coefficient de Darcy en régime laminaire : λ = 64/Re.

    Ce n'est pas une corrélation empirique mais la loi de Poiseuille réécrite :
    en injectant ΔP = 8ηLDv/(πR⁴) dans ΔP = λ(L/D)·½ρU², on tombe exactement sur 64/Re.
    """
    return 64.0 / re


def haaland_friction_factor(re: float | Array, rel_roughness: float | Array) -> float | Array:
    """
    Approximation **explicite** de Haaland — sert d'initialisation à Newton.

        1/√λ = −1,8 log₁₀[ (ε/3,7D)^1,11 + 6,9/Re ]

    Précise à ~2 % près, ce qui en fait un point de départ idéal : Newton converge
    alors en deux ou trois itérations.
    """
    inv_sqrt = -1.8 * np.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / re)
    return 1.0 / inv_sqrt**2


def colebrook_friction_factor(
    re: float | Array,
    rel_roughness: float | Array,
    *,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> float | Array:
    """
    Résout l'équation implicite de Colebrook–White par la méthode de Newton.

    On pose x = 1/√λ, ce qui rend l'équation nettement mieux conditionnée :

        F(x) = x + 2 log₁₀( ε/3,7D + 2,51 x/Re ) = 0
        F'(x) = 1 + (2/ln10) · (2,51/Re) / ( ε/3,7D + 2,51 x/Re )

    F est croissante et convexe, donc Newton initialisé par Haaland converge de façon
    monotone. L'implémentation est vectorisée : on peut passer des tableaux de Re et
    de rugosité, ce qui permet de tracer tout le diagramme de Moody d'un coup.
    """
    re = np.asarray(re, dtype=float)
    rr = np.asarray(rel_roughness, dtype=float)
    re, rr = np.broadcast_arrays(re, rr)

    x = 1.0 / np.sqrt(haaland_friction_factor(re, rr))  # initialisation
    ln10 = np.log(10.0)

    for _ in range(max_iter):
        inner = rr / 3.7 + 2.51 * x / re
        f = x + 2.0 * np.log10(inner)
        df = 1.0 + (2.0 / ln10) * (2.51 / re) / inner
        step = f / df
        x = x - step
        if np.all(np.abs(step) < tol):
            break

    lam = 1.0 / x**2
    return float(lam) if lam.ndim == 0 else lam


def friction_factor(re: float | Array, rel_roughness: float | Array) -> float | Array:
    """
    Coefficient de frottement λ sur toute la gamme de Re.

    - Re < 2300  : laminaire, λ = 64/Re (exact).
    - Re > 4000  : turbulent, Colebrook–White.
    - entre les deux : zone critique, où l'écoulement est instable et λ mal défini ;
      on interpole linéairement pour ne pas laisser de trou dans le tracé, mais
      **aucune valeur n'y est fiable** — c'est une zone qu'un ingénieur évite.
    """
    scalar = np.ndim(re) == 0 and np.ndim(rel_roughness) == 0
    re = np.atleast_1d(np.asarray(re, dtype=float))
    rr = np.atleast_1d(np.asarray(rel_roughness, dtype=float))
    re, rr = np.broadcast_arrays(re, rr)

    lam = np.empty_like(re)
    lam_lam = laminar_friction_factor(np.clip(re, 1e-6, None))
    lam_turb = colebrook_friction_factor(np.clip(re, RE_LAMINAR, None), rr)

    laminar = re <= RE_LAMINAR
    turbulent = re >= RE_TURBULENT
    critical = ~laminar & ~turbulent

    lam[laminar] = np.atleast_1d(lam_lam)[laminar]
    lam[turbulent] = np.atleast_1d(lam_turb)[turbulent]
    if np.any(critical):
        # interpolation linéaire en log(Re) entre les deux bornes de la zone critique
        a = laminar_friction_factor(RE_LAMINAR)
        b = colebrook_friction_factor(RE_TURBULENT, rr[critical])
        w = (np.log(re[critical]) - np.log(RE_LAMINAR)) / (
            np.log(RE_TURBULENT) - np.log(RE_LAMINAR)
        )
        lam[critical] = (1.0 - w) * a + w * b

    return float(lam[0]) if scalar else lam


def regular_head_loss(
    *, lam: float, length: float, diameter: float, rho: float, velocity: float
) -> float:
    """Perte de charge régulière ΔP = λ (L/D) ½ρU² (Darcy–Weisbach)."""
    return lam * (length / diameter) * 0.5 * rho * velocity**2


def singular_head_loss(*, k: float, rho: float, velocity: float) -> float:
    """
    Perte de charge **singulière** ΔP = K ½ρU².

    Contrairement à la perte régulière, elle n'est pas proportionnelle à la longueur :
    elle est localisée sur un accident de la conduite (coude, vanne, élargissement),
    où l'écoulement décolle et dissipe de l'énergie dans des tourbillons.
    """
    return k * 0.5 * rho * velocity**2


def bend_loss_coefficient(*, diameter: float, curvature_radius: float, angle_deg: float) -> float:
    """
    Coefficient de perte singulière d'un coude arrondi (corrélation de MF3-9) :

        Λ = [ 0,13 + 1,85 (D/2R_c)^3,5 ] · (θ/90)

    Le terme en (D/2R_c)^3,5 s'effondre dès que le coude est large : c'est la
    justification quantitative de la règle du métier « un coude doux coûte
    beaucoup moins cher qu'un coude sec ».
    """
    return (0.13 + 1.85 * (diameter / 2.0 / curvature_radius) ** 3.5) * (angle_deg / 90.0)


def velocity_from_head(
    *,
    head: float,
    length: float,
    diameter: float,
    rho: float,
    eta: float,
    roughness: float,
    include_kinetic: bool = True,
    tol: float = 1e-10,
    max_iter: int = 200,
) -> float:
    """
    Vitesse débitante réelle dans une conduite alimentée par une hauteur d'eau `head`.

    Le bilan de charge entre la surface libre et la sortie à l'air libre s'écrit

        ρ g H = λ (L/D) ½ρU²  [+ ½ρU² si on compte l'énergie cinétique emportée]

    λ dépendant lui-même de U via Re, c'est un point fixe : on itère U → Re → λ → U
    jusqu'à convergence (quelques dizaines d'itérations, amorties pour rester stable).
    """
    g = 9.81
    u = np.sqrt(2.0 * g * head)  # départ : estimation fluide parfait (Torricelli)

    for _ in range(max_iter):
        re = reynolds(rho=rho, velocity=u, diameter=diameter, eta=eta)
        lam = float(friction_factor(re, roughness / diameter))
        coeff = lam * length / diameter + (1.0 if include_kinetic else 0.0)
        u_new = np.sqrt(2.0 * g * head / coeff)
        u, delta = 0.5 * (u + u_new), abs(u_new - u)  # sous-relaxation
        if delta < tol:
            break

    return float(u)


def roughness_from_velocity(
    *,
    velocity: float,
    head: float,
    length: float,
    diameter: float,
    rho: float,
    eta: float,
    include_kinetic: bool = True,
) -> tuple[float, float, float]:
    """
    Problème **inverse** : on mesure U, on en déduit la rugosité ε de la conduite.

    Renvoie (ε, λ, Re). C'est la question 3 de MF3-6 : le bilan de charge donne
    directement λ, puis on retourne Colebrook — qui est explicite dans ce sens-là —
    pour extraire ε :

        ε = 3,7 D [ 10^(−1/(2√λ)) − 2,51/(Re√λ) ]

    Une mesure de débit suffit donc à caractériser l'état de surface d'un tuyau
    enterré qu'on ne peut pas inspecter.
    """
    g = 9.81
    coeff = 2.0 * g * head / velocity**2 - (1.0 if include_kinetic else 0.0)
    lam = coeff * diameter / length
    re = reynolds(rho=rho, velocity=velocity, diameter=diameter, eta=eta)

    inv_sqrt = 1.0 / np.sqrt(lam)
    roughness = 3.7 * diameter * (10.0 ** (-inv_sqrt / 2.0) - 2.51 / (re * np.sqrt(lam)))
    return float(roughness), float(lam), float(re)
