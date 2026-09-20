"""
Diffusion de quantité de mouvement : problème de Stokes et écoulements de Couette.

Un fluide visqueux au repos entre deux plaques ; à t = 0 on met la plaque du bas
en mouvement à la vitesse U. Le champ de vitesse v(z, t)·ex obéit à

    ∂v/∂t = ν ∂²v/∂z²        avec  ν = η/ρ  (viscosité cinématique)

C'est **exactement** une équation de diffusion : la viscosité transporte la
quantité de mouvement de proche en proche, comme la conduction transporte la
chaleur. D'où le vocabulaire commun (temps de diffusion L²/ν, épaisseur √(νt)).

Deux régimes, deux méthodes :
- aux temps courts la plaque du haut n'a « rien vu » : le fluide se comporte comme
  s'il était semi-infini et on a la solution analytique exacte en erfc
  (premier problème de Stokes) ;
- aux temps longs on tend vers le profil de Couette linéaire, et il faut un schéma
  numérique pour décrire tout le transitoire → Crank–Nicolson.

Exercices MF3-1, MF3-2 et MF3-4 du TD « Fluides réels en écoulement » (PSI).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import solve_banded
from scipy.special import erfc

Array = NDArray[np.float64]


# --------------------------------------------------------------------------- #
# 1. Premier problème de Stokes — solution analytique (milieu semi-infini)
# --------------------------------------------------------------------------- #
def stokes_first_problem(z: Array, t: float, *, plate_speed: float, nu: float) -> Array:
    """
    Profil de vitesse au-dessus d'une plaque infinie brusquement mise en mouvement.

        v(z, t) = U · erfc( z / (2√(νt)) )

    Le fluide est supposé semi-infini (z > 0). La solution ne dépend de z et t que
    par la variable de similitude ζ = z/(2√(νt)) : tous les profils se superposent
    une fois tracés en fonction de ζ — c'est la signature d'une diffusion pure.
    """
    if t <= 0.0:
        return np.zeros_like(z)
    return plate_speed * erfc(z / (2.0 * np.sqrt(nu * t)))


def similarity_variable(z: Array, t: float, *, nu: float) -> Array:
    """Variable de similitude ζ = z/(2√(νt)) du problème de Stokes."""
    return z / (2.0 * np.sqrt(nu * t))


def boundary_layer_thickness(t: float, *, nu: float) -> float:
    """
    Épaisseur de la couche limite visqueuse δ(t) ≈ √(νt).

    C'est la distance sur laquelle la plaque s'est « fait sentir » après un temps t.
    Convention usuelle : à z = 4√(νt) la vitesse est retombée à 0,5 % de U.
    """
    return float(np.sqrt(nu * t))


def diffusion_time(length: float, *, nu: float) -> float:
    """
    Temps caractéristique d'établissement du régime stationnaire : τ = L²/ν.

    Question 2 de MF3-1 : avec ν = 10⁻⁶ m²/s et L = 10 cm, τ = 10⁴ s ≈ 2 h 45.
    La quantité de mouvement diffuse *très* lentement — d'où l'intérêt de brasser.
    """
    return length**2 / nu


# --------------------------------------------------------------------------- #
# 2. Transitoire complet entre deux plaques — Crank–Nicolson
# --------------------------------------------------------------------------- #
def solve_startup(
    *,
    gap: float,
    n_points: int,
    nu: float,
    plate_speed: float,
    t_max: float,
    n_steps: int = 2000,
    n_rannacher: int = 8,
) -> tuple[Array, Array, Array]:
    """
    Intègre ∂v/∂t = ν ∂²v/∂z² sur [0, L] et renvoie (z, times, V).

    Conditions : v(0, t) = U (plaque du bas entraînée), v(L, t) = 0 (plaque du haut
    fixe), v(z, 0) = 0 (fluide au repos). `V[k]` est le profil complet à times[k].

    Schéma de **Crank–Nicolson** : moyenne du laplacien entre l'instant n et n+1.
    D'ordre 2 en temps comme en espace et *inconditionnellement stable*, alors qu'un
    schéma explicite imposerait ν·dt/dz² ≤ 1/2 — rédhibitoire ici, puisqu'il faut
    simuler jusqu'à τ = L²/ν. Le système est tridiagonal : O(n) via `solve_banded`.

    ⚠️ Subtilité : la condition initiale est **discontinue** (le fluide est au repos
    mais la plaque démarre à U). Crank–Nicolson, qui n'amortit pas les modes de haute
    fréquence, fait alors apparaître des oscillations parasites près de la paroi. On
    applique le remède classique de **Rannacher** : les `n_rannacher` premiers pas
    sont faits en Euler implicite, fortement dissipatif, qui tue ces modes ; on bascule
    ensuite sur Crank–Nicolson et on récupère l'ordre 2 pour tout le reste du calcul.
    """
    z = np.linspace(0.0, gap, n_points)
    dz = z[1] - z[0]
    times = np.linspace(0.0, t_max, n_steps + 1)
    dt = times[1] - times[0]

    def banded_operator(coeff: float) -> Array:
        """Matrice tridiagonale (I − coeff·Δ) en stockage bande, Dirichlet aux bords."""
        ab = np.zeros((3, n_points))
        ab[0, 2:] = -coeff  # sur-diagonale
        ab[1, :] = 1.0 + 2.0 * coeff  # diagonale
        ab[2, :-2] = -coeff  # sous-diagonale
        ab[1, 0] = ab[1, -1] = 1.0  # lignes des conditions aux limites
        ab[0, 1] = ab[2, -2] = 0.0
        return ab

    r_cn = nu * dt / (2.0 * dz**2)  # le « demi » vient de la moyenne de Crank–Nicolson
    r_be = nu * dt / dz**2  # Euler implicite : pas de moyenne
    ab_cn, ab_be = banded_operator(r_cn), banded_operator(r_be)

    V = np.zeros((n_steps + 1, n_points))
    v = np.zeros(n_points)
    v[0] = plate_speed  # la plaque démarre à t = 0⁺
    V[0] = v

    for k in range(1, n_steps + 1):
        if k <= n_rannacher:
            rhs = v.copy()  # Euler implicite : le second membre est v^n tel quel
            ab = ab_be
        else:
            rhs = v.copy()  # Crank–Nicolson : second membre explicite (I + r·Δ) v^n
            rhs[1:-1] = v[1:-1] + r_cn * (v[2:] - 2.0 * v[1:-1] + v[:-2])
            ab = ab_cn
        rhs[0], rhs[-1] = plate_speed, 0.0  # conditions aux limites
        v = solve_banded((1, 1), ab, rhs)
        V[k] = v

    return z, times, V


def couette_steady(z: Array, *, gap: float, plate_speed: float) -> Array:
    """
    Profil de Couette plan stationnaire : v(z) = U (1 − z/L).

    Question 3 de MF3-1 : en régime permanent ∂v/∂t = 0 donc v'' = 0, la vitesse est
    affine. C'est l'état vers lequel converge `solve_startup`.
    """
    return plate_speed * (1.0 - z / gap)


def wall_shear_stress(*, eta: float, plate_speed: float, gap: float) -> float:
    """
    Contrainte tangentielle sur la plaque en régime de Couette : τ = η U / L.

    Question 4 de MF3-1. Le fluide freine la plaque : la force surfacique vue par la
    plaque vaut −τ·ex. C'est aussi la base du viscosimètre plan (MF3-4).
    """
    return eta * plate_speed / gap


def viscosity_from_sliding_block(
    *, mass: float, area: float, angle: float, speed: float, film: float, g: float = 9.81
) -> float:
    """
    Exercice MF3-4 : un solide glisse à vitesse constante sur un film d'huile.

    À vitesse constante, la composante du poids le long de la pente équilibre
    exactement le frottement visqueux du film (profil de Couette d'épaisseur e) :

        M g sin α = η S v₀ / e      ⇒      η = M g sin α · e / (S v₀)

    `angle` est en radians.
    """
    return mass * g * np.sin(angle) * film / (area * speed)


# --------------------------------------------------------------------------- #
# 3. Viscosimètre de Couette cylindrique (MF3-2)
# --------------------------------------------------------------------------- #
def couette_cylindrical_constants(
    *, r_inner: float, r_outer: float, omega: float
) -> tuple[float, float]:
    """
    Constantes (A, B) du profil orthoradial v(r) = A r + B/r entre deux cylindres.

    Le cylindre intérieur (rayon R₁) est immobile, l'extérieur (R₂) tourne à ω.
    Les conditions d'adhérence v(R₁) = 0 et v(R₂) = ω R₂ donnent

        A = ω R₂² / (R₂² − R₁²),        B = −ω R₁² R₂² / (R₂² − R₁²)

    Le terme `A r` est une rotation en bloc (pas de cisaillement), le terme `B/r`
    est le tourbillon libre. B = 0 correspondrait au fluide tournant *solidement*
    avec les deux cylindres : aucun cisaillement, donc aucun couple mesurable.
    """
    denom = r_outer**2 - r_inner**2
    a = omega * r_outer**2 / denom
    b = -omega * r_inner**2 * r_outer**2 / denom
    return a, b


def couette_cylindrical_profile(
    r: Array, *, r_inner: float, r_outer: float, omega: float
) -> Array:
    """Profil de vitesse orthoradiale v(r) = A r + B/r entre les deux cylindres."""
    a, b = couette_cylindrical_constants(r_inner=r_inner, r_outer=r_outer, omega=omega)
    return a * r + b / r


def couette_torque(
    *, eta: float, height: float, r_inner: float, r_outer: float, omega: float
) -> float:
    """
    Couple visqueux transmis d'un cylindre à l'autre (question 3 de MF3-2).

    La force surfacique vaut fs = η (dv/dr − v/r) = −2ηB/r² : le terme de rotation
    en bloc disparaît, seul le cisaillement compte. Le couple sur un cylindre de
    rayon r vaut Γ = fs · (2πrH) · r = 4πηH|B|, **indépendant de r** — c'est la
    conservation du moment cinétique en régime stationnaire, et c'est ce qui rend la
    mesure fiable : le couple lu sur le fil de torsion est celui imposé à l'extérieur.

        Γ = 4πηH · ω R₁²R₂² / (R₂² − R₁²)
    """
    _, b = couette_cylindrical_constants(r_inner=r_inner, r_outer=r_outer, omega=omega)
    return 4.0 * np.pi * eta * height * abs(b)


def viscosity_from_torque(
    *, torque: float, height: float, r_inner: float, r_outer: float, omega: float
) -> float:
    """
    Inversion du viscosimètre (question 4 de MF3-2) : on remonte à η depuis Γ.

    C'est la raison d'être de l'appareil : Γ est proportionnel à η, avec un facteur
    purement géométrique connu. La torsion α du fil donne Γ = C α, d'où η.
    """
    _, b = couette_cylindrical_constants(r_inner=r_inner, r_outer=r_outer, omega=omega)
    return torque / (4.0 * np.pi * height * abs(b))
