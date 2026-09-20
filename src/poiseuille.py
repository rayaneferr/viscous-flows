"""
Écoulements internes établis : loi de Hagen–Poiseuille.

Un fluide visqueux poussé dans une conduite par un gradient de pression atteint,
loin de l'entrée, un profil stationnaire où la poussée de pression équilibre
exactement le cisaillement visqueux. En conduite cylindrique :

    v(r) = ΔP (R² − r²) / (4ηL)          (profil parabolique)
    Dv   = π R⁴ ΔP / (8ηL)               (loi de Hagen–Poiseuille)

La dépendance en **R⁴** est la grande affaire de ce chapitre : diviser le rayon
par deux divise le débit par seize. C'est elle qui rend une sténose artérielle si
coûteuse pour le cœur (MF3-7) et qui explique qu'une fissure de 10 µm laisse fuir
mille fois moins d'eau que ne le prédirait Torricelli (MF3-8).

Exercices MF3-7 et MF3-8 du TD « Fluides réels en écoulement » (PSI).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


# --------------------------------------------------------------------------- #
# 1. Conduite cylindrique (MF3-7 : l'artère)
# --------------------------------------------------------------------------- #
def poiseuille_profile(
    r: Array, *, radius: float, pressure_drop: float, length: float, eta: float
) -> Array:
    """
    Profil parabolique v(r) = ΔP (R² − r²) / (4ηL).

    Obtenu en intégrant l'équilibre  (1/r) d/dr( r dv/dr ) = −ΔP/(ηL)  avec la
    condition d'adhérence v(R) = 0 et la régularité sur l'axe (dv/dr = 0 en r = 0).
    """
    return pressure_drop * (radius**2 - r**2) / (4.0 * eta * length)


def poiseuille_flow_rate(
    *, radius: float, pressure_drop: float, length: float, eta: float
) -> float:
    """
    Débit volumique Dv = π R⁴ ΔP / (8ηL) (question 3 de MF3-7).

    L'intégrale ∫₀^R v(r) 2πr dr du profil parabolique. Le R⁴ vient de la
    combinaison « section en R² × vitesse moyenne en R² ».
    """
    return np.pi * radius**4 * pressure_drop / (8.0 * eta * length)


def pressure_drop_for_flow(
    *, flow_rate: float, radius: float, length: float, eta: float
) -> float:
    """Chute de pression nécessaire pour imposer un débit : ΔP = 8ηL Dv / (πR⁴)."""
    return 8.0 * eta * length * flow_rate / (np.pi * radius**4)


def hydraulic_resistance(*, radius: float, length: float, eta: float) -> float:
    """
    Résistance hydraulique Rh = 8ηL/(πR⁴), définie par ΔP = Rh · Dv.

    Analogie électrique directe : ΔP joue le rôle de la tension, Dv celui du
    courant. Les conduites en série ajoutent leurs résistances, comme des
    résistors — c'est la base du calcul de réseaux hydrauliques.
    """
    return 8.0 * eta * length / (np.pi * radius**4)


def mean_velocity(*, flow_rate: float, radius: float) -> float:
    """Vitesse débitante (moyenne de section) U = Dv/(πR²). Ici U = v_max/2."""
    return flow_rate / (np.pi * radius**2)


def stenosis_pressure_ratio(reduction: float | Array) -> float | Array:
    """
    Surcoût de pression dû à une sténose, **à débit constant** (question 5 de MF3-7).

    Si le rayon perd une fraction x, ΔP ∝ R⁻⁴ impose

        ΔP / ΔP₀ = (1 − x)⁻⁴

    Pour x = 5 %, cela fait +22,8 % de pression à fournir pour le même débit
    sanguin : une plaque d'athérome à peine visible coûte déjà très cher au cœur.
    """
    return (1.0 - reduction) ** -4


def mercury_pressure(height: float, *, rho_mercury: float = 13.8e3, g: float = 9.81) -> float:
    """
    Convertit une hauteur de colonne de mercure en pression (P = ρ g h).

    L'énoncé donne la tension « 12–8 » comme ΔP_cœur = 4 cm de mercure ; les
    tensiomètres se lisent encore en cmHg pour des raisons historiques.
    """
    return rho_mercury * g * height


# --------------------------------------------------------------------------- #
# 2. Fente plane (MF3-8 : la citerne fêlée)
# --------------------------------------------------------------------------- #
def slit_profile(
    x: Array, *, width: float, thickness: float, head: float, eta: float,
    rho: float = 1000.0, g: float = 9.81,
) -> Array:
    """
    Profil de vitesse dans une fissure plane d'épaisseur b et de longueur e.

        v(x) = ρgh/(2ηe) · (b²/4 − x²),      x ∈ [−b/2, b/2]

    Même équilibre que Poiseuille mais en géométrie plane : la pression motrice est
    ici la pression hydrostatique ρgh au fond de la citerne, et la longueur sur
    laquelle elle chute est l'épaisseur e de la paroi.
    """
    return rho * g * head / (2.0 * eta * thickness) * (width**2 / 4.0 - x**2)


def slit_flow_rate(
    *, width: float, span: float, thickness: float, head: float, eta: float,
    rho: float = 1000.0, g: float = 9.81,
) -> float:
    """
    Débit à travers la fissure : Dv = ρ g h b³ a / (12 η e).

    Intégrale du profil parabolique sur la largeur b, multipliée par la longueur a
    de la fissure. Noter le **b³** : la fuite est hypersensible à l'ouverture.
    """
    return rho * g * head * width**3 * span / (12.0 * eta * thickness)


def torricelli_flow_rate(*, width: float, span: float, head: float, g: float = 9.81) -> float:
    """
    Débit qu'aurait prédit un modèle de **fluide parfait** (Bernoulli/Torricelli).

        v = √(2gh),   Dv = a b √(2gh)

    Question 4 de MF3-8 : la comparaison est spectaculaire. Le fluide parfait ignore
    la dissipation visqueuse dans la fissure, qui est justement ce qui limite la
    fuite — il surestime le débit de plusieurs ordres de grandeur.
    """
    return span * width * np.sqrt(2.0 * g * head)
