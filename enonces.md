# Énoncés et résolution analytique

Ce projet part d'un TD de mécanique des fluides (PSI, *Fluides réels en écoulement*).
On reprend ici six exercices du sujet, avec l'énoncé puis la résolution « à la main ».
Le code du dépôt ne fait que mettre ces résultats en images — l'idée est que tu puisses
refaire les calculs toi-même avant de regarder la simulation.

Notation : `η` est la viscosité dynamique (Pa·s), `ν = η/ρ` la viscosité cinématique
(m²/s), `v` le champ de vitesse, `ΔP` une différence de pression motrice, `Dv` un débit
volumique.

---

## 0. Le point de départ : la viscosité diffuse la quantité de mouvement

Dans un fluide newtonien, deux couches voisines qui glissent l'une sur l'autre
s'échangent une force tangentielle proportionnelle au taux de cisaillement :

$$
f_s = \eta \, \frac{\partial v_x}{\partial z}
$$

C'est la définition même de η. Faisons alors un bilan de quantité de mouvement sur une
tranche de fluide comprise entre `z` et `z + dz`, de section `S`, dans un écoulement
unidirectionnel `v = v(z,t) e_x` (donc sans gradient de pression) :

$$
\rho\,S\,dz\,\frac{\partial v}{\partial t}
= \eta S \left[\frac{\partial v}{\partial z}\Big|_{z+dz} - \frac{\partial v}{\partial z}\Big|_{z}\right]
= \eta S\, \frac{\partial^2 v}{\partial z^2}\,dz
$$

d'où **l'équation de diffusion de la quantité de mouvement** :

$$
\boxed{\ \frac{\partial v}{\partial t} = \nu\,\frac{\partial^2 v}{\partial z^2},
\qquad \nu = \frac{\eta}{\rho}\ }
$$

C'est *mot pour mot* l'équation de la chaleur, avec ν à la place de la diffusivité
thermique. Tout le vocabulaire suit : temps de diffusion `τ = L²/ν`, épaisseur de
pénétration `δ ≈ √(νt)`, profil stationnaire affine. La viscosité n'est rien d'autre
qu'un **coefficient de diffusion de la quantité de mouvement** — c'est l'idée qu'il
faut retenir de tout le chapitre.

Quand il y a en plus un gradient de pression et que le régime est stationnaire, le
terme de gauche s'annule et il reste un équilibre pression ↔ viscosité : c'est
Poiseuille (partie 3).

---

## 1. Mise en mouvement d'une plaque — exercice MF3-1

> Un fluide incompressible (masse volumique µ, viscosité η, viscosité cinématique
> ν = 10⁻⁶ m²/s) est limité par deux plaques infinies, en `z = 0` et `z = L`. Tout est
> au repos pour `t < 0` ; à `t = 0` on met la plaque `z = 0` en mouvement à la vitesse
> constante `U`.
> 1. Montrer que l'accélération d'une particule de fluide s'écrit simplement `∂v/∂t`.
> 2. Montrer que `∂v/∂t = ν ∂²v/∂z²`. Ordre de grandeur de la durée d'établissement
>    du régime stationnaire pour `L = 10 cm` ?
> 3. En régime stationnaire, déterminer `v(z)`.
> 4. Déterminer la force surfacique exercée par le fluide sur la plaque inférieure.

**1. L'accélération.** L'accélération d'une particule est la dérivée *particulaire*

$$
\vec a = \frac{\partial \vec v}{\partial t} + (\vec v \cdot \vec\nabla)\vec v .
$$

Ici les lignes de champ sont des droites parallèles à `ex`, et `v` ne dépend que de `z`.
Le terme convectif s'écrit `v ∂v/∂x · ex` — or `v` ne dépend pas de `x`, donc il est
**identiquement nul**. Il reste `a = ∂v/∂t · ex`. C'est ce qui rend le problème linéaire
et exactement soluble : la non-linéarité de Navier–Stokes disparaît par géométrie.

**2. L'équation.** C'est le bilan de la partie 0 : `∂v/∂t = ν ∂²v/∂z²`. Par analyse
dimensionnelle, `[ν] = L²T⁻¹`, donc le seul temps qu'on puisse fabriquer avec ν et L est

$$
\tau \sim \frac{L^2}{\nu} = \frac{(0{,}1)^2}{10^{-6}} = 10^4\ \text{s} \approx 2\ \text{h}\ 45 .
$$

Ce résultat est frappant : il faut **presque trois heures** pour que le mouvement d'une
plaque se propage à travers 10 cm d'eau. La diffusion visqueuse est un mécanisme de
transport très lent — d'où la nécessité de brasser un liquide plutôt que d'attendre.

**3. Régime stationnaire.** `∂v/∂t = 0` donne `v'' = 0`, donc `v` est affine. Avec
`v(0) = U` et `v(L) = 0` :

$$
\boxed{\, v(z) = U\left(1 - \frac{z}{L}\right) \,}
$$

C'est **l'écoulement de Couette plan**. Le retrouver par un bilan de forces est
instructif : sur un parallélépipède de section `S` entre `z` et `z+dz`, à l'équilibre,
les contraintes visqueuses des deux faces se compensent, donc `η dv/dz` est uniforme,
donc `v` est affine. Même résultat, sans écrire d'EDP.

**4. Force sur la plaque.** La contrainte tangentielle vaut `η |dv/dz| = ηU/L`. Le
fluide **freine** la plaque :

$$
\boxed{\, \vec f_s = -\frac{\eta U}{L}\,\vec e_x \,}
$$

C'est la base du viscosimètre plan (exercice MF3-4 ci-dessous) : mesurer une force et
une vitesse donne η.

### Aux temps courts : le premier problème de Stokes

L'énoncé s'arrête au régime stationnaire, mais tout le transitoire est calculable.
Tant que `δ(t) = √(νt) ≪ L`, le fluide « ignore » la plaque du haut : le problème est
celui d'un milieu **semi-infini**. On cherche alors une solution auto-semblable
`v = U f(ζ)` avec `ζ = z/(2√(νt))` ; l'EDP se réduit à l'équation différentielle
ordinaire `f'' + 2ζ f' = 0`, qui s'intègre en

$$
\boxed{\, v(z,t) = U\,\mathrm{erfc}\!\left(\frac{z}{2\sqrt{\nu t}}\right) \,}
$$

C'est le **premier problème de Stokes**. Tous les profils se superposent une fois
tracés en fonction de ζ : c'est la signature d'une diffusion pure, et c'est exactement
ce que vérifie `figures/stokes_similarity.png`.

Simulation : `uv run main.py --mode stokes`

---

## 2. Viscosimètre de Couette — exercice MF3-2

> Un fluide visqueux est placé entre deux cylindres concentriques de hauteur `H` et de
> rayons `R₁` et `R₂`. Le cylindre intérieur est maintenu fixe par un fil de torsion,
> l'extérieur tourne à la vitesse angulaire constante `ω`. L'écoulement est
> orthoradial, `v(r) = A r + B/r`.
> 1. Déterminer `A` et `B`. À quoi correspondrait `B = 0` ?
> 2-3. Déterminer les couples subis par les deux cylindres. On admet
>    `f_s = η (dv/dr − v/r)`.
> 4. Comment déduire la viscosité du fluide de la mesure ?

**1. Les constantes.** Conditions d'adhérence aux parois : `v(R₁) = 0` (cylindre fixe)
et `v(R₂) = ω R₂` (entraîné). Le système

$$
A R_1 + \frac{B}{R_1} = 0, \qquad A R_2 + \frac{B}{R_2} = \omega R_2
$$

donne, en éliminant B = −A R₁² :

$$
\boxed{\, A = \frac{\omega R_2^2}{R_2^2 - R_1^2}, \qquad
B = -\frac{\omega R_1^2 R_2^2}{R_2^2 - R_1^2} \,}
$$

Les deux termes ont un sens physique distinct : `A r` est une **rotation en bloc** (le
fluide tourne comme un solide, sans se cisailler) et `B/r` est un **tourbillon libre**.
Si `B = 0`, il n'y a aucun cisaillement, donc aucune contrainte visqueuse, donc *aucun
couple mesurable* : ce serait le cas des deux cylindres tournant ensemble à la même
vitesse angulaire. Le viscosimètre ne fonctionne que grâce au terme en `B`.

**2-3. Le couple.** Avec `v = Ar + B/r`, on a `dv/dr = A − B/r²` et `v/r = A + B/r²`,
donc

$$
f_s = \eta\left(\frac{dv}{dr} - \frac{v}{r}\right) = -\frac{2\eta B}{r^2} .
$$

Le terme de rotation en bloc a bien disparu, comme annoncé. Le couple exercé sur un
cylindre de rayon `r`, de surface latérale `2πrH`, avec un bras de levier `r`, vaut

$$
\Gamma = f_s \cdot 2\pi r H \cdot r = 4\pi \eta H |B|
$$

**qui ne dépend pas de `r`.** Ce n'est pas un hasard : en régime stationnaire, le moment
cinétique du fluide ne varie pas, donc le couple qui entre par le cylindre extérieur
ressort intégralement par l'intérieur. C'est *exactement* ce qui rend l'appareil
utilisable — on mesure sur le cylindre fixe le couple imposé à l'autre.

$$
\boxed{\, \Gamma = 4\pi \eta H\,\frac{\omega R_1^2 R_2^2}{R_2^2 - R_1^2} \,}
$$

**4. La mesure.** Γ est **proportionnel à η**, avec un facteur purement géométrique
connu. Le fil de torsion de raideur `C` se tord d'un angle `α` tel que `Γ = Cα`, d'où

$$
\eta = \frac{C\,\alpha\,(R_2^2 - R_1^2)}{4\pi H \omega R_1^2 R_2^2} .
$$

Simulation : `uv run main.py --mode stokes` (deuxième figure)

---

## 3. Mesure d'une viscosité par glissement — exercice MF3-4

> Un solide parallélépipédique de masse `M = 10 kg`, d'aire de base `S = 0,1 m²`, glisse
> sur un plan incliné de `α = 25°` à la vitesse `v₀ = 80 cm/s`. Une couche d'huile
> (ρ = 900 kg/m³) d'épaisseur `e = 1 mm` le sépare du plan. Le profil de vitesse est
> celui d'un écoulement de Couette, `v(z) = v₀ z/e`. Quel est le coefficient de
> viscosité de l'huile ? Vérifier la validité de l'hypothèse d'écoulement laminaire et
> stationnaire.

Le solide glisse à **vitesse constante** : son accélération est nulle, donc la
composante du poids le long de la pente équilibre exactement le frottement visqueux du
film. Ce dernier vaut `η S v₀/e` (contrainte de Couette × surface) :

$$
M g \sin\alpha = \eta\,\frac{S v_0}{e}
\qquad\Longrightarrow\qquad
\boxed{\, \eta = \frac{M g \sin\alpha \cdot e}{S\, v_0} \,}
$$

A.N. : `η = 10 × 9,81 × sin(25°) × 10⁻³ / (0,1 × 0,8) ≈ 0,52 Pa·s`, soit 500 fois l'eau —
c'est bien l'ordre de grandeur d'une huile moteur.

**Les deux vérifications** (la partie de l'exercice qu'on a tendance à bâcler) :

- *Laminaire ?* `Re = ρ v₀ e/η = 900 × 0,8 × 10⁻³/0,52 ≈ 1,4`. On est à trois ordres de
  grandeur sous la transition (≈ 2300) : l'écoulement est non seulement laminaire, il
  est **rampant** (les effets visqueux dominent complètement l'inertie).
- *Stationnaire ?* Le temps d'établissement est `e²/ν` avec `ν = η/ρ ≈ 5,8·10⁻⁴ m²/s`,
  soit `≈ 1,7 ms`. Le profil de Couette s'installe instantanément à notre échelle : on a
  bien le droit de le supposer établi.

---

## 4. Poiseuille : l'artère — exercice MF3-7

> On étudie la circulation du sang (viscosité η, masse volumique ρ) dans une artère
> cylindrique d'axe `Oz`, de longueur `L` et de rayon `R`. On note `ΔP = p(0) − p(L) > 0`.
> On postule `v(M) = v(r,z) e_z` et `p(M) = p(r,z)`, l'écoulement étant stationnaire.
> 1. Pourquoi peut-on négliger la pesanteur ? Justifier `v(r,z) = v(r)` et `p(r,z) = p(z)`.
> 2. Montrer que `v(r) = ΔP(R² − r²)/4ηL`.
> 3. Établir l'expression du débit volumique.
> 4. A.N. : `L = 1 m`, `R = 0,5 cm`, `Dv = 80 cm³/s`, `η = 4·10⁻³ Pl`. La tension
>    « 12-8 » correspond à `ΔP_cœur = 4 cm` de mercure. Calculer `ΔP_artère` et conclure.
> 5. Un dépôt de cholestérol réduit le rayon de 5 %. Quelle est la variation relative de
>    pression correspondante, à débit sanguin constant ?

**1. Les simplifications.** L'artère est horizontale et fine : la variation de pression
hydrostatique sur son diamètre (`ρg·2R ≈ 100 Pa`) est petite devant ΔP. L'**incompressibilité**
`div v = ∂v/∂z = 0` impose ensuite `v = v(r)`. Enfin la projection de Navier–Stokes sur
`e_r` donne `∂p/∂r = 0`, donc `p = p(z)`.

**2. Le profil.** Projeté sur `e_z` en régime stationnaire, Navier–Stokes se réduit à un
équilibre pression ↔ viscosité :

$$
0 = -\frac{dp}{dz} + \frac{\eta}{r}\frac{d}{dr}\!\left(r\frac{dv}{dr}\right)
$$

Le membre de gauche ne dépend que de `z`, celui de droite que de `r` : les deux sont donc
**constants**, et `dp/dz = −ΔP/L`. Il reste à intégrer deux fois :

$$
\frac{1}{r}\frac{d}{dr}\!\left(r\frac{dv}{dr}\right) = -\frac{\Delta P}{\eta L}
\quad\Longrightarrow\quad
v(r) = -\frac{\Delta P}{4\eta L}r^2 + C_1 \ln r + C_2
$$

La régularité sur l'axe (`v` fini en `r = 0`) impose `C₁ = 0`, et l'adhérence à la paroi
`v(R) = 0` fixe `C₂` :

$$
\boxed{\, v(r) = \frac{\Delta P}{4\eta L}\left(R^2 - r^2\right) \,}
$$

**3. Le débit.** En intégrant le profil sur la section, `Dv = ∫₀^R v(r)\,2\pi r\,dr` :

$$
\boxed{\, D_v = \frac{\pi R^4 \Delta P}{8 \eta L} \,}
\qquad\text{(loi de Hagen–Poiseuille)}
$$

Le **R⁴** est le résultat central du chapitre. Il combine deux effets : la section
croît comme `R²`, et la vitesse moyenne elle-même croît comme `R²` (parce qu'une
conduite large cisaille moins son fluide).

**4. A.N.** `ΔP_artère = 8ηL Dv/(πR⁴) = 8 × 4·10⁻³ × 1 × 80·10⁻⁶ / (π × (5·10⁻³)⁴) ≈ 1{,}3·10^3` Pa.
Le cœur fournit `ΔP_cœur = ρ_Hg g h = 13,8·10³ × 9,81 × 0,04 ≈ 5{,}4·10^3` Pa. L'artère ne
consomme donc que **24 %** de la pression disponible : le cœur a de la marge.

*Remarque honnête :* avec ce débit, `Re = ρUD/η ≈ 2700`, soit pile la zone de transition.
Le profil de Poiseuille reste une bonne approximation, mais ce n'est plus un résultat
exact — le débit de 80 cm³/s de l'énoncé est celui du débit cardiaque total, généreux
pour une seule artère.

**5. La sténose.** À **débit constant** (l'organisme régule pour maintenir l'irrigation),
`ΔP ∝ R⁻⁴`. Si le rayon perd une fraction `x` :

$$
\boxed{\ \frac{\Delta P}{\Delta P_0} = (1-x)^{-4}\ }
\qquad\Longrightarrow\qquad
x = 5\ \%\ :\ (0{,}95)^{-4} = 1{,}228
$$

Soit **+22,8 % de pression** à fournir pour un rétrécissement de 5 % à peine visible à
l'imagerie. C'est toute la brutalité de la puissance quatrième : une sténose de 50 %
demanderait **seize fois** plus de pression. Le cœur compense en augmentant la tension
artérielle — d'où le lien direct entre athérosclérose et hypertension.

Simulation : `uv run main.py --mode poiseuille`

---

## 5. La citerne fêlée — exercice MF3-8

> Le fond d'une citerne remplie d'eau sur `h = 1 m` comporte une fissure de longueur
> `a = 5 cm` et de largeur `b = 10 µm`. Les parois ont une épaisseur `e = 2 cm`, la
> viscosité de l'eau vaut `η = 10⁻³ Pl`.
> 1. Proposer une forme concevable pour le champ des vitesses, et le déterminer.
> 2. Vérifier le caractère laminaire de l'écoulement.
> 3. Déterminer le débit et la quantité d'eau perdue chaque jour.
> 4. Comparer à ce qu'on obtiendrait en considérant l'eau comme un fluide parfait.

**1. Le profil.** Même physique que Poiseuille, mais en géométrie plane : la fissure est
beaucoup plus longue (`a`) que large (`b`), donc on cherche `v = v(x) e_z` où `x` est
la coordonnée en travers de la fente. La pression motrice est la pression hydrostatique
`ρgh` du fond de la citerne, qui chute sur l'épaisseur `e` de la paroi, d'où
`dp/dz = −ρgh/e`. L'équilibre devient `η d²v/dx² = −ρgh/e`, et avec l'adhérence
`v(±b/2) = 0` :

$$
\boxed{\, v(x) = \frac{\rho g h}{2\eta e}\left(\frac{b^2}{4} - x^2\right) \,}
$$

**3. Le débit.** En intégrant sur la largeur et en multipliant par la longueur `a` :

$$
\boxed{\, D_v = \frac{\rho g h\, b^3 a}{12\,\eta\, e} \,}
$$

Noter le **b³** : la fuite est hypersensible à l'ouverture de la fissure. A.N. :
`Dv ≈ 2,0·10⁻⁹ m³/s ≈ 2 mm³/s`, soit **0,18 L par jour**. Une fissure de 10 µm ne vide
pas une citerne.

**2. Laminaire ?** La vitesse débitante vaut `U = Dv/(ab) ≈ 4·10⁻³ m/s`, donc
`Re = ρUb/η ≈ 0,04`. On est en **écoulement rampant** : la viscosité écrase totalement
l'inertie.

**4. Et si l'eau était parfaite ?** Bernoulli entre la surface libre et la sortie donne
la formule de Torricelli `v = √(2gh) = 4,4 m/s`, soit `Dv = ab√(2gh) ≈ 2,2·10⁻⁶ m³/s`,
c'est-à-dire **191 L par jour** — plus de **mille fois** le débit réel.

C'est la leçon de l'exercice : dans un canal micrométrique, la dissipation visqueuse
n'est pas une petite correction, elle est *le* phénomène qui fixe le débit. Le modèle de
fluide parfait n'y a aucun sens, alors qu'il marche très bien pour un trou de 1 cm.

---

## 6. Le château d'eau — exercice MF3-6

> Un château d'eau de hauteur 60 m alimente un village. Une conduite de longueur 100 m
> et de rayon `R = 1,5 cm` part de son pied et débouche à l'air libre sur un robinet.
> 1. Quel débit peut-on attendre si l'écoulement est laminaire ?
> 2. Calculer alors le nombre de Reynolds. Conclure.
> 3. En réalité on mesure une vitesse débitante de 3,2 m/s. Déterminer le coefficient de
>    friction et estimer la rugosité du tuyau à l'aide du diagramme de Moody.

**1. La prédiction laminaire.** La pression motrice est celle de la colonne d'eau,
`ΔP = ρgH = 1000 × 9,81 × 60 ≈ 5,9·10⁵ Pa`. Poiseuille donne

$$
D_v = \frac{\pi R^4 \rho g H}{8\eta L} \approx 1{,}2 \cdot 10^{-1}\ \mathrm{m^3/s} = 117\ \mathrm{L/s}
$$

**2. Le verdict.** Cela correspond à une vitesse débitante `U = Dv/(πR²) ≈ 166 m/s` (!)
et donc

$$
Re = \frac{\rho U D}{\eta} \approx 5\cdot 10^6 \ggg 2300 .
$$

L'hypothèse laminaire **se contredit elle-même** : elle prédit une vitesse tellement
grande que l'écoulement ne peut pas être laminaire. C'est tout l'intérêt de l'exercice —
Poiseuille n'est pas « un peu imprécis » ici, il est faux d'un facteur 50.

**3. Le vrai calcul.** En turbulent, on abandonne le calcul exact pour la formulation de
**Darcy–Weisbach**, où toute l'ignorance est rangée dans un coefficient sans dimension λ :

$$
\Delta P = \lambda\,\frac{L}{D}\,\frac{1}{2}\rho U^2
$$

Le bilan de charge entre la surface libre et la sortie (en comptant l'énergie cinétique
emportée) s'écrit `ρgH = (λL/D + 1)·½ρU²`, d'où, avec `U = 3,2 m/s` mesuré :

$$
\lambda = \frac{D}{L}\left(\frac{2gH}{U^2} - 1\right) \approx 0{,}034,
\qquad Re = \frac{\rho U D}{\eta} \approx 9{,}6\cdot 10^4
$$

Il reste à remonter à la rugosité. L'équation de **Colebrook–White**

$$
\frac{1}{\sqrt\lambda} = -2\log_{10}\!\left(\frac{\varepsilon}{3{,}7 D} + \frac{2{,}51}{Re\sqrt\lambda}\right)
$$

est implicite en λ (c'est pour ça qu'on lit habituellement un abaque), mais elle est
**explicite en ε** ! Il suffit de l'inverser :

$$
\boxed{\, \varepsilon = 3{,}7\,D\left[10^{-1/(2\sqrt\lambda)} - \frac{2{,}51}{Re\sqrt\lambda}\right] \,}
$$

A.N. : `ε ≈ 0,20 mm`, à comparer aux **0,24 mm** de l'énoncé (obtenus par lecture
graphique sur l'abaque de Moody). L'écart de 15 % est celui qu'on attend entre un calcul
exact et une lecture à l'œil sur un diagramme log-log — le résultat physique est le même :
de l'acier commercial légèrement corrodé.

Ce qu'il faut retenir : **une simple mesure de débit permet de caractériser l'état de
surface d'une conduite enterrée** qu'on ne peut pas inspecter.

Simulation : `uv run main.py --mode moody`

---

## 7. L'oléoduc d'Alaska — exercice MF3-9

> L'oléoduc de l'Alaska a un diamètre de 1,2 m. Le pétrole brut (ρ = 800 kg/m³,
> η = 0,3 Pa·s) s'écoule avec un débit de 3400 L/s.
> 1. Estimer la chute de pression par mètre en supposant la loi de Poiseuille vérifiée.
>    L'hypothèse sous-jacente est-elle raisonnable ?
> 2. Estimer la chute de pression par mètre à l'aide du diagramme de Moody, sachant que
>    les aspérités des parois sont de l'ordre de 0,2 mm.
> 3. Le coefficient de perte de charge singulière d'un coude arrondi vaut
>    `Λ = [0,13 + 1,85 (D/2R_c)^3,5]·(θ/90)`. Quelle est la perte de charge due à un
>    coude de 10° de rayon de courbure 10 m ? Comment la minimiser ?

**1. Poiseuille, pour voir.** `U = Dv/(πD²/4) = 3,4/1,13 ≈ 3,0 m/s`, et

$$
\frac{\Delta P}{L} = \frac{8\eta D_v}{\pi R^4} \approx 20\ \mathrm{Pa/m} .
$$

Mais `Re = ρUD/η = 800 × 3,0 × 1,2/0,3 ≈ 9,6·10³ > 4000` : on est en régime turbulent,
l'hypothèse **n'est pas valable**. (Le pétrole est pourtant 300 fois plus visqueux que
l'eau — c'est le diamètre de 1,2 m qui fait basculer le Reynolds.)

**2. Moody.** Avec `ε/D = 2·10⁻⁴/1,2 ≈ 1,7·10⁻⁴` et `Re ≈ 9,6·10³`, Colebrook donne
`λ ≈ 0,031`, d'où

$$
\frac{\Delta P}{L} = \lambda\,\frac{1}{D}\,\frac{1}{2}\rho U^2 \approx 95\ \mathrm{Pa/m} ,
$$

soit **4,7 fois** la prédiction laminaire. Sur les 1300 km de l'oléoduc, cela représente
des dizaines de stations de pompage — l'erreur n'est pas académique, elle se compte en
mégawatts.

**3. Le coude.** `Λ = [0,13 + 1,85 (0,6/10)^3,5] × (10/90) ≈ 0,0145`, d'où

$$
\Delta P = \Lambda\,\frac{1}{2}\rho U^2 \approx 52\ \mathrm{Pa} .
$$

C'est une perte de charge **singulière** : contrairement à la perte régulière, elle n'est
pas proportionnelle à la longueur mais localisée sur l'accident de conduite, où
l'écoulement décolle et dissipe de l'énergie dans des tourbillons. Ces 52 Pa équivalent à
**55 cm de conduite droite** seulement — le coude est doux, donc peu coûteux. Pour la
minimiser : augmenter le rayon de courbure `R_c` (le terme en `(D/2R_c)^3,5` s'effondre
très vite) et réduire l'angle. C'est la justification quantitative de la règle du métier :
*un coude large coûte beaucoup moins cher qu'un coude sec.*

---

## Pour aller plus loin

Le TD contient d'autres situations qui reposent sur la même physique et qu'on pourrait
ajouter au projet : l'écoulement sur un plan incliné (MF3-3) et son profil semi-parabolique,
le viscosimètre à écoulement qui mesure η en chronométrant une vidange (MF3-5), la perfusion
et sa hauteur de flacon minimale (MF3-11), les lois de similitude en soufflerie (MF3-12), ou
le viscosimètre à bille et la traînée de Stokes `F = 6πηRv` (MF3-13). Toutes se ramènent au
même squelette : identifier la géométrie qui annule le terme convectif, écrire l'équilibre
pression ↔ viscosité, puis **vérifier le Reynolds** — c'est lui qui dit si le calcul qu'on
vient de faire a le droit d'exister.
