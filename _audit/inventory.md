# Inventory — Options-Greeks-Monte-Carlo

## Branche / repo
- Branch: `audit-2026-05` (tagguée `pre-audit-2026-05/main` = 7f751b9)
- Initial commit unique: `7f751b9` (Initial commit), puis `ca8863d` (audit phase 0).
- Auteur git: `Dan Allouche <dan.allouche@icloud.com>` (un seul auteur, OK).

## Fichiers trackés (git ls-files)
1. `.gitignore`
2. `Greeks_Monte_Carlo.ipynb` — notebook unique, 35 cells (13 code, 22 markdown), kernel python3 venv, python 3.12.4
3. `README.md` — 192 lignes, FR/EN mélangé, section Author = "Dan Allouche"
4. `_audit/.gitkeep`
5. `rapport_latex/rapport.tex` — 783 lignes, source LaTeX du rapport (correspond à Final_project.pdf)
6. `rapport_latex/figures/delta_comparaison.png`
7. `rapport_latex/figures/euler_convergence.png`
8. `rapport_latex/figures/gamma_comparaison.png`
9. `rapport_latex/figures/vega_comparaison.png`
10. `requirements.txt` — numpy/scipy/matplotlib/jupyter/ipykernel/notebook (versions souples >=)

## Fichiers physiques NON trackés (gitignored — voir .gitignore)
- `AlloucheDan.zip` (1.05 MB) — soumission ENSAE/Dauphine. Contenu: `Greeks_Monte_Carlo.ipynb` + `MonteCarlo.pdf` + macresource forks. Représente la livraison originale.
- `Final_project.pdf` (124 KB) — compilation du `rapport_latex/rapport.tex` (cohérent : titre, table des valeurs, structure 9 sections, paramètres S₀=100, K=100, etc.). C'est le rapport "officiel" court.
- `MonteCarlo.pdf` (725 KB) — autre version (probablement consigne / énoncé du projet, **inclus dans la soumission AlloucheDan.zip**). À examiner pour confirmer.
- `rapport_oral_complet.pdf` (347 KB) — version longue type "support de révision pour oral" (895 lignes tex), parallel au rapport court. Distinct.
- `rapport_oral_complet.tex` (30 KB) — source de ci-dessus, non tracké, untracked dans git status.
- `intro-calcul-sto.pdf` (957 KB), `Lecture notes.pdf` (777 KB), `polybouchard.pdf` (759 KB) — **trois PDFs de référence externe** (vraisemblables polycopiés/notes de cours, peut-être issus du cours Master MASEF). Pas signés "Dan Allouche". Risque licence-incertaine si copies de polys d'enseignants : conserver localement mais NE PAS publier (déjà gitignored — OK).
- `venv/` — environnement virtuel local, **NON tracké** (`git ls-files venv/` retourne 0 fichier). Source des 16 hits trufflehog raw (tous dans `venv/lib/python3.12/site-packages/...` tornado/pandas/pip, verified=0, exemples génériques). Hors zone de risque réelle.
- `.DS_Store` — non tracké (gitignored).

## Rôle des 5 PDFs externes
- `intro-calcul-sto.pdf`, `Lecture notes.pdf`, `polybouchard.pdf` : matériel pédagogique extérieur (probable polycopié L3/M1 calcul stochastique + cours Bouchard). Pas du code/contenu de l'auteur. Pas publiés (gitignored).
- `MonteCarlo.pdf` : présent ET dans le zip de soumission — possible énoncé/consigne du projet (à confirmer si signé/non).
- `Final_project.pdf` : rapport officiel court.
- `rapport_oral_complet.pdf` : doc de prep oral (895 lignes vs 783 du rapport.tex). Doublon partiel mais autre angle (avec code embarqué).
- Doublons potentiels : aucun strict (rapport.tex ≠ rapport_oral_complet.tex après diff).

## Cohérence Final_project.pdf vs rapport_latex/rapport.tex
- `rapport_latex/rapport.tex` est le source qui compile vers `Final_project.pdf`: titre identique, auteur identique, abstract identique. Structure 9 sections, paramètres S₀=100 K=100 T=1 r=5% σ=20% N=50,000 identiques. Tables Black-Scholes (Prix=10.4506, Δ=0.6368, Γ=0.0188, V=37.52) cohérentes. Donc rapport.tex = source de Final_project.pdf. **Confirmé**.

## Notebook (Greeks_Monte_Carlo.ipynb)
- 35 cells, structure identique à README: Intro → Toolbox → FD → Pathwise → Malliavin → Comparatif → Euler → EDP.
- Toutes les fonctions Monte Carlo utilisent `seed=42` par défaut + `np.random.default_rng(seed)`. **Déterministe**.
- BS analytiques calculés (cell 3): Prix=10.450584, Δ=0.636831, Γ=0.018762, V=37.524035. Cohérent rapport.
- Sorties intégrées dans le `.ipynb` (cells exécutées).
- Pas d'authors dans `metadata.authors` (absent). Kernelspec=python3.
- Méthodes implémentées: 3 méthodes (FD, Pathwise, Malliavin) pour Delta/Gamma/Vega + EDP Crank-Nicolson + Euler convergence + analyse epsilon biais-variance + histogramme distribution estimateurs.

## Sécurité / portabilité
- Secrets: 16 hits trufflehog raw, tous dans `venv/`, verified=0. **Pas de risque** (venv non tracké et déjà gitignored).
- Chemins absolus: aucun (abs_paths.txt vide).
- Hook pre-push installé (exit 1). OK.

## Identité
- Cas B (mon nom déjà présent normalisé) :
  - `rapport_latex/rapport.tex:29` → `Dan Allouche` ✓
  - `rapport_oral_complet.tex:81` → `Dan Allouche` ✓
  - `README.md:176` → `Dan Allouche` ✓
- Aucun autre nom humain. Aucune mention IA détectée.
