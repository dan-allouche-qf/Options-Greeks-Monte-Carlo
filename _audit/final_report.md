# Final report — Options-Greeks-Monte-Carlo (audit 2026-05)

## Statut

**Tier A confirmé** (post-régénération rapport 2026-05-24).

DECISIONS_REQUIRED [13] : **RÉSOLU**.

## Synthèse

Repo Greeks Monte Carlo (Delta / Gamma / Vega) avec 3 méthodes (Différences Finies, Pathwise, Malliavin) + schéma d'Euler + validation EDP Crank-Nicolson sur Black-Scholes ATM (S₀=100, K=100, T=1, r=5%, σ=20%, N=50 000, seed=42).

### Pipeline source unique de vérité

- `_audit/extract_numbers.py` — re-exécute toutes les fonctions du notebook `Greeks_Monte_Carlo.ipynb` (compute_delta, compute_gamma, compute_vega, simulate_euler, solve_bs_pde) avec mêmes seeds, écrit 58 `\newcommand{...}` dans `_audit/numbers.tex`.
- `rapport_latex/rapport.tex` — chaque chiffre numérique des tables Delta/Gamma/Vega/Euler/EDP est désormais inséré via `\input{../_audit/numbers.tex}` + macro nommée.
- Idempotent. Re-exécutable via `.venv-audit/bin/python _audit/extract_numbers.py`.

### 8 P0 corrigés (cf. DECISIONS [13])

| # | Quantité | Rapport (FAUX) | Notebook (recalculé) | Newcommand |
|---|---|---|---|---|
| 1 | Delta SE FD | 0.0021 | **0.0006** | `\deltaFdSe` |
| 2 | Delta SE PW | 0.0021 | **0.0006** | `\deltaPwSe` |
| 3 | Delta SE MAL | 0.0152 | **0.0042** | `\deltaMalSe` |
| 4 | Gamma FD | 0.01882 | **0.01828** | `\gammaFdEst` |
| 5 | Gamma MAL | 0.01871 | **0.01923** | `\gammaMalEst` |
| 6 | Vega FD | 37.52 | **37.61** | `\vegaFdEst` |
| 7 | Vega PW | 37.53 | **37.61** | `\vegaPwEst` |
| 8 | Vega MAL | 37.51 | **38.45** | `\vegaMalEst` |
| bonus | Synthèse Variance MAL | ~46× | **43×** (ratio variance) ou 6.5× (ratio écart-type) | `\varRatioMalPw` / `\deltaRatioMalPw` |
| bonus | Euler table | n=200 | **n=500** (vrai n du notebook) | `\eulerErrFiveHundred` |

### Convention statistique adoptée

- "Estimation" = `mean(estimateur)`
- "SE" (libellé "Écart-type" dans rapport, renommé "SE" pour clarté) = `sigma_sample / sqrt(N)` = erreur standard de la moyenne
- "IC 95%" = `±1.96 × SE`
- Pour FD, IC obtenu par bootstrap (200 ré-échantillons sur `W_T`), SE déduit par `IC/1.96`

Cette convention est celle du notebook (`ci = 1.96 * std(samples, ddof=1) / sqrt(n)`). La colonne "Écart-type" du rapport original confondait SE et σ_sample, ce qui expliquait le facteur 3.5× de divergence sur Delta.

### EDP (C25-C32) : NON TOUCHÉ

8 lignes EDP du rapport étaient déjà reproductibles parfaitement (cf. ancienne claims.csv). Conservées telles quelles, simplement passées par newcommand pour cohérence pipeline.

### Compilation

- `pdflatex rapport.tex` × 3 passes : 0 erreur, 0 référence cassée, warnings pré-existants (hyperref math chars dans titres section).
- PDF : `rapport_latex/rapport.pdf` (17 pages, 591 KB).
- `Final_project.pdf` racine régénéré identique (591 KB).
- PDFs non commités (convention repo : `*.pdf` gitignored).

### Cohérence claims.csv

- 32 claims : **32/32 MATCH** (post-régénération).
- Auparavant 8 NO-MATCH (C05–C08, C10, C12–C18, C20–C24 — soit 13 P0 dont 8 catégorisés P0 critiques + 5 P0 mineurs).

### Identité

Cas B : « Dan Allouche » partout (notebook, rapport, rapport_oral_complet, README). Inchangé.

### Non modifié

- `Greeks_Monte_Carlo.ipynb` — notebook correct, ré-exécuté en lecture seule.
- Notations math, équations, dérivations théoriques du rapport.
- Tables Black-Scholes (paramètres et valeurs analytiques).
- Bibliographie.

### Modifications mineures de présentation

- Colonne "Écart-type" → "SE" dans les 3 tables Delta/Gamma/Vega (clarification convention).
- Captions des tables enrichies : précision SE/IC convention + bootstrap pour FD.
- Phrase synthèse Delta reformulée : « ~7× plus grande » → « SE ~6.5× plus grande, donc variance ~43× plus élevée » (cohérent avec mesure réelle).
- Phrase synthèse Euler reformulée : à n≥50, erreur Euler du même ordre que SE MC, plus de gain net à augmenter n (saturation par bruit MC).

### Reproductibilité dynamique

Phase 2 dynamique exécutée : **OUI**.
- `.venv-audit/` créé from scratch via `python3 -m venv` (Python 3.12.x).
- `pip install -r requirements.txt` (numpy, scipy, matplotlib, jupyter, ipykernel, notebook).
- `_audit/extract_numbers.py` ré-exécuté avec succès, sorties bit-identiques au notebook committed.
- 3 passes `pdflatex` (TeX Live 2025 / TinyTeX).

## Commit

`audit(Options-Greeks-Monte-Carlo): regenerate rapport tables from notebook — fix 8 P0`

Fichiers modifiés :
- `rapport_latex/rapport.tex` (789 lignes, +6 vs initial)
- `_audit/extract_numbers.py` (NEW, 350 lignes)
- `_audit/numbers.tex` (NEW, 58 newcommands)
- `_audit/claims.csv` (refait, 32/32 MATCH)
- `_audit/final_report.md` (CE FICHIER, NEW)
- `.gitignore` (+1 ligne : `.venv-audit/`)

PDF non commité (convention repo `*.pdf` gitignored).
