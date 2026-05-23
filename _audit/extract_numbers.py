#!/usr/bin/env python3
"""
extract_numbers.py
==================
Source unique de verite pour le rapport LaTeX `rapport_latex/rapport.tex`.

Re-execute les fonctions cles du notebook `Greeks_Monte_Carlo.ipynb`
(simulation antithetique seed=42, 50 000 paires, methodes FD / Pathwise /
Malliavin pour Delta, Gamma, Vega, schema d'Euler 4 valeurs n, et EDP
Crank-Nicolson) et ecrit `_audit/numbers.tex` (un \newcommand par chiffre
publie dans le rapport).

Convention statistique adoptee (alignee sur le notebook) :
- "Estimation" = mean de l'estimateur Monte Carlo (point estimate).
- "Erreur standard" SE = sigma_sample / sqrt(N)  (= ecart-type de la moyenne).
- "IC 95%" = ±1.96 * SE.
- Pour FD, le IC est calcule par bootstrap (1.96 * std des replicas, B=200) ;
  on en deduit SE_FD = IC_FD / 1.96 pour homogeneite avec PW/MAL.

Le rapport publie colonnes "Ecart-type" et "IC 95%". Le terme "Ecart-type"
designe en fait l'erreur standard SE (convention finance quantitative
courante : ce que rapporte un IC ±1.96 SE). On garde le libelle existant.

Usage :
    .venv-audit/bin/python _audit/extract_numbers.py

Idempotent : peut etre relance autant de fois que necessaire.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import solve_banded
from scipy.stats import norm

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parent.parent
AUDIT_DIR = REPO_ROOT / "_audit"

numbers: dict[str, str] = {}


def reg(key: str, value, fmt: str = "{}") -> None:
    """Enregistre un chiffre dans le dict numbers (formatte si besoin)."""
    numbers[key] = fmt.format(value)


# =============================================================================
# Parametres (identiques au notebook cell 2)
# =============================================================================

S0 = 100.0
K = 100.0
T = 1.0
r = 0.05
sigma = 0.20
N = 50_000  # nombre de paires (donc 100 000 evaluations apres antithetique)
SEED = 42

# =============================================================================
# Formules Black-Scholes (notebook cell 3) - utilisees comme benchmark
# =============================================================================


def bs_d1(S0, K, T, r, sigma):
    return (np.log(S0 / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))


def bs_price(S0, K, T, r, sigma):
    d1 = bs_d1(S0, K, T, r, sigma)
    d2 = d1 - sigma * np.sqrt(T)
    return S0 * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)


def bs_delta(S0, K, T, r, sigma):
    return norm.cdf(bs_d1(S0, K, T, r, sigma))


def bs_gamma(S0, K, T, r, sigma):
    return norm.pdf(bs_d1(S0, K, T, r, sigma)) / (S0 * sigma * np.sqrt(T))


def bs_vega(S0, K, T, r, sigma):
    return S0 * norm.pdf(bs_d1(S0, K, T, r, sigma)) * np.sqrt(T)


# Valeurs BS de reference
price_bs = bs_price(S0, K, T, r, sigma)
delta_bs = bs_delta(S0, K, T, r, sigma)
gamma_bs = bs_gamma(S0, K, T, r, sigma)
vega_bs = bs_vega(S0, K, T, r, sigma)

# Parametre N reutilise dans les captions
reg("paramN", "50\\,000")

# BS reference (4 decimales comme rapport)
reg("bsPrice", price_bs, "{:.4f}")
reg("bsDelta", delta_bs, "{:.4f}")
reg("bsGamma", gamma_bs, "{:.4f}")
reg("bsVega", vega_bs, "{:.2f}")
# Versions a plus de chiffres pour les tables Gamma et Vega
reg("bsGammaFive", gamma_bs, "{:.5f}")  # 0.01876 dans le rapport
reg("bsVegaThree", vega_bs, "{:.3f}")   # 37.524 dans le rapport

# =============================================================================
# Simulation principale (notebook cell 5) - seed=42, antithetique
# =============================================================================


def simulate(S0, T, r, sigma, n_paths, seed=SEED):
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal(n_paths)
    W_T = np.sqrt(T) * Z
    drift = (r - 0.5 * sigma**2) * T
    S_plus = S0 * np.exp(drift + sigma * W_T)
    S_minus = S0 * np.exp(drift - sigma * W_T)
    return S_plus, S_minus, W_T


S_plus, S_minus, W_T = simulate(S0, T, r, sigma, N)

# =============================================================================
# Delta - 3 methodes (notebook cell 7)
# =============================================================================


def compute_delta(S_plus, S_minus, W_T, S0, K, T, r, sigma, eps=1.0, n_boot=200, seed=SEED):
    df = np.exp(-r * T)
    n = len(W_T)
    rng = np.random.default_rng(seed)

    def mc_price(S0_val, W):
        S_p = S0_val * np.exp((r - 0.5 * sigma**2) * T + sigma * W)
        S_m = S0_val * np.exp((r - 0.5 * sigma**2) * T - sigma * W)
        payoff = 0.5 * (np.maximum(S_p - K, 0) + np.maximum(S_m - K, 0))
        return df * payoff.mean()

    delta_fd = (mc_price(S0 + eps, W_T) - mc_price(S0 - eps, W_T)) / (2 * eps)

    fd_boots = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        W_boot = W_T[idx]
        d = (mc_price(S0 + eps, W_boot) - mc_price(S0 - eps, W_boot)) / (2 * eps)
        fd_boots.append(d)
    ci_fd = 1.96 * np.std(fd_boots)

    pw_plus = df * (S_plus > K) * S_plus / S0
    pw_minus = df * (S_minus > K) * S_minus / S0
    pw_samples = 0.5 * (pw_plus + pw_minus)
    delta_pw = pw_samples.mean()
    ci_pw = 1.96 * pw_samples.std(ddof=1) / np.sqrt(n)

    payoff_plus = np.maximum(S_plus - K, 0)
    payoff_minus = np.maximum(S_minus - K, 0)
    pi_plus = W_T / (S0 * sigma * T)
    pi_minus = -W_T / (S0 * sigma * T)
    mal_samples = 0.5 * (df * payoff_plus * pi_plus + df * payoff_minus * pi_minus)
    delta_mal = mal_samples.mean()
    ci_mal = 1.96 * mal_samples.std(ddof=1) / np.sqrt(n)

    return {
        "fd": (delta_fd, ci_fd),
        "pathwise": (delta_pw, ci_pw),
        "malliavin": (delta_mal, ci_mal),
    }


res_d = compute_delta(S_plus, S_minus, W_T, S0, K, T, r, sigma)
d_fd, ic_d_fd = res_d["fd"]
d_pw, ic_d_pw = res_d["pathwise"]
d_mal, ic_d_mal = res_d["malliavin"]

# SE = IC / 1.96 (convention notebook)
se_d_fd = ic_d_fd / 1.96
se_d_pw = ic_d_pw / 1.96
se_d_mal = ic_d_mal / 1.96

# Estimation - 4 decimales
reg("deltaFdEst", d_fd, "{:.4f}")
reg("deltaPwEst", d_pw, "{:.4f}")
reg("deltaMalEst", d_mal, "{:.4f}")

# Erreur standard (4 decimales pour Delta - meme nombre de chiffres significatifs que le rapport original)
reg("deltaFdSe", se_d_fd, "{:.4f}")
reg("deltaPwSe", se_d_pw, "{:.4f}")
reg("deltaMalSe", se_d_mal, "{:.4f}")

# Intervalle de confiance 95% (4 decimales avec ±)
reg("deltaFdIc", ic_d_fd, "{:.4f}")
reg("deltaPwIc", ic_d_pw, "{:.4f}")
reg("deltaMalIc", ic_d_mal, "{:.4f}")

# Erreur relative (en %)
err_d_fd = 100.0 * abs(d_fd - delta_bs) / delta_bs
err_d_pw = 100.0 * abs(d_pw - delta_bs) / delta_bs
err_d_mal = 100.0 * abs(d_mal - delta_bs) / delta_bs
reg("deltaFdErrPct", err_d_fd, "{:.2f}")
reg("deltaPwErrPct", err_d_pw, "{:.2f}")
reg("deltaMalErrPct", err_d_mal, "{:.2f}")

# Ratio variance MAL / PW (en termes de IC, donc d'ecart-type)
ratio_d_mal_pw = ic_d_mal / ic_d_pw
reg("deltaRatioMalPw", ratio_d_mal_pw, "{:.1f}")
# Version variance (ratio des SE au carre)
ratio_var_d_mal_pw = (se_d_mal / se_d_pw) ** 2
reg("deltaVarRatioMalPw", ratio_var_d_mal_pw, "{:.0f}")

# =============================================================================
# Gamma - FD et Malliavin (notebook cell 10)
# =============================================================================


def compute_gamma(S_plus, S_minus, W_T, S0, K, T, r, sigma, eps=1.0, n_boot=200, seed=SEED):
    df = np.exp(-r * T)
    n = len(W_T)
    rng = np.random.default_rng(seed)

    def mc_price(S0_val, W):
        S_p = S0_val * np.exp((r - 0.5 * sigma**2) * T + sigma * W)
        S_m = S0_val * np.exp((r - 0.5 * sigma**2) * T - sigma * W)
        payoff = 0.5 * (np.maximum(S_p - K, 0) + np.maximum(S_m - K, 0))
        return df * payoff.mean()

    C_up = mc_price(S0 + eps, W_T)
    C_mid = mc_price(S0, W_T)
    C_down = mc_price(S0 - eps, W_T)
    gamma_fd = (C_up - 2 * C_mid + C_down) / (eps**2)

    fd_boots = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        W_boot = W_T[idx]
        C_up_b = mc_price(S0 + eps, W_boot)
        C_mid_b = mc_price(S0, W_boot)
        C_down_b = mc_price(S0 - eps, W_boot)
        fd_boots.append((C_up_b - 2 * C_mid_b + C_down_b) / (eps**2))
    ci_fd = 1.96 * np.std(fd_boots)

    payoff_plus = np.maximum(S_plus - K, 0)
    payoff_minus = np.maximum(S_minus - K, 0)
    pi_plus = (W_T**2 - sigma * T * W_T - T) / (S0**2 * sigma**2 * T**2)
    pi_minus = (W_T**2 + sigma * T * W_T - T) / (S0**2 * sigma**2 * T**2)
    mal_samples = 0.5 * (df * payoff_plus * pi_plus + df * payoff_minus * pi_minus)
    gamma_mal = mal_samples.mean()
    ci_mal = 1.96 * mal_samples.std(ddof=1) / np.sqrt(n)

    return {"fd": (gamma_fd, ci_fd), "malliavin": (gamma_mal, ci_mal)}


res_g = compute_gamma(S_plus, S_minus, W_T, S0, K, T, r, sigma)
g_fd, ic_g_fd = res_g["fd"]
g_mal, ic_g_mal = res_g["malliavin"]

se_g_fd = ic_g_fd / 1.96
se_g_mal = ic_g_mal / 1.96

# Estimation - 5 decimales (rapport publie 0.01882, 0.01871, 0.01876)
reg("gammaFdEst", g_fd, "{:.5f}")
reg("gammaMalEst", g_mal, "{:.5f}")

# SE - 5 decimales pour coherence (rapport publie 0.00052, 0.00178)
reg("gammaFdSe", se_g_fd, "{:.5f}")
reg("gammaMalSe", se_g_mal, "{:.5f}")

reg("gammaFdIc", ic_g_fd, "{:.5f}")
reg("gammaMalIc", ic_g_mal, "{:.5f}")

err_g_fd = 100.0 * abs(g_fd - gamma_bs) / gamma_bs
err_g_mal = 100.0 * abs(g_mal - gamma_bs) / gamma_bs
reg("gammaFdErrPct", err_g_fd, "{:.2f}")
reg("gammaMalErrPct", err_g_mal, "{:.2f}")

# =============================================================================
# Vega - 3 methodes (notebook cell 13)
# =============================================================================


def compute_vega(S_plus, S_minus, W_T, S0, K, T, r, sigma, eps=0.01, n_boot=200, seed=SEED):
    df = np.exp(-r * T)
    n = len(W_T)
    rng = np.random.default_rng(seed)

    def mc_price(sigma_val, W):
        S_p = S0 * np.exp((r - 0.5 * sigma_val**2) * T + sigma_val * W)
        S_m = S0 * np.exp((r - 0.5 * sigma_val**2) * T - sigma_val * W)
        payoff = 0.5 * (np.maximum(S_p - K, 0) + np.maximum(S_m - K, 0))
        return df * payoff.mean()

    vega_fd = (mc_price(sigma + eps, W_T) - mc_price(sigma - eps, W_T)) / (2 * eps)

    fd_boots = []
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        W_boot = W_T[idx]
        v = (mc_price(sigma + eps, W_boot) - mc_price(sigma - eps, W_boot)) / (2 * eps)
        fd_boots.append(v)
    ci_fd = 1.96 * np.std(fd_boots)

    pw_plus = df * (S_plus > K) * S_plus * (W_T - sigma * T)
    pw_minus = df * (S_minus > K) * S_minus * (-W_T - sigma * T)
    pw_samples = 0.5 * (pw_plus + pw_minus)
    vega_pw = pw_samples.mean()
    ci_pw = 1.96 * pw_samples.std(ddof=1) / np.sqrt(n)

    payoff_plus = np.maximum(S_plus - K, 0)
    payoff_minus = np.maximum(S_minus - K, 0)
    pi_plus = (W_T**2 - sigma * T * W_T - T) / (sigma * T)
    pi_minus = (W_T**2 + sigma * T * W_T - T) / (sigma * T)
    mal_samples = 0.5 * (df * payoff_plus * pi_plus + df * payoff_minus * pi_minus)
    vega_mal = mal_samples.mean()
    ci_mal = 1.96 * mal_samples.std(ddof=1) / np.sqrt(n)

    return {
        "fd": (vega_fd, ci_fd),
        "pathwise": (vega_pw, ci_pw),
        "malliavin": (vega_mal, ci_mal),
    }


res_v = compute_vega(S_plus, S_minus, W_T, S0, K, T, r, sigma)
v_fd, ic_v_fd = res_v["fd"]
v_pw, ic_v_pw = res_v["pathwise"]
v_mal, ic_v_mal = res_v["malliavin"]

se_v_fd = ic_v_fd / 1.96
se_v_pw = ic_v_pw / 1.96
se_v_mal = ic_v_mal / 1.96

# Estimation - 2 decimales (rapport publie 37.52, 37.53, 37.51)
reg("vegaFdEst", v_fd, "{:.2f}")
reg("vegaPwEst", v_pw, "{:.2f}")
reg("vegaMalEst", v_mal, "{:.2f}")

# SE - 2 decimales (rapport publie 0.12, 0.11, 0.85)
reg("vegaFdSe", se_v_fd, "{:.2f}")
reg("vegaPwSe", se_v_pw, "{:.2f}")
reg("vegaMalSe", se_v_mal, "{:.2f}")

reg("vegaFdIc", ic_v_fd, "{:.2f}")
reg("vegaPwIc", ic_v_pw, "{:.2f}")
reg("vegaMalIc", ic_v_mal, "{:.2f}")

err_v_fd = 100.0 * abs(v_fd - vega_bs) / vega_bs
err_v_pw = 100.0 * abs(v_pw - vega_bs) / vega_bs
err_v_mal = 100.0 * abs(v_mal - vega_bs) / vega_bs
reg("vegaFdErrPct", err_v_fd, "{:.2f}")
reg("vegaPwErrPct", err_v_pw, "{:.2f}")
reg("vegaMalErrPct", err_v_mal, "{:.2f}")

# Ratio variance MAL / PW pour Vega
ratio_v_mal_pw = ic_v_mal / ic_v_pw
reg("vegaRatioMalPw", ratio_v_mal_pw, "{:.1f}")
ratio_var_v_mal_pw = (se_v_mal / se_v_pw) ** 2
reg("vegaVarRatioMalPw", ratio_var_v_mal_pw, "{:.0f}")

# =============================================================================
# Synthese : Variance relative MAL / Pathwise
# =============================================================================
# On utilise le ratio des variances (SE^2) sur Delta (referent) car c'est
# l'estimateur le plus discriminant et Pathwise est applicable.
# Note : la ligne "Pathwise ~0.85x" du rapport vient du notebook cell 26
# qui compare la variance FD vs Pathwise (FD ~0.85x celle de Pathwise).
ratio_var_pw_pw = 1.0
# FD vs PW (rapport variance) sur Delta
ratio_var_fd_pw = (se_d_fd / se_d_pw) ** 2
reg("varRatioFdPw", ratio_var_fd_pw, "{:.2f}")
reg("varRatioPwPw", ratio_var_pw_pw, "{:.0f}")
reg("varRatioMalPw", ratio_var_d_mal_pw, "{:.0f}")

# =============================================================================
# Schema d'Euler vs Exact (notebook cell 17) - n=10, 50, 100, 500
# =============================================================================


def simulate_euler(S0, T, r, sigma, n_paths, n_steps=100, seed=SEED):
    rng = np.random.default_rng(seed)
    dt = T / n_steps
    sqrt_dt = np.sqrt(dt)
    S_plus_e = np.full(n_paths, S0, dtype=np.float64)
    S_minus_e = np.full(n_paths, S0, dtype=np.float64)
    W_T_e = np.zeros(n_paths)
    for _ in range(n_steps):
        Z = rng.standard_normal(n_paths)
        dW = sqrt_dt * Z
        W_T_e += dW
        S_plus_e = S_plus_e * (1 + r * dt + sigma * dW)
        S_minus_e = S_minus_e * (1 + r * dt - sigma * dW)
    S_plus_e = np.maximum(S_plus_e, 1e-10)
    S_minus_e = np.maximum(S_minus_e, 1e-10)
    return S_plus_e, S_minus_e, W_T_e


# Reference : Exact MC (Pathwise Delta)
delta_exact_pw = res_d["pathwise"][0]

n_steps_list = [10, 50, 100, 500]
euler_delta = {}
for n_steps in n_steps_list:
    Sp_e, Sm_e, W_e = simulate_euler(S0, T, r, sigma, N, n_steps)
    res_e = compute_delta(Sp_e, Sm_e, W_e, S0, K, T, r, sigma)
    delta_euler_n = res_e["pathwise"][0]
    euler_delta[n_steps] = abs(delta_euler_n - delta_exact_pw)

# Erreur Euler sur Delta vs Exact (notebook ref) - 4 decimales
reg("eulerErrTen", euler_delta[10], "{:.4f}")
reg("eulerErrFifty", euler_delta[50], "{:.4f}")
reg("eulerErrHundred", euler_delta[100], "{:.4f}")
reg("eulerErrFiveHundred", euler_delta[500], "{:.4f}")

# =============================================================================
# EDP Crank-Nicolson (notebook cell 31)
# =============================================================================


def solve_bs_pde(S0, K, T, r, sigma, S_max=300, N_S=200, N_t=200):
    dS = S_max / N_S
    dt = T / N_t
    S = np.linspace(0, S_max, N_S + 1)
    C = np.maximum(S - K, 0)
    i = np.arange(1, N_S)
    Si = S[i]
    a = 0.5 * dt * (sigma**2 * Si**2 / dS**2 - r * Si / dS)
    b = 1 + dt * (sigma**2 * Si**2 / dS**2 + r)
    c = 0.5 * dt * (sigma**2 * Si**2 / dS**2 + r * Si / dS)

    for j in range(N_t):
        C[0] = 0
        C[N_S] = S_max - K * np.exp(-r * (N_t - j) * dt)
        d = C[1:N_S].copy()
        d[0] -= a[0] * C[0]
        d[-1] -= c[-1] * C[N_S]
        ab = np.zeros((3, N_S - 1))
        ab[0, 1:] = -c[:-1]
        ab[1, :] = b
        ab[2, :-1] = -a[1:]
        C[1:N_S] = solve_banded((1, 1), ab, d)

    idx = int(S0 / dS)
    if idx >= N_S - 1:
        idx = N_S - 2
    C_m = C[idx - 1] if idx > 0 else C[0]
    C_0 = C[idx]
    C_p = C[idx + 1]
    alpha = (S0 - S[idx]) / dS
    price = (1 - alpha) * C_0 + alpha * C_p
    delta = (C_p - C_m) / (2 * dS)
    gamma = (C_p - 2 * C_0 + C_m) / (dS**2)
    return price, delta, gamma


def pde_vega(S0, K, T, r, sigma, eps=0.01, **kwargs):
    C_up, _, _ = solve_bs_pde(S0, K, T, r, sigma + eps, **kwargs)
    C_down, _, _ = solve_bs_pde(S0, K, T, r, sigma - eps, **kwargs)
    return (C_up - C_down) / (2 * eps)


price_pde, delta_pde, gamma_pde = solve_bs_pde(S0, K, T, r, sigma)
vega_pde = pde_vega(S0, K, T, r, sigma)

reg("pdePrice", price_pde, "{:.4f}")
reg("pdeDelta", delta_pde, "{:.4f}")
reg("pdeGamma", gamma_pde, "{:.4f}")
reg("pdeVega", vega_pde, "{:.2f}")

reg("pdeErrPrice", abs(price_pde - price_bs), "{:.4f}")
reg("pdeErrDelta", abs(delta_pde - delta_bs), "{:.4f}")
reg("pdeErrGamma", abs(gamma_pde - gamma_bs), "{:.4f}")
reg("pdeErrVega", abs(vega_pde - vega_bs), "{:.2f}")

# =============================================================================
# Ecriture du fichier numbers.tex
# =============================================================================

out_lines = [
    "% Auto-generated by _audit/extract_numbers.py - DO NOT EDIT MANUALLY",
    "% Source unique de verite pour rapport_latex/rapport.tex.",
    "% Regenerer via: .venv-audit/bin/python _audit/extract_numbers.py",
    "%",
    "% Convention statistique : 'Ecart-type' designe l'erreur standard",
    "% SE = sigma_sample / sqrt(N) (= IC95 / 1.96). 'IC 95%' = +/-1.96 * SE.",
    "% Pour FD, IC obtenu par bootstrap 200 replicas (notebook cell 7).",
    "",
]
for k, v in numbers.items():
    out_lines.append(f"\\newcommand{{\\{k}}}{{{v}}}")

(AUDIT_DIR / "numbers.tex").write_text("\n".join(out_lines) + "\n")

print(f"OK - wrote {len(numbers)} numeric facts to _audit/numbers.tex")
print()
print("Key numbers extracted:")
for key in [
    "bsPrice", "bsDelta", "bsGamma", "bsVega",
    "deltaFdEst", "deltaFdSe", "deltaFdIc",
    "deltaPwEst", "deltaPwSe", "deltaPwIc",
    "deltaMalEst", "deltaMalSe", "deltaMalIc",
    "deltaRatioMalPw", "deltaVarRatioMalPw",
    "gammaFdEst", "gammaFdSe", "gammaMalEst", "gammaMalSe",
    "vegaFdEst", "vegaFdSe", "vegaPwEst", "vegaPwSe", "vegaMalEst", "vegaMalSe",
    "eulerErrTen", "eulerErrFifty", "eulerErrHundred", "eulerErrFiveHundred",
    "pdePrice", "pdeDelta", "pdeGamma", "pdeVega",
    "pdeErrPrice", "pdeErrDelta", "pdeErrGamma", "pdeErrVega",
]:
    print(f"  {key} = {numbers[key]}")
