# Options Greeks Monte Carlo

Monte Carlo computation of European option Greeks (Delta, Gamma, Vega) under Black-Scholes, comparing three approaches against analytical and PDE benchmarks.

## Overview

European call option pricing and Greek sensitivities under Black-Scholes by Monte Carlo. Three methods are implemented and compared:

- **Finite Differences**: numerical differentiation with common random numbers
- **Pathwise Derivative**: differentiation of the payoff (Lipschitz case)
- **Likelihood Ratio / Malliavin**: differentiation of the density via score functions

Antithetic variables are used for variance reduction. Results are validated against the closed-form Black-Scholes Greeks and against a Crank-Nicolson PDE solver.

## Methods Implemented

### 1. Finite Differences
- Approximates derivatives using perturbed simulations
- Works for all Greeks including Gamma
- Higher computational cost but always applicable

### 2. Pathwise Derivative
- Differentiates directly through the payoff function
- Lowest variance for Delta and Vega
- Cannot compute Gamma (discontinuous second derivative)

### 3. Likelihood Ratio Method
- Uses score function approach
- Applicable to all Greeks
- Moderate variance

## Getting Started

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Installation

1. Clone the repository:
```bash
git clone https://github.com/dan-allouche-qf/Options-Greeks-Monte-Carlo.git
cd Options-Greeks-Monte-Carlo
```

2. Create and activate a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install required packages:
```bash
pip install -r requirements.txt
```

### Usage

Open and run the Jupyter notebook:

```bash
jupyter notebook Greeks_Monte_Carlo.ipynb
```

The notebook is organized into sections:
1. **Introduction**: Theoretical background and problem setup
2. **Toolbox**: Implementation of simulation and Black-Scholes formulas
3. **Finite Differences**: Implementation and analysis
4. **Pathwise Derivative**: Implementation and analysis
5. **Likelihood Ratio**: Implementation and analysis
6. **Comparative Analysis**: Performance comparison across methods
7. **Conclusion**: Summary and recommendations

## Results

The project demonstrates:

- **Convergence**: All methods converge to Black-Scholes values with rate O(1/√n)
- **Variance**: Pathwise method shows lowest variance for Delta and Vega
- **Gamma Computation**: Only Finite Differences and Likelihood Ratio are applicable
- **Practical Recommendations**: Method selection based on specific use cases

### Sample Output

For a European call option with parameters:
- S₀ = 100 (initial stock price)
- K = 100 (strike price)
- T = 1 year (maturity)
- r = 5% (risk-free rate)
- σ = 20% (volatility)

The Monte Carlo estimates (50,000 paths) closely match Black-Scholes values:
- Price: ~10.45
- Delta: ~0.637
- Gamma: ~0.019
- Vega: ~37.52

## Technical Details

### Black-Scholes Model

The underlying asset follows a Geometric Brownian Motion:

```
dS_t = r S_t dt + σ S_t dW_t
```

With exact solution:

```
S_T = S₀ exp((r - σ²/2)T + σW_T)
```

### Variance Reduction

Antithetic variables technique:
- For each Z ~ N(0,1), also use -Z
- Reduces variance for monotonic payoffs
- Approximately 2x efficiency gain

## Project Structure

```
Options_Greeks_Monte_Carlo/
│
├── Greeks_Monte_Carlo.ipynb    # Main Jupyter notebook with all implementations
├── requirements.txt             # Python dependencies
├── README.md                    # This file
├── .gitignore                   # Git ignore rules
│
└── rapport_latex/               # LaTeX report (optional)
    ├── rapport.tex
    └── figures/
```

## Dependencies

- `numpy`: Numerical computations and random number generation
- `scipy`: Statistical functions (normal distribution)
- `matplotlib`: Plotting and visualizations
- `jupyter`: Interactive notebook environment

## References

The implementation is based on standard computational finance literature:
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*
- Jäckel, P. (2002). *Monte Carlo Methods in Finance*
- Black-Scholes-Merton option pricing model

## Author

**Dan Allouche** — Master MASEF, Université Paris-Dauphine

- GitHub: [@dan-allouche-qf](https://github.com/dan-allouche-qf)
