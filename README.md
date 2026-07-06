# Fluctuating environments favor extreme reproductive delays and penalize intermediate ones

**Authors:** Jorge Hidalgo, Lorenzo Fant, Rafael Rubio de Casas, and
Miguel A. Muñoz

This repository contains the code and data used in the paper\
**"Fluctuating environments favor extreme reproductive delays and penalize intermediate ones"**, arXiv:2512.05856 (2025).
It includes two main components:

1.  Numerical simulations of a **delayed stochastic differential equation** describing population dynamics under environmental fluctuations.
2.  An **agent-based evolutionary model** where reproductive delay is an evolving trait subject to mutation and natural selection.



# 1. `fixed-delay-stochastic-eq`: Delayed stochastic differential equation model

This folder contains C++ simulations of a population governed by a delayed logistic equation with environmental stochasticity. Results are analyzed in `figures.ipynb`.

## Model summary

The population density evolves as:

\[ `\dot{x}`{=tex}(t) = (b +
`\sigma `{=tex}`\xi`{=tex}\_`\tau`{=tex}(t)), x(t-`\alpha`{=tex})
`\left`{=tex}(1 - `\frac{x(t)}{K}`{=tex} `\right`{=tex}) - d, x(t). \]

Environmental noise 
is represented by a dichotomous Markov process with
autocorrelation:

\[ `\langle `{=tex}`\xi`{=tex}*`\tau`{=tex}(t),
`\xi`{=tex}*`\tau`{=tex}(s)`\rangle `{=tex}=
e\^{-2\|t-s\|/`\tau`{=tex}}. \]

### Compilation

Requires **GSL**. Compile with:

    g++ file.cpp -lm -lgsl -lgslcblas -o simulation

### Code structure

All simulation logic is implemented in `delayed_functions.h`. Each `.cpp` file is a
driver that calls functions from that header to perform a specific task:
| File | Task | Figures |
|------|------|---------|
| `example.cpp` | Minimal working example (linear growth rate) | — |
| `data-timeseries/timeseries.cpp` | Population time series *x(t)* and stationary PDF | Fig. 2 |
| `data-mean_population_diagram/mean_population_diagram_sigma.cpp` | Stationary mean population *x\** vs. α and σ | Fig. 3 |
| `data-linear_growth_rate/linear_growth_rate.cpp` | Mean growth rate *G* vs. delay α and noise amplitude σ | Figs. 4–5 |
| `data-mean_ext_times/mean_ext_times.cpp` | Mean extinction time *T* vs. system size *N* | Fig. 6 |
| `data-trade_off/trade_off.cpp` | Growth–extinction trade-off (*G* vs. *T*) | Fig. 7 |
| `data-appendix/b/linear_growth_rate_b.cpp` | Appendix: *G* vs. α for varying mean birth rate *b* | Appendix - Fig. 9|
| `data-mean_population_diagram/mean_population_diagram_tau.cpp` | Stationary mean population *x\** vs. α and τ | Appendix - Fig. 10 |
| `data-appendix/ou/linear_growth_rate_OU.cpp` | Appendix: *G* vs. α with Ornstein–Uhlenbeck noise | Appendix - Fig. 11 |
| `data-appendix/stochastic_d/linear_growth_rate_stochastic_d.cpp` | Appendix: *G* vs. α with stochastic death rate | Appendix - Fig. 12 |

### Usage

Compile and run the provided example:

    g++ example.cpp -lm -lgsl -lgslcblas -o simulation
    ./simulation example_params.dat

Output is written to `output_example_params.dat` in the same directory.

# 2. `evolutionary-algorithm`: Agent-based evolutionary model

This folder contains the Julia implementation of an agent-based evolutionary model in which the reproductive delay is a heritable trait subject to mutation and natural selection.
Results are plot with plot_figure.py (Fig. 8 in the manuscript).

## Model summary

Agents carry a heritable activation rate α (inverse of reproductive delay). Resources fluctuate between "good" and "bad" states as a dichotomous Markov process with correlation
time τ. Active agents consume resources, reproduce  roportionally to intake, and generate dormant propagules; mortality affects only active agents. Mutations shift α by a small increment at reproduction. All events are simulated in continuous time via a Gillespie algorithm.
## Code structure

- `resources_space.jl` — Core model definition: agent types 
  model parameters, initialization, and stepping functions (built on the
  [Agents.jl](https://juliadynamics.github.io/Agents.jl/) framework).
- `parallel_time_correlation.jl` — Driver script: runs independent evolutionary
  simulations in parallel across a grid of initial α values and saves results to CSV.

## Usage

Requires Julia with the packages listed in the project environment. To run:

    julia parallel_time_correlation.jl

This launches 16 parallel workers, sweeps over initial activation rates, and writes
output CSV files to the `data/` folder with a timestamp.

# Citation

If you use this repository, please cite:

> J. Hidalgo, L. Fant, R. Rubio de Casas, and M. A. Muñoz,\
> *Fluctuating environments favor extreme reproductive delays and penalize intermediate ones*, arXiv:2512.05856 (2025).
