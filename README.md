# Fluctuating environments favor extreme dormancy strategies and penalize intermediate ones

**Authors:** Jorge Hidalgo, Lorenzo Fant, Rafael Rubio de Casas, and
Miguel A. Muñoz

This repository contains the code and data used in the paper\
**"Fluctuating environments favor extreme dormancy strategies and
penalize intermediate ones"**.\
It includes two main components:

1.  Numerical simulations of a **delayed stochastic differential equation** describing population dynamics under environmental fluctuations.
2.  An **agent-based evolutionary model** where dormancy is an evolving trait subject to mutation and natural selection.

Both approaches support the central result of the paper: **variable environments select for extreme dormancy strategies while penalizing intermediate ones**.


# 1. `fixed-delay-stochastic-eq`: Delayed stochastic differential equation model

This folder contains C++ simulations of a population governed by a delayed logistic equation with environmental stochasticity. The trajectories are analyzed in `figures.ipynb`.


## Model summary

The population density evolves as:

$$
\dot{x}(t) = \left(b + \sigma \, \xi_\tau(t)\right) x(t-\alpha)
\left(1 - \frac{x(t)}{K}\right) - d \, x(t)
$$

Environmental noise is represented by a dichotomous Markov process with
autocorrelation:

$$
\langle \xi_\tau(t), \xi_\tau(s) \rangle =
e^{-2|t-s|/\tau}
$$

### Compilation

Requires **GSL**. Compile with:

    g++ file.cpp -lm -lgsl -lgslcblas -o simulation

# 2. `evolutionary-algorithm`: Agent-based evolutionary model

This folder contains simulations of an explicit ABM,implemented in Julia, where dormancy duration evolves through mutation and natural selection under a fluctuating environment. All ecological and evolutionary events are implemented via a **Gillespie algorithm**.

# Citation

If you use this repository, please cite:

> J. Hidalgo, L. Fant, R. Rubio de Casas, and M. A. Muñoz,\
> *Fluctuating environments favor extreme dormancy strategies and penalize intermediate ones* (2025).
