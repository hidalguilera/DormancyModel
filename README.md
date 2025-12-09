# Fluctuating environments favor extreme dormancy strategies and penalize intermediate ones

**Authors:** Jorge Hidalgo, Lorenzo Fant, Rafael Rubio de Casas, and
Miguel A. Muñoz

This repository contains the code and data used in the paper\
**"Fluctuating environments favor extreme dormancy strategies and penalize intermediate ones"**, [arXiv:2512.05856 (2025)](https://arxiv.org/abs/2512.05856)


It includes two main components:

1.  Numerical simulations of a **delayed stochastic differential equation** describing population dynamics under environmental fluctuations.
2.  An **agent-based evolutionary model** where dormancy is an evolving trait subject to mutation and natural selection.



# 1. stochastic-equation

This folder contains C++ simulations of a population governed by a delayed logistic equation with environmental stochasticity. The trajectories are analyzed in `figures.ipynb`.


In this model, the population density evolves as

$$
\frac{dx}{dt} = \left(b + \sigma  \xi_\tau(t)\right) x(t-\alpha)
\left(1 - \frac{x(t)}{K}\right) - d  x(t).
$$

Environmental noise is represented by a dichotomous Markov process that alternates between $\xi_\tau=\pm 1$ at constant rate $\tau^{-1}$.

## Compilation

Requires **GSL**. Compile with:

    g++ file.cpp -lm -lgsl -lgslcblas -o simulation

# 2. evolutionary-algorithm

This folder contains simulations of an explicit agent based model,implemented in Julia, where dormancy duration evolves through mutation and natural selection under a fluctuating environment. All ecological and evolutionary events are implemented via a **Gillespie algorithm**.

# Citation

If you use this repository, please cite:

> J. Hidalgo, L. Fant, R. Rubio de Casas, and M. A. Muñoz,  
> *Fluctuating environments favor extreme dormancy strategies and penalize intermediate ones*,  
> [arXiv:2512.05856 (2025)](https://arxiv.org/abs/2512.05856)

