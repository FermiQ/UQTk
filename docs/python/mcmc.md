# mcmc.py

## Overview

This file, `mcmc.py`, is part of the UQ Toolkit (UQTk) and provides implementations of Markov Chain Monte Carlo (MCMC) algorithms. These algorithms are used for sampling from probability distributions, particularly posterior distributions in Bayesian inference. The file includes a Hamiltonian MCMC (HMCMC) sampler and a more advanced Delayed Rejection Adaptive MCMC (DRAM) sampler. It also contains helper functions for a specific "banana-shaped" example posterior, used for testing and demonstrating the MCMC samplers.

## Key Components

*   **`HMCMC(U, grad_U, dt, nT, q)`**: Implements the Hamiltonian Monte Carlo (or Hybrid Monte Carlo) algorithm.
    *   `U`: A function that computes the potential energy (typically -log(posterior)).
    *   `grad_U`: A function that computes the gradient of the potential energy.
    *   `dt`: Time step for the leapfrog integrator.
    *   `nT`: Number of time steps in the leapfrog trajectory.
    *   `q`: The initial state (position vector) of the chain.
    *   *Description*: HMCMC uses Hamiltonian dynamics to propose new states, potentially leading to more efficient exploration of the state space compared to simpler MCMC methods, especially in higher dimensions. It uses a leapfrog integrator for simulating the dynamics.

*   **`dram(opts, cini, likTpr, lpinfo)`**: Implements the Delayed Rejection Adaptive MCMC (DRAM) algorithm.
    *   `opts`: A dictionary containing various options to control the sampler's behavior (see "Important Variables/Constants" for `opts` dictionary keys).
    *   `cini`: The initial state (1D NumPy array) of the MCMC chain.
    *   `likTpr`: A user-supplied function that computes the log-likelihood and log-prior. It takes the current sample and `lpinfo` as input and should return a list or tuple `[log_likelihood, log_prior]`.
    *   `lpinfo`: A user-supplied object containing any additional information or parameters needed by the `likTpr` function.
    *   *Description*: DRAM combines adaptive MCMC (where the proposal distribution is adapted based on the chain's history) with delayed rejection (which gives a second chance to proposals that would have otherwise been rejected). This can improve sampling efficiency and robustness.

*   **`norm_pdf_multivariate(x, mu, sigma)`**: Helper function to calculate the probability density function (PDF) of a multivariate normal distribution.
*   **`tranB(x1, x2, a)` / `invTranB(x1, x2, a)`**: Helper functions for coordinate transformation used in the banana-shaped PDF example.
*   **`plotBanana()`**: Helper function to plot the banana-shaped PDF.
*   **`postBanana(spl, postinfo)`**: Helper function to compute the log-posterior for the banana-shaped PDF example, designed to be used with the MCMC samplers.
*   **`dram_ex(method, nsteps)`**: An example function demonstrating how to use the `dram` sampler with the `postBanana` posterior.
*   **`logPropRatio(iq, spls)`**: Internal helper for `dram` to compute the log proposal ratio for delayed rejection stages.
*   **`logPostRatio(p1, p2)`**: Internal helper for `dram` to compute the log posterior ratio.
*   **`getAlpha(spls, post)`**: Internal helper for `dram` to compute the acceptance probability for delayed rejection stages.
*   **`ucov(spl, splmean, cov, lastup)`**: Internal helper for `dram` to update the covariance matrix during adaptation.

## Important Variables/Constants

*   **`opts` (Dictionary for `dram` function)**: This dictionary controls the behavior of the DRAM sampler. Key parameters include:
    *   `'method'`: `'am'` (Adaptive Metropolis) or `'dram'` (Delayed Rejection Adaptive MCMC).
    *   `'nsteps'`: Total number of MCMC steps.
    *   `'nburn'`: Number of burn-in steps (proposal covariance is fixed).
    *   `'nadapt'`: Adapt proposal covariance every `nadapt` steps after burn-in.
    *   `'nfinal'`: Stop adapting covariance after `nfinal` steps.
    *   `'inicov'`: Initial proposal covariance matrix.
    *   `'coveps'`: Small epsilon added to diagonal of covariance for numerical stability.
    *   `'burnsc'`: Factor to scale proposal during burn-in if acceptance is too high/low.
    *   `'gamma'`: Factor to multiply proposed jump size after burn-in (default 1.0).
    *   `'ndr'`: Number of delayed rejection stages (if `method='dram'`).
    *   `'drscale'`: List of scale factors for proposal covariance at each DR stage.
    *   `'spllo'`, `'splhi'`: Lower and upper bounds for chain samples.
    *   `'rnseed'`: Optional seed for the random number generator.
    *   `'tmpchn'`: Optional filename for saving intermediate chain states.
    *   `'ofreq'`: Frequency of saving intermediate chain states if `tmpchn` is specified.

*   **`Rmat`, `invRmat` (Global variables for `dram`)**: These store the Cholesky decomposition of proposal covariances and their inverses for different delayed rejection stages. They are managed internally by the `dram` function.

## Usage Examples

```python
# Conceptual Example for HMCMC (details depend on U and grad_U)

# Define potential energy function U(q) -> float
# def potential_energy(q):
#     # -log(posterior(q))
#     return -(-(q[0]**2)/2.0 - (q[1]**2)/2.0) # Example: standard normal

# Define gradient of potential energy grad_U(q) -> np.array
# def gradient_potential_energy(q):
#     return np.array([q[0], q[1]]) # Example: standard normal

# initial_q = np.array([0.0, 0.0])
# dt_hmcmc = 0.1
# nT_hmcmc = 10
# num_samples_hmcmc = 1000
# samples_hmcmc = [initial_q]

# for _ in range(num_samples_hmcmc - 1):
#     next_q = HMCMC(potential_energy, gradient_potential_energy, dt_hmcmc, nT_hmcmc, samples_hmcmc[-1])
#     samples_hmcmc.append(next_q)

# print(f"Generated {len(samples_hmcmc)} samples using HMCMC.")

# Example for DRAM (using the built-in banana example)

# Number of steps for the DRAM sampler
# n_dram_steps = 5000

# Run the DRAM example with 'dram' method
# dram_results = dram_ex(method='dram', nsteps=n_dram_steps)

# Accessing results from dram_results dictionary:
# chain_samples = dram_results['chain']
# map_estimate = dram_results['cmap']
# acceptance_rate = dram_results['accr']
# print(f"Generated {chain_samples.shape[0]} samples using DRAM.")
# print(f"MAP estimate: {map_estimate}")
# print(f"Acceptance rate: {acceptance_rate}")

# To use dram with a custom model, one would define:
# 1. cini_custom = np.array([...]) # Initial guess
# 2. opts_custom = {
#        'method': 'dram', 'nsteps': 10000, 'nburn': 2000, 
#        'nadapt': 200, 'inicov': np.eye(num_dimensions),
#        'spllo': np.array([...]), 'splhi': np.array([...]), ... 
#    } # Options
# 3. def custom_log_posterior(params, model_info):
#        # ... calculate log_likelihood and log_prior based on 'params' and 'model_info'
#        # return [log_likelihood, log_prior]
#    model_specific_info = { ... } # Any info needed by custom_log_posterior
# 4. custom_results = dram(opts_custom, cini_custom, custom_log_posterior, model_specific_info)
```

## Dependencies and Interactions

*   **`numpy` (as `npy`)**: Essential for numerical operations, array manipulations (samples, covariance matrices).
*   **`scipy.stats`**: Used by the example `norm_pdf_multivariate`, though not directly by the core MCMC algorithms.
*   **`scipy.linalg`**: Used for Cholesky decomposition (`scipy.linalg.cholesky`) and matrix inversion (`scipy.linalg.inv`) in the `dram` sampler for handling covariance matrices.
*   **`math`**: Standard math functions.
*   **`uuid`**: Used by `dram` to generate unique filenames if temporary chain saving is enabled without a specific name.
*   **`matplotlib.pyplot` (as `plt`)**: Used by the `plotBanana` example function, not by the MCMC algorithms themselves.

**Interactions:**
*   The MCMC algorithms (`HMCMC`, `dram`) require user-defined functions for the (log) posterior probability (or potential energy and its gradient for HMCMC).
*   `dram` internally uses helper functions like `ucov` for adaptive covariance updates and `logPropRatio`, `logPostRatio`, `getAlpha` for delayed rejection logic.
*   The example functions (`postBanana`, `dram_ex`, `plotBanana`, etc.) demonstrate how to set up and use the `dram` sampler for a specific problem.
```
