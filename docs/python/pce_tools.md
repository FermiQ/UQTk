# pce_tools.py

## Overview

This file, `pce_tools.py`, is part of the UQ Toolkit (UQTk). It provides a suite of Python functions for performing Polynomial Chaos Expansion (PCE) related computations. These tools facilitate tasks such as mapping random variables to PCE representations, evaluating PCE models, drawing samples from PCEs, performing Galerkin projections, regression, and Bayesian Compressive Sensing (BCS) for PCE coefficient determination. It also includes utilities for sensitivity analysis (Sobol indices) and Kernel Density Estimation (KDE).

## Key Components

*   **`UQTkMap2PCE(pc_model, rvs_in, verbose=0)`**: Obtains a PC representation for random variables described by samples. It uses a Rosenblatt transformation to map input RVs to the PC germ space.
*   **`UQTkEvalPC(pce_model, pce_coeffs, germ_sample)`**: **Deprecated.** Users are advised to use `UQTkEvaluatePCE` instead.
*   **`UQTkDrawSamplesPCE(pc_model, pc_coeffs, n_samples)`**: Draws samples of the underlying germ of a PC model and evaluates a given PCE for those samples.
*   **`UQTkEvaluatePCE(pc_model, pc_coeffs, samples)`**: Evaluates a PCE at a given set of samples of this PCE.
*   **`UQTkGalerkinProjection(pc_model, f_evaluations)`**: Obtains PC coefficients using Galerkin Projection via UQTk.
*   **`UQTkRegression(pc_model, f_evaluations, samplepts)`**: Obtains PC coefficients by regression.
*   **`UQTkBCS(pc_begin, xdata, ydata, eta=1.e-3, ...)`**: Performs Bayesian Compressive Sensing to obtain PC coefficients, potentially with basis growth.
*   **`UQTkOptimizeEta(pc_start, y, x, etas, niter, nfolds, ...)`**: Helper function for `UQTkBCS` to choose the optimum `eta` value via cross-validation.
*   **`UQTkEvalBCS(pc_model, f_evaluations, samplepts, sigma2, eta, ...)`**: Performs one iteration of Bayesian Compressive Sensing. Helper function for `UQTkBCS`.
*   **`UQTkCallBCSDirect(vdm_np, rhs_np, sigma2, eta=1.e-8, ...)`**: Calls C++ BCS routines directly with a Vandermonde matrix and right-hand side.
*   **`multidim_intersect(arr1, arr2)`**: Finds the intersection of two multi-dimensional NumPy arrays.
*   **`ind_split(ns, split_method, split_params)`**: Splits indices for cross-validation purposes.
*   **`UQTkGetQuadPoints(pc_model)`**: Generates quadrature points using UQTk and returns them as a NumPy array.
*   **`UQTkStDv(pc_model, pc_coeffs)`**: Computes the Standard Deviation of a PCE using UQTk.
*   **`UQTkGSA(pc_model, pc_coeffs)`**: Computes Sobol' sensitivity indices (main, total, and joint).
*   **`UQTkKDE(fcn_evals)`**: Performs Kernel Density Estimation on a set of function evaluations.
*   **`UQTkGetMultiIndex(pc_model, ndim)`**: Returns a 2D NumPy array of the PC multi-index.
*   **`UQTkPlotMiDims(pc_model, c_k, ndim, nord, type)`**: Creates a plot showing the magnitude of PC coefficients for each order.
*   **`kfold_split(nsamples, nfolds, seed=13)`**: Returns a dictionary of training and testing indices for k-fold cross-validation.
*   **`kfoldCV(x, y, nfolds=3, seed=13)`**: Splits data (x and y) into training/testing pairs for k-fold cross-validation.

## Important Variables/Constants

While the file primarily defines functions, some functions have important default parameters or internal constants that affect their behavior. For example:
*   **`UQTkBCS` parameters**: `eta` (stopping threshold), `niter` (iterations for order growth), `mindex_growth` (basis growth method), `sigma2` (initial noise variance).
*   **`UQTkCallBCSDirect` parameters**: `eta` (stopping threshold).
*   **Internal constants in `UQTkEvalBCS`**: `adaptive`, `optimal`, `scale` for BCS configuration.

## Usage Examples

```python
# Conceptual Example for UQTkEvaluatePCE
# Assuming 'pc_model' is a UQTk PC object, 'coeffs' is a NumPy array of PC coefficients,
# and 'samples_xi' is a NumPy array of samples in the germ space.

# pce_evaluations = UQTkEvaluatePCE(pc_model, coeffs, samples_xi)
# print(f"PCE evaluated at samples: {pce_evaluations}")

# Conceptual Example for UQTkGalerkinProjection
# Assuming 'pc_model' is a UQTk PC object and 'func_at_quad_pts' is a NumPy
# array of function evaluations at the PC model's quadrature points.

# pc_coefficients = UQTkGalerkinProjection(pc_model, func_at_quad_pts)
# print(f"Galerkin projected PC coefficients: {pc_coefficients}")
```

## Dependencies and Interactions

This file has several dependencies and interactions:

*   **`sys`**: Used to modify Python path for local UQTk module imports.
*   **`uqtkarray`**: Essential UQTk module for custom array types used throughout the functions.
*   **`quad` (as `uqtkquad`)**: UQTk module for quadrature-related functionalities.
*   **`pce` (as `uqtkpce`)**: Core UQTk module for PCE class definitions and basic operations.
*   **`tools` (as `uqtktools`)**: UQTk module providing various utility tools (e.g., Rosenblatt transformation).
*   **`bcs`**: UQTk module for Bayesian Compressive Sensing C++ backend.
*   **`utils.multiindex` (as `uqtkmi`)**: UQTk utility for multi-index operations.
*   **`numpy` (as `np`)**: Heavily used for numerical operations and array manipulations.
*   **`scipy.stats` and `math`**: Used for statistical operations (like KDE) and mathematical functions.
*   **`matplotlib.pyplot` and `matplotlib.rc`**: Used for plotting capabilities, specifically in `UQTkPlotMiDims` and `UQTkOptimizeEta`.
*   **`functools.reduce`**: Used in `UQTkBCS` for intersecting multi-indices.

The functions in this file often interact with each other. For instance, `UQTkBCS` uses `UQTkOptimizeEta` and `UQTkEvalBCS` as helper functions. Many functions rely on a `pc_model` object, which is typically an instance of a PCE class defined in `uqtkpce`.
```
