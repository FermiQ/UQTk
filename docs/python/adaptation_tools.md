# adaptation_tools.py

## Overview

This file, `adaptation_tools.py`, is a component of the UQ Toolkit (UQTk). It focuses on providing tools for dimensionality reduction and adaptation in Polynomial Chaos Expansions (PCE). The main idea is to find a rotation in the input parameter space such that most of the variance of the Quantity of Interest (QoI) is captured by a lower-dimensional manifold. This allows for more efficient PCE representations, especially for problems where the QoI is primarily sensitive to a few linear combinations of the original input variables.

## Key Components

*   **`gauss_adaptation(c_k, ndim, method = 0)`**: Computes an isometry (rotation matrix) based on first-order PCE coefficients (`c_k`) of a Gaussian system. This matrix is used to rotate the input space.
    *   `c_k`: 1D NumPy array of first-order PCE coefficients.
    *   `ndim`: Dimension of the problem.
    *   `method`: Integer (0-3) specifying the numerical method to compute the isometry (0: Gram-Schmidt on A with Gaussian coeffs and identity, 1: Orthogonal decomposition of `a*a.T`, 2: Orthogonal decomposition of Householder matrix, 3: Gram-Schmidt with sorted coeffs - recommended).
*   **`eta_to_xi_mapping(eta, A, zeta = None)`**: Maps points from a lower-dimensional space (`eta`) to the original higher-dimensional input space (`xi`) using the rotation matrix `A`.
    *   `eta`: N x d0 NumPy array of points in the reduced-dimension space.
    *   `A`: d x d NumPy array, the rotation matrix (isometry) from `gauss_adaptation`.
    *   `zeta`: (Optional) N x (d-d0) NumPy array to augment `eta` if `d0 < d`. Defaults to zeros.
*   **`mi_terms_loc(d1, d2, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0)`**: Locates the basis terms of a lower-dimensional PCE (`d1`) within the multi-index set of a higher-dimensional PCE (`d2`). This is useful for projecting coefficients between spaces of different dimensionality.
    *   `d1`: Lower dimension.
    *   `d2`: Higher dimension.
    *   `nord`: PC order.
    *   `pc_type`: Polynomial type (e.g., "HG" for Hermite).
    *   `param`: Quadrature level or parameter.
    *   `sf`: Sparsity flag ("sparse" or "full").
    *   `pc_alpha`, `pc_beta`: Parameters for the polynomial basis (e.g., for Jacobi).
*   **`l2_error_eta(c_1, c_2, d1, d2, nord, pc_type, param, sf, pc_alpha=0.0, pc_beta=1.0)`**: Calculates the relative L2-norm error between PCE coefficients from a lower-dimensional expansion (`c_1`) and a higher-dimensional expansion (`c_2`), after projecting `c_1` into the higher-dimensional space.
    *   `c_1`: Coefficients of the lower-dimensional PCE.
    *   `c_2`: Coefficients of the higher-dimensional PCE.
    *   Other parameters are similar to `mi_terms_loc`.
*   **`transf_coeffs_xi(coeffs, nord, ndim, pc_type, param, R, sf="sparse", pc_alpha=0.0, pc_beta=1.0)`**: Transforms PCE coefficients from the reduced (`eta`) space back to the original (`xi`) space. This is primarily for assessing the accuracy of the adaptation method.
    *   `coeffs`: PCE coefficients in the `eta` space.
    *   `R`: The rotation matrix.
    *   Other parameters are similar to `mi_terms_loc`.

## Important Variables/Constants

The functions in this file are primarily driven by their input arguments. There are no standalone global constants defined in this file that dictate general behavior outside of the function calls themselves. Key parameters within functions include:
*   **`method` in `gauss_adaptation`**: Controls the algorithm for computing the rotation matrix. `method = 3` is noted as the recommended approach.
*   **`pc_type`, `param`, `sf` in `mi_terms_loc`, `l2_error_eta`, `transf_coeffs_xi`**: These define the characteristics of the Polynomial Chaos Expansion being used (e.g., Hermite polynomials, quadrature rule parameters, sparse/full grid).

## Usage Examples

```python
# Conceptual Example for using gauss_adaptation and eta_to_xi_mapping

# Assume 'first_order_coeffs' are the 1st order PCE coefficients for a 3D problem
# first_order_coeffs = np.array([0.8, 0.1, 0.05]) 
# ndim = 3

# 1. Compute the rotation matrix
# R = gauss_adaptation(first_order_coeffs, ndim, method=3)

# 2. Define points in a reduced 1D eta space (e.g., quadrature points)
# eta_points = np.array([[0.5], [1.0], [1.5]]) # N x d0, here d0=1

# 3. Map eta points to the original xi space
# xi_points = eta_to_xi_mapping(eta_points, R)
# print("Mapped xi points:", xi_points)

# Conceptual Example for transf_coeffs_xi
# Assume 'eta_coeffs' are PCE coeffs in a 2D eta-space derived from a 3D xi-space
# nord = 3
# ndim_eta = 2 
# ndim_xi = 3
# pc_type = "HG" # Hermite polynomials
# param = 4 # Quadrature level
# R_matrix = # ... rotation matrix used for adaptation ...

# To transform these eta_coeffs to equivalent xi_coeffs (for comparison/validation):
# This function internally performs a change of basis.
# Note: For transf_coeffs_xi, the 'ndim' argument should be the dimension of the space
# the 'coeffs' are currently in, and 'R' facilitates the transformation to a space of the
# same dimension but different orientation. If the intention is to project from a lower-dim
# eta-space to a higher-dim xi-space, a combination of mi_terms_loc and direct coefficient
# assignment would be used. The `transf_coeffs_xi` function as written appears to assume
# `coeffs` are for an `ndim` dimensional space that is being rotated.

# For checking convergence (comparing coefficients from a reduced model to a full model):
# d1_coeffs = # ... coeffs from a d1-dimensional adapted model ...
# d_full_coeffs = # ... coeffs from the original d-dimensional model ...
# error, projected_d1_coeffs = l2_error_eta(d1_coeffs, d_full_coeffs, d1, d, nord, ...)
# print(f"Relative L2 error: {error}")
```

## Dependencies and Interactions

*   **`sys`**: Used to modify the Python path for local UQTk module imports.
*   **`uqtkarray`**: UQTk module for custom array types, essential for passing data to and from other UQTk C++ routines.
*   **`quad` (as `uqtkquad`)**: UQTk module for quadrature rules, likely used by `PCSet` internally.
*   **`pce` (as `uqtkpce`)**: Core UQTk module providing the `PCSet` class, which defines the PCE basis, multi-indices, and quadrature rules. Functions in `adaptation_tools.py` extensively use `PCSet` objects.
*   **`tools` (as `uqtktools`)**: General UQTk tools.
*   **`pce_tools` (from `PyUQTk.PyPCE` or local path)**: Provides helper functions for PCE operations, such as `UQTkGetQuadPoints`. `transf_coeffs_xi` explicitly uses `pce_tools.UQTkGetQuadPoints`.
*   **`numpy` (as `np`)**: Foundational library for numerical computation in Python. All data arrays (coefficients, points, matrices) are handled as NumPy arrays. Functions like `np.linalg.qr`, `np.linalg.eigh`, `np.dot`, `np.hstack`, `np.where`, `np.all`, `np.linalg.norm` are used.

**Interactions:**
*   The workflow typically starts with `gauss_adaptation` to find a rotation.
*   This rotation is then used in `eta_to_xi_mapping` to transform sample points between the original and reduced spaces.
*   Functions like `mi_terms_loc` and `l2_error_eta` are used to analyze and compare PCEs constructed in these different spaces (e.g., a full PCE in `xi` space vs. an adapted PCE in `eta` space).
*   `transf_coeffs_xi` is a specialized tool for transforming coefficients under rotation, useful for validating the adaptation process.
*   All functions rely on `PCSet` objects (from `uqtkpce`) to define the PCE context (polynomial type, order, dimension, quadrature).
```
