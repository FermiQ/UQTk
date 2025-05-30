# quad.h / quad.cpp

## Overview

The `quad.h` file and its implementation `quad.cpp` define the `Quad` class within the UQ Toolkit (UQTk). This class is responsible for generating various types of quadrature rules, which are sets of points and corresponding weights used for numerical integration. It supports 1D and multi-dimensional quadrature rules, including full tensor products and sparse grids. A variety of rule types are available, such as Gauss-Legendre, Gauss-Hermite, Clenshaw-Curtis, Newton-Cotes, and rules based on custom probability density functions (via recursion coefficients).

## Key Components

*   **`Quad` (Class)**:
    *   **Constructors**:
        *   `Quad(char *grid_type, char *fs_type, int ndim, int param, double alpha=0.0, double betta=1.0)`: For isotropic quadrature rules (same rule type and parameters in all dimensions).
            *   `grid_type`: String for rule type (e.g., "LU", "HG", "CC").
            *   `fs_type`: String for "full" tensor product or "sparse" grid.
            *   `ndim`: Number of dimensions.
            *   `param`: Number of points per dimension (for "full") or level (for "sparse").
            *   `alpha`, `betta`: Parameters for certain PC types (e.g., Laguerre, Jacobi).
        *   `Quad(Array1D<string>& grid_types, char *fs_type, Array1D<int>& param, Array1D<double>& alphas, Array1D<double>& betas)`: For anisotropic quadrature rules (dimension-specific types and parameters).
            *   `grid_types`: Array of strings for rule types per dimension.
            *   `fs_type`: "full" or "sparse".
            *   `param`: Array of points per dimension or level parameters.
            *   `alphas`, `bettas`: Arrays of alpha/beta parameters per dimension.
        *   `Quad()`: Default empty constructor.
    *   **`SetRule()`**: The primary method that computes and sets the quadrature points and weights based on the parameters provided in the constructor or subsequent setter methods.
    *   **`GetRule(Array2D<double>& q, Array1D<double>& w)`**: Retrieves the computed quadrature points (`q`) and weights (`w`).
    *   **`SetDomain(Array1D<double>& aa, Array1D<double>& bb)` / `SetDomain(Array1D<double>& aa)`**: Sets the integration domain endpoints. For compact domains (e.g., Legendre), `aa` and `bb` define lower and upper bounds. For semi-infinite domains (e.g., Laguerre), `aa` might define a lower bound.
    *   **`GetQdpts(Array2D<double>& q)` / `GetWghts(Array1D<double>& w)`**: Retrieve quadrature points and weights separately.
    *   **`SetQdpts(Array2D<double>& q)` / `SetWghts(Array1D<double>& w)`**: Allows externally setting quadrature points and weights.
    *   **`GetNQ()`**: Returns the total number of quadrature points.
    *   **`SetVerbosity(int verbosity)`**: Sets the verbosity level for output during rule generation.
    *   **`nextLevel()`**: (Primarily for sparse grids) Computes the points and weights for the next hierarchical level.
    *   **`SetAlpha(double alpha)` / `SetBeta(double betta)`**: Set alpha/beta parameters for the rule.
    *   **`SetLevel(int param)`**: Sets the level parameter (relevant for sparse grids).
    *   **`init()`**: Internal initialization function called by constructors.

## Important Variables/Constants

*   **`QD_MAX`**: A preprocessor macro defining a cap on the maximum number of quadrature points for full tensor-product rules to prevent excessive memory allocation (approx. 20,000,000).
*   **`QuadRule rule_` (Private Member)**: A struct holding the `Array2D<double> qdpts` (quadrature points) and `Array1D<double> wghts` (quadrature weights). This is the primary storage for the generated rule.
*   **`string grid_type_` / `Array1D<string> grid_types_`**: Stores the type(s) of quadrature rule (e.g., "LU", "HG").
*   **`string fs_type_`**: Stores the sparseness type ("full" or "sparse").
*   **`int ndim_`**: Stores the number of dimensions.
*   **`int maxlevel_` / `Array1D<int> param_`**: Stores the level or points-per-dimension parameter(s).
*   **`double alpha_`, `double beta_` / `Array1D<double> alphas_`, `Array1D<double> betas_`**: Stores the alpha and beta parameters for relevant PC/quadrature types.
*   **`Array1D<double> aa_`, `Array1D<double> bb_`**: Stores the domain endpoints.
*   **`int quadverbose_`**: Controls verbosity.
*   **`Array1D<int> growth_rules_`**: Internal variable determining how the number of points grows with level for different 1D rule types in sparse grids.

## Usage Examples

```cpp
#include "quad.h"
#include "Array1D.h"
#include "Array2D.h"
#include <iostream>
#include <string>
#include <vector> // For std::vector in example

int main() {
    // Example 1: 1D Gauss-Legendre Quadrature, 5 points
    char grid_type_lu[] = "LU";
    char fs_type_full[] = "full";
    int ndim_1d = 1;
    int npoints_1d = 5;
    Quad quad1D(grid_type_lu, fs_type_full, ndim_1d, npoints_1d);
    quad1D.SetRule(); // Compute the rule

    Array2D<double> q_pts_1d;
    Array1D<double> q_wts_1d;
    quad1D.GetRule(q_pts_1d, q_wts_1d);

    std::cout << "1D Gauss-Legendre Quadrature (" << quad1D.GetNQ() << " points):" << std::endl;
    for (int i = 0; i < quad1D.GetNQ(); ++i) {
        std::cout << "Point: " << q_pts_1d(i, 0) << ", Weight: " << q_wts_1d(i) << std::endl;
    }

    // Example 2: 2D Full Tensor Product of Gauss-Hermite (3 points) and Gauss-Legendre (2 points)
    Array1D<std::string> grid_types_2d(2);
    grid_types_2d(0) = "HG";
    grid_types_2d(1) = "LU";

    Array1D<int> params_2d(2);
    params_2d(0) = 3; // 3 points for HG
    params_2d(1) = 2; // 2 points for LU
    
    Array1D<double> alphas_2d(2, 0.0); // Default alpha
    Array1D<double> betas_2d(2, 1.0);  // Default beta

    Quad quad2D_aniso(grid_types_2d, fs_type_full, params_2d, alphas_2d, betas_2d);
    quad2D_aniso.SetRule();

    Array2D<double> q_pts_2d;
    Array1D<double> q_wts_2d;
    quad2D_aniso.GetRule(q_pts_2d, q_wts_2d);

    std::cout << "\n2D Anisotropic Full Tensor Quadrature (" << quad2D_aniso.GetNQ() << " points):" << std::endl;
    // for (int i = 0; i < quad2D_aniso.GetNQ(); ++i) {
    //     std::cout << "Point: (" << q_pts_2d(i, 0) << ", " << q_pts_2d(i, 1) 
    //               << "), Weight: " << q_wts_2d(i) << std::endl;
    // }

    // Example 3: 2D Sparse Grid Clenshaw-Curtis, level 2
    char grid_type_cc[] = "CC";
    char fs_type_sparse[] = "sparse";
    int ndim_2d_sparse = 2;
    int level_sparse = 2; 
    Quad quadSparse(grid_type_cc, fs_type_sparse, ndim_2d_sparse, level_sparse);
    // For sparse grids, SetDomain might be important if not default [-1,1]
    // Array1D<double> domain_a(ndim_2d_sparse, -1.0);
    // Array1D<double> domain_b(ndim_2d_sparse, 1.0);
    // quadSparse.SetDomain(domain_a, domain_b);
    quadSparse.SetRule();

    Array2D<double> q_pts_sparse;
    Array1D<double> q_wts_sparse;
    quadSparse.GetRule(q_pts_sparse, q_wts_sparse);
    std::cout << "\n2D Sparse Clenshaw-Curtis Level " << level_sparse << " (" << quadSparse.GetNQ() << " points):" << std::endl;
    // Outputting all points can be lengthy for sparse grids.
    
    return 0;
}

```

## Dependencies and Interactions

*   **`Array1D.h`, `Array2D.h`**: Uses UQTk's `Array1D` and `Array2D` classes extensively for storing quadrature points, weights, domain boundaries, parameters, and internal working arrays.
*   **`error_handlers.h`**: For UQTk's exception handling mechanism (e.g., `Tantrum`).
*   **`combin.h`**: Contains combinatorial functions like `choose`, used in `getMultiIndexLevel` for sparse grid construction.
*   **`multiindex.h`**: (Likely, though not explicitly shown in `quad.cpp` direct includes, it's related to `getMultiIndexLevel` logic or used by `combin.h`). Used for managing multi-indices in sparse grid construction.
*   **`gq.h`**: Contains implementations of 1D Gauss quadrature rules (e.g., `gq`, `gq_gen`, `vandermonde_gq`) which are called by the `create1DRule_*` methods. This is a core dependency for generating the fundamental 1D rules.
*   **`arrayio.h`**: Used by `create1DRule_pdf` to read recursion coefficients from "ab.dat" and by `create1DRule_GP3` to read Gauss-Patterson rule data.
*   **`arraytools.h`**: Contains utility functions for array manipulations (e.g., `array1Dto2D`, `merge`, `paddMatCol`, `getRow`, `Trans`, `subMatrix_row`, `subMatrix_col`, `getCol`, `is_equal`), used in rule combination and compression.
*   **Standard C++ Libraries**: `<math.h>` (or `<cmath>`), `<assert.h>`, `<cfloat>`, `<iostream>`, `<string.h>` (or `<cstring>`), `<stdio.h>`, `<sstream>`.

**Interactions:**
*   The `Quad` class is often used by other UQTk components that require numerical integration, such as those performing Galerkin projection in Polynomial Chaos Expansions (e.g., `PCSet` might use `Quad` to define its integration rule).
*   The `SetRule()` method orchestrates calls to various private helper functions to first create 1D rules (`create1DRule_*`) and then combine them for multi-dimensional rules (e.g., `MultiplyTwoRules`, `AddTwoRules`, `SubtractTwoRules` for sparse grids).
*   The `compressRule()` method is important for sparse grids to merge duplicate quadrature points that arise from the combination formula and sum their weights.
```
