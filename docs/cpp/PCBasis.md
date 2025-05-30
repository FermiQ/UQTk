# PCBasis.h / PCBasis.cpp

## Overview

The `PCBasis.h` file and its corresponding implementation in `PCBasis.cpp` define the `PCBasis` class. This class is a fundamental component of the UQ Toolkit (UQTk) for Polynomial Chaos Expansions (PCE). It encapsulates the properties and operations related to univariate orthogonal polynomial bases. The class handles the initialization of different polynomial types (e.g., Hermite-Gaussian, Legendre-Uniform), calculation of quadrature points and weights, evaluation of basis polynomials and their derivatives at specified points, computation of basis polynomial norms, and generation of random samples according to the underlying probability distribution of the basis.

## Key Components

*   **`PCBasis` (Class)**:
    *   **Constructor `PCBasis(const string type, const double alpha, const double betta, const int maxord)`**: Initializes a univariate polynomial chaos basis.
        *   `type`: A string specifying the basis type (e.g., "HG" for Hermite-Gaussian, "LU" for Legendre-Uniform, "LG" for Laguerre-Gamma, "JB" for Jacobi-Beta, "SW" for Stieltjes-Wigert, "pdf" for custom from `ab.dat`).
        *   `alpha`, `betta`: Parameters for certain basis types (e.g., Laguerre, Jacobi).
        *   `maxord`: The maximum polynomial order for which basis functions and related quantities (like quadrature rules) are pre-computed.
    *   **`Init1dQuadPoints(int qdpts)`**: Initializes 1D quadrature points and weights appropriate for the selected basis type and order. The number of points `qdpts` is typically related to `maxord`.
    *   **`Eval1dBasisAtQuadPoints()`**: Evaluates the 1D basis polynomials at the pre-computed quadrature points and stores them.
    *   **`Eval1dBasisAtCustPoints(Array2D<double>& psi, int kord, const Array1D<double>& custPoints)`**: Evaluates basis polynomials up to order `kord` at user-specified custom points.
    *   **`EvalBasis(const double &xi, Array1D<double> &basisEvals) const`**: Evaluates all basis polynomials up to `maxord_` at a single point `xi`.
    *   **`EvalBasis(const double &xi, const int kord, double *basisEvals) const`**: Evaluates basis polynomials up to order `kord` at a single point `xi`, storing results in a raw pointer array.
    *   **`EvalDerivBasis(const double& xi, Array1D<double>& basisDEvals)`**: Evaluates the first derivatives of Legendre basis polynomials at `xi`. (Note: Currently implemented only for "LU" type).
    *   **`Eval1dDerivBasisAtCustPoints(Array2D<double>& dpsi, int kord, const Array1D<double>& custPoints)`**: Evaluates first derivatives of Legendre basis polynomials at custom points. (Note: Currently implemented only for "LU" type).
    *   **`Eval2ndDerivBasis(const double& xi, Array1D<double>& ddP)`**: Evaluates the second derivatives of Legendre basis polynomials at `xi`. (Note: Currently implemented only for "LU" type).
    *   **`Eval2ndDerivCustPoints(Array2D<double>& psi, int kord, Array1D<double>& custPoints)`**: Evaluates second derivatives of Legendre basis polynomials at custom points. (Note: Currently implemented only for "LU" type).
    *   **`Eval1dNormSq_Exact(int kord)`**: Computes the exact squared norms of the 1D basis polynomials up to order `kord`.
    *   **`Get1dNormsSq(Array1D<double>& psi1dSq) const`**: Returns the numerically computed squared norms of the basis functions (via quadrature).
    *   **`Get1dNormsSqExact(Array1D<double>& psi1dSqExact) const`**: Returns the analytically computed exact squared norms of the basis functions.
    *   **`GetRandSample(Array1D<double>& randSamples)` / `GetRandSample(double* randSamples, const int& nSamp)`**: Generates random samples from the probability distribution associated with the PC basis.
    *   **`SeedRandNumGen(const int& seed)`**: Seeds the internal random number generator.
    *   **`GetQuadRule(Array2D<double>& qPoints, Array1D<double>& qWeights, Array2D<int>& qIndices)`**: Retrieves the computed quadrature points, weights, and indices.
    *   **`GetQuadPoints(Array2D<double>& quadPoints) const`**: Retrieves the quadrature points.
    *   **`GetQuadWeights(Array1D<double>& quadWeights) const`**: Retrieves the quadrature weights.
    *   **`GetBasisAtQuadPoints(Array2D<double>& psi1d) const`**: Retrieves the basis polynomials evaluated at quadrature points.
    *   **`GetPCType() const`**: Returns the string type of the PC basis (e.g., "HG", "LU").
    *   **`GetAlpha() const` / `GetBeta() const`**: Return the alpha/beta parameters of the basis.
    *   **`GetSeed() const`**: Returns the current seed of the random number generator.

## Important Variables/Constants

*   **`string type_`**: Private member storing the type of the polynomial basis (e.g., "HG", "LU").
*   **`Array2D<double> quadPoints_`**: Private member storing the 1D quadrature points.
*   **`Array1D<double> quadWeights_`**: Private member storing the 1D quadrature weights.
*   **`Array2D<double> psi1d_`**: Private member storing the 1D basis functions evaluated at `quadPoints_`. `psi1d_(iqp, iord)` is the value of order `iord` polynomial at quadrature point `iqp`.
*   **`Array1D<double> psi1dSq_`**: Private member storing the squared norms of the 1D basis functions, computed via quadrature.
*   **`Array1D<double> psi1dSqExact_`**: Private member storing the analytically exact squared norms of the 1D basis functions.
*   **`int maxord_`**: Private member storing the maximum order of polynomials considered for this basis instance.
*   **`double alpha_`, `double beta_`**: Private members storing the parameters for certain polynomial families (e.g., Jacobi, Laguerre).
*   **`dsfmt_t rnstate_`**: Private member storing the state of the Mersenne Twister random number generator.
*   **`int rSeed_`**: Private member storing the seed used for the random number generator.

## Usage Examples

```cpp
#include "PCBasis.h"
#include "Array1D.h"
#include "Array2D.h"
#include <iostream>
#include <string>

int main() {
    // Example: Legendre-Uniform basis up to order 5
    std::string basisType = "LU";
    double alpha = 0.0; // Not used for LU
    double beta = 0.0;  // Not used for LU
    int maxOrder = 5;

    PCBasis legendreBasis(basisType, alpha, beta, maxOrder);

    std::cout << "Basis Type: " << legendreBasis.GetPCType() << std::endl;

    // Get Quadrature Rule
    Array2D<double> qPoints;
    Array1D<double> qWeights;
    Array2D<int> qIndices; // May not be populated for all rule types
    legendreBasis.GetQuadRule(qPoints, qWeights, qIndices);
    std::cout << "Number of quadrature points: " << qPoints.XSize() << std::endl;
    // for (int i = 0; i < qPoints.XSize(); ++i) {
    //     std::cout << "Point " << i << ": " << qPoints(i,0) << ", Weight: " << qWeights(i) << std::endl;
    // }

    // Evaluate basis functions at a custom point
    Array1D<double> custom_point(1);
    custom_point(0) = 0.5;
    Array1D<double> basis_evals_at_custom_point(maxOrder + 1);
    legendreBasis.EvalBasis(custom_point(0), basis_evals_at_custom_point);
    // std::cout << "\nBasis functions evaluated at " << custom_point(0) << ":" << std::endl;
    // for (int i = 0; i <= maxOrder; ++i) {
    //     std::cout << "P_" << i << "(" << custom_point(0) << ") = " << basis_evals_at_custom_point(i) << std::endl;
    // }

    // Get exact squared norms
    Array1D<double> norms_sq_exact;
    legendreBasis.Get1dNormsSqExact(norms_sq_exact); // Will compute if not already
    // std::cout << "\nExact squared norms:" << std::endl;
    // for (int i = 0; i <= maxOrder; ++i) {
    //     std::cout << "||P_" << i << "||^2 = " << norms_sq_exact(i) << std::endl;
    // }
    
    // Generate random samples
    Array1D<double> samples(10); // Generate 10 samples
    legendreBasis.GetRandSample(samples);
    // std::cout << "\nRandom samples from Uniform(-1,1):" << std::endl;
    // for (int i = 0; i < samples.XSize(); ++i) {
    //     std::cout << samples(i) << " ";
    // }
    // std::cout << std::endl;

    return 0;
}
```

## Dependencies and Interactions

*   **`Array1D.h`, `Array2D.h`**: Uses UQTk's `Array1D` and `Array2D` classes for storing quadrature points, weights, basis evaluations, norms, and samples.
*   **`ftndefs.h`**: Likely contains Fortran interoperability definitions, though not directly evident in the `PCBasis` public interface itself.
*   **`dsfmt_add.h`**: Header for the double-precision SIMD-oriented Fast Mersenne Twister (dSFMT) random number generator, used in `GetRandSample` and `SeedRandNumGen`.
*   **`error_handlers.h`**: For UQTk's exception handling (e.g., `Tantrum`).
*   **`uqtkconfig.h`**: Potentially for UQTk library configuration options.
*   **`quad.h`**: Contains the `Quad` class, which is used internally by `Init1dQuadPoints` to generate quadrature rules.
*   **`arrayio.h`**: Used internally, for instance, by `EvalBasis` when `type_` is "pdf" to read `ab.dat`.
*   **`pcmaps.h`**: Used internally by `GetRandSample` for certain basis types (e.g., "LG", "JB") to map samples from a uniform distribution to the target distribution.
*   **`combin.h`**: Not directly evident in the provided code snippets for `PCBasis.cpp` but might be used by other parts of the PCE module.
*   **Standard C++ Libraries**: `<iostream>`, `<string.h>` (or `<cstring>`), `<math.h>` (or `<cmath>`), `<stdio.h>`, `<sstream>`.
```
