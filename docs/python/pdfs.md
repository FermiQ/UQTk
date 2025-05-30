# pdfs.py

## Overview

This file, `pdfs.py`, is part of the PyUQTk plotting utilities within the UQ Toolkit (UQTk). It provides functions for visualizing one-dimensional and two-dimensional probability density functions (PDFs) from sample data. The visualizations can be rendered as kernel density estimates (KDE), histograms, or scatter plots of the samples themselves. The main (but incomplete in the provided snippet) function `plot_pdfs` seems intended to create a matrix of PDF plots, often used in MCMC analysis (a "triangle plot" or "corner plot").

## Key Components

*   **`plot_pdf1d(sams, pltype='hist', color='b', lw=1, nom_height_factor=10., ax=None)`**: Plots a 1D PDF from samples.
    *   `sams`: 1D NumPy array of samples.
    *   `pltype`: Type of plot. Options are:
        *   `'kde'`: Kernel Density Estimate.
        *   `'hist'`: Histogram (default).
        *   `'sam'`: Scatter plot of the samples (plotted along the x-axis at y=0).
        *   `'nom'`: Vertical lines at sample locations, typically used to show nominal parameter values or discrete samples.
    *   `color`: Color of the plot.
    *   `lw`: Line width.
    *   `nom_height_factor`: Factor to control the height of lines when `pltype='nom'`.
    *   `ax`: Matplotlib axes object to plot on. If `None`, uses the current axes (`plt.gca()`).
*   **`plot_pdf2d(samsx, samsy, pltype='kde', ncont=10, color=None, lwidth=1, mstyle='o', ax=None)`**: Plots a 2D PDF from samples.
    *   `samsx`: 1D NumPy array of samples for the x-axis.
    *   `samsy`: 1D NumPy array of samples for the y-axis.
    *   `pltype`: Type of plot. Options are:
        *   `'kde'`: Kernel Density Estimate, shown as contour lines (default).
        *   `'sam'`: Scatter plot of the (x, y) samples.
    *   `ncont`: Number of contour levels for KDE plot.
    *   `color`: Color of the plot. For KDE, can be a single color string or a list/tuple for multiple contour colors.
    *   `lwidth`: Line width for KDE contours.
    *   `mstyle`: Marker style for scatter plot (`pltype='sam'`).
    *   `ax`: Matplotlib axes object to plot on. If `None`, uses the current axes (`plt.gca()`).
*   **`plot_pdfs(ind_show=[], samples_file='pchain.dat', plot_type='tri', names_file=None, burnin=0, every=1, nominal_file=None, nsam_xi=1, prangeshow_file=None)`**:
    *   *(Incomplete in the provided script)* This function is intended to generate a matrix of 1D and 2D PDF plots. Typically, diagonal plots show 1D marginal PDFs, and off-diagonal plots show 2D joint PDFs.
    *   `ind_show`: Indices of parameters/variables to show.
    *   `samples_file`: Path to the file containing sample data (e.g., MCMC chain).
    *   `plot_type`: Type of plot matrix, e.g., `'tri'` for a triangle plot (lower or upper).
    *   `names_file`: Path to a file containing names for the parameters.
    *   `burnin`: Number of initial samples to discard (burn-in period).
    *   `every`: Subsampling rate (take every Nth sample).
    *   `nominal_file`: Path to a file with nominal parameter values to overlay.
    *   `nsam_xi`: (Purpose unclear from snippet, possibly related to plotting specific sample indices).
    *   `prangeshow_file`: Path to a file specifying plot ranges for parameters.

## Important Variables/Constants

This file does not define significant standalone variables or constants that affect general behavior outside of function arguments. The default arguments in the functions act as configurable constants for individual plot calls.

## Usage Examples

```python
import numpy as np
import matplotlib.pyplot as plt
# Assuming pdfs.py is in the Python path or PyUQTk is installed
# from PyUQTk.plotting import pdfs # Or appropriate import based on project structure

# Generate some sample data
samples_1d = np.random.normal(loc=0, scale=1, size=1000)
samples_2d_x = np.random.normal(loc=0, scale=1, size=1000)
samples_2d_y = 0.5 * samples_2d_x + np.random.normal(loc=0, scale=0.5, size=1000)

# Example for plot_pdf1d
plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
# pdfs.plot_pdf1d(samples_1d, pltype='hist', color='skyblue') # Replace with actual import
# plt.title("1D PDF (Histogram)")

plt.subplot(1, 2, 2)
# pdfs.plot_pdf1d(samples_1d, pltype='kde', color='salmon', lw=2) # Replace with actual import
# plt.title("1D PDF (KDE)")

# plt.tight_layout()
# plt.show()


# Example for plot_pdf2d
plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
# pdfs.plot_pdf2d(samples_2d_x, samples_2d_y, pltype='sam', color='green', mstyle='.') # Replace with actual import
# plt.title("2D PDF (Samples)")

plt.subplot(1, 2, 2)
# pdfs.plot_pdf2d(samples_2d_x, samples_2d_y, pltype='kde', ncont=15, color='purple', lwidth=1.5) # Replace with actual import
# plt.title("2D PDF (KDE Contours)")

# plt.tight_layout()
# plt.show()

# Conceptual usage for the (incomplete) plot_pdfs function
# This is speculative based on typical "triangle plot" functionality.
# figs, axarr = pdfs.plot_pdfs(
#     samples_file='path/to/mcmc_chain.dat',
#     plot_type='tri', # For a lower triangle plot
#     burnin=1000,
#     every=10,
#     names_file='path/to/parameter_names.txt'
# )
# if figs:
#     plt.show()
```
*Note: The usage examples for `plot_pdf1d` and `plot_pdf2d` would require uncommenting and ensuring `PyUQTk.plotting.pdfs` is correctly imported. The example for `plot_pdfs` is highly conceptual due to its incomplete state in the script.*

## Dependencies and Interactions

*   **`os`**: Used for operating system interactions (path manipulation, though not explicitly shown in the function bodies provided).
*   **`sys`**: Used for system-specific parameters and functions (e.g., `sys.exit()` in error cases).
*   **`PyUQTk.utils.pdf_kde`**: This module is imported and used by `plot_pdf1d` and `plot_pdf2d` when `pltype='kde'` to perform the kernel density estimation.
*   **`numpy` (as `np`)**: Essential for numerical operations, especially array manipulation for sample data, grid generation for KDE, etc.
*   **`matplotlib.pyplot` (as `plt`)**: The core library used for generating all plots. Functions interact with `plt.gca()` to get current axes or can accept an `ax` object.

**Interactions:**
*   The plotting functions (`plot_pdf1d`, `plot_pdf2d`) are designed to be potentially used together, possibly by a higher-level function like the incomplete `plot_pdfs`, to create comprehensive visualizations of multidimensional sample data.
*   They rely on `pdf_kde.get_pdf` for the KDE calculations.
```
