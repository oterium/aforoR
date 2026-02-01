<div align="center">
  <div align="right" style="float: right; margin-left: 20px;">
    <img src="man/figures/logo.png" width="150">
  </div>
  <div align="left">
    <h1>aforoR</h1>
    <p><strong>Shape Analysis of Fish Otoliths</strong><br>
    <em>Anàlisi de FORmes d'Otòlits in R</em></p>
  </div>
</div>

<br clear="right"/>

[![Version](https://img.shields.io/badge/version-0.1.0-6b94ff.svg?style=flat-square)](https://github.com/oterium/aforoR)
[![License: MIT](https://img.shields.io/badge/license-MIT-6b94ff.svg?style=flat-square)](https://github.com/oterium/aforoR/blob/main/LICENSE)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-6b94ff.svg?style=flat-square)](https://cran.r-project.org/)
[![Documentation](https://img.shields.io/badge/docs-vignettes-6b94ff.svg?style=flat-square)](https://oterium.github.io/aforoR)

---

### Overview

**aforoR** is an R package designed for high-precision **Shape Analysis of Fish Otoliths**. It implements the methodology established by the *AFORO* team (Parisi-Baradad et al., 2005), providing a complete workflow from image preprocessing to advanced morphometric descriptors.

While specifically optimized for otoliths, the package is equally powerful for analyzing other closed biological structures such as leaves, statoliths, beaks, and shells.

### Key Features

*   ✨ **Automatic Scale Detection**: Automatically identifies 1mm scale bars in images to convert pixel measurements to physical units (mm).
*   📐 **Advanced Morphometrics**: Comprehensive calculation of geometric indices:
    *   *Basic*: Area, Perimeter, Length, Width.
    *   *Indices*: Roundness, Form Factor, Circularity, Rectangularity, Ellipticity, and Aspect Ratio.
*   🌊 **Multi-scale Wavelet Analysis**: Decomposition of contours using discrete wavelet transforms for detailed frequency-based shape characterization.
*   🎯 **Elliptic Fourier Descriptors (EFDs)**: Standard Fourier-based shape analysis.
*   📈 **Curvature Scale Space (CSS)**: Geometric approach for detecting inflection points and lobes at multiple smoothing scales.
*   📍 **Optimal Point Selection**: Implementation of the Hall & Bathia (2012) algorithm for identifying the most discriminative points in functional data.
*   📦 **Batch Processing**: Efficiently process entire directories of images with automated output generation.

---

### Installation

```r
# Install development version from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("oterium/aforoR")

# Load the package
library(aforoR)
```

### Quick Start

The easiest way to use **aforoR** is the `process_images` function, which automates the entire analysis pipeline for a folder containing `.jpg` images.

```r
library(aforoR)

# 1. Define your working directory with images
my_images <- "path/to/your/images"

# 2. Run the automated pipeline
# This will perform: Preprocessing -> Scale Detection -> Contour -> Wavelets -> Morphometrics
process_images(
  folder = my_images,
  detect_scale = TRUE,   # Automatically find 1mm scale bars
  wavelets = TRUE,       # Compute wavelet coefficients
  ef = TRUE,             # Compute Elliptic Fourier Descriptors
  save = TRUE,           # Export results as CSV files
  testing = TRUE         # Generate diagnostic visualizations
)
```

---

### Learning Resources (Vignettes)

The package includes detailed tutorials (vignettes) to guide you through different analytical scenarios:

| Tutorial | Description |
| :--- | :--- |
| [**Obtaining Contour**](https://oterium.github.io/aforoR/articles/Obtaining_Contour.html) | Introduction to the core image processing and contour extraction workflow. |
| [**Otolith Morphometry**](https://oterium.github.io/aforoR/articles/Otolith_Morphometry.html) | Guide to calculating shape indices and modeling allometric relationships (fish-otolith). |
| [**Working Wavelets**](https://oterium.github.io/aforoR/articles/Working_Wavelets.html) | Deep dive into wavelet-based population comparison, PCA, and classification models. |
| [**Analyzing Lobes**](https://oterium.github.io/aforoR/articles/Analyzing_Lobes.html) | Specialized analysis of specific morphological features (lobes). |
| [**Shape Descriptors**](https://oterium.github.io/aforoR/articles/Shape_Descriptors.html) | Comparison between EFD, Wavelets, and CSS for otolith analysis. |
| [**Point Selection**](https://oterium.github.io/aforoR/articles/Point_Selection.html) | Sequential point selection for classification using the Hall & Bathia algorithm. |

---

### Citation

If you use **aforoR** in your research, please cite it as:

Otero-Ferrer, J.L., Tuset, V.M., Manjabacas, A., Lombarte, A. (2025). *aforoR*: Anàlisi de FORmes d'Otòlits / Shape Analysis of Fish Otoliths. R package version 0.1.0. [https://github.com/oterium/aforoR](https://github.com/oterium/aforoR)

### Acknowledgments

Developed by **JLOF, VMT, AM, and AL**. Based on the [AFORO](http://aforo.cmima.csic.es/) methodology for marine science and morphological research.

<div align="center">
  <br>
  <p>🐟 <b>Happy Otolith Analysis!</b> 🐟</p>
  <sub>Made with ❤️ for the fisheries science community</sub>
  <br><br>
  [![GitHub stars](https://img.shields.io/github/stars/oterium/aforoR?style=social)](https://github.com/oterium/aforoR)
</div>
