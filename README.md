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
[![Roxygen](https://img.shields.io/badge/roxygen-7.3.2-6b94ff.svg?style=flat-square)](https://github.com/oterium/aforoR)

---

### Overview

Developed for **Shape Analysis of Fish Otoliths** following the methodology established by the *AFORO* team (Parisi-Baradad et al., 2005). The package provides tools to extract coordinates, elliptic Fourier descriptors, and wavelets at multiple scales.

While inspired by the [AFORO website](http://aforo.cmima.csic.es/), this package allows for powerful batch processing and advanced morphometric analysis. It is also applicable to other closed biological structures like leaves, statoliths, and beaks.

#### Key Applications
*   **Species Identification**: Taxonomy, diversity, and evolution.
*   **Stock Differentiation**: Critical for fisheries management.
*   **Ontogenetic Modeling**: Detecting ecomorphological shifts.
*   **Phenotype Detection**: Advanced shape-based classification.

---

### Integration & Features (v0.1.0)

*   ✨ **Automatic 1mm Scale Detection**: Automatically detects the scale bar in images to convert pixels to actual physical units (mm).
*   📐 **Advanced Morphometrics**: Full calculation of Area, Perimeter, and Shape Indices (Roundness, Form Factor, Circularity, etc.) based on Tuset et al.
*   🌊 **Wavelet Decomposition**: Multi-scale analysis of otolith contours.

---

### Installation & Quick Start

```r
# Install the development version from GitHub:
install.packages("devtools")
devtools::install_github("oterium/aforoR")

# Load the package
library(aforoR)
```

#### Analytical Workflow
1.  **Image Preprocessing**: Cleaning and binarizer.
2.  **Contour Extraction**: High-precision boundary detection.
3.  **Distance Measurements**: Radial and perimeter metrics.
4.  **Morphometric Analysis**: Automated physical measurements.
5.  **Wavelet Decomposition**: Frequency-based shape descriptors.
6.  **Visualization & Export**: Publication-ready results.

---

### Citation

Otero-Ferrer, J.L., Tuset, V.M., Manjabacas, A., Lombarte, A. (2025). *AFOROR*: Anàlisi de FORmes d'Otòlits/Shape Analysis of Fish Otoliths. R package version 0.1.0. [GitHub Repository](https://github.com/oterium/aforoR).

### Contributions

Developed and optimized by **JLOF, VMT, AM, and AL**. Leading the way in digital otolith analysis and tutorial development.

<div align="center">
  <br>
  <p>🐟 <b>Happy Otolith Analysis!</b> 🐟</p>
  <sub>Made with ❤️ for the fisheries science community</sub>
  <br><br>
  [![GitHub](https://img.shields.io/github/stars/oterium/aforoR?style=social)](https://github.com/oterium/aforoR)
</div>
