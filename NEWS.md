# aforoR 0.3.0

## 🚀 AFORO Harmonization & Enhanced Workflows

* **AFORO Canonical References**: Integrated foundational and latest AFORO direct references across package documentation, Overview, and bibliography (Lombarte et al., 2006 and Manjabacas et al., 2025).
* **Clockwise Contour Decomposition**: Standardized contour extraction and multi-scale decomposition to follow the AFORO clockwise convention (Rostrum -> Ventral -> Posterior -> Dorsal -> Rostrum) via `.coo_ensure_clockwise()` in `extract_contour()` and `.resolve_contour_and_start()`, providing 100% phase and coefficient parity with AFORO web and MATLAB databases.
* **Modernized Wavelet Architecture**: Refactored stationary wavelet convolution into `wavelet_atrous_1d()` and multi-scale spline decomposition into `wavelet_spline_transform()`, eliminating legacy function names and variable structures while guaranteeing exact zero-difference numerical identity.
* **Piera et al. (2005) Coordinate Reference**: Added citation to Piera et al. (2005) in coordinate system documentation (*Knowing Wavelets* tutorial) as the pioneering work defining Cartesian vs. Polar representations for otolith shapes.
* **Quick Start Basic Morphometrics**: Enhanced README Quick Start with immediate inspection of `head(results$morphometrics)` and a self-contained demonstration using bundled sample otolith imagery.
* **Modular Architecture Refactoring**: Split monolithic image processing routines into specialized, high-cohesion modules (`R/preprocess.R`, `R/scale_and_interactive.R`, `R/visualization.R`, `R/distances.R`, `R/wavelets.R`, and `R/shape_tools.R`), reducing code duplication and improving maintainability.
* **Interactive Manual Scale Bar Calibration**: Added `detect_scale = "manual"` via `select_scale_bar()`, allowing users to click physical scale bar endpoints with Euclidean distance calculation.
* **Interactive Contour Starting Point Selection**: Added `start_points = "manual"` via `select_starting_point()`, snapping interactive mouse clicks to the nearest contour vertex while supporting automatic fallback.
* **Adaptive Specimen Size Threshold**: Upgraded `min_area` to `"auto"` by default, dynamically calculating 0.5% of total image area to prevent false rejection of small specimens (e.g. *Ostium*) while retaining explicit numeric overrides.
* **Wavelet Signal Detail Correction**: Clarified and corrected `detail = TRUE` (default) to return the detail/wavelet coefficients signal ($w_j = a_{j-1} - a_j$, centered around 0) and `detail = FALSE` to return the approximation signal ($a_j$).
* **Harmonized Homologous Starting Point**: Standardized `start_point = "auto" | "manual" | c(x, y)` across `calculate_css()` and `calculate_distances()`, removing non-homologous `"none"`.
* **Updated Authorship & Copyright**: Designate Jose Luis Otero Ferrer as author, creator, and copyright holder (`aut`, `cre`, `cph`) and all co-authors as copyright holders (`cph`).
* **Comprehensive Test Suite & Documentation**: 598 unit and integration tests passing (`FAIL 0 | WARN 0`), fully refreshed documentation and vignettes.

# aforoR 0.2.1

## 🚀 Autonomous & CRAN-Ready Architecture

* **Pure-R Native Geometric Algorithms**: Implemented zero-dependency, native pure-R geometric functions in `R/shape_tools.R` (`.coo_area`, `.coo_perim`, `.coo_align`, `.coo_length`, `.coo_width`, `.coo_calliper`, and `.coo_interpolate`), completely replacing legacy external dependencies.
* **Native Normalized Elliptic Fourier Analysis (EFA)**: Implemented `.efourier()` following Kuhl & Giardina (1982), computing harmonic Fourier coefficients ($A_n, B_n, C_n, D_n$) with full invariant normalization for size, orientation, and starting phase.
* **Native Generalized Procrustes Analysis (GPA)**: Implemented `.procrustes_align()` for multi-contour iterative consensus shape registration using singular value decomposition (SVD).
* **Elimination of `Momocs` Dependency**: Removed `Momocs` and `Remotes:` from `DESCRIPTION`, achieving full CRAN policy compliance and eliminating build/check issues associated with unmaintained external packages.
* **Package Weight & Size Optimization**: Optimized vignette raster graphics down from 16.4 MB to 3.1 MB and purged unreferenced image copies from `man/figures/`, reducing the installed package size by 72% (from 46.0 MB to 12.9 MB) and clearing all CRAN help-size notes.
* **HPC & Parallel Safety**: Added automatic core count limits (`n_cores` clamped to $\le 2$ in testing/CRAN environments) in `select_points_hall()` and `process_images()`.

# aforoR 0.2.0

## 🚀 Major Features & Algorithms

* **Intelligent Scale Detection & Visual Overlay**: Enhanced `find_scale_bar_length()` with dynamic area scaling to reliably detect physical scale bars (e.g., 2mm bars) in high-resolution and 4K photographs. Added automated visual bounding boxes, technical metric annotations, and scale bar entries to diagnostic validation plots (`save_visualization` and `save_visualization_perimeter`) for quality control.
* **Intelligent Adaptive Tri-Class Threshold**: Implemented `calculate_adaptive_threshold()` utilizing a 3-class Otsu variance maximization to accurately delineate translucent edges and shaded margins in otolith photographs without clipping.
* **Generalized Procrustes Analysis (GPA)**: Added native support for GPA alignment (`procrustes = TRUE` in `process_images()`). It mathematically registers and scales a collection of otolith contours to a global "consensus shape," minimizing non-biological shape variance before extracting descriptors.
* **Advanced Morphometric Metrics**: Upgraded `calculate_morphometrics()` to extract cutting-edge geometric shape descriptors:
  * **Feret_Max** (Maximum Feret/Caliper Diameter): The maximum distance between any two contour points, serving as a highly robust measure of maximum otolith length.
  * **Feret_Min** (Minimum Feret/Caliper Diameter): Calculated via an angular sweep (360 degrees) of parallel tangents, providing a highly robust and standard measure of otolith width.
  * **PCA_Angle**: Computes the main orientation angle of the major axis of inertia (in degrees, scale-invariant, and normalized to $[0, 180)$) using Singular Value Decomposition (SVD).
* **Optimal Point Selection (Hall & Bathia)**: Fully integrated the advanced **Hall & Bathia sequential point selection** machine learning algorithm for optimized shape feature classification and contour point selection.
* **Curvature Descriptors**: Integrated **Curvature Scale Space (CSS)** shape descriptors alongside standard Wavelet and EFD methods.
* **Packaged Aphanopus Datasets**: Bundled real-world biological datasets for *Aphanopus carbo* and *Aphanopus intermedius* (`Aphanopus`, `Aphanopus_EF`, `Aphanopus_W5` datasets) to support reproducible examples and scientific workflows.

## 🎨 Design, Styling & Vignette Upgrades

* **Responsive Vignettes**: Injected responsive, high-performance styling blocks across all vignettes. Images (`<img>` tags) now occupy **100% of the page width** with correct aspect ratio (`height: auto`) and clean centring, offering a premium reading experience on the documentation website and offline HTML viewers.
* **High-Definition Graphics**: Replaced legacy, low-resolution raster plots with crisp, publication-grade high-resolution PNG graphics (`VIS_plot.png` and `ROS-ROL.png`) in the main vignettes, ensuring native HTML compatibility and responsiveness.
* **Navigation Restructuring**: Renamed the tutorial **"Basic Morphometry"** to **"Applications of Otolith Size"** and moved it to the **AFORO Experimental** section of the package documentation (`_pkgdown.yml`) to segregate standard and advanced biological size workflows.
* **Asset Cleanup & Weight Reduction**: Audited the assets directory and deleted **21 unlinked/orphaned images** from `vignettes/images/`, dramatically optimizing the package bundle weight and accelerating CRAN building workflows.

## 🛠️ Code Refactoring, Bug Fixes & Compliance

* **CRAN Namespace Compliance**: Qualified all calls to `var()` as `stats::var()` inside `R/morphometrics.R`, resolving strict R CMD check undefined global function notes.
* **Translation & Professionalization**: Translated the entire package interface, console logs, error messages, helper scripts, and 9 vignettes into publication-grade English.
