# aforoR 0.2.0

## 🚀 Major Features & Algorithms

* **Generalized Procrustes Analysis (GPA)**: Added native support for GPA alignment (`procrustes = TRUE` in `process_images()`). It mathematically registers and scales a collection of otolith contours to a global "consensus shape," minimizing non-biological shape variance before extracting descriptors.
* **Advanced Morphometric Metrics**: Upgraded [calculate_morphometrics()](file:///c:/aforoR/R/morphometrics.R) to extract cutting-edge geometric shape descriptors:
  * **Feret_Max** (Maximum Feret/Caliper Diameter): The maximum distance between any two contour points, serving as a highly robust measure of maximum otolith length.
  * **Feret_Min** (Minimum Feret/Caliper Diameter): Calculated via an angular sweep (360 degrees) of parallel tangents, providing a highly robust and standard measure of otolith width.
  * **PCA_Angle**: Computes the main orientation angle of the major axis of inertia (in degrees, scale-invariant, and normalized to $[0, 180)$) using Singular Value Decomposition (SVD).
* **Optimal Point Selection (Hall & Bathia)**: Fully integrated the advanced **Hall & Bathia sequential point selection** machine learning algorithm for optimized shape feature classification and contour point selection.
* **Curvature Descriptors**: Integrated **Centroid Distance Fourier (CDF)** and **Curvature Scale Space (CSS)** shape descriptors alongside standard EFD and Wavelet methods.
* **Packaged Aphanopus Datasets**: Bundled real-world biological datasets for *Aphanopus carbo* and *Aphanopus intermedius* (`Aphanopus`, `Aphanopus_EF`, `Aphanopus_W5` datasets) to support reproducible examples and scientific workflows.

## 🎨 Design, Styling & Vignette Upgrades

* **Responsive Vignettes**: Injected responsive, high-performance styling blocks across all vignettes. Images (`<img>` tags) now occupy **100% of the page width** with correct aspect ratio (`height: auto`) and clean centring, offering a premium reading experience on the documentation website and offline HTML viewers.
* **High-Definition Graphics**: Replaced legacy, low-resolution raster plots with crisp, publication-grade high-resolution PNG graphics (`VIS_plot.png` and `ROS-ROL.png`) in the main vignettes, ensuring native HTML compatibility and responsiveness.
* **Navigation Restructuring**: Renamed the tutorial **"Basic Morphometry"** to **"Applications of Otolith Size"** and moved it to the **AFORO Experimental** section of the package documentation (`_pkgdown.yml`) to segregate standard and advanced biological size workflows.
* **Asset Cleanup & Weight Reduction**: Audited the assets directory and deleted **21 unlinked/orphaned images** from `vignettes/images/`, dramatically optimizing the package bundle weight and accelerating CRAN building workflows.

## 🛠️ Code Refactoring, Bug Fixes & Compliance

* **CRAN Namespace Compliance**: Qualified all calls to `var()` as `stats::var()` inside [R/morphometrics.R](file:///c:/aforoR/R/morphometrics.R#L83), resolving strict R CMD check undefined global function notes.
* **Translation & Professionalization**: Translated the entire package interface, console logs, error messages, helper scripts, and 9 vignettes into publication-grade English.
* **Robust Code Quality**:
  * Refactored the monolithic `process_images` function into highly modular, testable, and reusable functions.
  * Reorganized and cleaned up output folder generations, splitting outputs into dedicated `Polar/` and `Cartesian/` directories.
* **Continuous Integration (CI/CD)**: Established a robust automated GitHub Actions pipeline with `R-CMD-check` workflows across operating systems and automated `pkgdown` website building and deployment to GitHub Pages.
