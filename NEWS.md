# aforoR 0.2.0

## Minor Changes & Documentation Enhancements
* **Responsive Vignettes**: Updated all vignettes to use dynamic, responsive image formatting. All images now automatically scale to 100% of the page width with preserved aspect ratio for a premium viewing experience across all devices.
* **Vector Graphics Upgrade**: Replaced low-resolution raster images with high-definition vector graphics (`VIS_plot.pdf` and `ROS-ROL.pdf`) in the main vignettes.
* **Navigation Restructuring**: Renamed and moved the tutorial **"Basic Morphometry"** to **"Applications of Otolith Size"** under the **AFORO Experimental** section of the package website, separating experimental and basic workflows.
* **Asset Optimization**: Conducted a thorough audit of vignette assets, deleting 21 unlinked/orphaned images in `vignettes/images/` to drastically reduce package size and optimize standard building workflows.

## Bug Fixes & Refactoring
* **CRAN Namespace Compliance**: Qualified all calls to `var()` as `stats::var()` inside `R/morphometrics.R`, resolving R CMD static checking notes regarding undefined global functions.
