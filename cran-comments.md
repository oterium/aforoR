## CRAN Submission Comments

This is a new submission of 'aforoR' (version 0.3.0).

### Test environments
* Windows 10 x64, R 4.4.2 / 4.4.3 (local)
* Ubuntu Linux (via GitHub Actions: R-release, R-devel)
* macOS (via GitHub Actions: R-release)
* win-builder (R-devel, R-release)

### R CMD check results
There were 0 ERRORS and 0 WARNINGS.

### Notes
* `New submission`: This is the first submission of aforoR to CRAN.
* `installed size is 13.9Mb`: The package contains 9 comprehensive HTML documentation vignettes (doc/ 11.0Mb) and high-resolution biological otolith reference datasets (extdata/ 1.3Mb) demonstrating multi-scale wavelet and morphological otolith shape analysis.
* `Additional_repositories: https://bioconductor.org/packages/release/bioc`: Required to declare the repository source for the Bioconductor package `EBImage` dependency.
* `checking for future file timestamps`: Local build artifact on Windows machine during check.

### Method References
The methodology implemented in this package is described in the following peer-reviewed publications cited in the DESCRIPTION and package documentation:
* Parisi-Baradad et al. (2005) <doi:10.1071/MF04162>
* Lombarte et al. (2006) <doi:10.3989/scimar.2006.70n1147>
* Manjabacas et al. (2025) <doi:10.1016/j.fishres.2024.107242>
* Piera et al. (2005) <doi:10.1071/MF04163>

### Downstream dependencies
There are currently no downstream dependencies on CRAN.
