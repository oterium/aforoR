# aforoR Test Infrastructure & Architecture Guide

## 1. Overview
This document specifies the testing architecture, feature inventory, coverage thresholds, execution protocols, and safety guarantees for the `aforoR` R package.

The testing suite provides end-to-end (E2E), opaque-box, and requirements-driven verification across all functional requirements (R1–R3) and empirical multi-species biological datasets (R4).

---

## 2. Testing Hierarchy: 4-Tier Test Suite

The test infrastructure is structured into 4 distinct verification tiers:

```
+-------------------------------------------------------------------------------+
|                       aforoR 4-Tier Testing Architecture                      |
+-------------------------------------------------------------------------------+
|  TIER 1: Feature Coverage (tests/testthat/test-e2e_acceptance.R)              |
|  - Requirement R1: Scale bar span, endpoints, red mask, plaque, visualization |
|  - Requirement R2: Float clamp, /255 guard, Otsu breaks, 2-ch TIFF, border 85%|
|  - Requirement R3: Flexible plot.css() S3 method with user ... overrides     |
|  - Threshold: >= 5 independent test cases per feature (>= 60 tests total)     |
+-------------------------------------------------------------------------------+
|  TIER 2: Boundary & Corner Cases (tests/testthat/test-e2e_acceptance.R)       |
|  - Empty / zero inputs (NULL, 0x0, empty matrix)                              |
|  - Pure red chromatic mask vs subtle red / orange / pink                      |
|  - Plaque vs thin line contrast inversion and eccentricity                    |
|  - Vertical vs horizontal scale bar orientation axis                          |
|  - Float overshoot near 1.0 (1.0 + 1e-16 vs 1.05)                             |
|  - 2-channel TIFF alpha strip array dimension and mode reduction              |
|  - Area >= 85% image frame filtering exact boundary                           |
|  - Small image min_area=1000 boundary (999 px vs 1001 px)                     |
+-------------------------------------------------------------------------------+
|  TIER 3: Cross-Feature Interactions (tests/testthat/test-e2e_acceptance.R)    |
|  - 2-channel TIFF + scale bar: alpha strip + span detection + diagnostic plot |
|  - Textured background + plaque scale bar: adaptive thresholding + plaque bar |
|  - Small image + custom plot.css: end-to-end processing pipeline              |
+-------------------------------------------------------------------------------+
|  TIER 4: Real-World Multi-Species Benchmark (tests/test_e2e_runner.R)        |
|  - 38 verified raw biological images across all 19 collections                |
|  - Read-Only Drive D: Safety Protocol: staged strictly in tempdir()           |
|  - Scale bar error <= 2 px (ratio 1.00 +- 0.02) on 13 scale collections       |
|  - Contour area < 85% total image area on 100% of tested samples              |
|  - Zero channel-tinting on 2-channel TIFFs                                    |
|  - Exports machine-readable CSV and Markdown summary reports                  |
+-------------------------------------------------------------------------------+
```

---

## 3. Feature Inventory & Coverage Mapping

| Feature # | Feature Name | Requirement | Target File | Interface Contract | Min Tests |
|---|---|---|---|---|---|
| 1 | Physical Scale Bar Span | R1 | `R/process_images.R` | `pixels_per_mm = bar_length_px / scale_mm`; true span $\max(\Delta x, \Delta y)$ | 5 |
| 2 | Scale Bar Endpoint Calculation | R1 | `R/process_images.R` | `attr(val, "scale_info")$endpoints` contains integer `x1, y1, x2, y2` | 5 |
| 3 | Red Scale Bar Detection | R1 | `R/process_images.R` | Chromatic red mask $(R > 0.45 \land R > 1.4\max(G, B))$ | 5 |
| 4 | Plaque & Vertical Scale Bars | R1 | `R/process_images.R` | Detect vertical span ($\Delta y > \Delta x$) and inner bar in plaques | 5 |
| 5 | Scale Bar Endpoint Visualization | R1 | `R/process_images.R` | High-contrast markers and span line rendered in diagnostic plots | 5 |
| 6 | FFT Float Clamping | R2 | `R/process_images.R` | `foto2` clamped strictly to $[0, 1]$ post `EBImage::filter2()` | 5 |
| 7 | Threshold Rescaling Guard | R2 | `R/process_images.R` | Check `max(img) > 1.05` before `/ 255` in `calculate_adaptive_threshold()` | 5 |
| 8 | Otsu Range Safeguard | R2 | `R/process_images.R` | Image values bounded strictly to prevent `breaks do not span range` error | 5 |
| 9 | 2-Channel TIFF Handling | R2 | `R/process_images.R` | Strip alpha channel (`img[,,1]`) and set `colorMode <- Grayscale` | 5 |
| 10 | Outer Frame Filtering | R2 | `R/process_images.R` | Contours with area $\ge 85\%$ of total image area discarded | 5 |
| 11 | Small Image `min_area` Support | R2 | `R/process_images.R` | Support `min_area = 1000` for 480x480 `Ostium` specimen contours | 5 |
| 12 | Flexible `plot.css()` S3 Method | R3 | `R/curvature.R` | Merge user `...` parameters using `utils::modifyList()` | 5 |
| 13 | Multi-Species Benchmark (38 samples) | R4 | `tests/test_e2e_runner.R` | 38 images from 19 collections tested from temporary staging | 38 |

---

## 4. Safety Protocol: Read-Only Biological Datasets

Source images reside in `D:/Back-up_Uvigo/mecanico/Dropbox 2017/dropbox antiguo`.
To guarantee preservation of the original scientific collections:
1. All test scripts and runners treat `D:/` strictly as **read-only**.
2. No temporary files, modified images, or derivative folders are written to `D:/`.
3. The benchmark runner (`tests/test_e2e_runner.R`) copies raw test images into a scratch directory under `tempdir()` before processing.
4. Temporary staging directories are completely unlinked on exit via `on.exit(unlink(staging_root, recursive = TRUE))`.

---

## 5. Execution Commands

### Run Unit and Acceptance Tests (Tiers 1, 2, 3)
```powershell
& "C:\Program Files\R\R-4.2.1\bin\Rscript.exe" -e "devtools::test(filter = 'e2e_acceptance')"
```

### Run Full Test Suite Across All Packages
```powershell
& "C:\Program Files\R\R-4.2.1\bin\Rscript.exe" -e "devtools::test()"
```

### Run Tier 4 Real-World 19-Collection Multi-Species Benchmark
```powershell
& "C:\Program Files\R\R-4.2.1\bin\Rscript.exe" tests/test_e2e_runner.R
```

### Run R CMD check
```powershell
& "C:\Program Files\R\R-4.2.1\bin\Rcmd.exe" check --no-manual .
```
