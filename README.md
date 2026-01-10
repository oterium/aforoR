# aforoR: Shape Analysis of Fish Otoliths
### *Anàlisi de FORmes d'Otòlits*

[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/oterium/aforoR)
[![License: MIT](https://img.shields.io/badge/license-MIT%20%2B%20file%20LICENSE-blue.svg)](https://github.com/oterium/aforoR/blob/main/LICENSE)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://cran.r-project.org/)
[![Roxygen: 7.3.2](https://img.shields.io/badge/roxygen-7.3.2-blue.svg)](https://github.com/oterium/aforoR)

**aforoR** is a specialized R package for automated morphometric analysis of fish otoliths (ear stones). Otoliths are calcium carbonate structures in the inner ear of fish that are crucial for taxonomic identification, population studies, and fisheries research. This package provides a comprehensive toolkit for extracting, analyzing, and comparing otolith shape descriptors using advanced image processing and statistical methods.

## 🐟 What are Otoliths?

Otoliths are biomineralized structures found in the inner ear of bony fish that function in hearing and balance. Their species-specific shapes make them invaluable for:

- **Species identification** and taxonomic studies
- **Population dynamics** research in fisheries science  
- **Paleontological** reconstructions
- **Stock assessment** and management
- **Age determination** and growth studies

## ✨ Key Features

- **Automated image preprocessing** with adaptive thresholding
- **Contour detection and extraction** from otolith images
- **Polar and perimeter distance measurements** 
- **Wavelet analysis** for shape characterization
- **Elliptic Fourier descriptors** for shape quantification
- **Zonal analysis** with customizable regions
- **High-quality visualizations** and diagnostic plots
- **Batch processing** of image collections
- **CSV export** for statistical analysis

## 📋 System Requirements

### R Environment
- **R version**: ≥ 4.0.0
- **Operating System**: Windows, macOS, or Linux

### Required Dependencies
The package depends on the following R packages:
- `EBImage` - For image processing and analysis
- `Momocs` - For morphometric analysis and shape descriptors  
- `utils` - For utility functions (base R)

### Suggested Dependencies
- `devtools` - For development and installation from GitHub

## 🚀 Installation

### From GitHub (Development Version)

The package is currently available as a development version from GitHub. You'll need to install the required dependencies first:

```r
# Install required packages from Bioconductor and CRAN
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install EBImage from Bioconductor
BiocManager::install("EBImage")

# Install Momocs from CRAN
install.packages("Momocs")

# Install devtools if not already installed
install.packages("devtools")

# Install aforoR from GitHub
devtools::install_github("oterium/aforoR")
```

### Troubleshooting Installation

If you encounter issues installing `EBImage`:
```r
# Alternative Bioconductor installation
BiocManager::install("EBImage", force = TRUE)
```

For `Momocs` installation issues on Linux, you may need system dependencies:
```bash
# Ubuntu/Debian
sudo apt-get install libfftw3-dev

# CentOS/RHEL
sudo yum install fftw-devel
```
## 💻 Usage

### Quick Start

```r
# Load the package
library(aforoR)

# Set your working directory to the folder containing otolith images
setwd("path/to/your/otolith/images")

# Process all JPG images in the current folder
process_images(
  folder = ".",
  wavelets = TRUE,
  ef = TRUE,
  testing = TRUE,
  save = TRUE
)
```

### Detailed Workflow

The package follows a systematic workflow for otolith shape analysis:

#### 1. Image Preprocessing
```r
# Preprocess a single image
processed_img <- preprocess_image("otolith_sample.jpg", threshold = NULL)

# The function applies:
# - Grayscale conversion
# - Gaussian filtering (31px brush, σ=9)  
# - Automatic thresholding (Otsu method) or manual threshold
```

#### 2. Contour Extraction
```r
# Extract the main otolith contour
contour <- extract_contour(processed_img)

# Filters contours by:
# - Minimum 1000 points
# - Minimum area of 10^5 pixels
```

#### 3. Distance Measurements
```r
# Calculate polar and perimeter distances
distances <- calculate_distances(contour, n_points = 512)

# Returns:
# - polar$radii: distances from centroid to boundary
# - perimeter$dist: distances along perimeter
# - normalized versions of both measurements
```

#### 4. Wavelet Analysis
```r
# Perform wavelet decomposition
wavelets <- calculate_wavelets_analysis(distances, n_scales = 9)

# Generates 9 scales of wavelet coefficients for:
# - Polar distance data
# - Perimeter distance data
```

#### 5. Visualization and Export
```r
# Save diagnostic visualizations
save_visualization(processed_img, distances, wavelets, 
                  "sample_image", "./output/")

# Export results to CSV files  
save_analysis_results(results_list, "./output/", "polar")
```

### Advanced Usage

#### Batch Processing with Custom Parameters
```r
# Process images in subfolders with custom settings
process_images(
  folder = "otolith_dataset",
  subfolder = TRUE,           # Process subfolders
  threshold = 0.5,            # Manual threshold
  wavelets = TRUE,            # Include wavelet analysis
  ef = TRUE,                  # Include Elliptic Fourier descriptors
  testing = TRUE,             # Save diagnostic images
  pseudolandmarks = "both",   # Use both polar and perimeter
  save = TRUE                 # Export CSV files
)
```

#### Working with Individual Functions
```r
# Step-by-step analysis of a single image
image_path <- "my_otolith.jpg"

# 1. Preprocess
binary_img <- preprocess_image(image_path)

# 2. Extract contour
contour <- extract_contour(binary_img)

# 3. Calculate distances  
distances <- calculate_distances(contour, n_points = 512)

# 4. Wavelet analysis
wavelets <- calculate_wavelets_analysis(distances, n_scales = 9)

# 5. Create structured results
results <- create_result_structure(distances$polar, wavelets$polar, 
                                  basename(image_path), result_type = "polar")
```

### Output Structure

The package generates several types of output:

#### CSV Files
- `DistanciaEN.csv` - Raw distance measurements
- `Distancia_NormEN.csv` - Normalized distance measurements  
- `Wavelet_1EN.csv` to `Wavelet_9EN.csv` - Wavelet coefficients by scale
- `EllipticCoeEN.csv` - Elliptic Fourier coefficients
- `Coords.csv` - Coordinate data

#### Visualization Files
- `zonal [filename].jpg` - Zonal analysis plots showing colored regions
- `wavelet [filename].jpg` - Wavelet visualization plots

#### Analysis Directories
- `Procesamiento/` - Polar distance analysis results
- `Procesamiento2/` - Perimeter distance analysis results

## 🔧 Function Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `process_images()` | Main function for batch processing of otolith images |
| `preprocess_image()` | Image preprocessing with filtering and binarization |
| `extract_contour()` | Contour detection and extraction from binary images |
| `calculate_distances()` | Compute polar and perimeter distance measurements |
| `calculate_wavelets_analysis()` | Wavelet decomposition of distance data |

### Analysis Functions

| Function | Description |
|----------|-------------|
| `regularradius()` | Calculate distances from centroid to boundary points |
| `dper()` | Compute perimeter-based distance measurements |
| `ild()` | Calculate Euclidean distance between two points |
| `fwaveletspl_3()` | Specialized wavelet transform function |
| `fatrous1d()` | Additional mathematical analysis function |

### Output Functions

| Function | Description |
|----------|-------------|
| `save_visualization()` | Generate and save diagnostic visualizations |
| `save_analysis_results()` | Export analysis results to CSV files |
| `create_result_structure()` | Create standardized result objects |

## 📊 Data Requirements

### Input Images
- **Format**: JPEG (.jpg) files
- **Content**: High-contrast otolith images with clear boundaries
- **Background**: Preferably uniform, contrasting with otolith
- **Resolution**: Higher resolution images provide better results
- **Quality**: Images should be well-focused and properly lit

### Recommended Image Preparation
1. **Standardized photography setup** with consistent lighting
2. **Clean otoliths** free from debris or tissue
3. **Contrasting background** (dark otolith on light background or vice versa)
4. **Consistent orientation** for comparative studies
5. **High resolution** (minimum 1000x1000 pixels recommended)

## 🚨 Troubleshooting

### Common Issues

**Error: "No contours found in image"**
- Check image quality and contrast
- Try adjusting the threshold parameter
- Ensure otolith is clearly visible against background

**Error: "No contour found meeting size criteria"**  
- Image resolution may be too low
- Otolith may be too small in the image
- Adjust size criteria in `extract_contour()` if needed

**Installation issues with EBImage**
- Ensure Bioconductor is properly installed
- Try installing from source: `BiocManager::install("EBImage", type="source")`

**Memory issues with large image sets**
- Process images in smaller batches
- Reduce image resolution if acceptable for your analysis
- Increase R memory limit: `memory.limit(size = 8000)` (Windows)

### Performance Tips
- Use SSD storage for faster I/O operations
- Process images in parallel using `parallel` package if needed
- Consider image resizing for very large datasets

## 🤝 Contributing

Contributions to aforoR are welcome! We encourage:

### Types of Contributions
- **Bug reports** and feature requests via GitHub Issues
- **Code contributions** through Pull Requests
- **Documentation improvements**
- **Example datasets** and use cases
- **Testing** with different otolith species and image types

### Development Workflow
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes following R package best practices
4. Add tests for new functionality
5. Update documentation using roxygen2
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Code Style
- Follow standard R package conventions
- Use descriptive function and variable names
- Include comprehensive documentation with examples
- Add input validation for all functions
- Write unit tests for new functionality

## 📝 Citation

If you use aforoR in your research, please cite the package as:

```bibtex
@software{otero_ferrer_2024_aforoR,
  author = {Otero-Ferrer, Jose Luis},
  title = {aforoR: Shape Analysis of Fish Otoliths},
  subtitle = {Anàlisi de FORmes d'Otòlits},
  year = {2024},
  version = {0.1.0},
  url = {https://github.com/oterium/aforoR},
  note = {R package version 0.1.0}
}
```

**APA Style:**
Otero-Ferrer, J. L. (2024). *aforoR: Shape Analysis of Fish Otoliths* (Version 0.1.0) [Computer software]. https://github.com/oterium/aforoR

## 📄 License

aforoR is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

This license allows for:
- ✅ Commercial use
- ✅ Modification  
- ✅ Distribution
- ✅ Private use

With the following conditions:
- 📋 License and copyright notice must be included
- ❌ No warranty or liability

## 👤 Contact

**Jose Luis Otero-Ferrer**
- 📧 Email: [joseluisotero@gmail.com](mailto:joseluisotero@gmail.com)
- 🔬 Alternative: [jootero@uvigo.es](mailto:jootero@uvigo.es)  
- 🐙 GitHub: [@oterium](https://github.com/oterium)

For questions about:
- **Technical issues**: Open a GitHub Issue
- **Research collaboration**: Contact via email
- **Feature requests**: GitHub Discussions or Issues

## 🙏 Acknowledgements

We extend our gratitude to:
- **Contributors** to the aforoR package development
- **R Core Team** for the R statistical computing environment
- **Bioconductor team** for the EBImage package
- **Momocs developers** for morphometric analysis tools
- **Scientific community** using otolith analysis in their research

Special thanks to researchers and institutions that have provided feedback, testing, and validation of the methods implemented in this package.

## 📚 Related Resources

### Scientific Background
- [Otolith morphology in fisheries science](https://example.com) (general reference)
- [Morphometric analysis methods](https://example.com) (methodological background)
- [Wavelet analysis in biological shape analysis](https://example.com)

### Related R Packages
- [`Momocs`](https://cran.r-project.org/package=Momocs) - Morphometric analysis
- [`EBImage`](https://bioconductor.org/packages/EBImage/) - Image processing
- [`shapes`](https://cran.r-project.org/package=shapes) - Statistical shape analysis

### Datasets and Examples
- Example otolith image datasets (to be added)
- Tutorial notebooks and workflows (to be added)

---

<div align="center">

**🐟 Happy Otolith Analysis! 🐟**

*Made with ❤️ for the fisheries science community*

[![GitHub](https://img.shields.io/github/stars/oterium/aforoR?style=social)](https://github.com/oterium/aforoR)

</div>
