# Hyperspectral Imaging (HSI) Automation Tools in R

<p align="justify">
This repository contains R scripts designed to streamline preprocessing of Near-Infrared Hyperspectral Imaging (HSI) data. These tools automate extraction, compilation, and cleaning of spectral information, enabling efficient workflows for chemometrics, food authenticity research, and machine-learning applications using ENVI-formatted hyperspectral datasets.
</p>

The repository includes:

- **roi_extraction_envi_hsi.R** – Extracts elliptical Regions of Interest (ROI) from ENVI `.hdr/.dat` hyperspectral cubes and computes average spectra.  
- **compile_ground_coffee_HSI_txt_to_csv.R** – Compiles instrument-generated `.txt` spectral files into a structured wavelength × sample matrix.

---

## Installation

Ensure you have R (≥ 4.0) installed.

### Install Required Packages

```r
install.packages(c("terra", "ggplot2"))
```

---

## 1. ROI Extraction from ENVI Hyperspectral Cubes

**Script:** `roi_extraction_envi_hsi.R`

<p align="justify">
This script is designed for ENVI `.hdr/.dat` hyperspectral cubes. It allows you to define an elliptical ROI, extract mean reflectance spectra, visualize ROI placement, and batch-process entire folders or nested directory structures.
</p>

### Features

- Reads ENVI hyperspectral files (`.hdr` + `.dat`)  
- Defines elliptical ROI using center coordinates + axis lengths  
- Extracts mean spectra across ROI pixels  
- Batch directory processing  
- Recursive subdirectory processing  
- ROI mask visualization  
- Spectrum plotting  
- Saves extracted spectra to `.csv`

### Example Usage

#### Single File Extraction

```r
source("roi_extraction_envi_hsi.R")

result <- process_hyperspectral_image(
  file_path     = "path/to/image.dat",
  center_x      = 100,
  center_y      = 150,
  x_axis_length = 15,
  y_axis_length = 15
)

visualize_roi(result)
plot_spectrum(result$spectrum)

write.csv(result$spectrum, "sample_spectrum.csv", row.names = FALSE)
```

#### Batch Processing (Non-Recursive)

```r
all_spectra <- batch_process_directory(
  directory     = "path/to/folder",
  x_axis_length = 8,
  y_axis_length = 4
)

write.csv(all_spectra, "batch_extracted_spectra.csv", row.names = FALSE)
```

#### Batch Processing (Recursive)

```r
all_spectra <- batch_process_directory_recursive(
  root_directory = "path/to/root",
  x_axis_length  = 8,
  y_axis_length  = 4
)

write.csv(all_spectra, "recursive_extracted_spectra.csv", row.names = FALSE)
```

---

## 2. Compile HSI-Derived TXT Files Into a Spectral Matrix

**Script:** `compile_ground_coffee_HSI_txt_to_csv.R`

<p align="justify">
This script compiles multiple `.txt` spectral output files (exported from HSI software) into a clean wavelength × sample matrix. It automatically removes metadata rows, extracts the "mean" reflectance column, and aligns all samples against a common wavelength axis.
</p>

### Expected TXT Format

Each `.txt` file should contain:

1. Wavelength  
2. Min  
3. Mean − SD  
4. **Mean reflectance (used)**  
5. Mean + SD  
6. Max  
7–9. Extra columns (ignored)

### Example Usage

```r
setwd("path/to/txt/files")
source("compile_ground_coffee_HSI_txt_to_csv.R")
```

Output generated:

- `ground_coffee_compiled.csv` – Matrix of samples × wavelengths

---

## Repository Structure

```
HSI-automation/
├── README.md
├── roi_extraction_envi_hsi.R
├── compile_ground_coffee_HSI_txt_to_csv.R
└── example_data/
```

---

## Best Practices

- Use short ROI radii (5–20 px) when isolating single kernels, droplets, or ground particles.  
- Use larger ROIs (20–50 px+) for blended samples or aggregated material.  
- Always verify ROI placement using `visualize_roi()` before batch processing.  
- Keep your `.hdr` and `.dat` files in the same folder with matching names.

---

## Issues and Support

If you encounter:

- File loading errors  
- Misaligned wavelengths  
- ROI visualization problems  
- Unexpected pixel counts  

You can open a GitHub issue or contact the maintainer.

---

## Contributing

Contributions are welcome!  
You can:

- Add support for additional HSI formats  
- Improve visualization  
- Add new ROI shapes (rectangles, polygons)  
- Add PCA or preprocessing modules (SNV, MSC, SG)  
- Submit bug fixes  
- Propose new features

---

## License

This project is released under the **MIT License**.  
You may use, modify, and distribute the code with attribution.

