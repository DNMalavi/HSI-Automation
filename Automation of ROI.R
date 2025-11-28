# ======================================================================
# ROI Extraction for Hyperspectral ENVI Images (HDR + DAT)
# Author: Derick Malavi
#
# Description:
#   - Creates an elliptical ROI on an ENVI hyperspectral cube
#   - Extracts the average spectrum within the ROI
#   - Supports batch and recursive directory processing
#   - Outputs combined spectra as a CSV-ready data frame
#
# Dependencies: terra, ggplot2
# ======================================================================

# Install necessary packages if not already installed:
# install.packages(c("terra", "ggplot2"))

library(terra)
library(ggplot2)

# ----------------------------------------------------------------------
# Helper: read wavelengths from ENVI .hdr file (if available)
# ----------------------------------------------------------------------
read_envi_wavelengths <- function(hdr_file, expected_bands = NULL) {
  if (!file.exists(hdr_file)) {
    return(NULL)
  }
  
  wavelengths <- NULL
  
  try({
    hdr_text <- readLines(hdr_file, warn = FALSE)
    
    wl_start <- grep("wavelength\\s*=\\s*\\{", hdr_text, ignore.case = TRUE)
    if (length(wl_start) == 0) {
      return(NULL)
    }
    
    wl_end <- grep("\\}", hdr_text[wl_start:length(hdr_text)])[1]
    wl_end <- wl_start + wl_end - 1
    
    wl_text  <- paste(hdr_text[(wl_start + 1):(wl_end - 1)], collapse = " ")
    wl_clean <- gsub("\\s+", "", wl_text)
    wavelengths <- as.numeric(strsplit(wl_clean, ",")[[1]])
    
    if (!is.null(expected_bands) && length(wavelengths) != expected_bands) {
      wavelengths <- NULL
    }
  }, silent = TRUE)
  
  wavelengths
}

# ----------------------------------------------------------------------
# Core function: create elliptical ROI and extract average spectrum
# ----------------------------------------------------------------------
process_hyperspectral_image <- function(file_path,
                                        center_x = NULL,
                                        center_y = NULL,
                                        x_axis_length = 10,
                                        y_axis_length = 10,
                                        sample_name = NULL,
                                        verbose = TRUE) {
  stopifnot(is.character(file_path), length(file_path) == 1)
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  if (x_axis_length <= 0 || y_axis_length <= 0) {
    stop("x_axis_length and y_axis_length must be positive.")
  }
  
  file_ext <- tolower(tools::file_ext(file_path))
  
  if (file_ext == "hdr") {
    hdr_file  <- file_path
    data_file <- sub("\\.hdr$|\\.HDR$", ".dat", file_path, ignore.case = TRUE)
  } else if (file_ext == "dat") {
    data_file <- file_path
    hdr_file  <- sub("\\.dat$|\\.DAT$", ".hdr", file_path, ignore.case = TRUE)
  } else {
    stop("File must be either .hdr or .dat format.")
  }
  
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file)
  }
  if (!file.exists(hdr_file)) {
    stop("Header file not found: ", hdr_file)
  }
  
  if (is.null(sample_name)) {
    sample_name <- tools::file_path_sans_ext(basename(data_file))
  }
  
  hs_data <- tryCatch({
    rast(data_file)
  }, error = function(e) {
    message("  terra::rast() failed; trying raster::brick() fallback...")
    if (!requireNamespace("raster", quietly = TRUE)) {
      install.packages("raster")
    }
    suppressPackageStartupMessages(library(raster))
    
    temp_brick <- tryCatch({
      brick(data_file)
    }, error = function(e2) {
      stop("  Both terra::rast() and raster::brick() failed for file: ", data_file)
    })
    
    rast(temp_brick)
  })
  
  img_dims  <- dim(hs_data)
  num_rows  <- img_dims[1]
  num_cols  <- img_dims[2]
  num_bands <- img_dims[3]
  
  if (verbose) {
    message("Loaded image: ", basename(data_file),
            " [rows: ", num_rows,
            ", cols: ", num_cols,
            ", bands: ", num_bands, "]")
  }
  
  r         <- hs_data[[1]]
  all_cells <- 1:ncell(r)
  xy        <- xyFromCell(r, all_cells)
  
  if (is.null(center_x)) {
    center_x <- mean(range(xy[, 1]))
  }
  if (is.null(center_y)) {
    center_y <- mean(range(xy[, 2]))
  }
  
  dx        <- (xy[, 1] - center_x) / x_axis_length
  dy        <- (xy[, 2] - center_y) / y_axis_length
  dist_norm <- sqrt(dx^2 + dy^2)
  
  mask_values <- dist_norm <= 1
  
  mask_raster <- r
  values(mask_raster) <- mask_values
  
  roi_cells <- all_cells[mask_values]
  
  if (length(roi_cells) == 0) {
    stop("Ellipse ROI does not include any pixels. Check center and axis lengths.")
  }
  
  extracted <- terra::extract(hs_data, roi_cells, cells = FALSE)
  pixel_values <- as.matrix(extracted[, -1, drop = FALSE])
  
  avg_spectrum <- colMeans(pixel_values, na.rm = TRUE)
  
  wavelengths <- read_envi_wavelengths(hdr_file, expected_bands = num_bands)
  if (is.null(wavelengths)) {
    wavelengths <- seq_len(num_bands)
  }
  
  result_df <- data.frame(
    sample_name = sample_name,
    center_x    = center_x,
    center_y    = center_y
  )
  
  spectral_df <- as.data.frame(t(avg_spectrum))
  colnames(spectral_df) <- paste0("Band_", wavelengths)
  
  result_df <- cbind(result_df, spectral_df)
  
  list(
    spectrum    = result_df,
    mask        = mask_raster,
    center      = c(center_x, center_y),
    axes        = c(x_axis_length, y_axis_length),
    num_pixels  = length(roi_cells)
  )
}

# ----------------------------------------------------------------------
# Batch processing: non-recursive
# ----------------------------------------------------------------------
batch_process_directory <- function(directory,
                                    center_x = NULL,
                                    center_y = NULL,
                                    x_axis_length = 10,
                                    y_axis_length = 10,
                                    file_pattern = "\\.dat$|\\.DAT$") {
  if (!dir.exists(directory)) {
    stop("Directory does not exist: ", directory)
  }
  
  data_files <- list.files(
    directory,
    pattern     = file_pattern,
    full.names  = TRUE,
    ignore.case = TRUE
  )
  
  if (length(data_files) == 0) {
    stop("No matching files found in the specified directory: ", directory)
  }
  
  all_spectra      <- vector("list", length(data_files))
  successful_count <- 0L
  
  for (i in seq_along(data_files)) {
    cat(sprintf("Processing file %d of %d: %s\n",
                i, length(data_files), basename(data_files[i])))
    
    result <- try(
      process_hyperspectral_image(
        file_path      = data_files[i],
        center_x       = center_x,
        center_y       = center_y,
        x_axis_length  = x_axis_length,
        y_axis_length  = y_axis_length,
        verbose        = FALSE
      ),
      silent = TRUE
    )
    
    if (!inherits(result, "try-error")) {
      successful_count <- successful_count + 1L
      all_spectra[[successful_count]] <- result$spectrum
      cat(sprintf("  Extracted ROI with %d pixels\n", result$num_pixels))
    } else {
      warning("Error processing ", basename(data_files[i]))
    }
  }
  
  if (successful_count == 0L) {
    stop("No files were successfully processed.")
  }
  
  all_spectra <- all_spectra[seq_len(successful_count)]
  combined_spectra <- do.call(rbind, all_spectra)
  
  combined_spectra
}

# ----------------------------------------------------------------------
# Batch processing: recursive search
# ----------------------------------------------------------------------
batch_process_directory_recursive <- function(root_directory,
                                              center_x = NULL,
                                              center_y = NULL,
                                              x_axis_length = 10,
                                              y_axis_length = 10,
                                              file_pattern = "\\.dat$|\\.DAT$") {
  if (!dir.exists(root_directory)) {
    stop("Root directory does not exist: ", root_directory)
  }
  
  cat("Searching for matching files in all subdirectories...\n")
  
  data_files <- list.files(
    root_directory,
    pattern     = file_pattern,
    full.names  = TRUE,
    recursive   = TRUE,
    ignore.case = TRUE
  )
  
  cat(sprintf("Found %d matching files in total.\n", length(data_files)))
  
  if (length(data_files) == 0) {
    stop("No matching files found in the directory tree.")
  }
  
  root_path <- normalizePath(root_directory, winslash = "/", mustWork = TRUE)
  
  all_spectra      <- vector("list", length(data_files))
  successful_count <- 0L
  
  for (i in seq_along(data_files)) {
    file_path <- normalizePath(data_files[i], winslash = "/", mustWork = TRUE)
    
    rel_path <- if (startsWith(file_path, root_path)) {
      substring(file_path, nchar(root_path) + 2L)
    } else {
      basename(file_path)
    }
    
    cat(sprintf("Processing file %d of %d: %s\n",
                i, length(data_files), rel_path))
    
    result <- try(
      process_hyperspectral_image(
        file_path      = file_path,
        center_x       = center_x,
        center_y       = center_y,
        x_axis_length  = x_axis_length,
        y_axis_length  = y_axis_length,
        verbose        = FALSE
      ),
      silent = TRUE
    )
    
    if (!inherits(result, "try-error")) {
      successful_count <- successful_count + 1L
      result$spectrum$folder_path <- dirname(rel_path)
      all_spectra[[successful_count]] <- result$spectrum
      cat(sprintf("  Extracted ROI with %d pixels\n", result$num_pixels))
    } else {
      warning("Error processing ", rel_path)
    }
  }
  
  if (successful_count == 0L) {
    stop("No files were successfully processed.")
  }
  
  all_spectra <- all_spectra[seq_len(successful_count)]
  combined_spectra <- do.call(rbind, all_spectra)
  
  combined_spectra
}

# ----------------------------------------------------------------------
# Visualization helpers
# ----------------------------------------------------------------------
visualize_roi <- function(result) {
  mask     <- result$mask
  spectrum <- result$spectrum
  center   <- result$center
  
  mask_df <- as.data.frame(mask, xy = TRUE)
  names(mask_df) <- c("x", "y", "mask_value")
  
  ggplot(mask_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = as.factor(mask_value))) +
    scale_fill_manual(
      values = c("0" = NA, "1" = "red"),
      labels = c("Outside ROI", "Inside ROI"),
      na.value = NA,
      guide = guide_legend(title = "ROI")
    ) +
    coord_fixed() +
    geom_point(
      data = data.frame(x = center[1], y = center[2]),
      aes(x = x, y = y),
      color = "blue",
      size  = 3
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(
      paste0(
        "ROI for sample: ", spectrum$sample_name[1],
        "\nCenter: (", round(spectrum$center_x[1], 2), ", ",
        round(spectrum$center_y[1], 2), ")",
        "\nNumber of pixels in ROI: ", result$num_pixels
      )
    )
}

plot_spectrum <- function(spectrum) {
  sample_name <- spectrum$sample_name[1]
  center_x    <- spectrum$center_x[1]
  center_y    <- spectrum$center_y[1]
  
  spectrum_data <- spectrum[, !names(spectrum) %in%
                              c("sample_name", "center_x", "center_y", "folder_path"),
                            drop = FALSE]
  
  wavelengths <- as.numeric(gsub("Band_", "", colnames(spectrum_data)))
  
  plot_data <- data.frame(
    Wavelength = wavelengths,
    Reflectance = as.numeric(spectrum_data[1, ])
  )
  
  ggplot(plot_data, aes(x = Wavelength, y = Reflectance)) +
    geom_line() +
    theme_minimal() +
    xlab("Wavelength") +
    ylab("Reflectance") +
    ggtitle(
      paste0(
        "Average Spectrum for Sample: ", sample_name,
        "\nCenter: (", round(center_x, 2), ", ", round(center_y, 2), ")"
      )
    )
}

# ----------------------------------------------------------------------
# Example usage (commented – adjust paths and uncomment to run)
# ----------------------------------------------------------------------
# setwd("Y:/Derick/2025_Coffee")
#
# # 1. Single file
# result <- process_hyperspectral_image(
#   "test_sam_2025-04-11_07-50-13_refl.dat",
#   center_x      = 100,
#   center_y      = 150,
#   x_axis_length = 15,
#   y_axis_length = 15
# )
# spectrum <- result$spectrum
# visualize_roi(result)
# plot_spectrum(spectrum)
#
# # 2. Non-recursive batch
# all_spectra <- batch_process_directory(
#   directory      = "Y:/Derick/2025_Coffee",
#   x_axis_length  = 8,
#   y_axis_length  = 4
# )
# write.csv(all_spectra, "Instant_Coffee_extracted_spectra_non_recursive.csv", row.names = FALSE)
#
# # 3. Recursive batch
# all_spectra_recursive <- batch_process_directory_recursive(
#   root_directory = "Y:/Derick/2025_Coffee",
#   x_axis_length  = 8,
#   y_axis_length  = 4
# )
# write.csv(all_spectra_recursive, "Instant_Coffee_extracted_spectra_recursive.csv", row.names = FALSE)
