# MALDI-TOF MS Data Analysis Pipeline for Microbial spectra comparison
# Packages and Methods:
# - MALDIquant/MALDIquantForeign: Mass spectrometry data processing
# - ggplot2: Advanced visualization
# - cowplot: Professional plot formatting

rm(list = ls())
gc()

# 1. INSTALL AND LOAD PACKAGES -------------------------------------------------
required_packages <- c("MALDIquant", "MALDIquantForeign", "ggplot2", "cowplot")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(MALDIquant)
library(MALDIquantForeign)
library(ggplot2)
library(cowplot)

# 2. SETUP PARAMETERS AND PATHS -------------------------------------------------
base_path <- "your_path_where_are_strain_spectra"
strain_folders <- c("strain_1", "strain_2")
output_dir <- "path_to_dir_results"
dir.create(output_dir, showWarnings = FALSE)

# Scientific names with proper formatting
strain_names <- list(
  "strain_1" = expression(paste(italic("name_of_strain_1"), "collection_number_1")),
  "strain_2" = expression(paste(italic("name_of_strain_2"), "collection_number_2"))
)

# Processing parameters
SNR_threshold <- 2
halfWindowSize <- 10
analysis_range <- c(5000, 10000) # m/z range for analysis

# 3. SPECTRA PROCESSING FUNCTION -----------------------------------------------
process_spectra <- function(spectra) {
  tryCatch({
    # 1. Average technical replicates
    if (length(spectra) > 1) {
      spectra <- averageMassSpectra(spectra, method = "mean")
    }
    
    # 2. Trim to analysis range
    spectra <- trim(spectra, range = analysis_range)
    
    # 3. Transform intensities (variance stabilization)
    spectra <- transformIntensity(spectra, method = "sqrt")
    
    # 4. Smooth spectra (noise reduction)
    spectra <- smoothIntensity(spectra, method = "SavitzkyGolay", 
                               halfWindowSize = halfWindowSize)
    
    # 5. Remove baseline
    spectra <- removeBaseline(spectra, method = "SNIP", iterations = 100)
    
    # 6. Normalize by TIC
    spectra <- calibrateIntensity(spectra, method = "TIC")
    
    return(spectra)
  }, error = function(e) {
    message("Processing error: ", e$message)
    return(NULL)
  })
}

# 4. DATA PROCESSING PIPELINE --------------------------------------------------
results <- list()

for (strain in strain_folders) {
  message("Processing: ", strain)
  
  # Import raw spectra
  spectra <- import(file.path(base_path, strain), type = "auto")
  
  # Process spectra
  processed <- process_spectra(spectra)
  
  if (!is.null(processed)) {
    # Detect peaks (if needed)
    peaks <- detectPeaks(processed, method = "MAD", 
                         SNR = SNR_threshold, 
                         halfWindowSize = halfWindowSize)
    
    results[[strain]] <- list(
      spectrum = processed,
      peaks = peaks
    )
    
    saveRDS(results[[strain]], file.path(output_dir, paste0(strain, "_processed.rds")))
  }
}

# 5. VISUALIZATION -------------------------------------------------------------
if (length(results) == 2) {
  # Prepare plot data
  plot_data <- data.frame()
  
  for (strain in names(results)) {
    # Spectrum data
    df_spectrum <- data.frame(
      mz = mass(results[[strain]]$spectrum),
      intensity = intensity(results[[strain]]$spectrum),
      strain = strain
    )
    
    # Normalize intensity to percentage (relative to the maximum intensity)
    df_spectrum$intensity <- df_spectrum$intensity / max(df_spectrum$intensity) * 100
    
    plot_data <- rbind(plot_data, df_spectrum)
  }
  
  # Create spectra plot
  p_spectra <- ggplot(plot_data, aes(x = mz, y = intensity, color = strain)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    scale_color_manual(
      name = "",
      values = c("#1E88E5", "#D81B60"),
      labels = strain_names
    ) +
    labs(title = "MALDI-TOF MS Spectra Comparison",
         x = "m/z", 
         y = "Intensity (%)") +  # Updated y-axis label
    coord_cartesian(xlim = analysis_range) +
    theme_cowplot(14) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(
        size = 12, 
        face = "italic",
        hjust = 0.5,
        margin = margin(t = 5)
      ),
      legend.title = element_blank(),
      legend.justification = "center",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  # Save and display
  ggsave(file.path(output_dir, "spectra_comparison.png"), 
         p_spectra, width = 12, height = 6, dpi = 300, bg = "white")
  print(p_spectra)
} else {
  message("Visualization requires exactly 2 processed spectra")
}