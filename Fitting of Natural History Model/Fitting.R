#################################################################
#  Fitting_selected_sites.R
#
#  Simple fitting script for:
#    - lung
#    - colorectal
#    - pancreas
#
#  No SIR multipliers are applied.
#
#  Fits:
#    lung:        (OMST, DMST) = (3.5, 0.5), (3.5, 1)
#    colorectal:  (4, 0.5), (4, 1)
#    pancreas:    (2, 0.5), (2, 1)
#################################################################

library(here)
library(openxlsx)

base_path <- here("Fitting of Natural History Model")
# -------------------------------------------------------------
# Load fitting functions
# -------------------------------------------------------------
source(file.path(base_path, "natural_history_code.R"))

# -------------------------------------------------------------
# Load cancer incidence data
# -------------------------------------------------------------
filepath <- file.path(base_path, "MCED_Cancer_List_US_2011_2015.xlsx")
sheets <- openxlsx::getSheetNames(filepath)
cancer_list <- lapply(sheets, openxlsx::read.xlsx, xlsxFile = filepath)
names(cancer_list) <- tolower(sheets)

# -------------------------------------------------------------
# Output folder
# -------------------------------------------------------------
filepath2 <- file.path(base_path, "Fitted_SEER_Data_Selected")
if (!dir.exists(filepath2)) dir.create(filepath2, recursive = TRUE)

# -------------------------------------------------------------
# Read and process SEER data
# -------------------------------------------------------------
read_data <- function(the_data, multiplier = 1) {
  the_data$PY <- the_data$Pop * the_data$years
  the_data$rate <- multiplier * the_data$Count / the_data$PY * 1e5
  the_data$Count <- the_data$Count * multiplier
  the_data <- the_data[the_data$midage > 0 & the_data$midage < 80, , drop = FALSE]
  the_data
}

# -------------------------------------------------------------
# OMST-DMST pairs
# -------------------------------------------------------------
fit_specs <- list(
  lung = list(pairs = list(c(3.5, 0.5), c(3.5, 1))),
  colorectal = list(pairs = list(c(4, 0.5), c(4, 1))),
  pancreas = list(pairs = list(c(2, 0.5), c(2, 1)))
)

# -------------------------------------------------------------
# Fit models and save
# -------------------------------------------------------------
fit_selected_models <- function(fit_specs) {
  
  for (site_clean in names(fit_specs)) {
    
    pairs <- fit_specs[[site_clean]]$pairs
    OMST_site <- sapply(pairs, `[[`, 1)
    DMST_site <- sapply(pairs, `[[`, 2)
    
    the_data <- read_data(cancer_list[[site_clean]])
    
    metadata_var <- paste0("metadata_", site_clean)
    assign(
      metadata_var,
      data.frame(
        fitID = seq_along(OMST_site),
        site = site_clean,
        OMST = OMST_site,
        DMST = DMST_site
      ),
      envir = .GlobalEnv
    )
    
    output_var <- paste0(site_clean, "out")
    assign(
      output_var,
      mapply(
        FUN = "get_fit",
        OMST_site,
        DMST_site,
        MoreArgs = list(
          the_data = the_data,
          num_seeds = 1,
          k = 16,
          mean1 = 3
        ),
        SIMPLIFY = FALSE
      ),
      envir = .GlobalEnv
    )
    
    save(
      list = c(metadata_var, output_var),
      file = file.path(filepath2, paste0(site_clean, "_outfile.Rdata"))
    )
  }
}

# -------------------------------------------------------------
# Run fitting
# -------------------------------------------------------------
fit_selected_models(fit_specs)