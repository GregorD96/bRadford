#' Calculates protein concentration
#'
#' @name bradford_analysis
#'
#' @description All-in-one function for computing protein concentrations determined by a bradford protein assay.
#' Requires protein standards of known concentration to plot a standard curve used for the accurate determination
#' of sample concentrations.
#'
#' NOTE: Concentration of diluted sample should be in the rang of measured standards.
#'
#' @param file path to input xls-file
#' @param sheet_calibration sheet containing the standards
#' @param sheet_sample sheet containing the samples
#' @param dil_fctr numeric value specifying the sample dilution (e.g: 1:10 dilution -> dil_fctr = 10 )
#' @param output_dir directory to which the output table and standard curve is saved to as a csv and
#' pdf file, respectively
#'
#' @return A tibble containing calculated protein concentrations and a standard curve plot
#' @export
#'
#' @examples
#' bradford_analysis(file = "path/to/Bradford_test.xls", sheet_calibration = 4, sheet_sample = 5,
#' dil_fctr = 4, output_dir = "/your/favorite/output/directory")




bradford_analysis <- function(file,sheet_calibration,sheet_sample,
                     dil_fctr,output_dir){

  # Read in both files as individual dataframes ----

  dil_series <- readxl::read_xls(path = file, sheet = sheet_calibration)
  samples <- readxl::read_xls(path = file, sheet = sheet_sample)

  # Save additional info in vectors ----

  dil_fctr <- dil_fctr                      # times how much to reach your initial sample concentration

  # Take the mean absorbance and normalize by subtracting the blank ----

  dil_series_clean <- dil_series %>%
    stats::na.omit() %>%
    dplyr::group_by(dplyr::across(c(1, 2))) %>%
    dplyr::summarize(absorbance = mean(.data$absorbance_series))

  dil_series_clean <- dil_series_clean %>%
    dplyr::group_by(dplyr::across(c(1, 2))) %>%
    dplyr::summarize(absorbance = (.data$absorbance - dil_series_clean$absorbance[dil_series_clean[1] == "blank"])) %>%
    dplyr::arrange(.data$concentration)

  # Calculate the OLS for the standard curve ----

  bradford_fit <- stats::lm(absorbance ~ concentration, data = dil_series_clean)

  # Plot the standard curve including the most important information from the calculated model ----

  standard_curve <- ggplot2::ggplot(bradford_fit$model, ggplot2::aes_string(x = names(bradford_fit$model)[2], y = names(bradford_fit$model)[1])) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::labs(title = paste("R^2 = ",signif(summary(bradford_fit)$r.squared, 4),
                       "Intercept =",signif(bradford_fit$coef[[1]],4 ),
                       " Slope =",signif(bradford_fit$coef[[2]], 4),
                       " P =",signif(summary(bradford_fit)$coef[2,4], 4)),
         x = paste("concentration [\u00B5g/\u00B5l]"),y = "Absorbance [595nm]") +       # 00B5 is Unicode for µ-symbol and /u escapes it to be converted -> needed since R packages can only contain ASCII symbols
    ggplot2::theme_linedraw()

  print(standard_curve)

  # Save the plot in a previously specified directory ----

  ggplot2::ggsave(filename = "standard_curve_bradford.pdf",
         device = "pdf",
         plot = standard_curve,
         path = output_dir,
         height = 10,
         width = 15,
         units = "in")

  # Calculate the sample concentration ----

  samples_clean <- samples %>%
    stats::na.omit() %>%
    dplyr::group_by(names) %>%
    dplyr::summarize(absorbance = mean(.data$absorbance_sample))

  samples_clean <- samples %>%
    dplyr::mutate(absorbance_sample = .data$absorbance_sample - samples$absorbance_sample[samples[1] == "blank"])
  samples_clean <- samples_clean[samples_clean$absorbance_sample != 0,]

  samples_clean <- samples_clean %>%
    dplyr::mutate(concentration_mg_per_ml = (((.data$absorbance_sample - bradford_fit$coefficients[[1]]) / bradford_fit$coefficients[[2]]) * dil_fctr))

  utils::write.csv(samples_clean, file = paste0(output_dir,"/","results_bradford.csv"))

  print(samples_clean, n = nrow(samples_clean))

}
