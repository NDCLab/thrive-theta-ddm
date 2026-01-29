library(dplyr)
library(tibble)
library(flextable)

lmer_export_apa <- function(model, path) {
  
  output_model <- as.data.frame(summary(model)$coefficients)
  conf_model <- as.data.frame(confint(model))
  conf_model <- conf_model[3:nrow(conf_model), 1:2]
  rownames(conf_model) <- rownames(output_model)
  colnames(conf_model) <- c("CI_lower", "CI_upper")
  
  # Convert row names to a column in both data frames
  output_model$rowname <- rownames(output_model)
  # <- output_model %>% rownames_to_column('rowname')
  conf_model$rowname <- rownames(conf_model)
  
  # Perform left join by the rowname column
  result <- left_join(output_model, conf_model, by = "rowname")
  
  # Optionally, set row names back to the rowname column
  result <- result %>% column_to_rownames('rowname')
  colnames(result)[1] <- "Beta"
  colnames(result)[2] <- "SE"
  colnames(result)[3] <- "df"
  colnames(result)[4] <- "t"
  colnames(result)[5] <- "p"
  result <- add_column(result, Parameter = rename_parameters(output_model$rowname), .before=1)
  
  # Print the result
  final_table <- nice_table(result, italics = c(2, 3, 4, 5, 6), col.format.p = 6, note = "Significance codes: * p < .05, ** p < .01, *** p < .001. Degrees of freedom for the fixed effects were estimated using Satterthwaite's approximation.")

    
  final_table <- final_table %>%
    set_header_labels(Beta = "\u03B2") %>%  # Map 'Beta' to Greek symbol
    set_table_properties(width = 1, layout = "autofit") %>%
    line_spacing(space = 1, part = "all") %>%
    padding(padding.top = .5, padding.bottom = .5, part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    fontsize(size = 9, part = "footer")
    # padding(i = ~ grepl("\\*", Parameter), j = 1, padding.left = 20)
  flextable::save_as_docx(final_table, path = path)
  
}

rename_parameters <- function(parameter_vector) {
  
  # Define the mapping for single term renaming
  # This named vector stores the old base names (as names) and new names (as values)
  name_map <- c(
    "acc" = "Accuracy",
    "soc" = "Condition",
    "age_m" = "Age",
    "sex" = "Sex",
    "dp_inperson" = "Peer mode",
    "ICPS_OCC_diff_collapsed" = "ICPS posterolateral",
    "ICPS_DLPFC_diff_collapsed" = "ICPS frontolateral",
    "ICPS_MOTOR_diff_collapsed" = "ICPS midlateral",
    "ICPS_early_OCC_diff_collapsed" = "ICPS posterolateral",
    "ICPS_early_DLPFC_diff_collapsed" = "ICPS frontolateral",
    "ICPS_early_MOTOR_diff_collapsed" = "ICPS midlateral",
    "bfne_b_scrdTotal_s1_r1_e1" = "BFNE"
    # Add more mappings here if needed in the future
  )
  
# Helper function to check map with or without digits
  get_mapped_name <- function(x, map) {
    # 1. Check exact match first (fixes your BFNE issue)
    if (x %in% names(map)) {
      return(map[[x]])
    }
    
    # 2. Check match after removing trailing digits (for factor levels)
    x_base <- gsub("\\d+$", "", x)
    if (x_base %in% names(map)) {
      return(map[[x_base]])
    }
    
    # 3. No match found, return original
    return(x)
  }
  
  # Initialize an empty vector
  new_names <- character(length(parameter_vector))
  
  # Iterate over each column name
  for (i in seq_along(parameter_vector)) {
    original_name <- parameter_vector[i]
    processed_name <- original_name 
    
    # Check if interaction term
    if (grepl(":", original_name)) {
      parts <- strsplit(original_name, ":")[[1]]
      
      # Rename each part using the helper
      renamed_parts <- sapply(parts, get_mapped_name, map = name_map, USE.NAMES = FALSE)
      
      processed_name <- paste(renamed_parts, collapse = " * ")
      
    } else {
      # Single term processing
      processed_name <- get_mapped_name(original_name, name_map)
    }
    
    new_names[i] <- processed_name
  }
  
  return(new_names)
}