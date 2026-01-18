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
  colnames(result)[1] <- "β"
  colnames(result)[2] <- "SE"
  colnames(result)[3] <- "df"
  colnames(result)[4] <- "t"
  colnames(result)[5] <- "p"
  result <- add_column(result, Parameter = rename_parameters(output_model$rowname), .before=1)
  
  # Print the result
  final_table <- nice_table(result, italics = c(2, 3, 4, 5, 6), col.format.p = 6, note = "Significance codes: * p < .05, ** p < .01, *** p < .001. Degrees of freedom for the fixed effects were estimated using Satterthwaite's approximation.")
  final_table <- final_table %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    line_spacing(space = 1, part = "all") %>%
    padding(padding.top = .5, padding.bottom = .5, part = "all")
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
    "ICPS_early_MOTOR_diff_collapsed" = "ICPS midlateral"
    # Add more mappings here if needed in the future
  )
  
  # Initialize an empty vector to store the new names
  new_names <- character(length(parameter_vector))
  
  # Iterate over each column name in the input vector
  for (i in seq_along(parameter_vector)) {
    original_name <- parameter_vector[i]
    processed_name <- original_name # Default to original if no rules apply
    
    # Check if the name contains interaction terms (indicated by ":")
    if (grepl(":", original_name)) {
      # Split the name into individual components based on ":"
      parts <- strsplit(original_name, ":")[[1]]
      
      # Rename each part
      renamed_parts <- sapply(parts, function(part) {
        # Remove trailing digits for lookup
        part_base <- gsub("\\d+$", "", part) 
        if (part_base %in% names(name_map)) {
          # If the base part is in the map, use the mapped name
          return(name_map[part_base])
        } else {
          # Otherwise, keep the original part (with digits if any)
          return(part) 
        }
      }, USE.NAMES = FALSE)
      
      # Join the renamed parts with " * "
      processed_name <- paste(renamed_parts, collapse = " * ")
      
    } else {
      # If it's not an interaction term, process it as a single term
      # Remove trailing digits for lookup
      original_name_base <- gsub("\\d+$", "", original_name)
      if (original_name_base %in% names(name_map)) {
        # If the base name is in the map, use the mapped name
        processed_name <- name_map[original_name_base]
      }
      # If not in the map, it remains as 'original_name' (already set in processed_name)
    }
    
    # Store the processed name
    new_names[i] <- processed_name
  }
  
  # Return the vector of new column names
  return(new_names)
}