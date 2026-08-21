library(officer)

combine_docs <- function(path, order_pattern = c("power", "ITPS", "ERN_laplacian", "ERN_", "ICPS", "ratio", "a_diff", "pea_", "peri_rt_", "pes_")) {

    output_path = sprintf("%s/combined", path)
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    # Define the input files and the output file name
    input_files <- Sys.glob(sprintf("%s/*.docx", path))
    
    # Define your desired order keywords
    # order_pattern <- c("power", "ITPS", "ERN_laplacian", "ERN_", "ICPS", "ratio", "a_diff", "pea_", "peri_rt_", "pes_")
    
    # Reorder the file list
    # Use unlist to turn the results into a simple vector of indices
    indices <- unlist(lapply(order_pattern, function(x) grep(x, input_files)))
    
    # Reorder the file list based on those indices
    input_files <- unique(input_files[indices])
    output_file <- sprintf("%s/combined_tables.docx", output_path)
    
    # Start with the first document as the base
    # read_docx() loads the document
    combined_doc <- read_docx(path = input_files[1])
    
    # Loop through the remaining documents (starting from the second one)
    for (i in 2:length(input_files)) {
      # Add a page break before the next document's content
      combined_doc <- body_add_break(combined_doc)
      
      # Append the content of the next document
      # body_add_docx() uses a Microsoft Word feature to insert the file content
      combined_doc <- body_add_docx(combined_doc, src = input_files[i])
    }
    
    # Save the final combined document
    print(combined_doc, target = output_file)
}
