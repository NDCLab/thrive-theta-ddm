#' Find Newest File Utility
#'
#' Provides a function to locate the most recently modified file matching a specific pattern.

#' Find the Newest File Matching a Pattern
#'
#' Searches for files matching the given glob pattern and returns the path
#' to the file with the most recent modification time.
#'
#' @param path_pattern A string containing the glob pattern to match files (e.g., "path/to/*.csv").
#'
#' @return A string containing the path to the newest file, or NULL if no files are found.
#' @export
find_newest_file <- function(path_pattern) {
  # Sys.glob expands wildcards similarly to python's glob
  matching_files <- Sys.glob(path_pattern)
  
  # Check if any files were found
  if (length(matching_files) == 0) {
    print("No matching files found.")
    return(NULL)
  } else {
    # Get metadata for all matching files
    info <- file.info(matching_files)
    
    # Find the row index with the maximum modification time (mtime)
    newest_index <- which.max(info$mtime)
    
    # Extract the filename (row names of the info dataframe contain paths)
    new_file_path <- row.names(info)[newest_index]
    
    print(sprintf("The newest file is: %s", new_file_path))
    
    return(new_file_path)
  }
}

# Usage Example:
# newest <- find_newest_file("path/to/files/*.csv")