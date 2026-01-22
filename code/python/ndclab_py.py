import os
import pandas as pd
from glob import glob
import numpy as np
from datetime import datetime
import re

def replace_outliers_with_nan(df, sd_thresh=3, exclude_cols=None):
    """
    Iterates through ALL columns in the dataframe.
    Replaces values exceeding +/- sd_thresh with NaN.
    
    Parameters:
    - exclude_cols: List of column names to skip (e.g., ['SubjectID', 'Block'])
    """
    # Handle default mutable argument
    if exclude_cols is None:
        exclude_cols = []
        
    for column in df.columns:
        # 1. Skip excluded columns explicitly
        if column in exclude_cols:
            continue

        # 2. Skip non-numeric columns (like string types)
        if not pd.api.types.is_numeric_dtype(df[column]):
            continue

        # 3. Calculate Stats
        mean = df[column].mean()
        std = df[column].std()
        
        upper = mean + (sd_thresh * std)
        lower = mean - (sd_thresh * std)

        # 4. Vectorized Mask (True = Outlier)
        outlier_mask = (df[column] > upper) | (df[column] < lower)
        
        # 5. Logging
        n_outliers = outlier_mask.sum()
        if n_outliers > 0:
            print(f"[Log] {column}: Replaced {n_outliers} outliers ({n_outliers/len(df):.2%} of data).")
        
        # 6. In-place Replacement
        df.loc[outlier_mask, column] = np.nan
        
    return df

def replace_outliers_with_nan_cols(df, columns_to_check, sd_thresh=3):
    """
    Replaces outliers with NaN.
    Raises KeyError if a column is missing.
    Returns the modified dataframe if successful.
    """
    # 1. Validation: Check all columns exist first
    # If this fails, the error is raised and nothing is returned.
    missing_cols = [col for col in columns_to_check if col not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing columns in dataframe: {missing_cols}")

    # 2. Process columns
    for column in columns_to_check:
        if not pd.api.types.is_numeric_dtype(df[column]):
            print(f"Skipping non-numeric column: {column}")
            continue

        mean = df[column].mean()
        std = df[column].std()
        
        upper = mean + (sd_thresh * std)
        lower = mean - (sd_thresh * std)

        # Vectorized replacement
        outlier_mask = (df[column] > upper) | (df[column] < lower)
        
        # Logging
        n_outliers = outlier_mask.sum()
        if n_outliers > 0:
            print(f"[Log] {column}: Replaced {n_outliers} outliers ({n_outliers/len(df):.2%} of data).")
        else:
            print(f"[Log] {column}: No outliers found.")
        
        # In-place Replacement
        df.loc[outlier_mask, column] = np.nan

    return df


def find_newest_file(path): 
    matching_files = glob(path)
    
    # Check if any files were found
    if not matching_files:
        print("No matching files found.")
    else:
        # Find the newest file based on modification time
        new_file_path = max(matching_files, key=os.path.getmtime)
        print(f"The newest file is: {new_file_path}")

        return new_file_path
