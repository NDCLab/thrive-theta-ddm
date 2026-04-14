# Code

This folder contains the code used to process behavioral and EEG data, fit models, compute statistics, and generate figures for the project.


### Directory Structure

* **`behavior/`**: Scripts for analyzing and processing behavioral log files, and creating valid behavioral datasets.
* **`preprocessing-eeg/`**: Scripts and pipelines (e.g., MADE) for preprocessing raw EEG data. Also includes a `tf/` subfolder with time-frequency scripts.
* **`postprocessing/`**: Scripts for computing ERP means, generating time-frequency (TF) arrays, and combining preprocessed data into formats suitable for modeling and analysis (e.g., CSVs).
* **`ddm/`**: Code for running and fitting the Shrinking-Spotlight Drift Diffusion Model (SSP-DDM) on the processed data, as well as code for parameter recovery.
* **`statistics/`**: R scripts and Jupyter notebooks for performing statistical analyses (e.g., Linear Mixed Models) on the resulting datasets.
* **`figures/`**: Scripts for generating plots and visualizations (e.g., ERP plots).
* **`matlab/`**: Shared utility functions and scripts for MATLAB operations.
* **`python/`**: Shared utility modules and functions for Python operations.

### Expected Order of Execution

1. **Preprocessing**:
   - Run behavioral processing scripts in `behavior/`.
   - Run EEG preprocessing scripts (like the MADE pipeline) in `preprocessing-eeg/`.
2. **Postprocessing**:
   - Use scripts in `postprocessing/` to compute ERP means, TF arrays, and generate unified datasets (e.g., the DDM CSV).
3. **Modeling**:
   - Run the SSP-DDM model fitting on the processed behavioral and EEG data using scripts in `ddm/`.
4. **Statistics**:
   - Perform statistical analyses on the final processed datasets and model outputs using scripts in `statistics/`.
5. **Figures**:
   - Generate final visualizations using scripts in `figures/`.

### Project Notes
