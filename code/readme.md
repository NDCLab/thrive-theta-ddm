# Code

## Overview
This directory contains all source code for the project, organized by analysis stage and modality.

## Subdirectories

### `behavior/`
Python scripts for parsing, cleaning, and analyzing behavioral data from the Flanker task. Includes quality control checks (`check_status.py`, `check_subject_csv.py`).

### `ddm/`
Code for the Shrinking Spotlight Protocol (SSP) Drift Diffusion Model.
*   **R/C++**: `SSP_DDM_fitting.R` and `simSSP_model_GB_noScale.cpp` implement the model fitting and simulation.
*   **Python**: Scripts for batch submission and result aggregation.

### `preprocessing-eeg/`
The MADE (Maryland Analysis of Developmental EEG) preprocessing pipeline.
*   **MATLAB**: Core pipeline functions (`MADE_pipeline.m`, `preprocess_eeg_piece.m`).
*   **Python**: Batch management scripts (`run_MADE_batch.py`, `check_preprocessed_files.py`).

### `postprocessing/`
Scripts to extract measures from preprocessed EEG data.
*   **ERPs**: `compute_erp_means.py`
*   **Time-Frequency**: `create_tf_arrays.py`, `compute_means_TF.py`, `compute_means_ICPS.py`.

### `figures/`
Code to generate visualizations, such as Grand Average ERP plots (`plot_erp.py`, `plot_erp.m`).

### `statistics/`
R scripts for running Linear Mixed-Effects Models (LMM) and exporting results to APA-style tables.
