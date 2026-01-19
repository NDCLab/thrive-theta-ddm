# Thrive: Theta DDM Analysis

## Project Goal
This repository contains the data analysis pipeline for the "Thrive" project. The goal is to investigate the neural mechanisms underlying social influence on decision-making using the Shrinking Spotlight Protocol (SSP) Drift Diffusion Model (DDM) and EEG time-frequency analysis. We aim to understand how social observation affects cognitive control processes, specifically examining ERP components (ERN, Pe) and theta band oscillations.

## Background & Design
The study employs a Flanker task with social (Observed) and non-social (Alone) conditions.
Key analyses include:
1.  **Behavioral Analysis**: Analyzing accuracy and reaction times (RT) to model decision-making processes.
2.  **Computational Modeling**: Fitting the SSP-DDM to behavioral data to estimate parameters like boundary separation, non-decision time, and attentional focus.
3.  **EEG Analysis**:
    *   **Preprocessing**: Automated pipeline (MADE) including filtering, artifact rejection (FASTER/ADJUST), and ICA.
    *   **ERP Analysis**: Computing Event-Related Potentials (ERN, Pe) for error and correct trials.
    *   **Time-Frequency Analysis**: Examining power and phase synchrony (ITPS, ICPS) in the theta and delta bands.

## Directory Structure

```yml
project-name
├── code
│   ├── behavior            # Scripts for behavioral data analysis and cleaning
│   ├── ddm                 # Scripts for SSP-DDM model fitting (R, C++, Python)
│   ├── figures             # Scripts to generate plots (ERPs, Topographies)
│   ├── matlab              # Helper MATLAB functions
│   ├── postprocessing      # Scripts for aggregating EEG metrics (ERP, TF, ICPS)
│   ├── preprocessing-eeg   # EEG preprocessing pipeline (MADE)
│   └── statistics          # R scripts for statistical analysis (LMM)
├── derivatives             # Generated data (preprocessed EEG, summary CSVs)
├── sourcedata              # Raw input data (checked)
└── results                 # Figures and statistical outputs
```

## Setup & Usage

### Prerequisites
*   **Python**: `pandas`, `numpy`, `scipy`, `mne`, `matplotlib`, `h5py`
*   **R**: `DEoptim`, `Rcpp`, `dplyr`, `flextable`
*   **MATLAB**: EEGLAB with plugins (FASTER, ADJUST, firfilt)
*   **C++**: Compiler compatible with Rcpp

### Workflow
1.  **Behavioral Processing**:
    *   Run `code/behavior/behavior_analysis.py` to aggregate PsychoPy data.
    *   Run `code/behavior/create_valid_behav.py` to filter subjects and create summary CSVs.

2.  **DDM Fitting**:
    *   Compile the C++ model: `code/ddm/simSSP_model_GB_noScale.cpp`.
    *   Run `code/ddm/run_ddm_batch.py` to submit fitting jobs.
    *   Aggregate results with `code/ddm/fitted.py`.

3.  **EEG Preprocessing**:
    *   Use `code/preprocessing-eeg/run_MADE_batch.py` to submit preprocessing jobs (calls `MADE_pipeline.m`).
    *   Verify outputs with `code/preprocessing-eeg/check_preprocessed_files.py`.

4.  **Post-Processing & Statistics**:
    *   Compute ERP means: `code/postprocessing/compute_erp_means.py`.
    *   Compute TF metrics: `code/postprocessing/compute_means_TF.py`.
    *   Run statistical models using R scripts in `code/statistics/`.

## Contributors
*   NDCLab Team
