# Spatiotemporal_PCA_Data_Driven
An Attempt at a data-driven decomption of stpatiotemporal brain waves using EEG

# Spatiotemporal PCA for EEG (MATLAB + EEGLAB)

This repository provides a **single MATLAB wrapper script** for running pooled **spatiotemporal principal component analysis (PCA)** on EEG datasets (`.set` files).  
It is designed for reproducible, publication-ready analysis of event-related potentials (ERPs) with **temporal × spatial decomposition**, using [EEGLAB](https://sccn.ucsd.edu/eeglab/).

---

## Background

Event-related potentials (ERPs) are multichannel, time-locked neural signals. Standard analyses often average over electrodes or time windows, but this risks missing distributed or temporally overlapping components.  

**Spatiotemporal PCA** addresses this by:
- Performing **temporal PCA** across all channels × subjects to extract latent time factors.  
- Performing **spatial PCA** within each temporal factor to localize scalp distributions.  
- Using **parallel analysis** to decide the optimal number of factors.  
- Allowing multiple **rotation methods** (e.g., Varimax, Promax) for interpretability.  

This implementation provides:
- Automated preprocessing (epoching, baseline correction).  
- Pooled analysis across conditions.  
- Robust reproducibility (logs, factor cutoffs, scree plots).  
- Detailed figures of temporal loadings and spatial topoplots.

---

## Scientific Rationale: Data-Driven Decomposition

Traditional ERP analyses often rely on **a priori electrode/time window choices** (e.g., “take mean amplitude at Fz between 300–500 ms”). While simple, this risks:
- Ignoring distributed components with non-canonical topographies.
- Overlooking overlapping or latent processes in time.
- Imposing assumptions that bias interpretation.

This toolbox implements a **data-driven decomposition strategy** via **spatiotemporal PCA**:

1. **Temporal PCA**  
   - Collapses across subjects, electrodes, and trials to extract latent time courses (factors).  
   - Factors are retained using **parallel analysis**, a permutation-based, data-driven criterion.  

2. **Spatial PCA**  
   - Within each temporal factor, electrode covariance patterns are decomposed into spatial factors.  
   - Provides interpretable scalp maps aligned with each temporal mode.  

3. **Rotation Methods**  
   - Multiple rotations (`Varimax`, `Promax`, etc.) are supported to improve factor interpretability.  
   - Rotation does not change the variance explained, but redistributes loadings to yield simpler, more meaningful factors.  

Together, this provides a **data-driven alternative** to “hand-drawn” ERP windows, allowing the data itself to reveal **when** and **where** the brain organizes information.  

---

## What Makes This Pipeline Different
- **Data-driven retention**: Factor cutoffs are not arbitrary — they are determined by real vs. permuted eigenvalue distributions.  
- **Reproducibility**: All factor choices are logged (`parallel_analysis_report.txt`).  
- **Exploratory + Confirmatory**: Can be used both to discover latent factors and to confirm known ERP effects (e.g., N2, P3).  
- **Condition-Pooled**: Designed to load and jointly decompose across multiple groups (e.g., Placebo vs THC), ensuring consistent factors across datasets.  

This decomposition aligns with modern calls for **precision EEG/ERP**: letting statistical structure, not subjective ROI choices, guide the discovery of meaningful components.


---
## Requirements

- **MATLAB R2021a+**  
- [**EEGLAB**](https://github.com/sccn/eeglab) (tested with 2024)  
  Required functions: `pop_loadset`, `pop_epoch`, `pop_rmbase`, `topoplot`.  

⚠️ The script automatically removes conflicting PCA plugins (`AAR131130`) from EEGLAB’s path.

---

## Usage

### Basic Example
```matlab
run_spatiotemporal_pca_wrapper( ...
    'example_data/', ...        % folder with .set files
    'results', ...              % output directory
    {'S  6','S  7'}, ...        % event codes
    [-0.2 0.7], ...             % epoch window (seconds)
    [-200 0], ...               % baseline correction (ms)
    [300 500], ...              % ERP window (ms)
    500, ...                    % permutations for parallel analysis
    {{'promax','Power',4}});    % rotation method


| Argument          | Type        | Default                  | Description                           |
| ----------------- | ----------- | ------------------------ | ------------------------------------- |
| `folder`          | `string`    | *required*               | Path to `.set` files.     			 |
| `outputDir`       | `string`    | *required*               | Directory where results are saved.    |
| `events`          | `cell`      | — (must be set)          | Event codes to epoch around.          |
| `epochWindow`     | `[min max]` | `[-0.2 0.7]`             | Epoch window in seconds.              |
| `BaselineCorrect` | `[min max]` | `[-200 0]`               | Baseline correction in ms.            |
| `erpWindow`       | `[min max]` | `[300 500]`              | ERP search window in ms.              |
| `nPermPA`         | `int`       | `500`                    | # permutations for parallel analysis. |
| `rotationMethods` | `cell`      | `{{'promax','Power',4}}` | Factor rotation(s).                   |


📈 Outputs

All outputs are written to the chosen outputDir:
Reports
	missing_epochs_all.txt — list of files without valid events.
	mismatch_length_all.txt — list of files excluded due to epoch length mismatch.
	parallel_analysis_report.txt — temporal & spatial factor retention criteria.

Figures
	Temporal PCA loadings (raw + SPSS-scaled).
	Spatial PCA topoplots (raw + scaled).
	Scree plots for temporal & spatial PCA.
	Combined summary (variance explained, cutoffs, thresholds).

Each figure is labeled with:
	Temporal factor ID (TF#)
	Rotation method (VARIMAX, PROMAX, etc.)
	Variance explained (%) and peak latency.

Methods Notes
	Temporal PCA: performed across subjects × channels × timepoints.
	Parallel analysis: determines component number (real eigenvalues vs random).
	Spatial PCA: applied within each temporal factor.
	Rotations: EEGLAB’s rotatefactors is used with user-defined methods.
	Baseline correction: applied in ms relative to event onset.

Citation

If you use this code in publications, please cite:
Gunn, M. (2025). Spatiotemporal PCA Wrapper for EEG in MATLAB (EEGLAB). GitHub Repository.

Contact
Maintained by Matthew Gunn, Ph.D.
For questions, suggestions, or bug reports: please open a GitHub issue.