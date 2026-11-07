# Spatiotemporal PCA for ERP/EEG

A standalone MATLAB workflow using EEGLAB datasets. This repository does not
install an EEGLAB menu. It pools each dataset's selected-event ERP, performs
temporal PCA, and then spatial PCA on the temporal scores.

## Setup and use

Tested with MATLAB R2023b, Statistics and Machine Learning Toolbox, and EEGLAB
2024.2. Initialize EEGLAB, then add this repository's root to the path:

```matlab
addpath('/path/to/eeglab'); eeglab('nogui');
addpath('/path/to/Spatiotemporal_PCA_Data_Driven');
rng(42); % record this seed with the analysis
run_spatiotemporal_pca_wrapper( ...
    '/data/group1', '/data/group2', '/results/new_analysis', ...
    {'S  6','S  7'}, [-0.2 0.7], [-200 0], [300 500], ...
    500, {{'varimax'},{'promax','Power',4}});
```

The second folder is required positionally; use `''` for one group. At least
three usable datasets are required. Use continuous `.set` inputs containing
matching events and channel locations. Each input contributes one averaged
ERP over all selected event types; conditions are not analyzed separately.
Input time grids and channel labels/order must match. Incompatible files are
listed in `mismatch_length_all.txt`. Files with no matching events are listed
in `missing_epochs_all.txt`.

MATLAB's Statistics PCA is placed first on the path for this call to avoid
EEGLAB plugin name conflicts. The original path is restored afterward.

## Outputs and interpretation

A new output folder is recommended for each run. It receives:

- `parallel_analysis_report.txt`: factor counts and permutation settings.
- `analysis_metadata.mat`: input folders, groups and analysis options.
- `pca_<rotation>_TF<number>.mat`: loadings, scores, eigenvalues and locations.
- `.fig` and `.png` figures: scree plots, temporal curves and scalp maps.

Parallel analysis independently permutes each variable's observations,
retaining the original marginal scale. Temporal and global spatial counts
use real eigenvalues exceeding the mean permuted eigenvalues. The spatial
count is capped by the subject/channel rank and reused for each temporal
factor; it is **not** a separate per-factor retention analysis. No factors
are forced when parallel analysis retains zero.

Single-factor rotations are identity transforms. Multi-factor scores use the
inverse transpose of the loading rotation, preserving reconstruction for
orthogonal and oblique rotations. The displayed scaled curves are
covariance-weighted loadings, not a claim of exact SPSS equivalence. Rotated
score-variance percentages describe the retained scores; with oblique
rotations factors overlap and these are not additive explained-variance
contributions. Spatial variance labels refer to unrotated eigenvalues.

This remains an exploratory analysis workflow. The smoke tests verify file
handling, matching, decomposition, output and reconstruction, not clinical
validity or confirmatory inference. Temporal autocorrelation and the chosen
permutation null should be reviewed for each study.

## Verification

The readiness audit processed eight synthetic EEG datasets, rejected a
channel-order mismatch, exercised one- and two-factor solutions with Varimax
and Promax, checked reconstruction and path restoration, and saved results.
See `tests/smoke_spatiotemporal_pca.m` for a reproducible synthetic test.

## Downloadable synthetic example

[Run the bundled synthetic example](examples/synthetic/README.md) to check your setup against
included reference outputs. It contains actual input files, a one-command runner,
expected tables/plots, and numerical comparisons that report PASS or a mismatch.
The example guide links to the main-branch ZIP, including the readiness fixes.
