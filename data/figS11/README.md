# Group statistics for the empirical alpha-BOLD correlation map

Two-tailed one-sample sign-flip permutation test of the mean alpha-BOLD
correlation against zero, Benjamini-Hochberg FDR correction across voxels,
and a BrainNet Viewer rendering of the surviving voxels.

Produces the FDR-thresholded map reported in the manuscript.

## Run

```matlab
addpath('<path-to>/fieldtrip');            % ft_defaults is called by the script
addpath(genpath('<path-to>/BrainNetViewer'));
alpha_bold_group_stats
```

All data paths inside the script are relative to the script's own location,
so it can be run from anywhere.

## Dependencies

| | |
|---|---|
| [FieldTrip](https://www.fieldtriptoolbox.org) | `ft_defaults`, `ft_sourceinterpolate`, `ft_volumewrite` |
| [BrainNet Viewer](https://www.nitrc.org/projects/bnv) | `BrainNet_MapCfg` |

MATLAB with the Statistics and Machine Learning Toolbox (`pdist2`).
Tested on R2024a.

## Bundled data (`data/`)

| file | contents |
|---|---|
| `AlphaBOLDCorr.mat` | `AlphaBOLDCorr` — 4417 voxels x 72 subjects, Pearson r; plus `inside`, `pos`, `unit` describing the source grid |
| `mni_sourcemodel.mat` | 6 mm grey matter source grid (cm) |
| `mni_mri.mat` | MNI template, interpolation target |
| `mycmap.mat` | colour map |
| `BNopt1.mat` | BrainNet Viewer option struct |

## Method

At each of the 4417 grey matter voxels, the null distribution of the group
mean correlation is built from 10<sup>5</sup> random sign flips of the
individual subject maps. The two-tailed p-value is the proportion of
permutations whose absolute mean equals or exceeds the observed absolute
mean. P-values are then corrected across voxels with the
Benjamini-Hochberg procedure (Benjamini & Hochberg, 1995).

## Output

Written next to the script:

| file | contents |
|---|---|
| `alpha_bold_group_stats.mat` | `meanR`, `p`, `qval`, `mask` |
| `alpha_bold_fdr.nii` | thresholded map |
| `alpha_bold_fdr_q0p05.png` | rendered figure |

At q < 0.05, 2275 of 4417 voxels are significant (2087 positive, 188
negative).

Runtime is about 2 minutes for 10<sup>5</sup> permutations. `nPerm`, `qFDR`
and the random `seed` are set at the top of the script.

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
a practical and powerful approach to multiple testing. *Journal of the Royal
Statistical Society: Series B*, 57(1), 289–300.

Xia, M., Wang, J., & He, Y. (2013). BrainNet Viewer: a network visualization
tool for human brain connectomics. *PLoS ONE*, 8(7), e68910.
