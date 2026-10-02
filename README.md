# Code for "Laser ablation rapid evaporative ionisation mass spectrometry imaging with few-micrometer sampling pixels"

Data: Zenodo, https://doi.org/10.5281/zenodo.21981110. Licence: MIT.

## hdi_txt_to_imzml.py

Converts a Waters HDImaging (Maldichrom) peak-picked `.txt` export into a centroid imzML + `.ibd` pair.

```
python hdi_txt_to_imzml.py INPUT.txt OUTPUT_STEM --width W --pixel-size P --polarity negative [--scan-direction bottom_up]
```

The mouse brain acquisitions (Fig. 3, `2025_05_21_MB_Res_*`) scan bottom-up and need
`--scan-direction bottom_up`. The DRG (Fig. 2) and FFPE (Fig. 4) acquisitions use the default.
Raster widths: Fig. 2 400; Fig. 3 at 2, 3, 5, 10, 20 um: 550, 367, 220, 110, 55; Fig. 4 at 20 um 425, at 30 um 283.

Tested with Python 3.13.5, numpy 2.1.3, pyimzml 1.5.5.

## edge_fitting.py

Step-edge spatial-resolution fit (Supplementary Fig. 7). An error-function edge model is fitted
to an intensity profile across an edge, and the 16-84% rise distance is computed from its sigma.

```python
from edge_fitting import fit_edge, resolution_16_84
popt, perr = fit_edge(distance_um, intensity)   # popt = low, high, x0, sigma
width_um = resolution_16_84(popt[3])
```

The edge profiles (distance in um, intensity in counts, one block per edge line) are in
`s_fig7.zip` of the data record. Fitting the intensity columns reproduces the published widths.

Tested with Python 3.13.5, numpy 2.1.3, SciPy 1.15.3.

## ffpe_segmentation.R

Segmentation of the FFPE oesophagus section (Fig. 4h,i). Written by Duncan Roberts.

```
Rscript ffpe_segmentation.R <2024_05_10_Human_MS19_20um_10Hz.raw folder> <2024_05_10_Human_MS19_20um_10Hz_1000peaks.txt> [cluster_plot.pdf]
```

Both inputs are in `fig4.zip` of the data record. The `.raw` folder is read only for the raster
size (its `_extern.inf`, 425 x 400 pixels).

Run with R 4.6.1 (Windows 11 x64). The UMAP and clustering steps are sensitive to the R and
package versions, so use these:

| Package | Version | Source |
|---|---|---|
| rawsewage | 0.0.0.9000 | GitHub `drobertsicl/rawsewage`, commit `99513cc5f6036892cc5e92ffd87a13baa34a95c7` |
| irlba | 2.3.5.1 | CRAN |
| genieclust | 1.1.5-2 | CRAN |
| vegan | 2.7-1 | CRAN (the script calls the internal `vegan:::.calc_clr`) |
| uwot | 0.2.4 | CRAN |
| lumbermark | 0.9.0 | CRAN |
| deadwood | 0.9.2 | CRAN |
| quitefastmst | 0.9.2 | CRAN |
| imager | 1.0.8 | CRAN |
| ggplot2 | 4.0.3 | CRAN |
| dplyr | 1.2.1 | CRAN |
| scales | 1.4.0 | CRAN |
| Matrix | 1.7-5 | CRAN |

```r
install.packages("remotes")
remotes::install_github("drobertsicl/rawsewage@99513cc5f6036892cc5e92ffd87a13baa34a95c7")
remotes::install_version("irlba", "2.3.5.1")
remotes::install_version("genieclust", "1.1.5-2")
remotes::install_version("vegan", "2.7-1")
remotes::install_version("uwot", "0.2.4")
remotes::install_version("lumbermark", "0.9.0")
```
