# Stationary stochastic P-bifurcations in a toe-like nonlinear isolator

Research code associated with the article:

**Stationary Stochastic P-Bifurcations and Damping-Regulated Probability-Modal Structures in a Bio-Inspired Toe-Like Nonlinear Isolator**

## Programs

| File | Purpose |
| --- | --- |
| `fpk_omega_probability_modal_map.m` | Scan excitation frequency and identify stationary amplitude-density modes and candidate transition intervals. |
| `fpk_zhat_probability_modal_map.m` | Scan excitation amplitude and identify stationary amplitude-density modes and candidate transition intervals. |
| `fpk_selected_omega_marginals_fig6_fixed.m` | Compute and plot FPK marginals at selected excitation frequencies (Fig. 6). |
| `fpk_selected_seta_marginals_fig9.m` | Compute and plot FPK marginals at selected noise intensities (Fig. 9). |
| `fpk_selected_Zhat_marginals_fig10.m` | Compute and plot FPK marginals at selected excitation amplitudes (Fig. 10). |
| `plot_fig11_fig13_separate_v2.m` | Generate deterministic response branches and stability maps for Figs. 11 and 13. |
| `compare_FPK_MC_marginals_fixed.m` | Compare saved FPK and MC amplitude/phase marginals and compute discrepancy metrics. |
| `plot_FPK_MC_main_figure.m` | Generate FPK/MC marginal comparisons and joint-density contours. |
| `mc_rk4_pdf_original.vfproj` | Intel Fortran / Visual Studio project referencing the original-system MC source. |
| `程序名及作用说明_MC_FPK.docx` | Original bilingual program descriptions supplied by the author. |

## Requirements and usage

Use MATLAB for the `.m` files. A specific MATLAB release and toolbox requirements have not been validated. Set MATLAB's current folder to this directory and inspect each function's settings before running it.

Examples using the actual supplied function names:

```matlab
results = fpk_selected_omega_marginals_fig6_fixed([0.1534,0.1543,0.1545,0.1552],5e-3);
results = fpk_omega_probability_modal_map([],5e-3);
results = fpk_zhat_probability_modal_map([],5e-3);
results = fpk_selected_seta_marginals_fig9([0.0001,0.0005,0.001,0.005,0.01]);
results = fpk_selected_Zhat_marginals_fig10([0.0005,0.0008,0.0010,0.0012,0.0015]);
results = plot_fig11_fig13_separate_v2();
```

These are separate analysis tasks; check each program's output settings and avoid overwriting results from earlier runs.

The comparison functions require `selected_Omega_stationary_FPK.mat` and `MC_PDF_aphi.txt` in the MATLAB current folder. The main comparison plot expects `results.Omega`, `results.a`, `results.phi`, and a joint-density field (preferably `results.P`) in the MAT file. Verify parameter and phase conventions match before comparison.

## Completeness and numerical interpretation

- The supplied collection does **not** include `mc_rk4_pdf_original.f90`, which the Fortran project references, or its referenced `ReadMe.txt`. The project file alone cannot reproduce MC simulations.
- The original Word description also mentions `fpk_selected_omega_marginals_v2.m`; that filename is not included in this collection. Use the supplied files and their actual function names.
- Precomputed MC/FPK comparison input data are not included.
- Several FPK solvers default to a weak reset regularization intended for qualitative mode-count and mode-location inspection. As documented in their source comments, relative modal weights, means, and variances from quick mode require further validation before publication use. Validate representative transition points with strict FPK and/or original-system MC.
- Programs have been inspected for documentation purposes but have not been executed or numerically validated as part of repository preparation.

## Citation

Please cite the article named above when using this work. Author list, journal, publication year, and DOI should be added when confirmed.

## License

No license has been selected by the author for this collection.

