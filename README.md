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
| `mc_rk4_pdf_original.f90` | Perform Monte Carlo simulations of the original stochastic system with Poincaré sampling and output amplitude–phase probability density statistics. |

## data

[data](data/) contains the article's CSV data, including parameter provenance records and FPK/MC curve comparison metrics.
