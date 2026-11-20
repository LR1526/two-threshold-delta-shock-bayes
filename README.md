# Bayesian inference for the two-threshold generalized δ-shock model

R code accompanying the manuscript *"Bayesian Inference for Survival Analysis in Two-Threshold Generalized δ-shock Models"*
(L. Rajabi and S. Khazaei), submitted to *Communications in Statistics – Theory and Methods*.

Only system lifetimes are used for inference. The likelihood and the posterior of the failure probability θ are computed **exactly**:
the likelihood of a sample is a finite polynomial in θ and (1−θ), so the posterior under a Beta prior is a finite mixture of Beta
densities (no MCMC, no truncation error).

## Files
| File | Content |
|---|---|
| `01_simulation.R` | Simulation study (tables and figures of the simulation section). |
| `02_nasa.R` | NASA C-MAPSS FD001 application (tables of the application section). |
| `train_FD001.txt` | NASA C-MAPSS FD001 training data (NASA Ames Prognostics Data Repository, Saxena and Goebel, 2008). |
| `LICENSE` | MIT license. |

## Requirements
R ≥ 4.0, base R only (no packages to install).

## How to run
Run from the repository root:
```r
source("01_simulation.R")   # or: Rscript 01_simulation.R
source("02_nasa.R")         # or: Rscript 02_nasa.R
```
* Both scripts use `set.seed(2026)`. Outputs go to `outputs_simulation/` and `outputs_nasa/` (created automatically); re-running overwrites the CSV tables.
* `01_simulation.R` first checks that the lifetime pmf sums to 1 and has the closed-form mean of Poursaeed (2021); it prints `pmf check passed`.
* Run time grows with the replication counts at the top of each script; lower them for a quick test (results then differ from the paper).
* `02_nasa.R` needs `train_FD001.txt` in the same folder.

## Simulation study (`01_simulation.R`)
Settings: p = 0.3, δ1 = 1, δ2 = 3, n = 150, θ ∈ {0.1, 0.3, 0.5, 0.7, 0.9}; prior Beta(1,1) unless stated.

| Output (in `outputs_simulation/`) | Content |
|---|---|
| `table1.csv` | Bias and MSE of the Bayes estimator (posterior mean) and the MLE (R = 500), Monte Carlo standard errors, MSE at n = 200, Cramér–Rao bound. |
| `table2.csv` | Coverage and mean width of 95% intervals for θ (credible, likelihood-ratio, nonparametric bootstrap, weighted-likelihood bootstrap) and rejection rates of the split predictive check for correct and misspecified data (50 replications). |
| `table3.csv` | Prior sensitivity: MSE of the posterior mean under six Beta priors and of the MLE (50 replications). |
| `table4_survival_paper.csv`, `table4_survival_full.csv` | Survival function S(t) at t = 5, 10, 20: true value, Bayes (posterior mean of S) and MLE (plug-in) estimates, MSE, credible and LR intervals, coverage (50 replications). The paper reports t = 10. |
| `figure1_estimates_*.png`, `figure2_running_mean_*.png` | Estimates of θ over replications and their running means (true θ = 0.5). |
| `figure3_ET_vs_p_*.png`, `figure_survival_*.png` | E(T) as a function of p; survival curve S(t). |

Not all outputs are reported in the paper. Figure file names contain a time stamp, so earlier files are not overwritten.

## NASA application (`02_nasa.R`)
Each engine is mapped to the model through the standardised residual of sensor 12 (cubic trend removed): a cycle with |z| > 1 is a shock; a shock after an inter-shock time ≤ δ1 = 1 is fatal; a shock after a time in (1, 3] is fatal if |z| > 1.5 (so θ is the probability that such a shock is severe). The lifetime of an engine is the cycle of its first fatal shock; only these 100 lifetimes are used for inference. The shock probability p is set to the observed shock rate and treated as known.

| Output (in `outputs_nasa/`) | Content |
|---|---|
| `table4_theta.csv` | Estimates and 95% intervals of θ (Bayes and MLE), 100 and 20 engines. |
| `table5_survival.csv` | Estimates of S(t) at t = 5, 10, 20 and absolute errors against the reference survival function. |
| `table6_gof.csv` | Goodness of fit: Pearson chi-square and split predictive check (200 random splits). |
| `table6b_negative_control.csv`, `table7_mapping_sensitivity.csv`, `table8_subsamples.csv` | Additional checks (not in the paper): negative control, other sensors/thresholds in the mapping, repeated random 20-engine subsamples. |

The file numbering in `outputs_nasa/` follows the script and differs from the table numbers in the paper.

## Citation
Please cite the paper (reference to be added after publication) and this repository (Zenodo DOI: https://doi.org/10.5281/zenodo.24301196). R: R Core Team, *R: A Language and Environment for Statistical Computing*, R Foundation for Statistical Computing, Vienna, Austria, https://www.R-project.org/.

## License
MIT (see `LICENSE`).
