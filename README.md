# BIRDMAn Manuscript Analyses

Data and code for Figure 2 of the **BIRDMAn** manuscript. The package itself is
at [`biocore/birdman`](https://github.com/biocore/birdman).

## Figure 2

```bash
python scripts/figure2.py        # -> figures/figure2.pdf
```

Runs from a fresh clone against committed data — no cluster, no refit. Needs
`birdman`, `cmdstanpy`, `biom-format`, `arviz`, `pandas`, `numpy`, `matplotlib`
(install `cmdstanpy` from conda-forge). Nothing to install from this repo;
`src/` is imported directly, so run from the repo root.

## Refitting

`results/relman_abx/` holds the original manuscript fits. To regenerate them run
`src/relman_abx/run_birdman_chunked.py`, then the numbered scripts in
`scripts/relman_abx/` in order — each documents its own inputs and outputs.
Refits are close but inexact: every feature is an independent MCMC run, so
posterior means shift with the seed (pooled Spearman ~0.97).

## Citation

Rahman G, Patel L, Chen Y, *et al.* BIRDMAn: a flexible Bayesian framework for
differential abundance analysis of microbiome data.
[bioRxiv 2023.01.30.526328](https://www.biorxiv.org/content/10.1101/2023.01.30.526328v1)

## License and maintainers

Code BSD-3-Clause (see [`LICENSE`](LICENSE)); processed data CC0 1.0.
Gibraan Rahman (`gibsramen@gmail.com`), Lucas Patel (`lpatel@ucsd.edu`),
Yang Chen (`yac027@ucsd.edu`).
