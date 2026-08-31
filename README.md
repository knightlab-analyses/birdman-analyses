# BIRDMAn Manuscript Analyses

Analysis code and processed data behind Figure 2 of the **BIRDMAn** (*Bayesian
Inferential Regression for Differential Microbiome Analysis*) manuscript. The
package itself is at [`biocore/birdman`](https://github.com/biocore/birdman).

## Figure 2

```bash
python scripts/figure2.py        # -> figures/figure2.{pdf,png}
```

Runs from a fresh clone against committed data — no cluster, no model refit.

- **a** — vaginal microbiome (Ravel *et al.* 2011); posterior differentials with
  95% HDI, from `data/qadabra/outputs/`
- **b** — dual-course ciprofloxacin (Dethlefsen & Relman 2011); log-ratio
  trajectories and posterior derivatives, from `results/relman_abx/`

## Refitting the models

The tables in `results/relman_abx/` are the original manuscript fits. To
regenerate them, run `src/relman_abx/run_birdman_chunked.py`, then the numbered
scripts in `scripts/relman_abx/` in order — each documents its own inputs and
outputs. Expect ~2.3 GB of per-feature posteriors and a close but inexact match:
every feature is an independent MCMC run, so posterior means shift with the seed
(pooled Spearman ~0.97 against the committed tables).

## Dependencies

`birdman`, `cmdstanpy`, `biom-format`, `arviz`, `pandas`, `numpy`, `matplotlib`;
install `cmdstanpy` from conda-forge. There is nothing to install from this repo
— `src/` is imported directly, so run from the repo root or add it to the path:

```bash
export PYTHONPATH=$PWD:$PYTHONPATH
```

## Citation

Rahman G, Patel L, Chen Y, *et al.* BIRDMAn: a flexible Bayesian framework for
differential abundance analysis of microbiome data.
[bioRxiv 2023.01.30.526328](https://www.biorxiv.org/content/10.1101/2023.01.30.526328v1)

## Maintainers

- **Gibraan Rahman** — original author (`gibsramen@gmail.com`)
- **Lucas Patel** — co-maintainer (`lpatel@ucsd.edu`)
- **Yang Chen** — co-maintainer (`yac027@ucsd.edu`)

## License

- **Code** (`src/`, `scripts/`, `notebooks/`): BSD 3-Clause — see [`LICENSE`](LICENSE).
- **Processed data** (`data/`): CC0 1.0 Universal Public Domain Dedication.
