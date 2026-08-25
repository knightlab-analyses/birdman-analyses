# BIRDMAn Manuscript Analyses

Analysis code, processed data, and notebooks accompanying the **BIRDMAn** (*Bayesian Inferential Regression for Differential Microbiome Analysis*) manuscript.

## About

[BIRDMAn](https://github.com/biocore/birdman) is a framework for fitting Bayesian differential abundance models to microbiome count data, built on top of [CmdStanPy](https://github.com/stan-dev/cmdstanpy). This repository contains the end-to-end analyses used to evaluate BIRDMAn in the manuscript, including the scripts used to fit models on a SLURM cluster, the processed input data, and the Jupyter notebooks that generate the manuscript figures. BIRDMAn itself is available at [`biocore/birdman`](https://github.com/biocore/birdman).

## Installation

Runtime dependencies include [`birdman`](https://github.com/biocore/birdman), [`cmdstanpy`](https://github.com/stan-dev/cmdstanpy), [`biom-format`](https://github.com/biocore/biom-format), [`arviz`](https://github.com/arviz-devs/arviz), `pandas`, `numpy`, and `matplotlib`. We suggest installing `cmdstanpy` from conda-forge.

There is nothing to install from this repo — `src/` is imported directly. Run the
scripts from the repo root, or put it on the path:

```bash
export PYTHONPATH=/path/to/birdman-analyses:$PYTHONPATH
```

## Reproducing Figure 2b (`relman_abx`)

`notebooks/relman_abx_figures.ipynb` runs from a fresh clone against the committed
summary tables — no cluster and no model refit required.

The committed tables under `results/relman_abx/` are the original model outputs used
for the manuscript figure. Refitting reproduces them closely but not exactly: each
feature is an independent MCMC run, so posterior means shift slightly with the seed
(pooled Spearman ~0.97 against the originals, ~30-37 of the top/bottom 40 features
per contrast unchanged).

To regenerate those tables from scratch, run the numbered scripts in order. All of
them resolve paths relative to the repo root; set `BIRDMAN_INFERENCE_DIR` to wherever
you want the per-feature posteriors written (~2.3 GB for 822 features).

```bash
export BIRDMAN_INFERENCE_DIR=/path/to/inferences

# fit one model per feature (SLURM array on a cluster, or --num-chunks 1 locally)
python src/relman_abx/run_birdman_chunked.py \
    --inference-dir $BIRDMAN_INFERENCE_DIR \
    --num-chunks 1 --chunk-num 1 \
    --chains 4 --num-iter 500 --num-warmup 1000 \
    --beta-prior 10.0 --disp-scale 0.5 --re-prior 3.0 \
    --logfile results/relman_abx/log/run.log

python scripts/relman_abx/2.01-summarize_birdman.py       # birdman_results.beta_var.tsv
python scripts/relman_abx/2.02-get_subj_effects.py        # birdman_results.subj.tsv
python scripts/relman_abx/2.03-concat_infs.py             # beta_var.nc (~79 MB, not committed)
python scripts/relman_abx/2.04-get_derivative_quantiles.py  # {fcp,scp}_derivative_quantiles.tsv
```

`beta_var.nc` is deliberately not committed — it is too large for plain git, and the
figure only consumes 6 covariates x 3 quantiles from it. `2.04` collapses it to the two
small tables the notebook reads, so step 3 is only needed if you are regenerating them.

`data/relman_abx/ref/97_otu_taxonomy.txt` is the Greengenes 13_8 97% OTU taxonomy,
subset to the 822 features in the processed table.

## Citation

If you use this code or the processed data, please cite the [BIRDMAn manuscript](https://www.biorxiv.org/content/10.1101/2023.01.30.526328v1):

## Maintainers

- **Gibraan Rahman** — original author (`gibsramen@gmail.com`)
- **Lucas Patel** — co-maintainer (`lpatel@ucsd.edu`)
- **Yang Chen** — co-maintainer (`yac027@ucsd.edu`)

## License

- **Code** (`src/`, `scripts/`, `notebooks/`): BSD 3-Clause — see [`LICENSE`](LICENSE).
- **Processed data** (`data/`): CC0 1.0 Universal Public Domain Dedication.
