# BIRDMAn Manuscript Analyses

Analysis code, processed data, and notebooks accompanying the **BIRDMAn** (*Bayesian Inferential Regression for Differential Microbiome Analysis*) manuscript.

## About

[BIRDMAn](https://github.com/biocore/birdman) is a framework for fitting Bayesian differential abundance models to microbiome count data, built on top of [CmdStanPy](https://github.com/stan-dev/cmdstanpy). This repository contains the end-to-end analyses used to evaluate BIRDMAn in the manuscript, including the scripts used to fit models on a SLURM cluster, the processed input data, and the Jupyter notebooks that generate the manuscript figures. BIRDMAn itself is available at [`biocore/birdman`](https://github.com/biocore/birdman).

## Installation

Runtime dependencies include [`birdman`](https://github.com/biocore/birdman), [`cmdstanpy`](https://github.com/stan-dev/cmdstanpy), [`biom-format`](https://github.com/biocore/biom-format), [`arviz`](https://github.com/arviz-devs/arviz), `pandas`, `numpy`, and `matplotlib`. We suggest installing `cmdstanpy` from conda-forge.


## Citation

If you use this code or the processed data, please cite the [BIRDMAn manuscript](https://www.biorxiv.org/content/10.1101/2023.01.30.526328v1):

## Maintainers

- **Gibraan Rahman** — original author (`gibsramen@gmail.com`)
- **Lucas Patel** — co-maintainer (`lpatel@ucsd.edu`)
- **Yang Chen** — co-maintainer (`yac027@ucsd.edu`)

## License

- **Code** (`src/`, `scripts/`, `notebooks/`): BSD 3-Clause — see [`LICENSE`](LICENSE).
- **Processed data** (`data/`): CC0 1.0 Universal Public Domain Dedication — see [`data/LICENSE`](data/LICENSE).
