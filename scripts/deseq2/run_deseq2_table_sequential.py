#!/usr/bin/env python
"""Replicate DESeq2 within BIRDMAn"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import biom
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent / "BIRDMAn"))
from birdman.model_base import TableModel

STAN_PATH = str(Path(__file__).parent / "BIRDMAn/birdman/templates/deseq2_table_sequential.stan")


def compute_deseq2_size_factors(table):
    """Median-of-ratios normalization (DESeq2-style)."""
    counts = table.matrix_data.toarray()
    n_samples = counts.shape[1]

    prevalence = np.sum(counts > 0, axis=1) / n_samples
    prevalent = prevalence >= 0.5

    if np.sum(prevalent) < 5:
        lib_sizes = counts.sum(axis=0)
        return lib_sizes / np.median(lib_sizes)

    subset = counts[prevalent, :]
    geo_means = np.zeros(subset.shape[0])
    for i in range(subset.shape[0]):
        nonzero = subset[i, :][subset[i, :] > 0]
        geo_means[i] = np.exp(np.mean(np.log(nonzero))) if len(nonzero) > 0 else np.nan

    valid = ~np.isnan(geo_means) & (geo_means > 0)
    ratios = np.where(subset[valid, :] > 0, subset[valid, :] / geo_means[valid, np.newaxis], np.nan)
    size_factors = np.nanmedian(ratios, axis=0)
    return np.where((size_factors > 0) & np.isfinite(size_factors), size_factors, 1.0)


def estimate_dispersions(table, size_factors):
    """Empirical Bayes dispersion estimation with mean-dispersion trend shrinkage."""
    counts = table.matrix_data.toarray()
    norm_counts = counts / size_factors[np.newaxis, :]
    feature_means = np.mean(norm_counts, axis=1)
    feature_vars = np.var(norm_counts, axis=1, ddof=1)

    # raw method-of-moments dispersions
    raw_dispersions = np.maximum(0.01, (feature_vars - feature_means) / (feature_means**2 + 1e-8))

    # fit mean-dispersion trend in log-log space
    valid = (feature_means > 0.01) & (raw_dispersions > 0)
    coeffs = np.polyfit(np.log(feature_means[valid]), np.log(raw_dispersions[valid]), deg=2)

    trend_dispersions = np.exp(np.polyval(coeffs, np.log(feature_means + 1e-8)))

    # adaptive shrinkage toward trend
    shrinkage_weights = np.clip(0.7 / (1 + feature_means), 0.2, 0.8)
    return (1 - shrinkage_weights) * raw_dispersions + shrinkage_weights * trend_dispersions


class DESeq2SequentialModel(TableModel):
    def __init__(self, table, formula, metadata, size_factors, dispersions,
                 num_iter=1000, num_warmup=500, chains=4, seed=42, beta_prior=5.0):
        super().__init__(table=table, model_path=STAN_PATH,
                         num_iter=num_iter, num_warmup=num_warmup, chains=chains, seed=seed)
        self.create_regression(formula=formula, metadata=metadata)

        ordered_sf = np.array([
            size_factors[list(table.ids(axis="sample")).index(s)]
            for s in self.sample_names
        ])
        fid_to_idx = {fid: i for i, fid in enumerate(table.ids(axis="observation"))}
        ordered_disp = np.array([dispersions[fid_to_idx[f]] for f in self.feature_names])

        self.add_parameters({
            "log_size_factor": np.log(ordered_sf),
            "dispersion": ordered_disp,
            "B_p": beta_prior,
        })
        self.specify_model(
            params=["beta_var"],
            dims={"beta_var": ["covariate", "feature"], "log_lhood": ["tbl_sample", "feature"]},
            coords={"covariate": self.colnames, "feature": self.feature_names, "tbl_sample": self.sample_names},
            include_observed_data=True, log_likelihood="log_lhood",
        )


def main():
    data_dir = Path(__file__).parent / "ravel_deseq2"
    output_dir = "deseq2_table"
    os.makedirs(output_dir, exist_ok=True)

    table = biom.load_table(str(data_dir / "table_filtered_collapsed.biom"))
    metadata = pd.read_csv(data_dir / "sample_metadata.tsv", sep="\t", index_col=0)
    deseq2_results = pd.read_csv(data_dir / "differentials.tsv", sep="\t")
    deseq2_results['feature_id'] = deseq2_results['row'].str.replace('F_', '')
    deseq2_results = deseq2_results.set_index('feature_id')

    # filter to healthy / bv samples
    valid_conditions = ["healthy", "bacterial_vaginosis"]
    filtered_metadata = metadata[metadata["study_condition"].isin(valid_conditions)].copy()
    common_samples = list(set(table.ids(axis="sample")) & set(filtered_metadata.index))
    table = table.filter(common_samples, axis="sample", inplace=False)
    filtered_metadata = filtered_metadata.loc[common_samples]
    print(f"Samples: {len(common_samples)}")

    size_factors = compute_deseq2_size_factors(table)

    # keep only features deseq2 kept (non-na padj)
    deseq2_kept = deseq2_results[deseq2_results['padj'].notna()].index.tolist()
    features_to_analyze = list(set(table.ids(axis="observation")) & set(deseq2_kept))
    table = table.filter(features_to_analyze, axis="observation", inplace=False)
    print(f"Features: {len(features_to_analyze)}")

    # pre-estimate dispersions
    print("Estimating dispersions (EB)...")
    dispersions = estimate_dispersions(table, size_factors)
    print(f"  median: {np.median(dispersions):.3f}, range: [{np.min(dispersions):.3f}, {np.max(dispersions):.3f}]")

    # compile and fit
    print("Compiling stan model...")
    from cmdstanpy import CmdStanModel
    sm = CmdStanModel(stan_file=STAN_PATH)

    model = DESeq2SequentialModel(
        table=table, formula="~ study_condition", metadata=filtered_metadata,
        size_factors=size_factors, dispersions=dispersions,
        num_iter=1000, num_warmup=500, chains=4, seed=42, beta_prior=5.0,
    )
    model.sm = sm

    print("Fitting model...")
    model.fit_model()
    inference = model.to_inference()

    # extract results
    beta_posterior = inference.posterior["beta_var"]
    covariates = list(beta_posterior.coords["covariate"].values)
    features = list(beta_posterior.coords["feature"].values)
    condition_coef = covariates[1] if len(covariates) > 1 else covariates[0]

    fid_to_idx = {fid: i for i, fid in enumerate(table.ids(axis="observation"))}
    results = []
    for fid in features:
        samples = beta_posterior.sel(covariate=condition_coef, feature=fid)
        deseq2_lfc = deseq2_results.loc[fid, 'log2FoldChange'] if fid in deseq2_results.index else np.nan
        base_mean = deseq2_results.loc[fid, 'baseMean'] if fid in deseq2_results.index else np.nan
        results.append({
            "feature_id": fid,
            "birdman_lfc": float(samples.mean()),
            "birdman_lfc_std": float(samples.std()),
            "birdman_lfc_q025": float(samples.quantile(0.025)),
            "birdman_lfc_q975": float(samples.quantile(0.975)),
            "birdman_disp": dispersions[fid_to_idx[fid]],
            "deseq2_lfc": deseq2_lfc,
            "baseMean": base_mean,
        })

    results_df = pd.DataFrame(results).set_index("feature_id")
    results_df = results_df.dropna(subset=["birdman_lfc", "deseq2_lfc"])

    # flip sign if needed to match deseq2 reference level
    r, _ = stats.spearmanr(results_df["birdman_lfc"], results_df["deseq2_lfc"])
    if r < 0:
        print("Flipping signs to match deseq2 reference level")
        results_df["birdman_lfc"] = -results_df["birdman_lfc"]
        results_df["birdman_lfc_q025"], results_df["birdman_lfc_q975"] = (
            -results_df["birdman_lfc_q975"], -results_df["birdman_lfc_q025"])

    spearman_r, _ = stats.spearmanr(results_df["birdman_lfc"], results_df["deseq2_lfc"])
    pearson_r, _ = stats.pearsonr(results_df["birdman_lfc"], results_df["deseq2_lfc"])
    results_df['ci_width'] = results_df['birdman_lfc_q975'] - results_df['birdman_lfc_q025']

    print(f"\nResults: Spearman r = {spearman_r:.4f}, Pearson r = {pearson_r:.4f}")
    print(f"ci widths: mean = {results_df['ci_width'].mean():.2f}, "
          f"range = [{results_df['ci_width'].min():.2f}, {results_df['ci_width'].max():.2f}]")

    results_df.to_csv(os.path.join(output_dir, "comparison_results.tsv"), sep="\t")
    print(f"Saved to {output_dir}/comparison_results.tsv")


if __name__ == "__main__":
    main()
