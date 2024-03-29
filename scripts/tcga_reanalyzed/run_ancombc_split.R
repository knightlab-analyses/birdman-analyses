#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

library(ANCOMBC)
library(biomformat)
library(phyloseq)

md_file <- "data/tcga_reanalyzed/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)

split_dir <- "data/tcga_reanalyzed/processed"
splits <- c(0:4)

design.formula <- "investigation + race + gender + data_submitting_center_label"

for (split in splits) {
    print(split)

    tbl_file <- paste0(
        "data/tcga_reanalyzed/processed/split_", split,
        "/merged_tbl.train.", split, ".biom"
    )
    tbl <- biomformat::read_biom(tbl_file)
    tbl <- as.matrix(biomformat::biom_data(tbl))

    samples <- colnames(tbl)
    subset_md <- subset(md, rownames(md) %in% samples)
    sample_order <- row.names(subset_md)
    tbl <- tbl[, sample_order]

    taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
    meta <- phyloseq::sample_data(subset_md)
    physeq <- phyloseq::phyloseq(taxa, meta)

    ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design.formula,
                                        zero_cut=1.0)
    results_beta <- as.data.frame(ancombc.results$res$beta)
    results_qval <- as.data.frame(ancombc.results$res$q_val)

    outdir <- paste0(
        "results/tcga_reanalyzed/ancombc_split/"
    )
    dir.create(outdir, showWarnings=FALSE)
    outfile_beta <- paste0(
        outdir, "split_", split, ".results.beta.tsv"
    )
    write.table(results_beta, file=outfile_beta, sep="\t")

    outfile_qval <- paste0(
        outdir, "split_", split, ".results.qval.tsv"
    )
    write.table(results_qval, file=outfile_qval, sep="\t")
}
