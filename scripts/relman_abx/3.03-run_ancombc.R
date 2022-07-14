#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(ANCOMBC)
library(biomformat)
library(phyloseq)

tbl_file <- "data/relman_abx/processed/processed_tbl.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/relman_abx/processed/processed_md.tsv"
md <- read.table(md_file, sep="\t", row.names=1, header=T)
md$antibiotic <- relevel(as.factor(md$antibiotic), "preCp")
md$host_subject_id <- relevel(as.factor(md$host_subject_id), "494.D")

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)
tbl <- tbl[, sample_order]

taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
meta <- phyloseq::sample_data(md)
physeq <- phyloseq::phyloseq(taxa, meta)

design_formula <- "antibiotic + host_subject_id"

ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design_formula,
                                    zero_cut=1.0)
results_beta <- as.data.frame(ancombc.results$res$beta)

outfile_beta <- "results/relman_abx/ancombc_results_beta.tsv"
write.table(results_beta, file=outfile_beta, sep="\t")
