#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

library(ANCOMBC)
library(biomformat)
library(phyloseq)

tbl_file <- "data/tcga/processed/species/processed_tbl.bacteria.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/tcga/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

tbl <- tbl[, sample_order]

taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
meta <- phyloseq::sample_data(md)
physeq <- phyloseq::phyloseq(taxa, meta)

design.formula <- "investigation + race + gender"

ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design.formula,
                                    zero_cut=1.0)
results_beta <- as.data.frame(ancombc.results$res$beta)
results_qval <- as.data.frame(ancombc.results$res$q_val)

outfile_beta <- "results/tcga/ancombc_results_beta.tsv"
write.table(results_beta, file=outfile_beta, sep="\t")

outfile_qval <- "results/tcga/ancombc_results_qval.tsv"
write.table(results_qval, file=outfile_qval, sep="\t")
