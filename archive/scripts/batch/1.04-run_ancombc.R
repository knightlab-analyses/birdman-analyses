#!/home/lpatel/software/miniconda3/envs/ancombc2/bin/Rscript
#SBATCH --chdir=/ddn_scratch/lpatel/projects/2024-07-17_q2-birdman/birdman-analyses-b2
#SBATCH --output=/ddn_scratch/lpatel/projects/2024-07-17_q2-birdman/birdman-analyses-b2/slurm_out/batch/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(ANCOMBC)
library(biomformat)
library(jsonlite)
library(phyloseq)

tbl_file <- "results/batch/sim/sim_counts.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "results/batch/sim/metadata.tsv"
md <- read.table(md_file, sep="\t", row.names=1, header=T)
md$case_ctrl <- relevel(as.factor(md$case_ctrl), "0")

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

params <- jsonlite::read_json("results/batch/sim/params.json")
log_depths <- unlist(params$log_depths)
md$log_depths <- log_depths

tbl <- tbl[, sample_order]

taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
meta <- phyloseq::sample_data(md)
physeq <- phyloseq::phyloseq(taxa, meta)

design_formula <- "case_ctrl + log_depths"

ancombc.results <- ANCOMBC::ancombc2(data = physeq, 
                                     fix_formula = design_formula,
                                     prv_cut = 0.0)

taxon_names <- ancombc.results$res$taxon

lfc_cols <- grep("^lfc_", names(ancombc.results$res), value = TRUE)
beta_res <- as.data.frame(ancombc.results$res[, lfc_cols])
colnames(beta_res) <- paste0(colnames(beta_res), "_beta")

q_cols <- grep("^q_", names(ancombc.results$res), value = TRUE)
q_res <- as.data.frame(ancombc.results$res[, q_cols])
colnames(q_res) <- paste0(colnames(q_res), "_qval")

results <- cbind(beta_res, q_res)

rownames(results) <- taxon_names

outfile <- "results/batch/ancombc_results.tsv"
write.table(results, file=outfile, sep="\t", row.names=TRUE)
