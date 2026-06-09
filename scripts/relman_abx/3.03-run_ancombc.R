#!/home/lpatel/software/miniconda3/envs/ancombc2/bin/Rscript
#SBATCH --chdir=/ddn_scratch/lpatel/projects/2024-07-17_q2-birdman/birdman-analyses-b2
#SBATCH --output=/ddn_scratch/lpatel/projects/2024-07-17_q2-birdman/birdman-analyses-b2/slurm_out/relman_abx/%x.out
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

rownames(beta_res) <- taxon_names
rownames(q_res) <- taxon_names

outfile_beta <- "results/relman_abx/ancombc_results_beta.tsv"
write.table(beta_res, file=outfile_beta, sep="\t")

outfile_qval <- "results/relman_abx/ancombc_results_qval.tsv"
write.table(q_res, file=outfile_qval, sep="\t")
