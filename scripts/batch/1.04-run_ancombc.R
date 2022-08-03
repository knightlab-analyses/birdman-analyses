#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/batch/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(ANCOMBC)
library(biomformat)
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

tbl <- tbl[, sample_order]

taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
meta <- phyloseq::sample_data(md)
physeq <- phyloseq::phyloseq(taxa, meta)

design_formula <- "case_ctrl"

ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design_formula,
                                    zero_cut=1.0)
column_names <- unlist(as.vector(attributes(ancombc.results$res)))
results <- as.data.frame(ancombc.results$res)
colnames(results) <- column_names

outfile <- "results/batch/ancombc_results.tsv"
write.table(results, file=outfile, sep="\t")
