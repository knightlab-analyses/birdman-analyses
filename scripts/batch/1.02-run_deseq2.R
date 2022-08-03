#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/batch/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(biomformat)
library(jsonlite)
library(DESeq2)

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

design_formula <- as.formula("~case_ctrl + log_depths")

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=tbl,
    colData=md,
    design=design_formula
)
dds_results <- DESeq2::DESeq(dds, sfType="poscounts")

results <- DESeq2::results(
    dds_results,
    format="DataFrame",
    tidy=TRUE,
    cooksCutoff=FALSE,
    contrast=c("case_ctrl", "1", "0")
)
row.names(results) <- results$row
outfile <- "results/batch/deseq2_results.tsv"
write.table(results, file=outfile, sep="\t")
