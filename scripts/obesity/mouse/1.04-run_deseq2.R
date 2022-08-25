#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/mouse/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(biomformat)
library(dplyr)
library(DESeq2)

tbl_file <- "data/obesity/processed/mouse/tbl_merged.genus.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/obesity/processed/mouse/metadata.merged.tsv"
md <- read.table(md_file, sep="\t", row.names=1, header=T)
md$diet <- relevel(as.factor(md$diet), "Standard")
md$instrument <- relevel(as.factor(md$instrument), "Illumina")

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

tbl <- tbl[, sample_order]

design_formula <- as.formula("~diet + instrument")
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
    contrast=c("diet", "HFD", "Standard")
)

row.names(results) <- results$row
results <- results %>% select(-c("row"))

outfile <- "results/obesity/mouse/deseq2_results.tsv"
write.table(results, file=outfile, sep="\t")
