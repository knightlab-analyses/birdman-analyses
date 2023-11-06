#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

library(biomformat)
library(DESeq2)

tbl_file <- "data/tcga_reanalyzed/processed/merged_tbl.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/tcga_reanalyzed/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

tbl <- tbl[, sample_order]

design.formula <- as.formula("~investigation + race + gender + data_submitting_center_label")

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=tbl,
    colData=md,
    design=design.formula
)
dds_results <- DESeq2::DESeq(dds, sfType="poscounts")

levels <- unique(md$investigation)

for (level in levels) {
    print(level)
    if (level == "BRCA") {
        print("Skipping")
        next
    }
    results <- DESeq2::results(
        dds_results,
        format="DataFrame",
        tidy=TRUE,
        cooksCutoff=FALSE,
        contrast=c("investigation", level, "BRCA")
    )
    row.names(results) <- results$row
    outfile <- paste0("results/tcga_reanalyzed/deseq2/", level, ".tsv")
    write.table(results, file=outfile, sep="\t")
}

int_results <- DESeq2::results(
    dds_results,
    format="DataFrame",
    tidy=TRUE,
    cooksCutoff=FALSE,
    name="Intercept"
)
row.names(int_results) <- int_results$row
outfile <- "results/tcga_reanalyzed/deseq2/Intercept.tsv"
write.table(int_results, file=outfile, sep="\t")
