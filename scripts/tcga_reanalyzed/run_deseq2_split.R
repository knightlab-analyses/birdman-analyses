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

md_file <- "data/tcga_reanalyzed/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)
levels <- unique(md$investigation)

split_dir <- "data/tcga_reanalyzed/processed"
splits <- c(0:4)

design.formula <- as.formula("~investigation + race + gender + data_submitting_center_label")

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

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData=tbl,
        colData=subset_md,
        design=design.formula
    )
    dds_results <- DESeq2::DESeq(dds, sfType="poscounts")

    outdir <- paste0(
        "results/tcga_reanalyzed/deseq2_split/split_", split
    )
    dir.create(outdir, showWarnings=FALSE)

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
        outfile <- paste0(outdir, "/", level, ".tsv")
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
    outfile <- paste0(outdir, "/Intercept.tsv")
    write.table(int_results, file=outfile, sep="\t")
}
