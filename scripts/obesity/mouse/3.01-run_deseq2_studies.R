#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(biomformat)
library(dplyr)
library(DESeq2)
library(stringr)

md_file <- "data/obesity/processed/processed_md.tsv"
md <- read.table(md_file, sep="\t", row.names=1, header=T)
md$diet <- relevel(as.factor(md$diet), "Standard")

run_deseq2 <- function(table, mdata) {
    samples <- colnames(table)
    mdata <- subset(mdata, rownames(mdata) %in% samples)
    sample_order <- row.names(mdata)

    table <- table[, sample_order]

    design_formula <- as.formula("~diet")
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData=table,
        colData=mdata,
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

    return(results)
}

study_dir <- "data/obesity/processed/study_tables"
study_tbl_files <- list.files(study_dir, full.name=T)

for (study_tbl_file in study_tbl_files) {
    study_name <- stringr::str_extract(study_tbl_file, "STUDY_\\d+")
    print(study_name)
    # Skip STUDY_10422 and STUDY_11548 since those two studies only
    # have HFD and no Standard
    if (study_name %in% c("STUDY_10422", "STUDY_11548")) {
        print("Skipping")
        next
    }

    tbl <- biomformat::read_biom(study_tbl_file)
    tbl <- as.matrix(biomformat::biom_data(tbl))
    results <- run_deseq2(tbl, md)

    outfile <- paste0(
        "results/obesity/study_tables/deseq2_results.",
        study_name,
        ".tsv"
    )
    write.table(results, file=outfile, sep="\t")
}


