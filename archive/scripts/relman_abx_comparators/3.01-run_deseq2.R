#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(biomformat)
library(DESeq2)

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

design.formula <- as.formula("~antibiotic + host_subject_id")

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=tbl,
    colData=md,
    design=design.formula
)
dds_results <- DESeq2::DESeq(dds, sfType="poscounts")

dir.create("results/relman_abx/deseq2")
levels <- c("preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp",
            "SecondWPC", "PostCp")
num_levels <- length(levels) - 1
for (i in c(1:num_levels)) {
    comp = paste0(levels[i+1], "_vs_", levels[i])

    print(comp)
    results <- DESeq2::results(
        dds_results,
        format="DataFrame",
        tidy=TRUE,
        cooksCutoff=FALSE,
        contrast=c("antibiotic", levels[i+1], levels[i])
    )
    row.names(results) <- results$row
    outfile <- paste0("results/relman_abx/deseq2/", comp, ".tsv")
    write.table(results, file=outfile, sep="\t")
}
