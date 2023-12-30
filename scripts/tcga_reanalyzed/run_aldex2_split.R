#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

library(biomformat)
library(ALDEx2)

md_file <- "data/tcga_reanalyzed/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)

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

    mm <- model.matrix(design.formula, subset_md)
    x <- ALDEx2::aldex.clr(tbl, mm)
    aldex2.results <- ALDEx2::aldex.glm(x)

    outdir <- paste0(
        "results/tcga_reanalyzed/aldex2_split/"
    )
    dir.create(outdir, showWarnings=FALSE)
    outfile <- paste0(
        outdir, "split_", split, ".results.tsv"
    )
    write.table(aldex2.results, file=outfile, sep="\t")
}
