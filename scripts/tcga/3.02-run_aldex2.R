#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

library(biomformat)
library(DESeq2)

tbl_file <- "data/tcga/processed/species/processed_tbl.bacteria.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/tcga/processed/processed_md.tsv"
md <- read.delim(md_file, sep="\t", row.names=1, header=T)
md$investigation <- gsub("TCGA-", "", md$investigation)
md$investigation <- relevel(as.factor(md$investigation), "BRCA")
md$race <- as.factor(md$race)
md$gender <- as.factor(md$gender)

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

tbl <- tbl[, sample_order]

design.formula <- as.formula("~investigation + race + gender")

mm <- model.matrix(design.formula, md)
x <- ALDEx2::aldex.clr(tbl, mm)
aldex2.results <- ALDEx2::aldex.glm(x)

outfile <- "results/tcga/aldex2_results.tsv"
write.table(aldex2.results, file=outfile, sep="\t")
