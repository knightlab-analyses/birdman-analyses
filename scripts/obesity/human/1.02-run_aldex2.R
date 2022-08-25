#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/human/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(biomformat)
library(ALDEx2)

tbl_file <- "data/obesity/processed/human/processed_tbl.genus.biom"
tbl <- biomformat::read_biom(tbl_file)
tbl <- as.matrix(biomformat::biom_data(tbl))

md_file <- "data/obesity/processed/human/processed_md.tsv"
md <- read.table(md_file, sep="\t", row.names=1, header=T)
md$diet <- relevel(as.factor(md$obese), "Lean")

samples <- colnames(tbl)
md <- subset(md, rownames(md) %in% samples)
sample_order <- row.names(md)

tbl <- tbl[, sample_order]

design_formula <- as.formula("~obese")

mm <- model.matrix(design_formula, md)
x <- ALDEx2::aldex.clr(tbl, mm)
aldex2.results <- ALDEx2::aldex.glm(x)

outfile <- "results/obesity/human/aldex2_results.tsv"
write.table(aldex2.results, file=outfile, sep="\t")
