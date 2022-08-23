#!/home/grahman/miniconda3/envs/da-R/bin/Rscript
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/mouse/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

library(ANCOMBC)
library(biomformat)
library(phyloseq)

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

taxa <- phyloseq::otu_table(tbl, taxa_are_rows=T)
meta <- phyloseq::sample_data(md)
physeq <- phyloseq::phyloseq(taxa, meta)

design_formula <- "diet + instrument"
ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design_formula,
                                    zero_cut=1.0)
results <- as.data.frame(ancombc.results$res)

outfile <- "results/obesity/mouse/ancombc_results.tsv"
write.table(results, file=outfile, sep="\t")
