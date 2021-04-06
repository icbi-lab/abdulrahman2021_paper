#!/usr/bin/env Rscript
'
run_de.R

Run a differential expression analysis using edgeR on single cell data given a predefined
SingleCellExperiment, design and contrast matrix.

Usage:
  run_de.R <input_data> <output_file> [--cpus=<n> --excel=<excel>]

Options:
  <input_data>      An .rda file containing three objects: adata, design, contrasts
  <output_file>     Results will be written to this file in tsv format.
  --excel=<excel>   Results will be (additionally) written to this file in excel format
  --cpus=<n>        Number of cpus to use. [default: 8].

' -> doc

library(docopt)
arguments = docopt(doc)

library(BiocParallel)
library(conflicted)
library(edgeR)
library(dplyr)
library(readr)
library(writexl)

print(arguments)

# set paralellism
options(mc.cores=as.numeric(arguments[['--cpus']]))
register(MulticoreParam(as.numeric(arguments[['--cpus']])))
print(registered())

# contains adata, design, constrasts
load(file=arguments[['<input_data>']])

var_data = data.frame(gene_symbol = rownames(tmp_counts))

dge = DGEList(counts=tmp_counts, samples=tmp_obs, genes=var_data)

message("Calculating NormFactors...")
dge <- calcNormFactors(dge)

message("Estimating Dispersion...")
dge <- estimateDisp(dge, design = design)

message("Fitting linear model...")
fit <- glmQLFit(dge, design = design)

message("Testing contrasts...")
qlfs = sapply(colnames(contrasts), function(cluster) {
    message(paste0("working on ", cluster))
    glmQLFTest(fit, contrast=contrasts[,cluster])
}, USE.NAMES=TRUE, simplify=FALSE)


message("Preparing results...")
tts = sapply(qlfs, function(qlf) {
    topTags(qlf, n=Inf, adjust.method="BH")
}, USE.NAMES=TRUE, simplify=FALSE)

all_results = bind_rows(lapply(names(tts), function(cluster) {
    tts[[cluster]]$table %>% mutate(cluster=cluster)
})) %>% as_tibble()

# save(design, dge, fit, qlfs, tts, all_results, file="../../../results/downstream_analysis/edger/results.rda", compress = FALSE)

message("Writing results...")
write_tsv(all_results, arguments[["<output_file>"]])

if(!is.null(arguments[['--excel']])) {
  write_xlsx(all_results, arguments[['--excel']])
}
