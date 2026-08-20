## ----setup1, include=FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup2,message=F,warning=F,echo=FALSE------------------------------------
# install.packages("devtools")
#devtools::install_github("boldsystems-central/BOLDconnectR")
library(BOLDconnectR)

## ----suggests,message=F,warning=F,echo=FALSE, quietly = TRUE------------------
seq_alignment_pkgs <- all(
  requireNamespace("msa", quietly = TRUE),
  requireNamespace("Biostrings", quietly = TRUE),
  requireNamespace("muscle", quietly = TRUE)
)

if (seq_alignment_pkgs) {
  library(msa)
  library(Biostrings)
  library(muscle)
}

## ----api-key, eval=T,include=FALSE--------------------------------------------
# bold.apikey("")

## ----fetch-data, message=F,warning=F,eval=T-----------------------------------
# Search
data_ids <- bold.public.search(taxonomy = list("Streptocephalus"))
# Fetch
# streptocephalus_data <- bold.fetch(
#   get_by = "processid",
#   identifiers = data_ids$processid,
#   filt_basecount = c(500, 670)
# )
streptocephalus_data <- test.data2

DT::datatable((head(streptocephalus_data, 10)))

## ----summarize-data, message=F,warning=F--------------------------------------
bcdm_summary <- bold.data.summarize(
  bold_df = streptocephalus_data,
  summary_type = "concise_summary"
)
DT::datatable(bcdm_summary$concise_summary)

## ----map-analysis, message=F,warning=F----------------------------------------
map_res <- bold.analyze.map(bold_df = streptocephalus_data)

## ----align-data, message=F,warning=F,eval=seq_alignment_pkgs------------------
seq_align <- bold.analyze.align(
  bold_df = streptocephalus_data,
  marker = "COI-5P",
  cols_for_seq_names = c("processid", "bin_uri"),
  align_method = "Muscle"
)
DT::datatable(
  seq_align[, c("aligned_seq", "msa.seq.name")],
  options = list(pageLength = 10, scrollX = TRUE)
)

## ----tree-analysis,message=F,warning=F,eval=seq_alignment_pkgs----------------
seq_tree <- bold.analyze.tree(
  bold_df = seq_align,
  dist_model = "K80",
  clus_method = "nj",
  tree_plot = TRUE,
  tree_plot_type = "p",
  save_dist_mat = TRUE,
  pairwise.deletion = TRUE
)

seq_tree$base_freq

## ----diversity-analysis_alpha, message=F,warning=F----------------------------
diversity_res <- bold.analyze.diversity(
  bold_df = streptocephalus_data,
  taxon_rank = "species",
  site_type = "locations",
  location_type = "country.ocean",
  diversity_profile = "richness"
)
DT::datatable(diversity_res$richness)

## ----diversity-analysis_beta, message=F,warning=F-----------------------------
beta_diversity_res <- bold.analyze.diversity(
  bold_df = streptocephalus_data,
  taxon_rank = "species",
  site_type = "locations",
  location_type = "country.ocean",
  diversity_profile = "beta",
  beta_index = "jaccard"
)

beta_diversity_res$total.beta

## ----export-data, eval=FALSE--------------------------------------------------
# bold.export(
#   bold_df = streptocephalus_data,
#   export_type = "fas",
#   cols_for_fas_names = c("bin_uri", "genus", "species"),
#   export = "~/Desktop/boldconnectr_sequences.fas"
# )

