### ## Run nichnet for our interest cluster at Pre; Post (Post vs Pre we run in `Fig5_Fibro_#C`, need some additive step)
### For multinichenet  v2.0.1

setwd("/gpfs/chenfeilab/Gaux/CRC_Article/Fig5/res/data/3c/multinichenet/")

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(multinichenetr))

# out dir
pre_out_dir <- "./res/PreR_vs_PreNR/PreR_vs_PreNR_multinichenet_lite_output.rds"
post_out_dir <- "./res/PostR_vs_PostNR/PostR_vs_PostNR_multinichenet_lite_output.rds"


cat(paste0("============================= Start run Multinichenet at ", Sys.time(), " ============================="), "\n")

# unified args
sample_id = "SampleID"
group_id = "Condition"
celltype_id = "split_cell_type"
covariates = NA
batches = NA
min_cells = 5
min_sample_prop = 0.50
fraction_cutoff = 0.05
empirical_pval = FALSE
logFC_threshold = 0.50
p_val_threshold = 0.05
p_val_adj = FALSE 
n.cores = 48
top_n_target = 250
scenario = "regular"
ligand_activity_down = FALSE

lr_network <- readRDS("./data/nichenet_lr_network_human.rds")
lt_matrix <- readRDS("./data/nichenet_lt_matrix_human.rds")

# oi
oi_coarse_cell <- c("CAF", "NLF", "CD4", "CD8", "MonoMph", "DC", "Epi_CNV_high", "Epi_CNV_low", "BC", "Endo")
senders_oi <- receivers_oi <- oi_coarse_cell

# get data

cat(paste0("============================= Get sce obj at ", Sys.time(), " ============================="), "\n")
merge_seurat <- readRDS("./data/merged_seurat_obj.rds")
sce <- Seurat::as.SingleCellExperiment(merge_seurat, assay = "RNA")
sce <- alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
sce_sub <- sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]
sce_sub <- sce_sub[, SummarizedExperiment::colData(sce_sub)[,group_id] %in% c("NR_BL", "NR_END", "R_BL", "R_END")]

### Pre ###
cat(paste0("============================= Run Multinichenet for Pre ", Sys.time(), " ============================="), "\n")

contrasts_oi_pre = c("'R_BL-NR_BL','NR_BL-R_BL'")
contrast_tbl_pre = tibble(contrast = c("R_BL-NR_BL","NR_BL-R_BL"), group = c("R_BL","NR_BL"))
sce_sub_pre <- sce_sub[, SummarizedExperiment::colData(sce_sub)[,group_id] %in% contrast_tbl_pre$group]

multinichenet_output_pre <- multi_nichenet_analysis(
  sce = sce_sub_pre, 
  celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
  batches = batches, covariates = covariates, 
  lr_network = lr_network, ligand_target_matrix = lt_matrix, 
  contrasts_oi = contrasts_oi_pre, contrast_tbl = contrast_tbl_pre, 
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
)
# save
multinichenet_output_pre <- make_lite_output(multinichenet_output_pre)
saveRDS(multinichenet_output_pre, pre_out_dir)

cat("\n","\n")

### Post ###
cat(paste0("============================= Run Multinichenet for Post ", Sys.time(), " ============================="), "\n")

contrasts_oi_post = c("'R_END-NR_END','NR_END-R_END'")
contrast_tbl_post = tibble(contrast = c("R_END-NR_END","NR_END-R_END"), group = c("R_END","NR_END"))
sce_sub_post <- sce_sub[, SummarizedExperiment::colData(sce_sub)[,group_id] %in% contrast_tbl_post$group]

multinichenet_output_post <- multi_nichenet_analysis(
  sce = sce_sub_post, 
  celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
  batches = batches, covariates = covariates, 
  lr_network = lr_network, ligand_target_matrix = lt_matrix, 
  contrasts_oi = contrasts_oi_post, contrast_tbl = contrast_tbl_post, 
  senders_oi = senders_oi, receivers_oi = receivers_oi,
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, 
  min_sample_prop = min_sample_prop,
  scenario = scenario, 
  ligand_activity_down = ligand_activity_down,
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj, 
  empirical_pval = empirical_pval, 
  top_n_target = top_n_target, 
  n.cores = n.cores, 
  verbose = TRUE
)

multinichenet_output_post <- make_lite_output(multinichenet_output_post)
saveRDS(multinichenet_output_post, post_out_dir)

cat("\n","\n")
cat(paste0("============================= Finish at ", Sys.time(), " ============================="), "\n")
