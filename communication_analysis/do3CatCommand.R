#### #### Get Cell-Cell-Communication (CCC) analysis data at command line #### ####
## Description: Get a merged and a list of CellChat object to do CCC for CRC cohort
## Here we do 3C analysis for our interested target subtype of cell


# args
args = commandArgs(T)
workDir = args[1]
annoDir = args[2]
outDir = args[3]
responseInfoDir = args[4]
# "main" or "sub": if main run cellchat for all main cluster
cellAnno = args[5]
# a character string split by `:`
cellType = args[6]

sub_cell <- unlist(stringr::str_split(cellType, pattern = ":"))

cat(paste0("============================= Start run do3CatCommand at ", Sys.time(), " ============================="), "\n")
cat(paste0("Command line args: "), "\n")
cat(paste("do3CatCommand.R", workDir, annoDir, outDir, responseInfoDir, cellAnno, paste(sub_cell, collapse=" "), sep = " "),"\n")

# # test
# workDir = "/chenfeilab/Gaux/onGoingProject/wjw/CRC_project/cellbenderFilter/07_CellCellCommun/"
# annoDir = "./data/seurat_anno_latest/"
# outDir = "./res_2310/"
# responseInfoDir = "./data/meta/SampleToResponse.csv"
# cellAnno = "sub"
# cellType = "BC:CD4:CD8:Mye"

# Step0: set our env and load package
setwd(workDir)
options(stringsAsFactors = FALSE)
Sys.setenv(PATH = paste("/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/bin", Sys.getenv("PATH"), sep=":"))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(CellChat, lib.loc = "/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/lib/R/library"))

cat(paste0(Sys.time(), " Start get CCC Data"),"\n")
cat(paste0("Using CellChat V", packageVersion("CellChat"), " to do CCC at ", getwd()),"\n")

# Step1: prepare Seurat data
seurat_dir = paste0(outDir, "seurat_obj/")
seurat_obj_dir = paste0(seurat_dir, "response_timepoint_seurat_lst.rds")

if (!file.exists(seurat_obj_dir)) {
  cat(paste0(Sys.time(), " No seurat annotation obj, start prepare seruat annotation data"), "\n")
  dir.create(seurat_dir, showWarnings = TRUE, recursive = TRUE)
  # load all seurat anno
  subcell_dir <- annoDir
  subcell_file <- list.files(subcell_dir, pattern = '.*_anno_ccc.rds', ignore.case = FALSE) %>% paste0(subcell_dir, .)
  subcell_name <- list.files(subcell_dir, pattern = '.*_anno_ccc.rds', ignore.case = FALSE) %>% stringr::str_split_i(., "_", 1)
  anno_lst <- mapply(function(x, y) {x = readRDS(y)}, subcell_name, subcell_file, SIMPLIFY = FALSE)
  # merge seurat
  merge_seurat <- merge(anno_lst[[1]], anno_lst[2:length(anno_lst)], add.cell.ids = subcell_name, project = "merged_anno")
  merge_seurat[['cell_cluster']] <- stringr::str_split_i(colnames(merge_seurat), "_", 1)
  Idents(merge_seurat) <- "cell_cluster"
  # update response info
  responseInfo <- read.csv(responseInfoDir)
  #meta <- merge_seurat@meta.data %>%
  #  select(-Response) %>%
  #  left_join(responseInfo, by = "PatientID")
  meta <- merge_seurat@meta.data %>% 
    left_join(responseInfo, by = "PatientID")
  merge_seurat$Response <- meta$Response
  # remove no Response info cell
  merge_seurat <- merge_seurat[,merge_seurat@meta.data$PatientID != 'I008']
  # add pre and post info
  GetCombineID <- function(
    metadata
  ){
    # BL: baseline timepoint
    BL_sampleID <- metadata %>% filter(Treatment == 'baseline') %>% pull(SampleID) %>% unique()
    # MID: middle time point
    MID_sampleID <- metadata %>% filter(Treatment %in% c('post_IO', 'post_SCRT')) %>% pull(SampleID) %>% unique()
    # END: END time point
    END_sampleID <- metadata %>% filter(Treatment %in% c('surgery', 'post_SCRT+IO', 'post_IO+SCRT')) %>%
      arrange(PatientID, desc(Treatment)) %>%
      distinct(PatientID, .keep_all = TRUE) %>%
      pull(SampleID)
    ID_lst <- list(BL_sampleID = BL_sampleID, MID_sampleID = MID_sampleID, END_sampleID = END_sampleID)
    return(ID_lst)
  }
  ID_lst <- GetCombineID(merge_seurat@meta.data)
  merge_seurat@meta.data <- merge_seurat@meta.data %>% 
    mutate(TimePoint = ifelse(SampleID %in% ID_lst$BL_sampleID, "BL", ifelse(SampleID %in% ID_lst$MID_sampleID, "MID", "END")))
  # split
  response_lst <- SplitObject(merge_seurat, split.by = "Response")
  response_timepoint_lst <- lapply(response_lst, function(x){SplitObject(x, split.by = "TimePoint")}) %>% unlist()
  # save
  cat(paste0(Sys.time(), " Start save seurat object at ", seurat_obj_dir), "\n")
  saveRDS(response_timepoint_lst, seurat_obj_dir)
  cat(paste0(Sys.time(), " Finish prepare seruat annotation data"), "\n")
}else{
  cat(paste0(Sys.time(), " Already exist seurat annotation file, load seruat annotation data at ", seurat_obj_dir), "\n")
  response_timepoint_lst <- readRDS(seurat_obj_dir)
}


# Step2: get cellchat data
cat(paste0(Sys.time(), " Start run cellchat"), "\n")
# CR_nCR_TimePoint <- readRDS(paste0(seurat_dir, "response_timepoint_seurat_lst.rds"))

#' Get CellChat object
#'
#' @param seurat.obj a seurat object to creat cellchat object
#' @param cell.group a character string to define cell groups, must in seurat obj meta.data
#' @param db.use CellChat function to choose LR database from `CellChatDB.human`, `CellChatDB.mouse`, `CellChatDB.zebrafish` all subset of them also create your database
#' @param cell.threshold int. the min number of cells required in each cell group for cell-cell communication
#' @param get.signature logical. Should net signature be returned?
#' @param get.signature.script a character string. signature.script dir
#' @param out.file.name a character string. set the output dir
#' @return A CellChat object
#' @usage GetCellChatObj
#'
GetCellChatObj <- function(seurat.obj, cell.group = "cell_cluster", db.use = CellChatDB.human,
                           cell.threshold = 10, out.file.name)
{
  # print message
  mes <- paste(unique(seurat.obj$Response), unique(seurat.obj$TimePoint), sep = "_")
  cat(paste0(Sys.time(), " Start run cellchat for"), mes, "\n")
  # create basis cellchat obj
  cellchat <- createCellChat(object = seurat.obj, group.by = cell.group)
  # set LRs db to infer CCC
  cellchat@DB <- db.use
  # get the LRs(singnaling) expression data
  cellchat <- subsetData(cellchat)
  # identify over-expressed LRs in each cell group (so too much group will spend too much time)
  cellchat <- identifyOverExpressedGenes(cellchat)
  # then identify over-expressed LR interactions if either L or R is over-expressed.
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI, but we not use, may introduce artifact by this diffusion process,
  # although we run this step, but author suggest set `raw.use = TRUE`, so actually we don't use this data
  cellchat <- projectData(cellchat, PPI.human)
  # inferred intercellular communication network of each ligand-receptor pair
  invisible(capture.output(cellchat <- computeCommunProb(cellchat, population.size = FALSE, raw.use = TRUE), type = "message"))
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = cell.threshold)
  # infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  # calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  # # do Net Analysis computeCentrality
  # BUG: may cause by R version
  # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name="netP")

  saveRDS(cellchat, out.file.name)

  return(cellchat)
}

# unique(response_timepoint_lst$R.BL$cell_cluster)
# unique(response_timepoint_lst$R.BL$cell_type)

group_lst <- list(R_BL = response_timepoint_lst$R.BL, R_END = response_timepoint_lst$R.END,
                  NR_BL = response_timepoint_lst$NR.BL, NR_END = response_timepoint_lst$NR.END)


if (cellAnno == "main") {
  target_subtype <- c("mainType", "eachOther")
} else if (cellAnno == "sub") {
  # target
  target_subtype <- sub_cell
} else {
  message("Invalid cellAnno arg, run cellchat for main cell type")
  target_subtype <- c("mainType", "eachOther")
}

target_dir_name <- paste(target_subtype, collapse = '_vs_')
cc_dir = paste0(outDir, target_dir_name, "/data/00_cellchat_obj/")

if (!dir.exists(cc_dir)) {
  dir.create(cc_dir, recursive=TRUE)
  cat(paste0(Sys.time(), " Start run cellchat for"), paste(target_subtype, collapse = " vs "), "\n")
  # get substype seurat obj
  if (cellAnno == "sub") {
    subtype_seurat <- lapply(group_lst, function(x){subset(x, subset = cell_cluster %in% target_subtype)})
    cell_group <- "cell_type"
  } else {
    subtype_seurat <- group_lst
    cell_group <- "cell_cluster"
  } 
  
  future::plan("multisession", workers = 3)
  options(future.globals.maxSize = 80000 * 1024^2)
  # run
  response_timepoint_CC_obj <- lapply(names(subtype_seurat),
                                    function(x){GetCellChatObj(seurat.obj = subtype_seurat[[x]],
                                                               cell.group = cell_group,
                                                               db.use = CellChatDB.human,
                                                               cell.threshold = 10,
                                                               out.file.name = paste0(cc_dir, x, "_cellchat_obj.rds")
                                    )})
}else{
  cat(paste0(Sys.time(), " Already exist cellchat object in"), cc_dir, ", load them", "\n")
}


# # check if each group has all of subtype 
# subtype_seurat$R_BL$cell_type %>% table() %>% length()
# subtype_seurat$R_END$cell_type %>% table() %>% length()
# subtype_seurat$NR_BL$cell_type %>% table() %>% length()
# subtype_seurat$NR_END$cell_type %>% table() %>% length()


# Step3: do net analysis
# https://github.com/sqjin/CellChat/issues/140
# https://www.r-bloggers.com/2020/09/future-1-19-1-making-sure-proper-random-numbers-are-produced-in-parallel-processing/
# this PR:
# https://github.com/sqjin/CellChat/pull/703
na_dir = paste0(outDir, target_dir_name, "/data/01_netAnalysis_obj/")
if ((dir.exists(cc_dir)) & (!dir.exists(na_dir))) {
  cc_file <- list.files(cc_dir, pattern = '.*_cellchat_obj.rds', ignore.case = FALSE) %>% paste0(cc_dir, .)
  cc_name <- list.files(cc_dir, pattern = '.*_cellchat_obj.rds', ignore.case = FALSE) %>% stringr::str_extract(., "(R|NR)_(BL|END)")
  response_timepoint_CC_obj <- mapply(function(x, y) {x = readRDS(y)}, cc_name, cc_file, SIMPLIFY = FALSE)
  
  cat(paste0(Sys.time(), " Start run NetAnalysis"), "\n")
  dir.create(na_dir, recursive=TRUE)
  for (obj in names(response_timepoint_CC_obj)) {
    cat(paste0(Sys.time(), " Run NetAnalysis for"), obj, "\n")
    do.computeCentrality.command <- paste("Rscript", "doNetAnalysis.R", obj, cc_dir, na_dir)
    invisible(capture.output(system(do.computeCentrality.command)))
  }
  cat(paste0(Sys.time(), " Finish run NetAnalysis"), "\n")
}else{
  cat(paste0(Sys.time(), " Already exist cellchat NetAnalysis object in"), na_dir, ", load them", "\n")
}

# Step4:  merge
merged_dir = paste0(outDir, target_dir_name, "/data/02_merged_cellchat_obj/")
if ((dir.exists(na_dir)) & (!dir.exists(merged_dir))) {
  na_file <- list.files(na_dir, pattern = '.*_withCentrality.rds', ignore.case = FALSE) %>% paste0(na_dir, .)
  na_name <- list.files(na_dir, pattern = '.*_withCentrality.rds', ignore.case = FALSE) %>% stringr::str_extract(., "(R|NR)_(BL|END)")
  dir.create(merged_dir, recursive=TRUE)
  # list
  cat(paste0(Sys.time(), " Start merge cellchat NetAnalysis object and save"), "\n")
  response_timepoint_NA_obj <- mapply(function(x, y) {x = readRDS(y)}, na_name, na_file, SIMPLIFY = FALSE)
  saveRDS(response_timepoint_NA_obj, paste0(merged_dir, target_dir_name, "_cellchat_obj_lst.rds"))
  # merged
  cellchat <- mergeCellChat(response_timepoint_NA_obj, add.names = names(response_timepoint_NA_obj))
  saveRDS(cellchat, paste0(merged_dir, target_dir_name, "_merged_cellchat_obj.rds"))
  cat(paste0(Sys.time(), " Finish save merged and list cellchat NetAnalysis object"), "\n")
  cat(paste0("============================= Finish run do3CatCommand at ", Sys.time(), " ============================="), "\n")
}else{
  cat(paste0(Sys.time(), " Already exist merged and list cellchat NetAnalysis object"), "\n")
  cat(paste0("============================= Finish run do3CatCommand at ", Sys.time(), " ============================="), "\n")
}

