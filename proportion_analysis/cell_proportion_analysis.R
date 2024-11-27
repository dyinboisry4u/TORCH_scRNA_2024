# ---- pkgs ----
# # # pkg
library(dplyr)
library(tidyr)
library(circlize)
library(stringr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tibble)
library(purrr)
library(ggbeeswarm)
library(ggh4x)
library(gtools)
library(rstatix)
library(ggforce)
library(ggalluvial)
library(ggradar)
library(ComplexHeatmap)
library(patchwork)
library(ggcorrplot, lib.loc = "/share/home/Grape/software/install_pkg/miniconda3/envs/mainenv/lib/R/library")
# library(ggcorrplot2)

options(bitmapType='cairo')

# # # Func
# # ############## 1. Preprocess ##############

#' Prepare data for TORCH trial cohort cell prop analysis
#' 
#' @param obs.dir a character string indicating **Scanpy** `.obs` csv output dir 
#' @param resp.dir a character string indicating `PatientID - Response` table dir
#' @param cell.column description
#' @param anno.ord a vector indicating the order of cell type
#' @param remove.no.resp a bool, if remove I008_BL
#' @note here we update response info, remove I008_BL sample(out of group) and set cell factor level
#' @return a data.frame of cell proportion with group
#' 
PrepareData <- function(
    obs.dir,
    resp.dir,
    cell.column = "CellType",
    anno.ord = NULL,
    remove.no.resp = TRUE
){
  
  cell.column = as.name(cell.column)
  
  response_info <- read.csv(resp.dir)
  # add response info and remove I008_BL
  cell_anno <- read.csv(obs.dir, row.names = 1) %>%
    dplyr::mutate(PatientID = stringr::str_split_i(SampleID, "_", 1)) %>%
    dplyr::left_join(response_info, by = 'PatientID') %>%
    dplyr::rename(SubCluster = {{cell.column}})
  
  if (is.character(anno.ord) & all(anno.ord %in% cell_anno$SubCluster)) {
    cell_anno <- cell_anno %>% 
      dplyr::mutate(SubCluster = factor(SubCluster, levels = anno.ord))
  } else {
    cell_anno <- cell_anno %>% 
      dplyr::mutate(SubCluster = factor(SubCluster, levels = unique(cell_anno$SubCluster)))
  }
  
  if (remove.no.resp) {
    cell_anno <- cell_anno %>% dplyr::filter(SampleID != 'I008_BL')
  }
  
  return(cell_anno)
}


#' Get cell proportion from metadata
#' 
#' @description
#' Get `sample level` or `condition level` cell proportion data.frame form metadata, return long or width data
#' 
#' @param meta.data a data.frame to get cell proportion, from **Seurat** `@meta.data` or **Scanpy** `.obs`
#' @param group.by a character string or **dplyr** style data-variable indicating the group/column of cell proportion
#' @param cell.column a character string or **dplyr** style data-variable indicating the column of cell annotation, this column must be a **factor**
#' @param add.meta a bool, if add remained meta data to result? if not group by sample id, pls set to FALSE
#' @param return.long a bool, if return long data?
#' 
#' @return a data.frame of cell proportion with group
#' @usage **sample level** GetCellProport(meta.data, group.by = "SampleID", cell.column = "SubCluster", add.meta = TRUE, return.long = TRUE...)
#' 
GetCellProportion <- function(
    meta.data,
    group.by = "SampleID",
    cell.column = "SubCluster",
    add.meta = TRUE,
    return.long = TRUE
) {
  
  # get all cell column
  cell_anno <- meta.data %>% magrittr::extract2(cell.column) %>% levels()
  
  group.by = as.name(group.by)
  cell.column = as.name(cell.column)
  
  # get cell prop
  cell_proportion <- meta.data %>% 
    group_by({{group.by}}, {{cell.column}}) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by({{group.by}}) %>%
    mutate(Freq = n / sum(n)) %>%
    select(-n) %>%
    spread(key = {{cell.column}}, value = Freq, fill = 0) %>%
    ungroup() 
  
  if (add.meta) {
    # get metadata
    sample_info <- meta.data %>%
      select(-{{cell.column}}) %>%
      unique()
    cell_proportion <- left_join(cell_proportion, sample_info, by = join_by({{group.by}}))
  }
  
  # only keep cell column
  if (return.long) {
    # exclude all meta
    if (add.meta) {
      cell_proportion <- cell_proportion %>% 
        gather(key = 'Cluster', value = 'Proportion', -colnames(cell_proportion)[! colnames(cell_proportion) %in% cell_anno])
    } 
    # exclude group
    else {
      cell_proportion <- cell_proportion %>% 
        gather(key = 'Cluster', value = 'Proportion', -{{group.by}})
    }
  }
  
  return(cell_proportion)
}

#' Get our experimental design group SampleID or PatientID for any timecourse
#' 
#'
#' @param cell.prop a data.frame of wide type cell proportion 
#' @param deal.dup if deal duplicate samples if in same time point
#' @note
#' 
#' @return a list of SampleID
#' @usage 
#' 
GetCombineID <- function(
    cell.prop,
    deal.dup = TRUE
){
  
  sample_num <- unique(cell.prop$SampleID) %>% length()
  if (sample_num < 58) {
    warning("Some samples do not have any cell in current subcluster")
  }
  
  # BL: baseline timepoint
  BL_sampleID <- cell.prop %>%
    filter(Treatment == 'baseline') %>%
    pull(SampleID)
  # MID: middle time point
  MID_sampleID <- cell.prop %>% 
    filter(Treatment %in% c('post_IO', 'post_SCRT')) %>%
    pull(SampleID)
  # END: END time point
  if (deal.dup) {
    END_sampleID <- cell.prop %>% 
      filter(Treatment %in% c('surgery', 'post_SCRT+IO', 'post_IO+SCRT')) %>%
      arrange(PatientID, desc(Treatment)) %>%
      distinct(PatientID, .keep_all = TRUE) %>%
      pull(SampleID)
  } else {
    END_sampleID <- cell.prop %>% 
      filter(Treatment %in% c('surgery', 'post_SCRT+IO', 'post_IO+SCRT')) %>%
      pull(SampleID)
  }
  
  # paired patient
  paired_patientID <- cell.prop %>%
    group_by(PatientID) %>%
    filter(n() > 1) %>%
    distinct(PatientID) %>%
    pull(PatientID)
  # all paired patient
  all_paired_patientID <- cell.prop %>%
    group_by(PatientID) %>%
    filter(n() >= 3) %>%
    distinct(PatientID) %>%
    pull(PatientID)
  
  ID_lst <- list(BL_sampleID = BL_sampleID, MID_sampleID = MID_sampleID, END_sampleID = END_sampleID, 
                 paired_patientID = paired_patientID, all_paired_patientID = all_paired_patientID)
  
  return(ID_lst)
}


#' Get our experimental design group SampleID or PatientID for any treatment
#' 
#'
#' @param cell.prop a data.frame of wide type cell proportion 
#' @param deal.dup if deal duplicate samples if in same time point
#' 
#' @return a list of SampleID
#' @usage 
#' 
GetTreatmentID <- function(
    cell.prop
){
  
  sample_num <- unique(cell.prop$SampleID) %>% length()
  if (sample_num < 58) {
    warning("Some samples do not have any cell in current subcluster")
  }
  
  # BL sample ID
  BL_sampleID <- cell.prop %>% filter(Treatment == 'baseline') %>% pull(SampleID)
  # IO sample ID
  IO_sampleID <- cell.prop %>% filter(Treatment == 'post_IO') %>% pull(SampleID)
  # RT sample ID
  RT_sampleID <- cell.prop %>% filter(Treatment == 'post_SCRT') %>% pull(SampleID)
  # IO+RT
  IO_RT_sampleID <- cell.prop %>% filter(Treatment == 'post_IO+SCRT') %>% pull(SampleID)
  # RT+IO
  RT_IO_sampleID <- cell.prop %>% filter(Treatment == 'post_SCRT+IO') %>% pull(SampleID)
  
  ID_lst <- list(BL = BL_sampleID, IO = IO_sampleID, RT = RT_sampleID, IORT = IO_RT_sampleID, RTIO = RT_IO_sampleID)
  
  return(ID_lst)
}

#' Get our experimental design group prop
#' 
#' @param cell.prop a data.frame of long type cell proportion 
#' @param id id list of sample, `GetCombineID` output
#' @note
#' 
#' @return a list of all group sample proportion
#' @usage GetGroupProportion(cell.prop, id)
#' 
GetConditionProportion <- function(
    cell.prop,
    id
){
  
  # Tag: time course with response
  BL_prop <- cell.prop %>% filter(SampleID %in% id$BL_sampleID) %>% mutate(Timecourse = "Pre", Tag = paste0("Pre_", Response))
  MID_prop <- cell.prop %>% filter(SampleID %in% id$MID_sampleID) %>% mutate(Timecourse = "Mid", Tag = paste0("Mid_", Response))
  END_prop <- cell.prop %>% filter(SampleID %in% id$END_sampleID) %>% mutate(Timecourse = "Post", Tag = paste0("Post_", Response))
  # C / I
  C_BL_prop <- BL_prop %>% filter(Group == 'consolidation')
  I_BL_prop <- BL_prop %>% filter(Group == 'induction')
  C_MID_prop <- MID_prop %>% filter(Group == 'consolidation')
  I_MID_prop <- MID_prop %>% filter(Group == 'induction')
  C_END_prop <- END_prop %>% filter(Group == 'consolidation')
  I_END_prop <- END_prop %>% filter(Group == 'induction')
  
  condition_prop_lst <- list(BL = BL_prop, MID = MID_prop, END = END_prop, 
                             C_BL = C_BL_prop, I_BL = I_BL_prop, C_MID = C_MID_prop, I_MID = I_MID_prop, C_END = C_END_prop, I_END = I_END_prop)
  
  return(condition_prop_lst)
}

# # ############## 2. Statistic ##############

#' Do Mann-Whitney-U test / Unpaired Wilcoxon test
#' 
#' @description
#' Use `ggpubr` or `rstatix` do wilcox test for some variable
#' 
#' @param object a cell proportion object you want to test
#' @param var a character string indicating corresponding test variable, must be a **factor** for multiple 1 vs 1 test
#' @param var.ord a character vector indicating `var` factor level
#' @param group group of test
#' @param comparisons a list of compare
#' @param scale set the `rstatix::add_xy_posion()` c("fixed", "free", "free_y")
#' @param p.pos.increse P-value positions increase step: fraction of total height, if `scale = "free_y"` calculate by each plot panels
#' @param stat.method use `ggpubr` or `rstatix` to do wilcox test
#' @param rm.outlier a bool, if remove outlier before do test
#' 
#' @return wilcox test result
#' @usage CellProportionWilcox(arg1, arg2 = default, ...)
#' 
CellProportionWilcox <- function(
    object,
    var = "Response",
    var.ord = c("nCR", "CR"),
    value = "Proportion",
    group = "Cluster",
    comparisons = NULL,
    scale = "fixed",
    p.pos.increse = 0.25,
    position.group = NULL,
    stat.method = c("ggpubr", "rstatix"),
    rm.outlier = FALSE
){
  
  stat.method <- match.arg(stat.method)
  
  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }
  
  wilcox_formula = as.formula(paste(value, var, sep = "~"))
  
  var_name = as.name(var)
  value_name  = as.name(value)
  group_name = as.name(group)
  
  # must set factor for the multiple 1 vs 1 test
  object <- object %>% mutate(!!enquo(var_name) := factor(!!enquo(var_name), var.ord))
  
  if (rm.outlier) {
    object <- object %>% 
      group_by(!!enquo(var_name), !!enquo(group_name)) %>%
      mutate(outlier = is_outlier(!!enquo(value_name))) %>% 
      filter(outlier == FALSE)
  }
  
  if (stat.method == 'ggpubr') {
    res <- compare_means(
      wilcox_formula, data = object, group.by = group,
      method = "wilcox.test", paired = FALSE
    )
  }
  # but rstatix::wilcox_test() have bugs 
  if (stat.method == 'rstatix') {
    res <- object %>% 
      group_by(!!enquo(group_name)) %>% 
      wilcox_test(wilcox_formula, comparisons = comparisons, paired = FALSE, alternative = "two.sided") %>%
      # fix no facet x location
      # note: step.increse is a fraction
      add_xy_position(x = group, step.increase = p.pos.increse, scales = scale, group = position.group,) 
  }
  
  return(res)
}

#' Do t-test for paired sample
#' 
#' @description
#' Use `rstatix` do paired t-test for paired sample
#' 
#' @param object a paired cell proportion object
#' 
#' 
CellProportionT <- function(
    object,
    var = "Response",
    var.ord = c("nCR", "CR"),
    value = "Proportion",
    group = "Cluster",
    comparisons = NULL,
    scale = "free_y",
    p.pos.increse = 0.25
){
  
  t_formula = as.formula(paste(value, var, sep = "~"))
  
  var_name = as.name(var)
  value_name  = as.name(value)
  group_name = as.name(group)
  
  object <- object %>% mutate(!!enquo(var_name) := factor(!!enquo(var_name), var.ord))
  
  # use rstatix
  res <- object %>% 
    group_by(!!enquo(group_name)) %>% 
    t_test(t_formula, comparisons = comparisons, paired = TRUE, alternative = "two.sided") %>%
    add_xy_position(x = group, step.increase = p.pos.increse, scales = scale) 
  
  return(res)
}


## Do R*C cross table chi-square test to get Ro/e or pearson residual for two group
DoRCchiTest <- function(
    cell_anno,
    group_col = "Response",
    res
){
  
  cross_tab <- table(cell_anno$SubCluster, cell_anno[[group_col]])
  chi_res <- chisq.test(cross_tab)
  
  if (res == 'Roe') {
    Roe_res <- chi_res$observed / chi_res$expected
    return(Roe_res)
  }
  
  if (res == 'residual') {
    # pearson_residual <- chi_res$residuals
    pearson_residual <- (chi_res$observed - chi_res$expected) / sqrt(chi_res$expected)
    return(pearson_residual)
  }
  
  # res_lst = list(Roe_table = Roe_res, residual_table = pearson_residual)
}

#' Get Wilcoxon test result
#' 
#' @description
#' package for `CellProportionWilcox`, to get wilcox result for all sample
#' 
#' @param data test data
#' @param data.name a character string indicating task name
#' @param ... same as `CellProportionWilcox`
#' @return wilcox test result
#' @usage GetWilcoxRes(arg1, arg2 = default, ...)
#' 
GetWilcoxRes <- function(
    data, 
    data.name,
    var = "Response",
    var.order = c("nCR", "CR"),
    value = "Proportion",
    group = "Cluster",
    stat.method = 'ggpubr',
    outlier = FALSE
){
  
  res_name = if(outlier) paste0(data.name, '_rm_outlier_wilcox_res') else paste0(data.name, '_wilcox_res')
  
  print(paste("Use", stat_method, "to get", res_name, "with formula:", value, "~", var))
  
  res <- CellProportionWilcox(object = data,
                              var = var,
                              var.ord = var.order,
                              value = value,
                              group = group,
                              stat.method = stat.method,
                              rm.outlier = outlier
  )
  
  return(res)
}

#' get CI with bootstrap
#' 
#' @param conf.level confidence level
#' @param n.resample counts of resamping
#' @param is.response if consider response info?
#' 
# use `Hmisc` to obtain CI for the population mean without assuming normality
getCIwithBoot <- function(
    data,
    boot.var="Proportion",
    conf.level = 0.90,
    n.resample = 5000,
    is.response = TRUE
){
  # recurrence
  set.seed(444)
  
  print(paste("Resampling", n.resample, "times to get", conf.level, "confidence level CI"))
  
  boot <- function(x, var = boot.var) {
    CI_res <- rbind(Hmisc::smean.cl.boot(x[[var]], 
                                         conf.int = conf.level,
                                         B = n.resample,
                                         na.rm=TRUE)) %>%
      data.frame %>%
      setNames(c("mean","lwr","upr"))
    return(CI_res)
  }
  
  if (is.response) {
    res <- data %>%
      group_by(Response, Timecourse, Cluster) %>%
      do(boot(.))
  }
  else{
    res <- data %>%
      group_by(Timecourse, Cluster) %>%
      do(boot(.))
  }
  
  return(res)
}


# # ############## 3. Plot ##############

#' my ggplot2 publish theme 
#' 
#' @param base.size base size
#' @param base.family font
#' @note fix theme_test
#' 
#' @return ggplot2 theme
#' 
PublicationTheme <- function(
    base.size = 15, 
    base.family = "Helvetica"
) {
  
  half_line <- base.size/2
  # fix theme_test color to black
  theme_test(base_size = base.size, base_family = base.family) %+replace% 
    theme(
      # plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      # text = element_text(colour = "black"),
      # panel.background = element_rect(colour = NA),
      # plot.background = element_rect(colour = NA),
      # panel.border = element_rect(colour = NA),
      panel.border = element_rect(fill = NA, colour = "black"),
      # axis.title = element_text(face = "bold",size = rel(1)),
      # axis.title.y = element_text(angle = 90, vjust = 2),
      # axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(size = rel(0.8), colour = "black"),
      # axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # legend.key = element_rect(colour = NA),
      # legend.position = "bottom",
      # legend.direction = "horizontal",
      # legend.key.size= unit(0.2, "cm"),
      # plot.margin = unit(c(10,5,5,5),"mm"),
      strip.text = element_text(colour = "black", 
                                size = rel(0.8), 
                                margin = margin(0.8 * half_line, 0.8 * half_line, 0.8 * half_line, 0.8 * half_line)
      )
      # strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      # strip.text = element_text(face="bold")
    )
}

#' Plot alluvial plot
#' 
#' @description
#' Use `ggalluvial` to plot a alluvial/stack bar plot to show cell proportion change during treatment
#' 
#' @param cell.prop a data.frame of cell proportion
#' @param cell.col a color vector indicating cell type color 
#' @param group.col a color vector indicating `Group` color
#' @param treat.col a color vector indicating `Treatment` color
#' 
#' @return a list of SampleID
#' @usage 
#' 
PlotAlluvialStackBarPlot <- function(
    cell.prop,
    x = "Response",
    x.level = c("NR", "R"),
    cell.col = c("Epithelial"="#28A5AB", "Immune"="#6AAFCE", "Stromal"="#469970"),
    alpha = 0.5,
    facet.by = "Timepoint",
    facet.level = c("Pre", "Post"),
    facet.col = c("Pre" = "#69af86", "Post" = "#fea477")
){
  # 
  p <- cell.prop %>% 
    ggplot(., aes(x = factor(get(x), levels = x.level), y = Proportion, fill = Cluster, 
                  stratum = Cluster, alluvium = Cluster)) + 
    geom_stratum(color = NA, width = 0.6, size = 0.5) +
    geom_flow(alpha = alpha)
  
  if (all(is.character(facet.by), is.character(facet.level), is.character(facet.col))) {
    p <- p + facet_wrap2(~factor(get(facet.by), levels = facet.level), nrow = 1, 
                         strip = strip_themed(background_x = elem_list_rect(fill = facet.col)))
  }
  
  p <- p +
    PublicationTheme() +
    # theme_test(base_size = 15) +
    labs(x=NULL, y="Cell Proportion", title = NULL) +
    scale_fill_manual(values = cell.col) +
    theme(legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.ticks.x = element_blank()
    )
  
  return(p)
}

#' Plot facet boxplot
#'
#' @param cell.prop a data.frame of cell proportion
#' 
#' @note `fill_box = TRUE, add_point = FALSE`, `fill_box = FALSE, add_point = TRUE`. `fill_box = TRUE, add_point = TRUE`
#' 
PlotFacetBoxPlot <- function(
    object,
    stat_object,
    x,
    y,
    facet,
    x_order,
    x_color,
    fill_box = TRUE,
    add_point = FALSE,
    point_size = 1,
    facet_order,
    facet_color,
    facet_nrow = 1,
    scale = "free_y",
    stat_method = 'ggpubr',
    show_outlier = TRUE
){
  
  outlier = if(show_outlier & !add_point) NULL else NA
  
  p <- ggplot(object, aes(x = factor(get(x), levels = x_order), y = get(y)))
  
  if (fill_box) {
    p <- p + stat_boxplot(geom = "errorbar", cex = 0.3,
                          width = 0.3, color = "black", 
                          position = position_dodge(0.8)
    ) + 
      # to hide the vertical lines behind the box
      geom_boxplot(fatten = NA, colour = "white", fill = "white", 
                   width = 0.8, cex = 0.3, 
                   position = position_dodge(0.8), 
                   outlier.colour = NA) +
      geom_boxplot(aes(fill = factor(get(x), levels = x_order)), 
                   width = 0.8, color = "black",
                   fatten = 1.2, cex = 0.3, alpha = 0.6,
                   outlier.shape = outlier, outlier.alpha = 1,
                   position = position_dodge(0.8)
      ) +
      scale_fill_manual(values = x_color)
  } else {
    p <- p + stat_boxplot(geom = "errorbar", cex = 0.3,
                          aes(color = factor(get(x), levels = x_order)), 
                          width = 0.3, position = position_dodge(0.8)
    ) + 
      geom_boxplot(aes(color = factor(get(x), levels = x_order)), 
                   width = 0.8, fill = NA,
                   fatten = 1.2, cex = 0.3,
                   outlier.shape = outlier,
                   position = position_dodge(0.8)
      ) +
      scale_color_manual(values = x_color)
  }
  
  if (add_point) {
    # samller stroke(dot line width) to fit small figure
    p <- p + geom_point(aes(fill = factor(get(x), levels = x_order)),
                        position = position_jitter(width = 0.15, seed = 444), shape = 21, stroke = 0.25,
                        size = point_size, color = "black") +
      scale_fill_manual(values = x_color)
  }
  
  x_strip = strip_themed(background_x = elem_list_rect(fill = facet_color))
  
  p <- p +
    PublicationTheme() +
    labs(x = NULL, y = 'Cell Proportion') +
    theme(legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_text(color = 'black'),
          # axis.text.x = element_text(angle = 45, hjust = 1)
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    ) +
    # separate y axis by set scales = "free_y"
    facet_wrap2(~factor(get(facet), levels = facet_order), scales = scale,
                nrow = facet_nrow, strip = x_strip) +
    theme(strip.switch.pad.grid = unit(1, "inch"))
  
  if (stat_method == 'ggpubr') {
    p <- p + scale_y_continuous(limits = c(0, 0.8)) + 
      stat_pvalue_manual(stat_object, 
                         y.position = 0.75,
                         label = "p.format")
  }
  
  if (stat_method == 'rstatix') {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      stat_pvalue_manual(stat_object,
                         # hide.ns = "p",
                         tip.length = 0,
                         # step.increase = -0.02,
                         # y.position = 0.6,
                         label = "p")
  }
  
  return(p)
}


#' Plot facet paired sample dot line plot
#'
#' @param pairwise_data a data.frame of pairwise cell proportion
#' @param x x variable
#' @param x_order x variable factoe level
#' @param line_group group to add line link tow dot
#' @param dot_shape_group dot shape 
#' @note 
#' 
#' 
#' 
PlotPairwiseLine <- function(
    pairwise_data,
    stat_object,
    x = "Tag",
    x_order = c("Pre", "Post"),
    y = "Proportion",
    dot_group = "Response",
    dot_shape_group = "Group",
    line_group = "PatientID",
    scale = "free_y",
    facet = "Cluster",
    facet_order,
    facet_color,
    dot_color,
    facet_nrow = 2,
    dot_shape = c(21, 22)
    
){
  
  x_strip = strip_themed(background_x = elem_list_rect(fill = facet_color))
  
  p <- ggplot(pairwise_data, aes(x = factor(get(x), levels = x_order), y = get(y), group = get(dot_group))) + 
    geom_line(aes(group = get(line_group)), color = "black", linetype = 2, linewidth = 0.3) +
    geom_point(aes(color = factor(get(x), levels = x_order), shape = get(dot_shape_group)), fill = NA, stroke = 0.6, size = 2.5) +
    PublicationTheme() +
    # scale_color_manual(values = c("#9dd3a8", "#69af86", "#fad6a6", "#fea477")) +
    labs(x=NULL, y='Cell Proportion', title = NULL) +
    scale_colour_manual(values = dot_color) +
    scale_shape_manual(values = dot_shape) +
    # scale_colour_manual(values = line_color) +
    facet_wrap2(~factor(get(facet), levels = facet_order), scales = scale,
                nrow = facet_nrow, strip = x_strip) +
    theme(legend.title = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
    )
  # add pval
  if(!is.null(stat_object)){
    p <- p + scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      stat_pvalue_manual(stat_object,
                         # hide.ns = "p",
                         tip.length = 0,
                         # step.increase = -0.02,
                         # y.position = 0.6,
                         label = "p")
    
  }
  return(p)
  
}

#' Plot CI region line plot
#'
#' @param name description
#' @param name description
#' 
PlotLineCI <- function(
    data,
    prop.data = NULL, 
    x = "Timecourse",
    x.order = c("Pre", "Mid", "Post"),
    is.group = FALSE,
    add_point = FALSE,
    scale = "free_y",
    facet_nrow = 2,
    facet.order, 
    facet.color
){
  
  x_strip = strip_themed(background_x = elem_list_rect(fill = facet.color))
  
  data = data %>% mutate(Cluster = factor(Cluster, levels=facet.order))
  
  if (is.group) {
    p <- ggplot(data, aes(x = factor(get(x), levels = x.order)))
    
    if (add_point & is.data.frame(prop.data)) {
      p <- p + geom_point(data = prop.data, mapping = aes(x = factor(get(x), levels = x.order), y = Proportion, color = Response), 
                          shape = 21,fill = NA, stroke = 0.5, size = 1, position = position_dodge(0.8))
    }
    
    p <- p +
      geom_errorbar(aes(ymin = lwr, ymax = upr, group = Response), width = 0.25, linewidth = 0.4, position = position_dodge(0.8)) +
      geom_point(aes(x = factor(get(x), levels = x.order), y = mean, group = Response, fill = Response), position = position_dodge(0.8), color = "black", size=2, shape = 21) +
      geom_line(aes(y = mean, color = Response, group = Response), linetype = 2, position = position_dodge(0.8), linewidth = 0.4) +
      scale_fill_manual(values=c("#AFC38E", "#E2AA81")) +
      scale_colour_manual(values=c("#AFC38E", "#E2AA81")) +
      facet_wrap2(~factor(Cluster, levels = facet.order), scales = scale,
                  nrow = facet_nrow, strip = x_strip) +
      PublicationTheme() +
      theme(legend.title = element_blank(),
            plot.title = element_text(size = 16, hjust = 0.5),
            axis.title = element_text(size = 12),
            # axis.text.x = element_blank(),
            # axis.ticks.x = element_blank(),
            # legend.position = "none"
      ) + 
      labs(x=NULL, y='Cell Proportion')
    
    return(p)
  }
  
  else {
    p <- ggplot(data, aes(x = factor(get(x), levels = x.order), group = Cluster)) +
      geom_point(aes(x = factor(get(x), levels = x.order), y = mean, color = Cluster), size=2) +
      geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Cluster), alpha = 0.3) +
      geom_line(aes(y = mean, color = Cluster), linetype = 'solid', linewidth = 0.5) +
      facet_wrap2(~factor(Cluster, levels = facet.order), scales = scale,
                  nrow = facet_nrow, strip = x_strip) +
      scale_fill_manual(values=facet.color) +
      scale_colour_manual(values=facet.color) +
      PublicationTheme() +
      theme(legend.title = element_blank(),
            plot.title = element_text(size = 16, hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none"
      ) + 
      labs(x=NULL, y='Cell Proportion')
    return(p)
  }
  
}

#' Plot mean+/-SEM bar plot with sample point
#' 
#' @description
#' Plot mean+/-SEM bar plot with sample point: filp 
#' 
#' @param cell.prop a data.frame of cell proportion
#' @param groups group of get mean +/- se
#' @param cell.ord cell order
#' @param x.order
#' @param coord.flip if flip x-y
#' @param rev.y if reverse y axis
#' @param y.lim ylim
#' @note
#' 
#' @return a list of SampleID
#' @usage 
#' 
PlotMeanBarPlot <- function(
    cell.prop,
    groups = c("Tag", "Cluster"),
    cell.ord,
    x.order = c("Pre", "Post_nCR", "Post_CR"),
    x.color = c("#69af86", "#fad6a6", "#fea477"),
    fill.by = "Tag",
    facet = NULL,
    stat.obj = NULL,
    color.by = NULL,
    coord.flip = FALSE,
    rev.y = FALSE,
    y.lim = c(0, 0.86),
    y.tick = c(0, 0.25, 0.5, 0.75),
    y.title = "Cell Proportion"
){
  
  x.name = as.name(groups[1])
  
  cell.prop <- cell.prop %>%
    mutate("{{x.name}}" := factor({{x.name}}, levels = x.order))
  
  # get mean and sem with group
  cell_proport_sem <- cell.prop %>%
    desc_statby(measure.var = "Proportion", grps = groups)
  
  if (!is.null(facet)) {
    cell.name = as.name(groups[2])
    cluster_to_coarse <- unique(cell.prop[c(facet, groups[2])]) %>% tibble::column_to_rownames(groups[2])
    cell_proport_sem <- cell_proport_sem %>%
      mutate(Facet = sapply({{cell.name}}, function(x){cluster_to_coarse[x, facet]}))
    cell.prop <- cell.prop %>% rename(Facet = as.name(facet))
    
    if (!is.null(stat.obj)) {
      stat.obj <- stat.obj %>%
        mutate(Facet = sapply({{cell.name}}, function(x){cluster_to_coarse[x, facet]}))
    }
  }
  
  darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }
  
  point_color <- sapply(x.color, darken, USE.NAMES = FALSE)
  
  # other
  if (!is.null(fill.by)) {
    p <- ggplot(cell_proport_sem, aes(x = factor(Cluster, levels = cell.ord), y = mean, group = get(fill.by))) + 
      geom_bar(aes(fill = get(fill.by)), width = 0.8, stat = "identity", position = position_dodge(0.8), alpha = 0.6) + 
      geom_point(data = cell.prop, aes(x = Cluster, y = Proportion, colour = get(fill.by)), position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),width = 0.25, position = position_dodge(0.8)) +
      scale_fill_manual(values = x.color) +
      scale_colour_manual(values = point_color)
    
    if (!is.null(facet)) {
      p <- p + facet_grid2(~ factor(Facet, levels = c("Tc", "Bc", "Mye")), scales = "free", space = "free", strip = strip_themed(background_x = elem_list_rect(fill = coarse_color)))
    }
    
    if (!is.null(stat.obj)) {
      p <- p + 
        # scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        stat_pvalue_manual(stat.obj,
                           # hide.ns = "p",
                           tip.length = 0,
                           label = "p") +
        ylim(c(0,1.2))
    }
    
  }
  
  # fig1
  if (!is.null(color.by)) {
    p <- ggplot(cell_proport_sem, aes(x = factor(Cluster, levels = cell.ord), y = mean, color = get(color.by))) + 
      geom_bar(width = 0.6, fill = NA, linewidth = 0.6, stat = "identity", position = position_dodge(0.7), alpha = 1) +
      geom_point(data = cell.prop, aes(x = Cluster, y = Proportion, fill = get(color.by)), 
                 size = 1.5, position = position_dodge(0.7), alpha = 0.8,
                 shape = 21, stroke = 0) +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                    width = 0.3, position = position_dodge(0.7)) +
      scale_fill_manual(values = point_color) +
      scale_colour_manual(values = x.color)
  }
  
  if (coord.flip) {
    if (rev.y) {
      p <- p + scale_y_reverse(limits = rev(y.lim), breaks = rev(y.tick)) +
        scale_x_discrete(position = "top")
    } else {
      p <- p + scale_y_continuous(limits = y.lim, breaks = y.tick)
    }
    p <- p + coord_flip()
  }
  
  p <- p + 
    PublicationTheme() +
    # theme_bw(base_size = 15) +
    labs(x = NULL, y = y.title) +
    theme(legend.title = element_blank(),
          axis.title = element_text(size = 12))
  
  return(p)
  
}


#' Plot nest donut plot
#'
#' @param cell.prop a data.frame of cell proportion
#' 
#' 
PlotNestDonutPlot <- function(
    cell.prop,
    cell.col,
    y.title = 'Pre Subcell Proportion'
){
  
  p <- ggplot(cell.prop, aes(x = Response, y = Proportion, fill = Cluster)) + 
    geom_col(position = 'fill', width = 0.9, color = 'white', linewidth = 0.5) +
    coord_polar(theta = "y") +
    geom_col(aes(x = -0.5, y = 0)) +
    facet_wrap(~ Major, dir = "v", ncol = 1)  +
    theme_void() + 
    labs(x = NULL, y = y.title, title = NULL) +
    scale_fill_manual(values = cell.col) +
    theme(strip.text = element_blank(), panel.spacing.y = unit(0, "pt"),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          panel.grid.major = element_blank()) +
    theme(legend.title = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          axis.text = element_blank(),
          axis.ticks.x = element_blank()
    )
  
}


# LollipopRoe
PlotLollipopRoe <- function(
    Roe_table,
    group_level = c("NR", "R"),
    group_color = c('#AFC38E','#E2AA81')
){
  
  # get plot data
  Roe_plot_data <- Roe_table %>% 
    as.data.frame() %>%
    dplyr::rename(SubCluster=Var1, Group=Var2, Roe=Freq) %>%
    dplyr::mutate(Group = factor(Group, levels = group_level),
                  SubCluster = stringr::str_split_i(SubCluster, "_", 1))
  
  # plot
  Roe_plot <- ggplot(Roe_plot_data, aes(x = SubCluster, y = Roe, group = Group)) +
    # https://stackoverflow.com/questions/17521438/geom-rect-and-alpha-does-this-work-with-hard-coded-values
    annotate("rect", xmin = 0, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = '#FEB941', alpha = 0.08) +
    annotate("rect", xmin = 4.5, xmax = 9, ymin = -Inf, ymax = Inf, fill = '#00712D', alpha = 0.08) +
    geom_hline(yintercept = 1, linetype = 2, linewidth=0.3, color = "black") +
    geom_linerange(aes(xmin = SubCluster, xmax = SubCluster, ymin=0, ymax = Roe), position = position_dodge(0.8), linewidth = 1.5, linetype=1, color="darkgrey") +
    geom_point(aes(color = Group), size = 4.5, alpha = 1, position = position_dodge(0.8)) +
    theme_light() +
    # PublicationTheme() +
    geom_text(aes(x = SubCluster, y = Roe* 1, label = round(Roe, 2)), position = position_dodge(0.8)) + 
    scale_colour_manual(values = group_color) +
    # scale_size_continuous(range = c(1, 3), breaks = c(-log10(0.01), -log10(0.05)), labels = c("p<0.01", "p>0.05"), name="pval") +
    theme(
      panel.grid.major = element_blank(),
      axis.text = element_text(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      panel.border = element_rect(fill = NA, colour = "black"),
      legend.position = 'none',
      axis.text.x = element_text(angle=0, size = 10, vjust = 0.5, hjust = 0.5)
    ) +  
    ylim(-0.1, max(Roe_plot_data$Roe)+0.1) + 
    ylab("Ro/e") + 
    xlab("") +
    ggtitle("")
  
  return(Roe_plot)
}