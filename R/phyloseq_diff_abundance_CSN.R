# Modified from Florentin Constancias' original function found here: https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R
# I have added five MaAsLin3 optional parameters to the function: 
# Line 45. maaslin3_summary_plot_first_n. Description from the official [MaAsLin3 tutorial](https://github.com/biobakery/maaslin3): summary_plot_first_n (default 25): Include the top summary_plot_first_n features with significant associations. 
# Line 46. maaslin3_small_random_effects. Description from the official [MaAsLin3 tutorial](https://github.com/biobakery/maaslin3): small_random_effects (default FALSE): Automatically replace random effects with fixed effects in the logistic prevalence model to handle low numbers of observations per group.
# Line 47. maaslin3_heatmap_vars. Description from the official [MaAsLin3 tutorial](https://github.com/biobakery/maaslin3): heatmap_vars: Vector of variable names to be used in the heatmap section of the summary plot. Continuous variables should match the metadata column name, and categorical variables should be of the form "[variable] [level]".
# Line 48. maaslin3_reference. Description from the official [MaAsLin3 tutorial](https://github.com/biobakery/maaslin3): reference: For a variable with more than two levels supplied with fixed_effects, the factor to use as a reference provided as a string of 'variable,reference' semi-colon delimited for multiple variables.
# Line 49. maaslin3_min_prevalence. Description from the official [MaAsLin3 tutorial](https://github.com/biobakery/maaslin3): min_abundance (default 0): Features with abundances of more than min_abundance in more than min_prevalence of the samples will be included for analysis. The threshold is applied after normalization and before transformation.
# OBS! For these parameters to be active, it is also necessary to source the "phyloseq_maaslin3_CSN.R" function!
phyloseq_diff_abundance <- function(ps_tmp = ps_up %>%  subset_samples(Sample == "Saliva"),
                                    approach = c("run_lefse", 
                                                 "run_ancom",
                                                 "ancombc2",
                                                 "maaslin3",
                                                 "trans_diff_rf",
                                                 "classifier_rf"),
                                    glom = "Species",
                                    unclassified_name = "UNCLASSIFIED",
                                    taxa_rank = "all", #"OTU"
                                    density = "Quant",
                                    comp_group = "Time",
                                    
                                    plot_top_sig = NULL,
                                    
                                    palette = time_pal,
                                    pvalue_cutoff = 0.05,
                                    p_adjust = "BH",
                                    
                                    
                                    lefse_taxa_rank = taxa_rank,
                                    
                                    lefse_prv_cut = 0.2,
                                    lefse_sample_min = 10,
                                    lefse_lda_cutoff = 2,
                                    lefse_multigrp_strat = FALSE,
                                    lefse_strict = "0",
                                    lefse_bootstrap_n = 30,
                                    
                                    maaslin3_formula =  formula,
                                    maaslin3_fixed_effects = NULL,
                                    maaslin3_strata_effects = NULL,
                                    maaslin3_random_effects = NULL, 
                                    masslin3_output_dir = "~/test_masslin3/",
                                    maaslin3_prv_cut = 0,
                                    maaslin3_lefse_plot = TRUE,
                                    maaslin3_summary_plot_first_n = 25,
                                    maaslin3_small_random_effects=FALSE, 
                                    maaslin3_heatmap_vars = NULL,
                                    maaslin3_reference = NULL, 
                                    maaslin3_min_prevalence = 0, 
                                    
                                    formula =  "~ Time +  (1|Subject)",
                                    
                                    ancom_prv_cut = 0.1,
                                    ancom_confounders = character(0), # Subject
                                    
                                    ancombc2_forumla = formula,
                                    ancombc2_fix_formula = "Time",
                                    ancombc2_rand_formula = NULL,
                                    ancombc2_group = "Time",
                                    ancombc2_pairwise = TRUE,
                                    ancombc2_dunnet = TRUE,
                                    ancombc2_trend = TRUE,
                                    ancombc2_global = TRUE,
                                    ancombc2_lib_cut = 1000, # ADD TO FLO'S FUNCTION ON GITHUB
                                    
                                    
                                    linda_prv_cut = 0.1,
                                    linda_formula = formula,
                                    linda_comp_group = comp_group,
                                    
                                    rf_prv_cut = 0.1,
                                    rf_prop_train = 3/4,
                                    
                                    ref_train_max_mtry=3,
                                    ref_train_ntree = c(100, 500, 1000),
                                    feature_imp_nrep = 99,
                                    
                                    trans_diff_rf_MeanDecreaseGini_cutoff = 0){
  # To improve: store in list output of test so we can plot those in a loop or laply manner
  #
  # - Lefse
  # - ANCOMBC (en count et pas %)
  # - linda (count? percent?)
  # - maaslin3 (avec lefse plot) - check maaslin3 et low replication
  # 
  # https://microbiome.github.io/course_2022_radboud/differential-abundance-analysis-demo.html
  # Boxplot of the significant taxa + heatmap + vulcano?
  # Ideallement un modele global Sample Type + Time + 1:Subject
  # puis lefse sur les significant pour les pairwise comp.
  
  
  ########## ----- Source function and load packages
  require(microbiomeMarker);require(tidyverse);require(magrittr);require(ANCOMBC);require(microeco)
  
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_microeco.R")  
  
  suppressMessages(suppressWarnings({
    
    ########## ----- 
    
    out <- list()
    sumfilter = 0
    
    ps_tmp %>% 
      filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp
    
    # ps_tmp %>%
    # microbiome::core(detection = 0, prevalence = prv_cut) #-> ps_filtered
    
    # ps_tmp %>%
    # filter_taxa(function(x){sum(x > sumfilter) >  prv_cut*nsamples(ps_tmp)}, prune = TRUE) -> ps_tmp
    
    ########## ----- Agglomerate Taxa
    
    if(!is.null(glom))
    {
      ps_tmp %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = glom, rename_ASV = glom) -> ps_tmp
    }
    
    ########## ----- Prepare count object
    
    if(!is.null(density))
    {
      ps_tmp %>%  
        phyloseq_density_normalize(value_idx = density) -> ps_count
    }
    
    ########## ----- 
    
    if ("run_lefse"  %in% approach)
    {
      ########## ----- microbiomeMarker::run_lefse see:  https://github.com/yiluheihei/microbiomeMarker/issues/95
      # increase permutation to 999
      
      ps_tmp %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) >  lefse_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        microbiomeMarker::run_lefse(group = comp_group, 
                                    norm = "CPM",
                                    wilcoxon_cutoff = pvalue_cutoff,
                                    kw_cutoff = pvalue_cutoff,
                                    taxa_rank = taxa_rank, 
                                    sample_min = lefse_sample_min,
                                    multigrp_strat = lefse_multigrp_strat,
                                    lda_cutoff = lefse_lda_cutoff,
                                    bootstrap_n =  lefse_bootstrap_n,
                                    strict = lefse_strict) -> out$mmlefse
      
      microbiomeMarker::plot_abundance(out$mmlefse, group = comp_group) -> out$mmlefse_abp
      
      # out$mmlefse_abp$data %>% 
      #   mutate(prop = abd / 1e+6) %>% 
      #   ggplot()
      
      # reoder plot so it is in the same order of lefse see: https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
      
      microbiomeMarker::plot_ef_bar(out$mmlefse) + 
        scale_fill_manual(values = palette, na.value = "grey10") -> out$mmlefse_p
      
      microbiomeMarker::marker_table(out$mmlefse) %>%
        data.frame() -> out$mmlefse_df
      
      # save plot, mmlesfse and mmlefse_df
    }
    
    if ("lefse"  %in% approach)
    {
      ########## ----- microeco lefse https://chiliubio.github.io/microeco_tutorial/model-based-class.html
      
      ps_tmp %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) >  lefse_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      t1 <- trans_diff$new(dataset = data, 
                           lefse_min_subsam = lefse_sample_min, 
                           lefse_sub_strict = lefse_multigrp_strat, 
                           lefse_sub_alpha = NULL, 
                           lefse_bootstrap_n= lefse_bootstrap_n, 
                           boots = lefse_bootstrap_n,
                           group_choose_paired = NULL,
                           method = "lefse",
                           filter_thres = 0, #default 0; the abundance threshold,
                           alpha = pvalue_cutoff,
                           group = comp_group, 
                           lefse_taxa_rank = taxa_rank, 
                           p_adjust_method = p_adjust,
                           plot_pal = palette,
                           n.cores = 4,
                           lefse_subgroup = NULL)
      
      out$lefse$res_diff  <-  t1$res_diff
      
      if(is.null(plot_top_sig)){
        t1$res_diff %>%  nrow() -> plot_top_sig
        # print("toto")
      }
      
      
      # out$lefse$res_diff <- t1$res_diff
      
      out$lefse$g1 <- t1$plot_diff_bar(use_number = 1:plot_top_sig, color_values = palette, threshold = lefse_lda_cutoff)
      # plot the abundance using same taxa in g
      out$lefse$g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa, plot_type = "barerrorbar",
                                         add_sig = TRUE, errorbar_addpoint = FALSE, 
                                         errorbar_color_black = TRUE, 
                                         color_values = palette) # barerrorbar  ggboxplot
      # now the y axis in g1 and g2 is same, so we can merge them
      # remove g1 legend; remove g2 y axis text and ticks
      out$lefse$g1 <- out$lefse$g1 + theme(legend.position = "none")
      out$lefse$g2 <- out$lefse$g2 + theme(axis.text.y = element_blank(), 
                                           axis.ticks.y = element_blank(), 
                                           panel.border = element_blank()) + scale_y_continuous(trans='sqrt')
      out$lefse$p <- out$lefse$g1 %>% aplot::insert_right(out$lefse$g2 )
      
      out$lefse$all <- clone(t1)
      
      
    }
    
    
    
    if ("run_ancom"  %in% approach)
    {
      
      ########## ----- microbiomeMarker::run_ancom 
      
      ps_count %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        filter_taxa(function(x){sum(x > 0) >  ancom_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>%   
        microbiomeMarker::run_ancom(group = comp_group, 
                                    pvalue_cutoff = pvalue_cutoff,
                                    taxa_rank = taxa_rank, 
                                    confounders = ancom_confounders,
                                    p_adjust = p_adjust) -> out$mmancom
      
      microbiomeMarker::plot_ef_bar(out$mmancom) + 
        scale_fill_manual(values = palette, na.value = "grey10") -> out$mmancom_p
      
      microbiomeMarker::marker_table(out$mmancom) %>%
        data.frame() -> out$mmancom_df
      
      microbiomeMarker::plot_abundance(out$mmancom, group = comp_group) -> out$mmancom_abp
      
      # out$mmlefse_abp$data %>% 
      #   mutate(prop = abd / 1e+6) %>% 
      #   ggplot()
      
      # reoder plot so it is in the same order of lefse see: https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
      
    }
    
    ########## ----- microbiomeMarker::run_ancombc see: https://github.com/yiluheihei/microbiomeMarker/issues/97
    
    # ps_count %>%
    #   physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
    #   subset_taxa(Kingdom != unclassified_name) %>% 
    #   # subset_taxa(Order != "unassigned") %>% 
    #   microbiomeMarker::run_ancombc(group = comp_group, pvalue_cutoff = pvalue_cutoff,
    #                               taxa_rank = taxa_rank, 
    #                               confounders = ancom_confounders,
    #                               p_adjust = p_adjust) -> mmancombc
    # 
    # microbiomeMarker::plot_ef_bar(mmancombc) + 
    #   scale_fill_manual(values = palette, na.value = "grey10") -> mmancombc_p
    # 
    # microbiomeMarker::marker_table(mmancombc) %>%
    #   data.frame() -> mmancombc_df
    
    
    ########## ----- microbiomeMarker::run_ancombc see: https://github.com/yiluheihei/microbiomeMarker/issues/97
    
    # ps_count %>%
    #   physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
    #   subset_taxa(Kingdom != unclassified_name) %>% 
    #   # subset_taxa(Order != "unassigned") %>% 
    #   microbiomeMarker::run_ancombc(group = comp_group, pvalue_cutoff = pvalue_cutoff,
    #                               taxa_rank = taxa_rank, 
    #                               confounders = ancom_confounders,
    #                               p_adjust = p_adjust) -> mmancombc
    # 
    # microbiomeMarker::plot_ef_bar(mmancombc) + 
    #   scale_fill_manual(values = palette, na.value = "grey10") -> mmancombc_p
    # 
    # microbiomeMarker::marker_table(mmancombc) %>%
    #   data.frame() -> mmancombc_df
    
    if ("ancombc2"  %in% approach)
    {
      ########## ----- ancombc2 https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
      
      ps_count %>% # accepts count & 
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = ifelse(taxa_rank == "all", "Species", taxa_rank), rename_ASV = ifelse(taxa_rank == "all", "Species", taxa_rank)) %>%
        ANCOMBC::ancombc2(data = ., tax_level = ifelse(taxa_rank == "all", "Species", taxa_rank), 
                          group = ancombc2_group , pairwise = ancombc2_pairwise, #group = "Time",
                          fix_formula = ancombc2_fix_formula, 
                          rand_formula = ancombc2_rand_formula, 
                          p_adj_method = p_adjust, prv_cut = 0, lib_cut = ancombc2_lib_cut, # ADD TO FLO'S FUNCTION ON GITHUB
                          struc_zero = TRUE, neg_lb = TRUE, alpha = pvalue_cutoff, 
                          n_cl = 2, verbose = TRUE,
                          dunnet = ancombc2_dunnet, trend = ancombc2_trend,global = ancombc2_global,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                          trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                      nrow = 2, 
                                                                      byrow = TRUE)),
                                               node = list(2),
                                               solver = "ECOS",
                                               B = 10)) -> out$ancombc2
      # ancombc2$res %>%
      #   mutate_if(is.numeric, function(x) round(x, 2))
      # 
      # out$ancombc2$res_global %>%
      # mutate_if(is.numeric, function(x) round(x, 2)) %>%
      # dplyr::filter(passed_ss == "TRUE" & diff_abn == "TRUE") -> out$ancombc2global_signif
      
      # out$ancombc2$res_pair %>%
      # mutate_if(is.numeric, function(x) round(x, 2))
      
      # out$ancombc2$res_trend %>%
      # mutate_if(is.numeric, function(x) round(x, 2)) %>%
      # filter(diff_abn == TRUE & passed_ss == TRUE) -> out$ancombc2trend_signif
      
    }
    
    if ("maaslin3"  %in% approach)
    {
      
      ########## ----- maaslin3 see https://github.com/biobakery/biobakery/wiki/MaAsLin3#44-level-contrasts
      
      # By default, the (up to) two metadata variables with the most significant associations will be plotted in the coefficient plot, and the rest will be plotted in the heatmap. Because predicting the output variable names can be tricky, it is recommended to first run maaslin3 without setting coef_plot_vars or heatmap_vars, look at the names of the variables in the summary plot, and then rerun with maaslin_plot_results_from_output after updating coef_plot_vars and heatmap_vars with the desired variables.
      
      ps_tmp %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > maaslin3_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = ifelse(taxa_rank == "all", "Species", taxa_rank), rename_ASV = ifelse(taxa_rank == "all", "Species", taxa_rank)) %>%
        #phyloseq_maaslin3_CSN(formula = maaslin3_formula,
        phyloseq_maaslin3(formula = maaslin3_formula,
                          fixed_effects =  maaslin3_fixed_effects, 
                          strata_effects = maaslin3_strata_effects,
                          random_effects = maaslin3_random_effects, 
                          correction = p_adjust,
                          max_significance = pvalue_cutoff,
                          output_dir = masslin3_output_dir,
                          min_prevalence = maaslin3_min_prevalence,
                          # coef_plot_vars = TRUE,
                          augment = TRUE,
                          plot_associations = TRUE, 
                          save_models = TRUE, 
                          plot_summary_plot = TRUE,
                          cores = 2, 
                          summary_plot_first_n = maaslin3_summary_plot_first_n, 
                          small_random_effects = maaslin3_small_random_effects, 
                          heatmap_vars = maaslin3_heatmap_vars, 
                          verbosity = "ERROR", 
                          reference = maaslin3_reference) -> out$maaslin3 # random_effects = , group_effects = )
      
      # Save results in LEfSe format
      # maaslin_write_results_lefse_format(masslin3_output_dir, maaslin3$fit_data_abundance, maaslin3$fit_data_prevalence)
      # maaslin3::maaslin_plot_results_from_output()
      
      maaslin3:::preprocess_merged_results(rbind(out$maaslin3$fit_data_prevalence$results,
                                                 out$maaslin3$fit_data_abundance$results)) -> out$maaslin3$merged_res
      if(maaslin3_lefse_plot){
        run_make_coef_plot(merged_results_sig = out$maaslin3$merged_res %>% filter(qval_joint <= pvalue_cutoff), 
                           max_significance = pvalue_cutoff, 
                           class = comp_group) -> out$maaslin3$maaslin3_lefse_plot 
      }
      # maaslin3::maaslin_contrast_test()
      
    }
    
    if ("linda"  %in% approach)
    {
      ########## ----- microeco::linda
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > linda_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      
      t1 <- trans_diff$new(dataset = data, 
                           method = "linda", 
                           # formula = linda_formula,
                           remove_unknown = TRUE,
                           alpha = pvalue_cutoff, 
                           taxa_level = ifelse(taxa_rank == "all", "Species",taxa_rank),
                           group = linda_comp_group,
                           filter_thres = 0,
                           p_adjust_method = p_adjust)
      
      t1$res_diff -> out$linda$res_diff
      
      t1$res_diff  %<>%   subset(P.adj <= 0.05) # subset(Significance %in% c("*","**","***"))
      
      out$linda$diff_bar  <- t1$plot_diff_bar(keep_full_name = FALSE, 
                                              heatmap_cell =  "P.adj",
                                              heatmap_sig = "Significance",
                                              heatmap_x = "Factors",
                                              heatmap_y = "Taxa",
                                              heatmap_lab_fill = "P.adj")
      
      
      out$linda$all <-  clone(t1)
      
      # uses CLR?
      # using directly the linda package?
      
    }
    
    if ("classifier_rf"  %in% approach)
    {
      ########## ----- microeco::rf
      # or https://readingradio.github.io/J.nigra.Rmds/RF.GMW.Jnigra.html
      # error Error in serialize(data, node$con) : connection is not open
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > rf_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      # initialize: use "genotype" as response variable
      # x.predictors parameter is used to select the taxa; here we use all the taxa data in d1$taxa_abund
      t1 <- trans_classifier$new(dataset = data, y.response = comp_group, x.predictors = taxa_rank,  n.cores = 4)
      
      # generate train and test set
      t1$cal_split(prop.train = rf_prop_train)
      
      # Before training the model, we run the set_trainControl to invoke the trainControl function of caret package to generate the parameters used for training. 
      #Here we use the default parameters in trainControl function.
      t1$set_trainControl(method = "repeatedcv",
                          classProbs = TRUE,
                          savePredictions = TRUE)
      
      t1$cal_feature_sel(
        boruta.maxRuns = 300,
        boruta.pValue = 0.05,
        boruta.repetitions = 4)
      
      # use default parameter method = "rf"
      # require(doParallel)
      # library(caret)
      # library(randomForest)
      n <- parallel::detectCores()/2 # experiment!
      cl <- parallel::makeCluster(n)
      doParallel::registerDoParallel(cl)
      
      t1$cal_train(method = "rf", max.mtry = ref_train_max_mtry, ntree = ref_train_ntree)
      
      # t1$
      
      # t1$cal_caretList()
      
      # t1$cal_caretList_resamples()
      
      t1$cal_predict()
      
      out$classifier_rf$res_train <- t1$res_train
      
      # plot the confusionMatrix to check out the performance
      
      out$classifier_rf$res_confusion_stats <- t1$res_confusion_stats
      out$classifier_rf$res_confusion_fit <- t1$res_confusion_fit
      
      t1$plot_confusionMatrix()
      # t1$plot_confusion()
      
      t1$cal_ROC()
      #Using cal_ROC and plot_ROC can get the ROC (Receiver Operator Characteristic) curve.
      # out$Specificitysensitivity() <- t1$res_ROC$res_roc
      # out$RecallPrecision() <- t1$res_ROC$res_pr
      
      out$classifier_rf$plotROC  <- t1$plot_ROC(plot_method = FALSE)
      
      out$classifier_rf$plotROC2  <- t1$plot_ROC(plot_method = FALSE, plot_type = "PR")
      
      # default all groups
      #t1$plot_ROC(size = 0.5, alpha = 0.7)
      
      
      # default method in caret package without significance
      # t1$cal_feature_imp()
      
      # out$res_feature_imp <-  t1$res_feature_imp()
      
      # generate significance with rfPermute package
      t1$cal_feature_imp(rf_feature_sig = TRUE,group = comp_group,  num.rep = feature_imp_nrep, n.cores = 4) #, 
      
      out$classifier_rf$res_feature_imp <- t1$res_feature_imp
      
      out$classifier_rf$plot_feature_imp1 <- t1$plot_feature_imp(colour = "red", fill = "red", width = 0.6)
      
      # add_sig = TRUE: add significance label
      
      out$classifier_rf$plot_feature_imp2 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "red", fill = "red", width = 0.6, add_sig = TRUE)
      
      out$classifier_rf$plot_feature_imp3 <- t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
      
      out$classifier_rf$all <- clone(t1)
      
    } 
    
    if ("trans_diff_rf"  %in% approach)
    {
      ########## ----- microeco::rf
      # https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_diff-class
      # The ‘rf’ method depends on the random forest(Beck and Foster 2014; Yatsunenko et al. 2012) and the non-parametric test. The current method implements random forest by bootstrapping like the operation in LEfSe and employs the significant features as input. MeanDecreaseGini is selected as the indicator value in the analysis.
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > rf_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      t1 <- trans_diff$new(dataset = data, 
                           method = "rf",
                           filter_thres = 0, #default 0; the abundance threshold,
                           alpha = pvalue_cutoff,
                           group = comp_group, 
                           taxa_level = taxa_rank, 
                           p_adjust_method = p_adjust,
                           plot_pal = palette,
                           nresam = rf_prop_train, #boots = 999, 
                           n.cores = 4)
      
      t1$res_diff %>% 
        dplyr::filter(MeanDecreaseGini >= trans_diff_rf_MeanDecreaseGini_cutoff) -> t1$res_diff
      
      
      if(is.null(plot_top_sig)){
        t1$res_diff %>%  nrow() -> plot_top_sig
      }
      
      
      out$trans_diff_rf$res_diff <- t1$res_diff
      # plot the MeanDecreaseGini bar
      # group_order is designed to sort the groups
      out$trans_diff_rf$g1 <- t1$plot_diff_bar(use_number = 1:plot_top_sig, color_values = palette)
      # plot the abundance using same taxa in g
      out$trans_diff_rf$g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = TRUE, errorbar_addpoint = FALSE, errorbar_color_black = TRUE, color_values = palette) # barerrorbar  ggboxplot
      # now the y axis in g1 and g2 is same, so we can merge them
      # remove g1 legend; remove g2 y axis text and ticks
      out$trans_diff_rf$g1 <- out$trans_diff_rf$g1 + theme(legend.position = "none")
      out$trans_diff_rf$g2 <- out$trans_diff_rf$g2 + theme(axis.text.y = element_blank(), 
                                                           axis.ticks.y = element_blank(), 
                                                           panel.border = element_blank()) + scale_y_continuous(trans='sqrt')
      out$trans_diff_rf$p <- out$trans_diff_rf$g1 %>% aplot::insert_right(out$trans_diff_rf$g2 )
      
      out$trans_diff_rf$all <- clone(t1)
      
    } 
    
    
  }))
  
  return(out)
}
