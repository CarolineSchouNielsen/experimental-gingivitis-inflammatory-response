phyloseq_maaslin3_CSN <- function(phyloseq,
                              formula = NULL,
                              fixed_effects = NULL,
                              reference = NULL,
                              random_effects = NULL,
                              group_effects = NULL,
                              ordered_effects = NULL,
                              strata_effects = NULL,
                              feature_specific_covariate = NULL,
                              feature_specific_covariate_name = NULL,
                              feature_specific_covariate_record = NULL,
                              min_abundance = 0,
                              min_prevalence = 0,
                              zero_threshold = 0,
                              min_variance = 0,
                              max_significance = 0.1,
                              normalization = 'TSS',
                              transform = 'LOG',
                              correction = 'BH',
                              standardize = TRUE,
                              unscaled_abundance = NULL,
                              median_comparison_abundance = TRUE,
                              median_comparison_prevalence = FALSE,
                              median_comparison_abundance_threshold = 0,
                              median_comparison_prevalence_threshold = 0,
                              subtract_median = FALSE,
                              warn_prevalence = TRUE,
                              augment = TRUE,
                              evaluate_only = NULL,
                              plot_summary_plot = TRUE,
                              summary_plot_first_n = 25,
                              coef_plot_vars = NULL,
                              heatmap_vars = NULL,
                              plot_associations = TRUE,
                              max_pngs = 30,
                              cores = 5,
                              save_models = TRUE,
                              verbosity = 'FINEST',
                              output_dir = "~/test_masslin3/",
                              add_ASV_taxonomy = TRUE, 
                              small_random_effects = FALSE){
  
  ##---------------------------------------------
  require(tidyverse); require(maaslin3); require(phyloseq)
  
  
  
  maaslin3(input_data = phyloseq %>% otu_table() %>%  t() %>% data.frame() ,
           input_metadata = phyloseq %>% sample_data()  %>% data.frame(),
           output = output_dir,
           formula = formula,
           fixed_effects = fixed_effects,
           reference = reference,
           random_effects = random_effects,
           group_effects = group_effects,
           ordered_effects = ordered_effects,
           strata_effects = strata_effects,
           feature_specific_covariate = feature_specific_covariate,
           feature_specific_covariate_name = feature_specific_covariate_name,
           feature_specific_covariate_record = feature_specific_covariate_record,
           min_abundance = min_abundance,
           min_prevalence = min_prevalence,
           zero_threshold = zero_threshold,
           min_variance = min_variance,
           max_significance = max_significance,
           normalization = normalization,
           transform = transform,
           correction = correction,
           standardize = standardize,
           unscaled_abundance = unscaled_abundance,
           median_comparison_abundance = median_comparison_abundance,
           median_comparison_prevalence = median_comparison_prevalence,
           median_comparison_abundance_threshold = median_comparison_abundance_threshold,
           median_comparison_prevalence_threshold = median_comparison_prevalence_threshold,
           subtract_median = subtract_median,
           warn_prevalence = warn_prevalence,
           augment = augment,
           evaluate_only = evaluate_only,
           plot_summary_plot = plot_summary_plot,
           summary_plot_first_n = summary_plot_first_n,
           coef_plot_vars = coef_plot_vars,
           heatmap_vars = heatmap_vars,
           plot_associations = plot_associations,
           max_pngs = max_pngs,
           cores = cores,
           save_models = save_models,
           verbosity = verbosity, 
           small_random_effects = small_random_effects) -> out
  # 
  # if(add_ASV_taxonomy == TRUE){
  #   
  #   phyloseq %>%
  #     tax_table() %>%
  #     as.data.frame()%>%
  #     rownames_to_column(var = "feature") -> tax_table
  #   
  #   left_join(out$results,
  #             tax_table, by="feature") -> out$results_ASV_tax
  #   
  #   write_tsv(out$results_ASV_tax, file=paste0(output_dir, "results_tax_info.tsv"))
  #   
  # }
  ##---------------------------------------------
  
  # gc()
  
  ##---------------------------------------------
  
  return(out)
}
