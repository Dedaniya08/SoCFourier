#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# METHYLATION ARRAY FOURIER SEASONALITY ANALYSIS
# 
# Analysis performed on two cohorts
# 1. ENID: discovery; 450k; aged 2y
# 2. EMPHASIS: replication EPIC; aged 7-9y
#
# model year-round methylation using fourier regression
# significance (goodness of fit) of seasonality effect is determined by 
# comparing model with one pair of fourier terms (FTs) to covariates-only 
# model using a likelihood ratio test (LRT). 
# See paper for further details.
#
# NOTE: analysis is restricted to loci intersecting both arrays only
#
# output saved to project directory ('proj.dir') 
# which should contain following subdirectories:
# ./cohort_data/
# ./plots/
# ./analysis/annotation_tables/
# ./analysis/enrichment_tests/
# ./analysis/plots/
# ./analysis/R_objects/
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

require(tidyverse)
require(GenomicRanges)

proj.dir = '~/projects/soc_enid_emphasis/' # output directory
src.dir = '~/src/soc_enid_emphasis/' # source code directory

source(paste0(src.dir,'SoCFourier_modelling_functions.R'))
source(paste0(src.dir,'SoCFourier_CpG_annotation_functions.R'))
source(paste0(src.dir,'SoCFourier_plot_functions.R'))
source(paste0(src.dir,'SoCFourier_stats_functions.R'))
source(paste0(src.dir,'gamete_embryo_plot_and_stats_functions.R'))

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# LOAD DATA
# load methylation data and MEs on 450k
# cgs intersecting both cohorts
cgs.array.intersect = readRDS(paste0(proj.dir,'cohort_data/all_intersecting_cgs.RDS')) 
enid.M = 
  readRDS(paste0(proj.dir,'cohort_data/ENID_Mfil_norm_filt_no_replicates.RDS'))[cgs.array.intersect,]
emph.M = 
  readRDS(paste0(proj.dir,'cohort_data/EMPH_Mfil_norm_filt_no_replicates.RDS'))[cgs.array.intersect,]
me.cgs = readRDS(paste0(proj.dir,'cohort_data/all.MEs.RDS')) # known 450k MEs


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# FIT FOURIER MODELS ----
# ENID / DISCOVERY
if (file.exists(paste0(proj.dir,'analysis/R_objects/ENID_lrt.results.PCs1-6.RDS'))){
  
  enid.lrt.results = readRDS(paste0(proj.dir,'analysis/R_objects/ENID_lrt.results.PCs1-6.RDS'))
  cat('ENID LRT results loaded\n')
  
} else {
  
  enid.M = enid.M[cgs.array.intersect,]
  enid.covs = as.data.frame(readRDS(paste0(proj.dir,'cohort_data/enid.design.covs.RDS')))

  # RUN FOURIER MODEL AND compute LRT pValues
  enid.lrt.results = generate_fourier_lrt_pvals(meth = enid.M,design = enid.covs)
  cat('Determining LRT p-values from Fourier regression models for ENID cohort...')
  saveRDS(enid.lrt.results,file = paste0(proj.dir,'analysis/R_objects/ENID_lrt.results.PCs1-6.RDS'))
  cat('done\n\n')

}

# EMPHASIS/REPLICATION
if (file.exists(paste0(proj.dir,'analysis/R_objects/EMPH_lrt.results.PCs1-6.RDS'))){
  
  emph.lrt.results = readRDS(paste0(proj.dir,'analysis/R_objects/EMPH_lrt.results.PCs1-6.RDS') )
  cat('EMPHASIS LRT results loaded\n')
  
} else {
  
  emph.M = emph.M[cgs.array.intersect,]
  emph.covs = as.data.frame(readRDS(paste0(proj.dir,'cohort_data/emph.design.covs.RDS')))

  # RUN FOURIER MODEL AND compute LRT pValues
  cat('Determining LRT p-values from Fourier regression models for EMPHASIS cohort...')
  emph.lrt.results = generate_fourier_lrt_pvals(meth = emph.M,design = emph.covs)
  saveRDS(emph.lrt.results,file = paste0(proj.dir,'analysis/R_objects/EMPH_lrt.results.PCs1-6.RDS'))
  cat('done\n\n')

}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# identify discovery CpGs and replicating CpGs with FDR<=10% ----
# DISCOVERY
discovery.cgs = filter(enid.lrt.results,lrt.fdr<=0.1)$cg
cat(length(discovery.cgs),'\'discovery cgs\' with LRT FDR <= 10% identified in ENID cohort\n')
# REPLICATION / SOC-CPGS
emph.lrt.results$disc.lrt.fdr = NA
emph.lrt.results[discovery.cgs,'disc.lrt.fdr'] = p.adjust(emph.lrt.results[discovery.cgs,'lrt.pval'],method = 'fdr')
soc.cgs = filter(emph.lrt.results,disc.lrt.fdr<=0.1)$cg
cat(length(soc.cgs),'\'SoC-CpGs\': discovery cgs with LRT FDR <= 10% identified in EMPHASIS cohort\n')

saveRDS(discovery.cgs,file = paste0(proj.dir,'analysis/R_objects/discovery_cgs.RDS'))
saveRDS(soc.cgs,file = paste0(proj.dir,'analysis/R_objects/soc_cgs.RDS'))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Generate matched and random controls ----
# 'Matched controls' have similar Beta distributions to SoC-CpGs

# 1. KS-matched controls
if (!file.exists(paste0(proj.dir,'analysis/R_objects/ks.matches.5perCG.thresh.1.RDS'))){
  
  cat('generating KS-matched controls with similar methylation Beta distributions to SoC-CpGs...')
  ks.matched.controls = generate_ks_matches(
    discovery.cgs = discovery.cgs,
    soc.cgs = soc.cgs,
    me.cgs = me.cgs,
    all.array.cgs = cgs.array.intersect,
    n.matches.per.query = 5,
    ks.thresh = 0.1
  )
  saveRDS(ks.matched.controls,file = paste0(proj.dir,'analysis/R_objects/ks.matches.5perCG.thresh.1.RDS'))
  ctrl.cgs = as.vector(unlist(ks.matched.controls))
  cat(length(ctrl.cgs),'KS-matched controls generated and saved to file\n')
  
} else {
  
  ctrl.cgs = as.vector(unlist(
    readRDS(paste0(proj.dir,'analysis/R_objects/ks.matches.5perCG.thresh.1.RDS'))
  ))
  cat(length(ctrl.cgs),'KS-matched controls loaded from file\n')
  
}

# 2. random controls
if (!file.exists(paste0(proj.dir,'analysis/R_objects/rand.cgs.RDS'))){
  
  cat('generating set of random control cgs...')
  rand.cgs = sample(
    x = setdiff(cgs.array.intersect,union(union(ctrl.cgs,discovery.cgs),me.cgs)),
    size = 625,
    replace = F
  )
  saveRDS(rand.cgs,file = paste0(proj.dir,'analysis/R_objects/rand.cgs.RDS'))
  cat('done\n\n')
  
} else {
  
  rand.cgs = readRDS(file = paste0(proj.dir,'analysis/R_objects/rand.cgs.RDS'))
  cat(length(rand.cgs),'random controls loaded from file\n')
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Generate fitted values and save modelling summary stats for SoC-CpGs and controls ----
if (!file.exists(paste0(proj.dir,'analysis/R_objects/enid.fourier.fit.RDS')) |
  !file.exists(paste0(proj.dir,'analysis/R_objects/emph.fourier.fit.RDS'))){

    
    cat('generating fitted values at loci of interest using SoC Fourier regression model... ')
    cgs.to.fit = Reduce(union,list(discovery.cgs,rand.cgs,ctrl.cgs))
    
    # ENID
    enid.design = as.data.frame(readRDS(paste0(proj.dir,'cohort_data/ENID_design.PCs1-6.RDS')))
    enid.fit = lm(t(enid.M[cgs.to.fit,])~.+0,data=enid.design,)
    
    # EMPHASIS
    emph.design = as.data.frame(readRDS(paste0(proj.dir,'cohort_data/EMPHASIS_design.PCs1-6.RDS')))
    emph.fit = lm(t(emph.M[cgs.to.fit,])~.+0,data=emph.design)
    
    enid.fitted.vals = get_fourier_fitted_vals(model.fit = enid.fit, design = enid.design)
    emph.fitted.vals = get_fourier_fitted_vals(model.fit = emph.fit, design = emph.design)
    
    cat('done\n\n')
    
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Generate summary stats from fitted values for SoC-CpGs and controls ----
# methylation maxima / minima; effect size (amplitude)
# LRT pValues and FDR

    cat('summarising SoC Fourier regression model stats... ')
    enid.model.stats = generate_model_fit_stats(
      fitted.vals = enid.fitted.vals,
      lrt.results = enid.lrt.results,
      discovery = discovery.cgs,
      ctrl = ctrl.cgs,
      rand.cgs = rand.cgs,
      soc = soc.cgs
    )
    emph.model.stats = generate_model_fit_stats(
      fitted.vals = emph.fitted.vals,
      lrt.results = emph.lrt.results,
      discovery = discovery.cgs,
      ctrl = ctrl.cgs,
      rand.cgs = rand.cgs,
      soc = soc.cgs
    )
    # add FDR computed for discovery loci
    emph.model.stats$disc.lrt.fdr = NA
    emph.model.stats[filter(emph.model.stats,locus=='discovery')$cg,'disc.lrt.fdr'] = 
      emph.lrt.results[filter(emph.model.stats,locus=='discovery')$cg,'disc.lrt.fdr']
    
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    # Save fitted values and model fit summary stats ----
    saveRDS(
      list(
        effects = enid.model.stats,
        fitted.vals = enid.fitted.vals
      ),
      file = paste0(proj.dir,'analysis/R_objects/enid.fourier.fit.RDS')
    )
    
    saveRDS(
      list(
        effects = emph.model.stats,
        fitted.vals = emph.fitted.vals
      ),
      file = paste0(proj.dir,'analysis/R_objects/emph.fourier.fit.RDS')
    )
    cat('done\n\n')

} else {
  
  enid.fit.stats = readRDS(paste0(proj.dir,'analysis/R_objects/enid.fourier.fit.RDS'))
  emph.fit.stats = readRDS(paste0(proj.dir,'analysis/R_objects/emph.fourier.fit.RDS'))
  cat('ENID and EMPHASIS fitted values and SoC Fourier regression model stats loaded from file\n')
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Create annotated array for downstream analysis ----
cat('creating array-wide annotatation object...')
annot.gr = annotate_array(
  array.cgs = cgs.array.intersect,
  enid.model.stats = enid.fit.stats$effects,
  emph.model.stats = emph.fit.stats$effects
)
saveRDS(annot.gr,file = paste0(proj.dir,'analysis/R_objects/annotated.array.gr.RDS'))
cat('done\n')

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Create annotated tables for publication ----
cat('creating annotatation tables for publication...')
annotation.tables = create_annotation_tables(annot.gr = annot.gr)
write.csv(annotation.tables$annot.table.discovery,file = paste0(proj.dir,'analysis/annotation_tables/discovery.annotation.table.csv'),row.names = F)
write.csv(annotation.tables$annot.table.replication,file = paste0(proj.dir,'analysis/annotation_tables/replication.annotation.table.csv'),row.names = F)
write.csv(annotation.tables$annot.table.replication.clustered,file = paste0(proj.dir,'analysis/annotation_tables/replication.clustered.annotation.table.csv'),row.names = F)
cat('done\n')


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Get cluster sizes ----
cat('cluster size stats:\n')
print(
  get.cluster.sizes(test.gr = annot.gr[annot.gr$discovery],test.name = 'discovery.cgs',cluster.gap = 1000)
)
print(
  get.cluster.sizes(test.gr = annot.gr[annot.gr$replicated],test.name = 'replicated.cgs',cluster.gap = 1000)
)
cat('\n')

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Summarise features overlapping discovery/replicated CpGs
summarise_overlapping_features(annot.gr = annot.gr)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# PLOTS ----

# plot size params and save directory
plt.regular = get_figure_params(
  plot.dir = paste0(proj.dir,'analysis/plots/'),
  width = 8,
  height = 6,
  dpi = 600
)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot some Beta distribution of KS-matched controls vs CpG to be matched ----
plot_ks_matches(matches = 
  readRDS(paste0(proj.dir,'analysis/R_objects/ks.matches.5perCG.thresh.1.RDS')),
  n.query.cgs.to.plot = 5,
  betas = m2beta(enid.M),
  plot.save = T,
  plot.params = plt.regular
)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot mean and median Beta distributions for SoC-CpGs vs random controls ----
plot_beta_summary_stat_vs_matched(annot.gr = annot.gr,summary.stat = 'mean',plot.save = TRUE,plot.params = plt.regular)
plot_beta_summary_stat_vs_matched(annot.gr = annot.gr,summary.stat = 'median',plot.save = TRUE,plot.params = plt.regular)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot genome-wide SoC-CpG locations on karyogram  ----
plot_cpg_locations_karyo(annot.gr = annot.gr,plot.save = T,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot fitted fourier curves ----
enid.fitted.vals = readRDS(paste0(proj.dir,'analysis/R_objects/enid.fourier.fit.RDS'))$fitted.vals
emph.fitted.vals = readRDS(paste0(proj.dir,'analysis/R_objects/emph.fourier.fit.RDS'))$fitted.vals
annot.gr = readRDS(paste0(proj.dir,'analysis/R_objects/annotated.array.gr.RDS'))

# ENID
plt = plot_SoC_Fourier_curves_discovery(
  annot.gr = annot.gr,fitted.vals = enid.fitted.vals,plot.legend = F,plot.save = TRUE,plot.filename = 'enid.fourier.pdf',plot.params = plt.regular)
print(plt)
# EMPHASIS
plt = plot_SoC_Fourier_curves_discovery(
  annot.gr = annot.gr,fitted.vals = emph.fitted.vals,plot.legend = F,plot.save = TRUE,plot.filename='emph.fourier.pdf',plot.params = plt.regular)
print(plt)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot date of concpetion at methylation maxima/minima across cohorts ----
# and compute Spearman correlation of doc between cohorts at maxima
plot_cross_cohort_doc_max(plot.save = T,plot.params = plt.regular)



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot date of concpetion at methylation maxima/minima across cohorts ----
# and compute Spearman correlation of doc between cohorts at maxima
plot_doc_maxmin(annot.gr = annot.gr,plot.max.doc.only = TRUE,plot.save = TRUE,plot.params = plt.regular)
# repeat for doc at max and min for each cohort
plot_doc_maxmin(annot.gr = annot.gr,plot.max.doc.only = FALSE,plot.save = TRUE,plot.params = plt.regular)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot seasonal amplitude (distance between methylatin max and min) ----
plot_amplitude(annot.gr = annot.gr,plot.save = TRUE,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot date of conception at meth max vs seasonal amplitude ----
plot_doc_vs_amplitude(annot.gr = annot.gr,plot.save = TRUE,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot summary methylation Beta distributions vs array background ----
plot_beta_summary_stat_vs_bgd(annot.gr = annot.gr,summary.stat = 'mean',plot.save = TRUE,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot trans-gastrulation intermediate methylation from Guo dataset ----
plot_Guo_intermediate_meth(annot.gr = annot.gr,plot.save = TRUE,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot sperm hypomethylation from Okae dataset ----
plot_Okae_sperm_hypomethylation(annot.gr = annot.gr,plot.save = TRUE,plot.params = plt.regular)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot genomic context (CGI and gene position) ----
plot_genomic_context(annot.gr = annot.gr,plot.save = TRUE,plot.params = plt.regular,verbose = T)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Plot pairwise cg correlations ----
plot_pairwise_cpg_correlations(
  annot.gr = annot.gr,
  disc.fitted.vals = enid.fitted.vals,
  repl.fitted.vals = emph.fitted.vals,
  plot.save = T,
  plot.params = plt.regular
)
# cluster adjusted
plot_pairwise_cpg_correlations(
  annot.gr = annot.gr,
  disc.fitted.vals = enid.fitted.vals,
  repl.fitted.vals = emph.fitted.vals,
  cluster.adj = T,
  cluster.size = 1000,
  plot.save = T,
  plot.params = plt.regular
)



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# STATS AND ENRICHMENT TESTS ----

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Intermediate methylation enrichment stats Supp Table 5 ----
intermediate_methylation_enrichment(annot.gr = annot.gr)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Seasonal amplitude stats Supp Table 6 ----
amplitude_stats(annot.gr = annot.gr)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Enrichment tests Supp Table 7  ----
# (MEs, Teh GxE, ZFP57, TRIM28, CTCF, ERV1, ERVK, ERVL-MaLR)
enrichment_tests(
  annot.gr = annot.gr,
  output.to.file = TRUE,
  sink.file = paste0(proj.dir,'analysis/enrichment_tests/enrichment_stats.txt'),
  verbose.output = TRUE
)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Enrichment tests Supp Table 8  ----
# (Zink PofOm, oo.gDMRs and sp.gDMRs)
enrichment_tests_PofO(
  annot.gr = annot.gr,
  output.to.file = T,
  sink.file = paste0(proj.dir,'analysis/enrichment_tests/enrichment_stats_PofO.txt'),
  verbose
)
