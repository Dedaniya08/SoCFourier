#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# mQTL and GxE analysis using GEM
#
# Analysis performed on data from EMPHASIS cohort
# using imputed GSA genotype data
#
# 1. find G1 and G2 SNPs using GEM for 
#   a) SoC-CpGs
#   b) 625 KS-matched control CpGs
#   c) 625 random background control CpGs
#
#   GEM models:
#     E environment only
#     G genetic (mQTL) main effects not including environment in model
#     GxE gene-environment interactions; genetic and envioronmental main effects 
#     included in model
# 
#   all GEM models use single, most significant E-term (cos or sin) only
# 
# 2. compare adjusted R^2 for all G1 and G2 models and find 'winning' model using AIC
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

proj.dir = '~/projects/soc_enid_emphasis/' # output directory
src.dir = '~/src/soc_enid_emphasis/' # source code directory

source(paste0(src.dr,'GEM_analysis_functions.R'))
source(paste0(src.dir,'GEM_plot_functions.R'))

#^^^^^^^^^^^^^^^^^^^^^^^^^^
# GEM ANALYSIS

# Configure PLINK files ----
configuire_PLINK_files_for_GEM()

# Prepare text files for GEM analysis ----
prepare_GEM_data(annot.gr = readRDS(paste0(proj.dir,'analysis/R_objects/annotated.array.gr.RDS')))

# run GEM analysis ----
run_GEM_analysis()

# create GEM summary results tables ----
create_GEM_summary_results_tables()

# compare E, G and GxE models with winning SNPs identified by GEM ----
compare_E_G_GxE_models()

# get G1 and G2 SNP locations ----
get_G1_G2_locations()

# generate G1 and G2 annotation tables ----
generate_G1_G2_annotation_tables()

# characterise G1 and G2 SNPs ----
characterise_G1_G2_SNPs(snp.clust.dist = 5e6)


#^^^^^^^^^^^^^^^^^^^^^^^^^^
# GEM PLOTS
setwd(paste0(proj.dir,'analysis/emphasis_GEM/'))
gem.soc.res = readRDS('./R_objects/gem.soc.res.RDS')
gem.ctrl.res = readRDS('./R_objects/gem.ctrl.res.RDS')
gem.rnd.res = readRDS('./R_objects/gem.rnd.res.RDS')

# get plot dimensions etc ----
plt.params = get_figure_params(plot.dir = paste0(proj.dir,'analysis/emphasis_GEM/plots/'))

# % winning model pie charts ----
winning_model_pie_charts(soc.res = gem.soc.res,ctrl.res = gem.ctrl.res,rnd.res = gem.rnd.res,plot.params = plt.params)

# delta adjusted R2 barplots ----
delta_adj_r2_barplots(soc.res = gem.soc.res,ctrl.res = gem.ctrl.res,rnd.res = gem.rnd.res,plot.params = plt.params)

# SoC Fourier plots stratfied by genotype ----
plot_soc_fourier_stratified_by_genotype(
  soc.res = gem.soc.res,plot.params = plt.params,log.to.file = T,
  sink.file = paste0(proj.dir,'~/projects/soc_enid_emphasis/analysis/emphasis_GEM/plots/stratified_fourier_plots/strat_plots_stats.log')
)

# plot G1 and G2 SNP locations in karyogram
plot_G1_G2_locations_karyogram()
