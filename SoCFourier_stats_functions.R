#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Stats and enrichment tests for loci identified in
# the SoC discovery-relication analysis.
# These are called from SoCFourier_main_analysis.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' METHYLATION AMPLITUDE STATS FOR PAPER
#'
#' @param annot.gr 
#'
#' @return NULL
amplitude_stats = function(annot.gr=NULL){
  
  require(tidyverse)
  
  get_amplitude_stats = function(data=NULL){
    
  stats.df = data.frame(
    cohort=character(),
    test.locus=character(),
    control.locus=character(),
    wilcox.p=numeric()
  )
    
  # 1. compare replicating loci and MEs vs controls in each cohort
  for (cohort in c('ENID','EMPHASIS')){
    
    for (control.locus in c('random\ncontrols','matched\ncontrols')){
      
      for (test.locus in c('SoC-CpGs\nMEs','SoC-CpGs\nnon-MEs','ME\nnot replicated')){
        
        test.data = data[data$cohort==cohort,]
        x = test.data[test.data$locus==test.locus,'amplitude']
        y = test.data[test.data$locus==control.locus,'amplitude']
        test = wilcox.test(
          x,y,
          alternative = 'greater'
        )
        stats.df = rbind(stats.df,
         data.frame(
           cohort=cohort,
           test.locus=test.locus,
           control.locus=control.locus,
           median.change = 100*(median(x) - median(y)),
           wilcox.p=test$p.value
         )
        )
      }
    }
  }
  
  # 2. analyse cross-cohort changes at each locus type
  for (test.locus in c('SoC-CpGs\nMEs','SoC-CpGs\nnon-MEs','ME\nnot replicated','random\ncontrols','matched\ncontrols')){
    
    test.data = data[data$locus==test.locus,]
    x = test.data[test.data$cohort=='ENID','amplitude']
    y = test.data[test.data$cohort=='EMPHASIS','amplitude']
    test = wilcox.test(
      x,y,
      alternative = 'greater'
    )
    stats.df = rbind(stats.df,
       data.frame(
         cohort='both',
         test.locus=test.locus,
         control.locus='inter-cohort change',
         median.change = 100*(median(x) - median(y)),
         wilcox.p=test$p.value
       )
    )
  }
  
  return(stats.df)
    
  }
  
  
  amplitudes.df = as.data.frame(mcols(annot.gr)) %>%
    filter(replicated | ctrl | rand | (ME & discovery)) %>% # replicated, controls and MEs in discovery set
    select(cpg,ME,discovery,replicated,rand,ctrl,enid.amplitude,emph.amplitude,replicated) %>%
    gather(enid.amplitude,emph.amplitude,
           key = 'cohort',
           value = 'amplitude'
    )
  amplitudes.df$cpg = as.character(amplitudes.df$cpg)
  amplitudes.df[grepl('enid',amplitudes.df$cohort),'cohort'] = 'ENID'
  amplitudes.df[grepl('emph',amplitudes.df$cohort),'cohort'] = 'EMPHASIS'
  amplitudes.df[amplitudes.df$discovery & amplitudes.df$ME,'locus'] = 'ME\nnot replicated'
  amplitudes.df[amplitudes.df$replicated & amplitudes.df$ME,'locus'] = 'SoC-CpGs\nMEs'
  amplitudes.df[amplitudes.df$replicated & !amplitudes.df$ME,'locus'] = 'SoC-CpGs\nnon-MEs'
  amplitudes.df[amplitudes.df$rand,'locus'] = 'random\ncontrols'
  amplitudes.df[amplitudes.df$ctrl,'locus'] = 'matched\ncontrols'
  # table(amplitudes.df$locus)/2
  
  cat('WILCOX TESTs FOR AMPLITUDE DIFFERENCES\n')
  stats.df = get_amplitude_stats(data = amplitudes.df)
  stats.df$test.locus = recode(stats.df$test.locus,
       'SoC-CpGs\nMEs' = 'SoC-CpGs_MEs',
       'SoC-CpGs\nnon-MEs' = 'SoC-CpGs_non-MEs',
       'ME\nnot replicated' = 'ME_not_replicated',
       'random\ncontrols' = 'random_controls',
       'matched\ncontrols' = 'matched_controls'                           
  )
  stats.df$control.locus = recode(stats.df$control.locus,
       'random\ncontrols' = 'random_controls',
       'matched\ncontrols' = 'matched_controls',
       'inter-cohort change' = 'inter-cohort_change'
  )
  stats.df$median.change = round(stats.df$median.change,digits = 2)
  
  print(stats.df)

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Enrichment tests replicated loci and controls Supp Table 7
#' 
#' Test for enrichment of MEs, Teh GxE TFBS (CTCF, ZFP57, TRIM28),
#' and TEs (ERV1, ERVK, ERVL, ERVL-MaLR)
#'
#' @param annot.gr 
#' @param sink.file string file name with path to save verbose output
#'
#' @return
enrichment_tests = function(annot.gr=annot.gr,output.to.file=FALSE,sink.file=NULL,verbose.output=FALSE){

  require(tidyverse)
  require(GenomicRanges)
  
  # PARAMS
  feat.max.gap = 0 # max distance between cg and feature of interest (ME/GxE)
  tfbs.max.gap = 10000 # max distance between cg and TFBS
  erv.max.gap = 10000 # max distance between cg and TE
  cluster.gap = 1000 # inter-cg distance to define a cluster (for cluster-adjusted tests)
  include.erv = TRUE # test ERV enrichment?
  verbose <<- verbose.output

  if (output.to.file) sink(
    file = sink.file
  )
  
  # TFBS: these are loci tested for overlap within tfbs.max.gap
  # CTCF
  ctcf.bed = read.table('/data/hg19_annots/TFBS/regions/ENCFF002CDS.bed')[,1:3]
  colnames(ctcf.bed) = c('chr','start','end')
  ctcf.gr = GenomicRanges::makeGRangesFromDataFrame(ctcf.bed)
  cat(length(ctcf.gr),'ctcf features identified\n')
  
  # ZFP57
  zfp57.bed = read.table('/data/hg19_annots/ZFP57/regions/GSE78099.ZFP57.both_reps_union.bed')[,1:3]
  colnames(zfp57.bed) = c('chr','start','end')
  zfp57.gr = GenomicRanges::makeGRangesFromDataFrame(zfp57.bed)
  cat(length(zfp57.gr),'zfp57 features identified\n')
  
  # TRIM28
  trim28.bed = read.table('/data/hg19_annots/ChIP-seq/regions/ENCFF002CRN.bed')[,1:3]
  colnames(trim28.bed) = c('chr','start','end')
  trim28.gr = GenomicRanges::makeGRangesFromDataFrame(trim28.bed)
  cat(length(trim28.gr),'trim28 features identified\n')
  
  if (include.erv){
    cat('\n')
    te.bed = vector()
    for (bed.file in list.files('/data/hg19_annots/RepeatMasker/regions/')){
      if(grepl('ERV',bed.file)){
        if (!grepl('\\?',bed.file) & !grepl('ERV\\.',bed.file)){
          cat(bed.file,':\n')
          bed = read.table(paste0('/data/hg19_annots/RepeatMasker/regions/',bed.file))[,1:3]
          colnames(bed) = c('chr','start','end')
          bed.gr = GenomicRanges::makeGRangesFromDataFrame(bed)
          te.bed = c(te.bed,strsplit(bed.file,'.bed')[[1]])
          cat(length(bed.gr),'features identified\n')
          assign(gsub('.bed','.gr',bed.file),bed.gr)
        }
      }
    }
  }
  
  cat('\n**********************************************\n')
  cat('SoC ENRICHMENT ANALYSIS\n')
  cat('**********************************************\n\n')
  
  cat('max.gap =',feat.max.gap,'\n')
  cat('tfbs.max.gap =',tfbs.max.gap,'\n')
  cat('erv.max.gap =',erv.max.gap,'\n\n')
  
  # annotation data incl amplitudes, DoCs etc
  annot.gr$array.background = TRUE # background control analysis
  
  # enrichment tests
  all.res = data.frame(
    feat.set = character(),
    test.set = character(),
    n.test.cgs = numeric(),
    OR.cg = numeric(),
    p.value.cg = numeric(),
    n.test.cgs.clust = numeric(),
    OR.clust = numeric(),
    p.value.clust = numeric(),
    stringsAsFactors = F
  )
  
  # different methods depending for directly overlapping and proximal overlaps
  # 'FEATURES' are tested for DIRECT overlap (MEs, gDMRs, PofOm, teh gxe)
  do_FET_features = function(test.set=NULL,feat.set=NULL,maxgap=NULL,featname=NULL,testname=NULL){
    return(
      FET_enrichment(
        test.gr = annot.gr[mcols(annot.gr)[[test.set]]],
        feature.gr = annot.gr[mcols(annot.gr)[[feat.set]]],
        array.gr = annot.gr,
        maxgap = maxgap,
        clustergap = cluster.gap,
        featname = feat.set, 
        testname = test.set
      )
    )
  }
  # TFBS and TEs are tested for overlap within tfbs.max.gap / erv.max.gap
  do_FET_features_TFBS_TEs = function(test.set=NULL,feat.set=NULL,maxgap=NULL,featname=NULL,testname=NULL){
    return(
      FET_enrichment(
        test.gr = annot.gr[mcols(annot.gr)[[test.set]]],
        feature.gr = get(paste0(feat.set,'.gr')),
        array.gr = annot.gr,
        maxgap = maxgap,
        clustergap = cluster.gap,
        featname = strsplit(feat.set,'.gr')[[1]],
        testname = test.set
      )
    )
  }
  
  # MEs and Teh GxE
  for (feat.set in c('ME','teh.gxe')){
    for (test.set in c('replicated','ctrl')){
      all.res = rbind(all.res,
        do_FET_features(
          test.set = test.set,feat.set = feat.set,maxgap = feat.max.gap,featname = feat.set,testname = test.set
          )
      )
    }
  }
  
  # TFs
  for (feat.set in c('ctcf','zfp57','trim28')){
    for (test.set in c('replicated','ctrl')){
      all.res = rbind(all.res,
        do_FET_features_TFBS_TEs(
          test.set = test.set,feat.set = feat.set,maxgap = tfbs.max.gap,featname = feat.set,testname = test.set
                      )
      )
    }
  }
  
  # TEs
  if (include.erv){
    for (feat.set in te.bed){
      feature.set = gsub('.bed','.gr',bed.file)
      for (test.set in c('replicated','ctrl')){
        all.res = rbind(all.res,
          do_FET_features_TFBS_TEs(
            test.set = test.set,feat.set = feat.set,maxgap = erv.max.gap,featname = feat.set,testname = test.set
                        )
        )
      }
    }
  }  
  
  cat('\n-----------------------\n\n')
  
  # format for supp tables
  all.res$OR.cg = format(all.res$OR.cg,digits = 2)
  all.res$p.value.cg = format(all.res$p.value.cg,digits = 2)
  all.res$OR.clust = format(all.res$OR.clust,digits = 2)
  all.res$p.value.clust = format(all.res$p.value.clust,digits = 2)
  all.res$test.set = as.character(all.res$test.set)
  all.res[all.res=='discovery.cgs'] = 'discovery'
  all.res[all.res=='replicated'] = 'SoC-CpGs'
  all.res[all.res=='ctrl'] = 'matched_controls'
  
  print(all.res)
  
  if (output.to.file){
    
    sink()
    print(all.res)
    
  } 
  
    
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' PofOm Enrichment tests replicated loci and controls
#' for Supp Table 8
#' 
#' Test for enrichment of Zink PofOm, oo.gDMRs and sp.gDMRs
#' gDMR data provided by Dave Monk;
#' Sanchez-Delgado et al PLoS Gen 2016
#' 
#' @param annot.gr 
#' @param sink.file string file name with path to save verbose output
#'
#' @return
enrichment_tests_PofO = function(annot.gr=annot.gr,output.to.file=FALSE,sink.file=NULL,verbose.output=FALSE){
  
  require(tidyverse)
  require(GenomicRanges)
  
  # PARAMS
  cluster.gap = 1000 # max inter-cg distance to form cluster in test set

  if (output.to.file) sink(sink.file)

  # count gDMRs
  cat('analysis following number of gDMRs mapping to array background:\n')
  print(
    colSums(as.data.frame(mcols(annot.gr))
            [,c('oo.gDMR','oo.blast.gDMR','oo.blast.plac.gDMR','sp.gDMR','sp.blast.gDMR','zink.PofOm')])
  )
  
  # enrichment tests
  all.res = data.frame(
    test.set = character(),
    feat.set = character(),
    n.test.cgs = numeric(),
    OR.cg = numeric(),
    p.value.cg = numeric(),
    n.test.cgs.clust = numeric(),
    OR.clust = numeric(),
    p.value.clust = numeric(),
    stringsAsFactors = F
  )
  
  # REPLICATION
  for (feature.set in c('zink.PofOm','oo.gDMR','oo.blast.gDMR','oo.blast.plac.gDMR','sp.gDMR','sp.blast.gDMR')){
    all.res = rbind(all.res,
      FET_enrichment(
        test.gr = annot.gr[annot.gr$replicated],
        feature.gr = annot.gr[mcols(annot.gr)[[feature.set]]],
        array.gr = annot.gr,
        maxgap = 0,
        clustergap = cluster.gap,
        featname = feature.set, 
        testname = 'replicated'
      )
    )
  }
  
  # ALL MEs
  for (feature.set in c('zink.PofOm','oo.gDMR','oo.blast.gDMR','oo.blast.plac.gDMR','sp.gDMR','sp.blast.gDMR')){
    all.res = rbind(all.res,
      FET_enrichment(
        test.gr = annot.gr[annot.gr$ME],
        feature.gr = annot.gr[mcols(annot.gr)[[feature.set]]],
        array.gr = annot.gr,
        maxgap = 0,
        clustergap = cluster.gap,
        featname = feature.set, 
        testname = 'ME'
      )
    )
  }
  
  
  # ctrl
  for (feature.set in c('oo.gDMR','oo.blast.gDMR','oo.blast.plac.gDMR','sp.gDMR','sp.blast.gDMR','zink.PofOm')){
    all.res = rbind(all.res,
      FET_enrichment(
        test.gr = annot.gr[annot.gr$ctrl],
        feature.gr = annot.gr[mcols(annot.gr)[[feature.set]]],
        array.gr = annot.gr,
        maxgap = 0,
        clustergap = cluster.gap,
        featname = feature.set, 
        testname = 'ctrl'
      )
    )
  }
  
  # ALL PofOm
  for (feature.set in c('oo.gDMR','oo.blast.gDMR','oo.blast.plac.gDMR','sp.gDMR','sp.blast.gDMR')){
    all.res = rbind(all.res,
      FET_enrichment(
        test.gr = annot.gr[annot.gr$zink.PofOm],
        feature.gr = annot.gr[mcols(annot.gr)[[feature.set]]],
        array.gr = annot.gr,
        maxgap = 0,
        clustergap = cluster.gap,
        featname = feature.set, 
        testname = 'PofOm'
      )
    )
  }
  
  # format for supp tables
  all.res$OR.cg = format(all.res$OR.cg,digits = 2)
  all.res$p.value.cg = format(all.res$p.value.cg,digits = 2)
  all.res$OR.clust = format(all.res$OR.clust,digits = 2)
  all.res$p.value.clust = format(all.res$p.value.clust,digits = 2)
  all.res$test.set = as.character(all.res$test.set)
  all.res[all.res$test.set=='replicated','test.set'] = 'SoC-CpGs'
  all.res[all.res$test.set=='ctrl','test.set'] = 'ctrl'
  
  print(all.res)
  
  
  if (output.to.file){
    
    sink()
    print(all.res)
    
  } 
  
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Enrichment for intermediate methylation states compared to 
#' random controls
#'
#' @param annot.gr 
#' @param output.to.file send output to file?
#' @param sink.file sink filename
#'
#' @return NULL
intermediate_methylation_enrichment = function(annot.gr=annot.gr,output.to.file=FALSE,sink.file=NULL){

  require(GenomicRanges)
  require(tidyverse)
  
  inter.meth.df = as.data.frame(mcols(annot.gr)) %>%
    select(replicated,ME,rand,enid.mean.beta,emph.mean.beta) %>%
    gather(
      enid.mean.beta,emph.mean.beta,key = cohort,value = mean.beta
    )
  
  # pool cohorts
  inter.meth.df$cohort = NULL
  
  res.df = data.frame(
    pc.inter.meth.test=numeric(),
    pc.inter.meth.null=numeric(),
    chi.sq.p=numeric()
  )
  
  inter.meth.df$inter.meth = inter.meth.df$mean.beta <=0.75 & inter.meth.df$mean.beta >= 0.25
  inter.meth.df$soc.ME = inter.meth.df$replicated & inter.meth.df$ME
  inter.meth.df$soc.non.ME = inter.meth.df$replicated & !inter.meth.df$ME
  
  message('H1: ME replicated increased intermediate_vs_random_controls')
  data = inter.meth.df %>%
    filter(soc.ME | rand) %>%
    select(soc.ME,inter.meth)
  test=prop.test(table(data)[c(2,1),c(2,1)],alternative = 'greater')
  print(test)
  res.df['ME_replicated_vs_random',] = c(test$estimate[1],test$estimate[2],test$p.value)
  
  message('H1: other-SoC increased intermediate_vs_random_controls')
  data = inter.meth.df %>%
    filter(soc.non.ME | rand) %>%
    select(soc.non.ME,inter.meth)
  test=prop.test(table(data)[c(2,1),c(2,1)],alternative = 'greater')
  print(test)
  res.df['other_replicated_vs_random',] = c(test$estimate[1],test$estimate[2],test$p.value)

  message('H1: ME replicated increased intermediate_vs_array_background')
  data = inter.meth.df %>%
    select(soc.ME,inter.meth)
  test=prop.test(table(data)[c(2,1),c(2,1)],alternative = 'greater')
  print(test)
  res.df['ME_replicated_vs_array_background',] = c(test$estimate[1],test$estimate[2],test$p.value)
  
  message('H1: other-SoC increased intermediate_vs_random controls')
  data = inter.meth.df %>%
    select(soc.non.ME,inter.meth)
  test=prop.test(table(data)[c(2,1),c(2,1)],alternative = 'greater')
  print(test)
  res.df['other_replicated_vs_array_background',] = c(test$estimate[1],test$estimate[2],test$p.value)

  # convert to percentages
  res.df$pc.inter.meth.test = format(100*res.df$pc.inter.meth.test,digits = 3)
  res.df$pc.inter.meth.null = format(100*res.df$pc.inter.meth.null,digits = 3)
  res.df$chi.sq.p = formatC(100*res.df$chi.sq.p,format='e',digits = 2)
  
  
  print(res.df)

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate KS-matched controls for a set of query CpGs
#' 
#' KS-matched controls have similar methylation Beta distribution
#' to SoC-CpGs in ENID cohort
#'
#' @param annot.gr GR experiment annotations
#' @param n.matches.per.query # number of KS-matches per query cpg
#' @param ks.thresh # ks-pvalue threshold
#'      the larger this is the more likely we are to accept the null hypothesis that 
#'      the tested locus is drawn from the same distribution as the query locus 
#'      i.e. the closer the matched control CpG distribution is to the SoC-CpG query
#' @return List of KS matches - each query CpG mapped to n.matches targets
generate_ks_matches <- function(discovery.cgs=NULL,soc.cgs=NULL,me.cgs=NULL,all.array.cgs=NULL,n.matches.per.query=NULL,ks.thresh=NULL) {
  
  
  require(matrixStats)
  require(GenomicRanges)
  source('~/src/soc_multi_cohort_analysis/disc_repl/SoCFourier_modelling_functions.R')
  
  ks.thresh = 0.1 
  n.matches.per.query = 5 
  verbose = T # verbose output
  
  all.array.cgs = annot.gr$cpg
  enid.betas = m2beta(readRDS('~/projects/soc_enid_emphasis/cohort_data/ENID_Mfil_norm_filt_no_replicates.RDS')[all.array.cgs,])
  
  soc.cgs = annot.gr[annot.gr$replicated]$cpg
  # generate background set from which to draw KS-matched controls
  # this excludes discovery cgs (and hence SoC-CpGs) and known MEs
  cgs.background = setdiff(all.array.cgs,union(discovery.cgs,me.cgs))
  
  enid.median.betas = rowMedians(enid.betas)
  names(enid.median.betas) = row.names(enid.betas)
  
  ks.test.pvalue = function(x, y){
    
    ks.p = ks.test(x, y, alternative = c("two.sided"), exact = NULL)$p.value
    
    return(ks.p)
    
  }
  
  # search for KS matches ----
  # note this method is faster than vectorizing due to stopping conditions etc
  enid.median.betas.to.bgd = enid.median.betas[cgs.background]
  loci.to.ks.match = soc.cgs
  ks.matches.all = list()
  
  for (query.cg in loci.to.ks.match){
    
    # pre-filter on median +/- 0.2 to optimise search
    cg.median = enid.median.betas[query.cg]
    cgs.to.test = names(
      enid.median.betas.to.bgd[enid.median.betas.to.bgd<=cg.median+0.2 & enid.median.betas.to.bgd>=cg.median-0.2]
    )
    
    if (verbose) cat('matching',query.cg,'with median Beta=',cg.median,'\n')
    
    n.matches=0
    ks.matches = vector()
    
    # main KS search loop
    while(n.matches<5){
      test.cg = cgs.to.test[1]
      if(ks.test.pvalue(enid.betas[query.cg,],enid.betas[test.cg,]) > ks.thresh){
        ks.matches = c(ks.matches,test.cg)
        n.matches = n.matches + 1
        if (verbose) cat('n.matches=',n.matches,'\n')
      } 
      cgs.to.test = cgs.to.test[cgs.to.test!=test.cg]
    }
    
    ks.matches.all[[query.cg]] = ks.matches # store matches
    
    # remove KS matches from search list
    enid.median.betas.to.bgd = enid.median.betas.to.bgd[setdiff(names(enid.median.betas.to.bgd),ks.matches)]  
    if (verbose) cat(length(enid.median.betas.to.bgd),'background cgs remain\n\n')
    
  }
  
  return(ks.matches.all)
 
}


#^^^^^^^^^^^^^^^^^^^^^^^^^
#' Compute boostrap proportion
#'
#' @param d LOGICAL VECTOR from which proportion is calculated
#' @param n integer number of bootstraps
#'
#' @return
#' @export
#'
#' @examples
bootstrap.prop <- function(d, n){
  
  if (any(!is.logical(d))) {
    cat('require vector of logicals as input to bootstrap.prop function!\n')
    stop()
  }
  replicate(n, sum(sample(d, length(d), replace=TRUE))/length(d))
} 

alpha.quantiles.lower <- function(x, a) {
  x <- sort(x)
  x[ceiling(length(x) * a / 2)]
}
alpha.quantiles.upper <- function(x, a) {
  x <- sort(x)
  x[floor(length(x) * (1 - a / 2))]
}




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
require(GenomicRanges)

#' De-cluster genomic ranges object
#'
#' De-clustering is achieved by randomly sampling
#' one cg from each cluster. Clusters are contiguous
#' cgs using inter-cg distance of clustergap bp
#'
#' @param test.gr GR object to decluster
#' @param clustergap integer inter-cg distance to define cluster
#' @param verbose logical verbose output?
#'
#' @return de-clustered GR object where all cgs are at least
#'   clustergap bp apart
de_cluster_gr <- function(test.gr=NULL, clustergap=NULL, verbose=F) {
  
  # select random cg from multi-cg test clusters
  if (verbose) cat('sampling from',length(test.gr),'ranges\n')
  test.clust.gr = GenomicRanges::reduce(test.gr,min.gapwidth=clustergap)
  if (verbose) cat(length(test.clust.gr),'clusters using cluster inter-cg distance of',clustergap,'bp identified\n')
  olap = findOverlaps(test.clust.gr,test.gr)
  test.clust.n.cgs = as.vector(table(queryHits(olap)))
  test.1.cgs = which(test.clust.n.cgs==1) 
  test.gt1.cgs = which(test.clust.n.cgs>1) 
  if (verbose) cat(length(test.gt1.cgs),'clusters with more than 1 cg identified\n')
  
  # select random cg from muliptle cg clusters
  cg.sample = vapply(
    test.gt1.cgs,
    function(x) sample(subjectHits(olap[queryHits(olap)==x]),size = 1),
    FUN.VALUE = 1
  )
  if (verbose) cat(length(cg.sample),'random single cgs selected from multi-cg clusters\n')
  # create final feature set of singletons
  test.reduced.gr = test.gr[c(cg.sample,subjectHits(olap[queryHits(olap) %in% test.1.cgs]))]
  if (!all(width(test.reduced.gr)==2)) stop('declustering failed!!\n')
  if (verbose) cat('declustered object containing',length(test.reduced.gr),'cgs at least',clustergap,'bp apart created\n')
  
  return(test.reduced.gr)
  
}












#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Fisher Exact Test for enrichment
#' 
#' Test cgs are tested for overlap with feature cgs
#' within maxgap bp. Both cg-wise and cluster-wise
#' tests are performed. The latter reduces cluster bias
#' effect by sampling a single cg from multi-cg
#' test clusters. Clusters are defined as cgs within
#' maxgap bp.
#' 
#' @param test.gr GR object of test cgs
#' @param feature.gr GR object of feature cgs
#' @param array.gr GR array background
#' @param maxgap integer test cg considered to overlap 
#'  feature if within maxgap bp
#' @param clustergap integer inter-cg distance to define cluster
#' @param featname string feature name
#' @param testname string test cg name
#'
#' @return NULL
FET_enrichment = function(
  test.gr=NULL,feature.gr=NULL,array.gr=NULL,maxgap = NULL,clustergap=NULL,
  featname=NULL,testname=NULL
){
  
  require(GenomicRanges)
  
  if(maxgap==0) maxgap = -1 # GR convention for strict (non-adjacent) overlaps
  
  argnames = sys.call()
  if (!is.null(testname)){
    test.set.name = testname
  } else {
    test.set.name = strsplit(as.character(argnames[2]),'.gr')[[1]]
  }
  if (!is.null(featname)){
    feat.set.name = featname    
  } else {
    feat.set.name = strsplit(as.character(argnames[3]),'.gr')[[1]]
  }
  
  if (verbose){
    cat('\n*******\n')
    cat('Feature set:',feat.set.name,'\n')
    cat('Test set:',test.set.name,'\n\n')
  }
  
  # filter feature cgs to include only those within maxgap of array cg
  # all others have no effect on the analysis
  n.feat.orig = length(feature.gr)
  feature.gr = subsetByOverlaps(feature.gr,array.gr,maxgap = maxgap)
  if (verbose){
    cat(length(feature.gr),'features overlapping array cg within',maxgap,'bp (out of original',n.feat.orig,'features) will be in this analysis\n')
  }
  
  # cpg-wise enrichment
  # array cgs overlapping test
  if (verbose) cat('\n*** cg-wise enrichment:\n')
  olaps = queryHits(findOverlaps(array.gr,test.gr))
  mcols(array.gr)[,test.set.name] = F
  mcols(array.gr)[olaps,test.set.name] = T
  # array cgs overlapping feature within maxgap
  olaps = queryHits(findOverlaps(array.gr,feature.gr,maxgap = maxgap))
  mcols(array.gr)[,feat.set.name] = F
  mcols(array.gr)[olaps,feat.set.name] = T

  cont.table = table(mcols(array.gr[,c(test.set.name,feat.set.name)]))
  fet.cg = fisher.test(cont.table,alternative = 'greater')
  if (verbose){
    print(cont.table)
    cat('FET: OR=',fet.cg$estimate,'; p=',fet.cg$p.value,'\n\n')
  }
  
  # cluster-adusted enrichment
  if (verbose) cat('*** cluster-adjusted enrichment:\n')
  
  # de-cluster GR to create test gr of singletons
  test.reduced.gr = de_cluster_gr(test.gr = test.gr,clustergap = clustergap,verbose = F)

  # array cgs overlapping cluster-adjusted test within maxgap
  olaps = queryHits(findOverlaps(array.gr,test.reduced.gr))
  mcols(array.gr)[,paste0(test.set.name,'_red')] = F
  mcols(array.gr)[olaps,paste0(test.set.name,'_red')] = T
  if (verbose){
    cat(sum(mcols(array.gr)[,paste0(test.set.name,'_red')]),'array cgs overlap reduced (singleton)',test.set.name,'within',maxgap,'bp\n')
  }
  
  cont.table = table(mcols(array.gr[,c(paste0(test.set.name,'_red'),feat.set.name)]))
  if (verbose) cat(cont.table[2,2],feat.set.name,'overlap',test.set.name,'within',maxgap,'bp\n')
  fet.clust = fisher.test(cont.table,alternative = 'greater')
  if (verbose) cat('FET: OR=',fet.clust$estimate,'; p=',fet.clust$p.value,'\n\n')
  
  return(
    data.frame(
      test.set=test.set.name,
      feat.set=feat.set.name,
      n.test.cgs = length(test.gr),
      OR.cg=as.vector(fet.cg$estimate),
      p.value.cg=fet.cg$p.value,
      n.test.cgs.clust = length(test.reduced.gr),
      OR.clust=as.vector(fet.clust$estimate),
      p.value.clust=fet.clust$p.value
    )
  )
  
}
