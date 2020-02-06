#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Functions for annotating CpGs and CpG clusters identified in
# the SoC discovery-relication analysis.
# These are called from SoCFourier_main_analysis.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Annotate array CpGs following Fourier regression analysis
#' 
#' Annotations include:
#' genomic location; 
#' discovery, soc (replicated), ctrl, rand, ME, PofOm, GxE etc labels
#' Fourier regression model fit stats (date of conception at methylation max/min, 
#' amplitude etc);
#' CGI and gene annotations from Illumina 450k manifest
#'
#' @param array.cgs character vector cgs to be annotated
#' @param enid.model.stats df ENID summary stats from Fourier regression fit
#' @param emph.fit.stats df EMPHASIS summary stats from Fourier regression fit
#'
#' @return annotated array GR object
#' 
annotate_array = function(array.cgs=NULL, enid.model.stats=NULL, emph.model.stats = NULL){
  
  require(GenomicRanges)
  require(tidyverse)
  require(matrixStats)
  
  proj.dir = '~/projects/soc_enid_emphasis/' # output directory
  src.dir = '~/src/soc_enid_emphasis/' # source code directory
  
  # create annotated array GR object
  array.annot.gr = readRDS('/data/450k/cpglocs.450k.gr.RDS')
  array.annot.gr = array.annot.gr[array.annot.gr$cpg %in% array.cgs]
  
  # annotate SoC-associated cgs (discovery and replicated) and controls
  array.annot.gr$discovery = FALSE
  array.annot.gr$replicated = FALSE
  array.annot.gr$ctrl = FALSE
  array.annot.gr$rand = FALSE
  array.annot.gr[array.annot.gr$cpg %in% filter(enid.model.stats,locus=='discovery')$cg]$discovery = TRUE
  array.annot.gr[array.annot.gr$cpg %in% filter(enid.model.stats,replicated)$cg]$replicated = TRUE
  array.annot.gr[array.annot.gr$cpg %in% filter(enid.model.stats,locus=='ctrl')$cg]$ctrl = TRUE
  array.annot.gr[array.annot.gr$cpg %in% filter(enid.model.stats,locus=='rand')$cg]$rand = TRUE

  # MEs (Van Baak + Sci Adv overlapping 450k)
  me.cgs = readRDS(paste0(proj.dir,'cohort_data/all.MEs.RDS'))
  me.cgs = intersect(me.cgs,cgs.array.intersect)
  me.gr = array.annot.gr[array.annot.gr$cpg %in% me.cgs]
  
  array.annot.gr$ME = FALSE
  array.annot.gr$ME.100bp = FALSE
  array.annot.gr[array.annot.gr$cpg %in% me.cgs]$ME = TRUE
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,me.gr,maxgap = 1000))]$ME.100bp = TRUE
  # table(as.data.frame(mcols(array.annot.gr[array.annot.gr$discovery]))[,c('discovery','ME')])
  # table(as.data.frame(mcols(array.annot.gr[array.annot.gr$replicated]))[,c('replicated','ME')])
  
  # Monk gDMRs
  oo_gdmrs = readRDS('/data/hg19_annots/Monk_gDMRs/oocyte_gDMRs.RDS')
  sp_gdmrs = readRDS('/data/hg19_annots/Monk_gDMRs/sperm_gDMRs.RDS')
  oo_blast_gdmrs = readRDS('/data/hg19_annots/Monk_gDMRs/oocyte_blastocyst_gDMRs.RDS')
  sp_blast_gdmrs = readRDS('/data/hg19_annots/Monk_gDMRs/sperm_blastocyst_gDMRs.RDS')
  oo_blast_plac_gdmrs = readRDS('/data/hg19_annots/Monk_gDMRs/oocyte_blastocyst_placenta_gDMRs.RDS')
  oo_gdmrs.gr = makeGRangesFromDataFrame(oo_gdmrs,keep.extra.columns = T)
  oo_blast_gdmrs.gr = makeGRangesFromDataFrame(oo_blast_gdmrs,keep.extra.columns = T)
  sp_gdmrs.gr = makeGRangesFromDataFrame(sp_gdmrs,keep.extra.columns = T)
  sp_blast_gdmrs.gr = makeGRangesFromDataFrame(sp_blast_gdmrs,keep.extra.columns = T)
  oo_blast_plac_gdmrs.gr = makeGRangesFromDataFrame(oo_blast_plac_gdmrs,keep.extra.columns = T)
  
  array.annot.gr$oo.gDMR = F
  array.annot.gr$oo.blast.gDMR = F
  array.annot.gr$oo.blast.plac.gDMR = F
  array.annot.gr$sp.gDMR = F
  array.annot.gr$sp.blast.gDMR = F
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,oo_gdmrs.gr))]$oo.gDMR = T
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,oo_blast_gdmrs.gr))]$oo.blast.gDMR = T
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,oo_blast_plac_gdmrs.gr))]$oo.blast.plac.gDMR = T
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,sp_gdmrs.gr))]$sp.gDMR = T
  array.annot.gr[queryHits(findOverlaps(array.annot.gr,sp_blast_gdmrs.gr))]$sp.blast.gDMR = T

  # LOAD Zink et al PofOm cgs (Nat Gen 2018)
  # (see ~/src/soc_multi_cohort_analysis/disc_repl/misc/Zink_PofOm_to_h19_gr.R)
  imp.zink.gr = readRDS('/data/cpg_sets/imprinted/Zink_SuppTab1_hg19.gr')
  PofOm.cgs = subsetByOverlaps(array.annot.gr,imp.zink.gr)$cpg
  array.annot.gr$zink.PofOm = F
  array.annot.gr[array.annot.gr$cpg %in% PofOm.cgs]$zink.PofOm = T

  # Teh GxE
  gxe.cgs = as.character(read.csv('/data/cpg_sets/gVMR/Teh_Supplemental_Tables_6.csv',skip = 1)$CpG)
  gxe.gr = array.annot.gr[array.annot.gr$cpg %in% gxe.cgs]
  array.annot.gr$teh.gxe = F
  array.annot.gr[array.annot.gr$cpg %in% gxe.gr$cpg]$teh.gxe = T
  
  # add ENID/EMPHASIS model fit stats
  # ENID
  array.annot.gr$enid.max.doc.theta = NA
  array.annot.gr$enid.min.doc.theta = NA
  array.annot.gr$enid.amplitude = NA
  array.annot.gr$enid.lrt.pval = NA
  array.annot.gr$enid.lrt.fdr = NA
  array.annot.gr[match(enid.model.stats$cg,array.annot.gr$cpg)]$enid.max.doc.theta = 
    as.numeric(enid.model.stats$max.doc.theta)
  array.annot.gr[match(enid.model.stats$cg,array.annot.gr$cpg)]$enid.min.doc.theta = 
    as.numeric(enid.model.stats$min.doc.theta)
  array.annot.gr[match(enid.model.stats$cg,array.annot.gr$cpg)]$enid.amplitude = 
    as.numeric(enid.model.stats$amplitude)
  array.annot.gr[match(enid.model.stats$cg,array.annot.gr$cpg)]$enid.lrt.pval = 
    as.numeric(enid.model.stats$lrt.pval)
  array.annot.gr[match(enid.model.stats$cg,array.annot.gr$cpg)]$enid.lrt.fdr = 
    as.numeric(enid.model.stats$lrt.fdr)
  
  # EMPHASIS
  array.annot.gr$emph.max.doc.theta = NA
  array.annot.gr$emph.min.doc.theta = NA
  array.annot.gr$emph.amplitude = NA
  array.annot.gr$emph.lrt.pval = NA
  array.annot.gr$emph.lrt.fdr = NA
  array.annot.gr$emph.disc.lrt.fdr = NA
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.max.doc.theta = 
    as.numeric(emph.model.stats$max.doc.theta)
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.min.doc.theta = 
    as.numeric(emph.model.stats$min.doc.theta)
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.amplitude = 
    as.numeric(emph.model.stats$amplitude)
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.lrt.pval = 
    as.numeric(emph.model.stats$lrt.pval)
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.lrt.fdr = 
    as.numeric(emph.model.stats$lrt.fdr)
  array.annot.gr[match(emph.model.stats$cg,array.annot.gr$cpg)]$emph.disc.lrt.fdr = 
    as.numeric(emph.model.stats$disc.lrt.fdr)
  
  # mean, median and variance of methylation Betas
  if (!exists('enid.M')){
    enid.M = 
      readRDS(paste0(proj.dir,'cohort_data/ENID_Mfil_norm_filt_no_replicates.RDS'))
  }
  if (!exists('emph.M')){
  emph.M = 
    readRDS(paste0(proj.dir,'cohort_data/EMPH_Mfil_norm_filt_no_replicates.RDS'))
  }
  source(paste0(src.dir,'SoCFourier_modelling_functions.R'))
  array.annot.gr$enid.mean.beta = rowMeans(m2beta(enid.M[array.annot.gr$cpg,]))
  array.annot.gr$emph.mean.beta = rowMeans(m2beta(emph.M[array.annot.gr$cpg,]))
  array.annot.gr$enid.median.beta = rowMedians(m2beta(enid.M[array.annot.gr$cpg,]))
  array.annot.gr$emph.median.beta = rowMedians(m2beta(emph.M[array.annot.gr$cpg,]))
  array.annot.gr$enid.var.beta = rowVars(m2beta(enid.M[array.annot.gr$cpg,]))
  array.annot.gr$emph.var.beta = rowVars(m2beta(emph.M[array.annot.gr$cpg,]))
  
  # add illumina gene and CGI annotations
  annots.450k.df = readRDS('/data/450k/HumanMethylation450_15017482_v1-2.RDS')
  annots.450k.df = filter(annots.450k.df,grepl('^cg',Name)) # remove control probes etc
  row.names(annots.450k.df) = as.character(annots.450k.df$Name)
  array.annot.gr$UCSC_RefGene_Name = as.character(annots.450k.df[array.annot.gr$cpg,'UCSC_RefGene_Name'])
  array.annot.gr$UCSC_RefGene_Group = as.character(annots.450k.df[array.annot.gr$cpg,'UCSC_RefGene_Group'])
  array.annot.gr$Relation_to_UCSC_CpG_Island = as.character(annots.450k.df[array.annot.gr$cpg,'Relation_to_UCSC_CpG_Island'])
  array.annot.gr$Enhancer = as.character(annots.450k.df[array.annot.gr$cpg,'Enhancer'])
  
  return(array.annot.gr)
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Create annotated tables of discovery and SoC-CpGs for publication
#' 
#' These include details on genomic location, overlapping features,
#' gene info, and model fit stats.
#'
#'
#' @param annot.gr 
#'
#' @return df annotation table
create_annotation_tables = function(annot.gr = NULL){

  require(GenomicRanges)
  require(tidyverse)
  
  annot.gr = annot.gr[annot.gr$discovery]
  
  # tidy up gene lists
  annot.gr$UCSC_RefGene_450k = 
    unlist(lapply(
      annot.gr$UCSC_RefGene_Name, 
      function(x) ifelse(x=="", "", paste0(unique(unlist(strsplit(as.character(x),';'))),collapse=' ')))
    )

  # 1. DISCOVERY
  annot.discovery.df = as.data.frame(annot.gr)
  annot.discovery.df$chr = as.numeric(str_remove_all(annot.discovery.df$seqnames,'chr')) # easier to order by chr
  # table to print
  annot.discovery.df = annot.discovery.df %>%
    select(cpg,chr,loc=start,
           ME,ME.100bp,
           oo.gDMR,oo.blast.gDMR,oo.blast.plac.gDMR,sp.gDMR,sp.blast.gDMR,
           zink.PofOm,teh.gxe,
           enid.max.doc.theta,enid.min.doc.theta,enid.amplitude,enid.lrt.pval,enid.lrt.fdr,
           emph.max.doc.theta,emph.min.doc.theta,emph.amplitude,emph.lrt.pval,emph.lrt.fdr,
           UCSC_RefGene_Group,Enhancer,Relation_to_UCSC_CpG_Island,UCSC_RefGene_450k) %>%
    arrange(chr,loc)

  # 2. REPLICATION
  annot.replication.gr = annot.gr[annot.gr$replicated]
  annot.replication.df = as.data.frame(annot.replication.gr)
  annot.replication.df$chr = as.numeric(str_remove_all(annot.replication.df$seqnames,'chr')) # easier to order by chr
  annot.replication.df = annot.replication.df %>%
    select(cpg,chr,loc=start,replicated,
           ME,ME.100bp,
           oo.gDMR,oo.blast.gDMR,oo.blast.plac.gDMR,sp.gDMR,sp.blast.gDMR,
           zink.PofOm,teh.gxe,
           enid.max.doc.theta,enid.min.doc.theta,enid.amplitude,enid.lrt.pval,enid.lrt.fdr,
           emph.max.doc.theta,emph.min.doc.theta,emph.amplitude,emph.lrt.pval,emph.lrt.fdr,
           UCSC_RefGene_Group,Enhancer,Relation_to_UCSC_CpG_Island,UCSC_RefGene_450k) %>%
    arrange(chr,loc)

  # 3. REPLICATION CLUSTERED
  annot.replication.clust.gr = GenomicRanges::reduce(annot.replication.gr,min.gapwidth=1000)
  # add overlapping genes
  olaps = findOverlaps(annot.replication.clust.gr,annot.replication.gr)
  olap.n.cpgs = table(queryHits(olaps))
  for (clust in seq(length(annot.replication.clust.gr))){
    annot.replication.clust.gr$n.cpgs[clust] = olap.n.cpgs[clust]
    annot.replication.clust.gr$UCSC_RefGene_450k[clust] = 
      paste(unique(as.character(annot.replication.gr[subjectHits(olaps[queryHits(olaps)==clust])]$UCSC_RefGene_450k)),
            collapse = ',')
  }
  annot.replication.clust.df = as.data.frame(annot.replication.clust.gr)
  annot.replication.clust.df$chr = as.numeric(str_remove_all(annot.replication.clust.df$seqnames,'chr')) # easier to order by chr
  annot.replication.clust.df = annot.replication.clust.df %>%
    select(chr,loc=start,n.cpgs,width,UCSC_RefGene_450k) %>%
    arrange(chr,loc)
  
  return(
    list(
      annot.table.discovery=annot.discovery.df,
      annot.table.replication=annot.replication.df,
      annot.table.replication.clustered=annot.replication.clust.df
    )
  )
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Determine cpg cluster sizes 
#'
#' @param test.gr GR object loci of interest
#' @param test.name name for test set
#' @param cluster.gap inter-cpg distance (bp) to define clusters
#'
#' @return NULL
get.cluster.sizes = function(test.gr=NULL,test.name=NULL,cluster.gap=NULL){
  test.clust.gr = GenomicRanges::reduce(test.gr,min.gapwidth=cluster.gap)
  message(test.name)
  cat(length(test.gr),'cgs map to',length(test.clust.gr),'cg clusters\n')
  
  # cluster sizes
  olaps = findOverlaps(query = test.clust.gr,subject = test.gr)
  print(table(table(queryHits(olaps))))
  cat('\n')
}
