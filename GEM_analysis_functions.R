#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# mQTL and GxE analysis using GEM
# functions called from SoCFourier_GEM_analysis.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' configure PLINK imputed GSA genotype files for GEM
#'
#' 1. reinstate correct chr in .bim file
#'  2. merge .bim files
#'  3. filter on MAF 10%
#'  4. recode to dosage format required by GEM
#' @return NULL
configuire_PLINK_files_for_GEM = function(){
  
  require(tidyverse)
  
  # 1. add correct chromosome label to 1st column of .bim file
  
  for (chr in seq(1,22)){
    
    cat('reformatting chr',chr,'.bim file\n')  
    bim = read_tsv(paste0('/data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/by_chromosome/fcgene_plink_chr',chr,'_info09.bim'),col_names = F)
    bim$X1 = chr
    write_tsv(
      bim,
      paste0(
        '/data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/',
        'fcgene_plink_chr',chr,'_info09.bim'
      ),
      col_names = F
    )
  }
  
  # 2. merge into single PLINK file
  system(
    paste(
      '/usr/bin/PLINK/plink',
      '--bfile /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/by_chromosome/fcgene_plink_chr1_info09',
      '--make-bed',
      '--merge-list /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/merge_chromosomes_list.txt',
      '--out /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/EMPHASIS_GSA_imputed_info_09_ALL',
      sep = ' '
    )
  )
  
  # 3. filter on MAF
  system(
    paste(
      '/usr/bin/PLINK/plink',
      '--bfile /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/EMPHASIS_GSA_imputed_info_09_ALL',
      '--make-bed',
      '--maf 0.1',
      '--out /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/working/EMPHASIS_GSA_imputed_info_09_ALL_maf_10',
      sep = ' '
    )
  )
  
  # 4. recode to dosage format and long format
  # split into chunks to avoid segmentation fault
  system(
    paste(
      '/usr/bin/plink2',
      '--bfile /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/working/EMPHASIS_GSA_imputed_info_09_ALL_maf_10',
      '--recode A-transpose',
      '--out /data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/working/EMPHASIS_GSA_imputed_info_09_ALL_maf_10_dose',
      sep = ' '
    )
  )
  
  # load to tibble
  snps = read_tsv('/data/EMPHASIS_GSA/EMPHASIS_GSA_imputed_info_09/merged/working/EMPHASIS_GSA_imputed_info_09_ALL_maf_10_dose.traw')
  
  # save RDS
  saveRDS(object = snps,file = '~/projects/soc_enid_emphasis_analysis/GxE/R_objects/GSA_SNPs_MAF10_imputed.rds')
  
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' create data files (txt) for GEM
#'
#' @param annot.gr 
#' @param output.to.file 
#' @param sink.file 
#'
#' @return NULL
prepare_GEM_data = function(annot.gr=NULL,output.to.file=FALSE,sink.file=NULL){
  
  require(tidyverse)
  require("GEM")
  require(GenomicRanges)
  require(lubridate)
  require(reshape2)
  
  # PARAMS
  proj.dir = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/'
  hz.filter.freq = 10 # minimum number of hz variants (for SNP filtering)
  
  if (output.to.file){
    sink(paste0(proj.dir,'results/GEM_GxE_data_prep_',timestamp(),'.log'))
  }
  
  # get SoC-CpGs and controls
  soc.cgs = annot.gr[annot.gr$replicated]$cpg
  cn.cgs = annot.gr[annot.gr$cn.var]$cpg
  rand.cgs = annot.gr[annot.gr$cn.rand]$cpg
  rm(annot.gr)
  
  # load methylation and genotype data 
  design <- readRDS("~/projects/soc_enid_emphasis/cohort_data/EMPHASIS_design.PCs1-6.RDS")
  EMPH_Mfil <- readRDS("~/projects/soc_enid_emphasis/cohort_data/EMPH_Mfil_norm_filt_no_replicates.RDS")
  # replace array IDs with EMPH IDs
  EMPH_pdata <- readRDS("~/projects/soc_enid_emphasis/cohort_data/EMPHASIS_pdata.RDS")
  rownames(design) = EMPH_pdata[rownames(design),'Subject_ID']
  colnames(EMPH_Mfil) = EMPH_pdata[colnames(EMPH_Mfil),'Subject_ID']

  # initialise SNP data
  GMB_SNPs <- readRDS(paste0(proj.dir,"R_objects/GSA_SNPs_MAF10_imputed.RDS"))
  # combine CHR and SNP names (chr missing from some SNP names)
  GMB_SNPs$SNP_ID = paste(paste0('chr',GMB_SNPs$CHR),GMB_SNPs$SNP,sep = '_')
  # drop superfluous columns
  GMB_SNPs = GMB_SNPs %>% select(SNP_ID,matches('^[0-9]'))
  # plink A-transpose option codes MAJOR allele dosage; swith to MINOR allele dosage for 
  # compatability / sanity checks with non-imputed data
  GMB_SNPs[,2:289] = 2 - GMB_SNPs[,2:289]
  # reformat sample IDs to match methylation data
  GMB_SNPs <- GMB_SNPs %>% setNames(gsub("\\d+_", "",names(.)))
  # check all IDs present in methylation data
  snp.samples.to.keep = names(GMB_SNPs[,2:289])[names(GMB_SNPs[,2:289]) %in% colnames(EMPH_Mfil)]
  GMB_SNPs = GMB_SNPs[,c('SNP_ID',snp.samples.to.keep)]
  # generate new filtered SNP list with hz variants > 10
  GMB_SNPs = GMB_SNPs[
    rowSums(GMB_SNPs[,2:ncol(GMB_SNPs)]==2,na.rm = T)>=hz.filter.freq,
  ]
  cat(prettyNum(nrow(GMB_SNPs),big.mark = ','),'SNPs with MAF>10% and',hz.filter.freq,
      'or more homozgous variant individuals will be considered for this analysis\n\n')
  # SNP IDs to rownames
  GMB_SNPs = GMB_SNPs %>% column_to_rownames('SNP_ID')

  #match sample order between design, imputed GSA, and Mfil 
  design <- design[snp.samples.to.keep,]
  EMPH_Mfil <- EMPH_Mfil[,snp.samples.to.keep]
  identical(rownames(design),colnames(EMPH_Mfil))
  identical(rownames(design),colnames(EMPH_Mfil))
  identical(rownames(design),colnames(GMB_SNPs))
  identical(colnames(EMPH_Mfil),colnames(GMB_SNPs))
  
  # SoC-CpGs
  soc_Mfil <- EMPH_Mfil[rownames(EMPH_Mfil) %in% soc.cgs,]
  identical(colnames(soc_Mfil),rownames(design))
  
  ctrl_Mfil <- EMPH_Mfil[rownames(EMPH_Mfil) %in% cn.cgs,]
  identical(colnames(ctrl_Mfil),rownames(design))
  
  rnd_Mfil <- EMPH_Mfil[rownames(EMPH_Mfil) %in% rand.cgs,]
  identical(colnames(rnd_Mfil),rownames(design))
  
  dim(GMB_SNPs)
  dim(soc_Mfil)
  dim(ctrl_Mfil)
  dim(rnd_Mfil)
  dim(design)
  
  identical(colnames(GMB_SNPs),rownames(design))
  identical(colnames(soc_Mfil),rownames(design))
  identical(colnames(GMB_SNPs),colnames(soc_Mfil))
  
  rm(EMPH_Mfil)
    
  #^^^^^^^^^^^^^^^^^^^^^^^^^
  # save R objects to file
  saveRDS(design,file = paste0(proj.dir,'R_objects/design_imputed_GEM.RDS'))
  saveRDS(GMB_SNPs,file = paste0(proj.dir,'R_objects/GMB_SNPs_imputed_GEM.RDS'))
  
  # save as text files for GEM
  write.table(GMB_SNPs,paste0(proj.dir,"data/GMB_SNPs.txt"),sep="\t")
  write.table(GMB_SNPs[1:100,],paste0(proj.dir,"data/GMB_SNPs_test.txt"),sep="\t") # for GEM testing
  saveRDS(soc_Mfil,paste0(proj.dir,"R_objects/soc_Mfil.RDS"))
  
  write.table(soc_Mfil,paste0(proj.dir,"data/soc_Mfil.txt"),sep="\t")
  saveRDS(ctrl_Mfil,paste0(proj.dir,"R_objects/ctrl_Mfil.RDS"))
  write.table(ctrl_Mfil,paste0(proj.dir,"data/ctrl_Mfil.txt"),sep="\t")
  saveRDS(rnd_Mfil,paste0(proj.dir,"R_objects/rnd_Mfil.RDS"))
  write.table(rnd_Mfil,paste0(proj.dir,"data/rnd_Mfil.txt"),sep="\t")
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # covariate data for GEM
  # 
  # generate following objects required for GEM
  # ('E' var of interest is always last column of cov file):
  # 
  # env_sin: sin(doc.theta) E variable of interest (for GEM E model)
  # env_cos: cos(doc.theta) E variable of interest (for GEM E model)
  # covs_sin: adjustment covariates with sin as E variable of interest (for GEM GxE model)
  # covs_cos: adjustment covariates with cos as E variable of interest (for GEM GxE model)
  # covs_only: adjustment covariates with no E variable (for GEM G model)
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
  env_sin <- t(design[,"sin(doc.theta)",drop=F])
  env_cos <- t(design[,"cos(doc.theta)",drop=F])
  
  covs_sin <- design[,c("SexM","Age","MasterGroupNo2",
                        "PC1","PC2","PC3","PC4","PC5","PC6",
                        "sin(doc.theta)")]
  covs_sin[,'SexM'] <- as.character(covs_sin[,'SexM'])
  covs_sin = t(covs_sin)
  
  covs_cos <- design[,c("SexM","Age","MasterGroupNo2",
                        "PC1","PC2","PC3","PC4","PC5","PC6",
                        "cos(doc.theta)")]
  covs_cos[,'SexM'] <- as.character(covs_cos[,'SexM'])
  covs_cos = t(covs_cos)
  
  covs_only <- design[,c("SexM","Age","MasterGroupNo2",
                         "PC1","PC2","PC3","PC4","PC5","PC6")]
  covs_only[,'SexM'] <- as.character(covs_only[,'SexM'])
  covs_only = t(covs_only)
  
  
  dim(env_sin)
  row.names(env_sin)
  dim(env_cos)
  row.names(env_cos)
  dim(covs_sin)
  rownames(covs_sin)
  covs_sin[,c(1,2)]
  dim(covs_cos)
  rownames(covs_cos)
  covs_cos[,c(1,2)]
  
  identical(colnames(env_sin),rownames(design))
  identical(colnames(env_sin),colnames(GMB_SNPs))
  identical(colnames(env_sin),colnames(soc_Mfil))
  identical(colnames(env_cos),rownames(design))
  identical(colnames(env_cos),colnames(GMB_SNPs))
  identical(colnames(env_cos),colnames(soc_Mfil))
  identical(colnames(covs_sin),rownames(design))
  identical(colnames(covs_sin),colnames(GMB_SNPs))
  identical(colnames(covs_sin),colnames(soc_Mfil))
  identical(colnames(covs_cos),rownames(design))
  identical(colnames(covs_cos),colnames(GMB_SNPs))
  identical(colnames(covs_cos),colnames(soc_Mfil))
  
  #save as text files
  write.table(env_sin,paste0(proj.dir,"data/env_sin.txt"),sep="\t",quote = F)
  write.table(env_cos,paste0(proj.dir,"data/env_cos.txt"),sep="\t",quote = F)
  write.table(covs_sin,paste0(proj.dir,"data/covs_sin.txt"),sep="\t",quote = F)
  write.table(covs_cos,paste0(proj.dir,"data/covs_cos.txt"),sep="\t",quote = F)
  write.table(covs_only,paste0(proj.dir,"data/covs_only.txt"),sep="\t",quote = F)

  if (output.to.file) sink()
  
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' GEM mQTL and GxE analysis using imputed genotypes from the 
#' EMPHASIS cohort
#'
#' 1. find G1 and G2 SNPs using GEM for 
#'   a) 125 SoC-CpGs
#'   b) 625 KS-matched controls
#'   c) 625 random background control
#'   
#'   GEM models:
#'    E environment only
#'    G genetic (mQTL) main effects not including environment in model
#'    GxE gene-environment interactions; genetic and envioronmental main effects included in model
#'    
#'   all GEM models use single, most significant E-term (cos or sin) only
#'   
#'   2. compare adjusted R^2 for all G1 and G2 models and find 'winning' model using AIC
#'   
#' @param write.to.file 
#'
#' @return NULL
run_GEM_analysis = function(write.to.file=FALSE){
  
  require(tidyverse)
  require("GEM")
  require(GenomicRanges)
  
  setwd('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/')
  
  if (write.to.file) {
    
    sink(paste0(proj.dir,'results/GEM_run',timestamp(),'.log'))
    cat('log file saved to:',paste0(proj.dir,'results/GEM_run',timestamp(),'.log'))
    
  }
  
  
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # Run E-models ----
  # these are run purely to identify the most significant Fourier term (sin/cos) for GxE analysis
  # do this using R lm combining sin, cos and cov terms

  env_sin = read.table(file = './data/env_sin.txt',header = T,sep = '\t')
  env_cos = read.table(file = './data/env_cos.txt',header = T,sep = '\t')
  covs_only = read.table(file = './data/covs_only.txt',header = T,sep = '\t')
  design = rbind(env_sin,env_cos,covs_only)
  
  # SoC-CpGs
  soc_Mfil = read.table(file = './data/soc_Mfil.txt',header = T,sep = '\t')
  lm.fit = lapply(
    row.names(soc_Mfil),
    function(x) summary(lm(formula(paste(x,'~.')),data = as.data.frame(t(rbind(soc_Mfil[x,],design)))))$coeff[2:3,4]
    # function(x) colnames(as.data.frame(x))
  )
  names(lm.fit) = row.names(soc_Mfil)
  lm.fit.pvals = t(as.data.frame(lm.fit))
  colnames(lm.fit.pvals) = c('sin','cos')
  
  soc_model_sel <- data.frame(which.E=character(),stringsAsFactors = F)
  for (cpg.id in row.names(lm.fit.pvals)){
    soc_model_sel[cpg.id,'which.E'] = 
      ifelse(
        lm.fit.pvals[cpg.id,'sin'] <= lm.fit.pvals[cpg.id,'cos'],
        'sin', 'cos'
      )
  }
  saveRDS(soc_model_sel,file = './R_objects/soc_model_sel.RDS')
  
  # generate SoC_M objects corresponding to dominant sin and cos terms
  soc_Mfil_sin = soc_Mfil[row.names(soc_model_sel[soc_model_sel$which.E=='sin',,drop=F]),]
  soc_Mfil_cos = soc_Mfil[row.names(soc_model_sel[soc_model_sel$which.E=='cos',,drop=F]),]
  write.table(soc_Mfil_sin,"./data/soc_Mfil_sin.txt",sep="\t")
  write.table(soc_Mfil_cos,"./data/soc_Mfil_cos.txt",sep="\t")
  
  # ctrl CpGs
  ctrl_Mfil = read.table(file = './data/ctrl_Mfil.txt',header = T,sep = '\t')
  lm.fit = lapply(
    row.names(ctrl_Mfil),
    function(x) summary(lm(formula(paste(x,'~.')),data = as.data.frame(t(rbind(ctrl_Mfil[x,],design)))))$coeff[2:3,4]
    # function(x) colnames(as.data.frame(x))
  )
  names(lm.fit) = row.names(ctrl_Mfil)
  lm.fit.pvals = t(as.data.frame(lm.fit))
  colnames(lm.fit.pvals) = c('sin','cos')
  
  ctrl_model_sel <- data.frame(which.E=character(),stringsAsFactors = F)
  for (cpg.id in row.names(lm.fit.pvals)){
    ctrl_model_sel[cpg.id,'which.E'] = 
      ifelse(
        lm.fit.pvals[cpg.id,'sin'] <= lm.fit.pvals[cpg.id,'cos'],
        'sin', 'cos'
      )
  }
  saveRDS(ctrl_model_sel,file = './R_objects/ctrl_model_sel.RDS')
  
  # generate ctrl_M objects corresponding to dominant sin and cos terms
  ctrl_Mfil_sin = ctrl_Mfil[row.names(ctrl_model_sel[ctrl_model_sel$which.E=='sin',,drop=F]),]
  ctrl_Mfil_cos = ctrl_Mfil[row.names(ctrl_model_sel[ctrl_model_sel$which.E=='cos',,drop=F]),]
  write.table(ctrl_Mfil_sin,"./data/ctrl_Mfil_sin.txt",sep="\t")
  write.table(ctrl_Mfil_cos,"./data/ctrl_Mfil_cos.txt",sep="\t")
  
  # random control CpGs
  rnd_Mfil = read.table(file = './data/rnd_Mfil.txt',header = T,sep = '\t')
  lm.fit = lapply(
    row.names(rnd_Mfil),
    function(x) summary(lm(formula(paste(x,'~.')),data = as.data.frame(t(rbind(rnd_Mfil[x,],design)))))$coeff[2:3,4]
    # function(x) colnames(as.data.frame(x))
  )
  names(lm.fit) = row.names(rnd_Mfil)
  lm.fit.pvals = t(as.data.frame(lm.fit))
  colnames(lm.fit.pvals) = c('sin','cos')
  
  rnd_model_sel <- data.frame(which.E=character(),stringsAsFactors = F)
  for (cpg.id in row.names(lm.fit.pvals)){
    rnd_model_sel[cpg.id,'which.E'] = 
      ifelse(
        lm.fit.pvals[cpg.id,'sin'] <= lm.fit.pvals[cpg.id,'cos'],
        'sin', 'cos'
      )
  }
  saveRDS(rnd_model_sel,file = './R_objects/rnd_model_sel.RDS')
  
  # generate rnd_M objects corresponding to dominant sin and cos terms
  rnd_Mfil_sin = rnd_Mfil[row.names(rnd_model_sel[rnd_model_sel$which.E=='sin',,drop=F]),]
  rnd_Mfil_cos = rnd_Mfil[row.names(rnd_model_sel[rnd_model_sel$which.E=='cos',,drop=F]),]
  write.table(rnd_Mfil_sin,"./data/rnd_Mfil_sin.txt",sep="\t")
  write.table(rnd_Mfil_cos,"./data/rnd_Mfil_cos.txt",sep="\t")


  # Run GEM G-models ----
  GEM_Gmodel("./data/GMB_SNPs.txt","./data/covs_only.txt","./data/soc_Mfil.txt",
             0.001, "./GEM_results/GEM_Gmodel_SoC_p001.txt")
  GEM_Gmodel("./data/GMB_SNPs.txt","./data/covs_only.txt","./data/ctrl_Mfil.txt",
             0.001, "./GEM_results/GEM_Gmodel_ctrl_p001.txt")
  GEM_Gmodel("./data/GMB_SNPs.txt","./data/covs_only.txt","./data/rnd_Mfil.txt",
             0.001, "./GEM_results/GEM_Gmodel_rnd_p001.txt")
  

  # check all cgs covered at this sig threshold
  gem.soc.g.res = read.table('GEM_results/GEM_Gmodel_SoC_p001.txt',header = T,sep = '\t')
  length(unique(gem.soc.g.res$cpg))
  gem.rnd.g.res = read.table('GEM_results/GEM_Gmodel_rnd_p001.txt',header = T,sep = '\t')
  length(unique(gem.rnd.g.res$cpg))
  gem.ctrl.g.res = read.table('GEM_results/GEM_Gmodel_ctrl_p001.txt',header = T,sep = '\t')
  length(unique(gem.ctrl.g.res$cpg))
  rm(gem.soc.g.res,gem.rnd.g.res,gem.ctrl.g.res)
  

  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  # Run GEM GxE-models ----
  # use E-term (sin or cos) determined in E-only models above
  # do this by running GEM separately for sin/cos-associated CpGs
  # note E-term must be last row of GEM covariate file

  # SoC-CpGs
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_sin.txt",
               "./data/soc_Mfil_sin.txt", 0.001,
               "./GEM_results/GEM_GEmodel_SoC_sindoctheta.txt",topKplot = 0, savePlot=F)
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_cos.txt",
               "./data/soc_Mfil_cos.txt", 0.001,
               "./GEM_results/GEM_GEmodel_SoC_cosdoctheta.txt",topKplot = 0, savePlot=F)
  
  # # check all cgs covered at this sig threshold
  # gem.soc.gxe.sin.res = read.table('./GEM_results/GEM_GEmodel_SoC_sindoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.soc.gxe.sin.res$cpg))
  # gem.soc.gxe.cos.res = read.table('./GEM_results/GEM_GEmodel_SoC_cosdoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.soc.gxe.cos.res$cpg))

  # # sanity check - GEM output should be same as equivalent lm interaction model
  # # check top CpG-SNP pair
  # gem.soc.gxe.sin.res = read.table('./GEM_results/GEM_GEmodel_SoC_sindoctheta.txt',header = T,sep = '\t')
  # covs_sin = read.table('./data/covs_sin.txt',header = T,sep = '\t')
  # gmb.snps = read.table('./data/GMB_SNPs.txt',header = T,sep = '\t')
  # 
  # top.cpg = as.character(gem.soc.gxe.sin.res[1,'cpg'])
  # top.snp = as.character(gem.soc.gxe.sin.res[1,'snp'])
  # # top.cpg = as.character(gem.soc.gxe.sin.res[2,'cpg'])
  # # top.snp = as.character(gem.soc.gxe.sin.res[2,'snp'])
  # snp = gmb.snps[top.snp,]
  # soc_Mfil = read.table(file = './data/soc_Mfil_sin.txt',header = T,sep = '\t')
  # m.cpg = soc_Mfil[top.cpg,]
  # 
  # design = t(rbind(m.cpg,covs_sin,snp))
  # colnames(design)[11] = 'sin.theta'
  # colnames(design)[12] = unlist(strsplit(top.snp,':'))[1]
  # 
  # lm.fit = lm(
  #   paste(
  #     top.cpg,'~',
  #     paste(
  #       colnames(design)[2:10],collapse = '+'
  #     ),'+',
  #     paste(
  #       colnames(design)[11:12],collapse = '*'
  #     )
  #   ),
  #   data = as.data.frame(design)
  # )
  # filter(gem.soc.gxe.sin.res,cpg==top.cpg & snp==top.snp)[,1:5]
  # summary(lm.fit)
  # # ** this confirms that GEM runs full interaction model including E and G terms as main effects **
  # 
  
  
  # ctrl CpGs
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_sin.txt",
               "./data/ctrl_Mfil_sin.txt", 0.001,
               "./GEM_results/GEM_GEmodel_ctrl_sindoctheta.txt",topKplot = 0, savePlot=F)
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_cos.txt",
               "./data/ctrl_Mfil_cos.txt", 0.001,
               "./GEM_results/GEM_GEmodel_ctrl_cosdoctheta.txt",topKplot = 0, savePlot=F)
  
  # # check all cgs covered at this sig threshold
  # gem.ctrl.gxe.sin.res = read.table('./GEM_results/GEM_GEmodel_ctrl_sindoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.ctrl.gxe.sin.res$cpg))
  # gem.ctrl.gxe.cos.res = read.table('./GEM_results/GEM_GEmodel_ctrl_cosdoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.ctrl.gxe.cos.res$cpg))

  # random CpGs
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_sin.txt",
               "./data/rnd_Mfil_sin.txt", 0.001,
               "./GEM_results/GEM_GEmodel_rnd_sindoctheta.txt",topKplot = 0, savePlot=F)
  GEM_GxEmodel("./data/GMB_SNPs.txt", "./data/covs_cos.txt",
               "./data/rnd_Mfil_cos.txt", 0.001,
               "./GEM_results/GEM_GEmodel_rnd_cosdoctheta.txt",topKplot = 0, savePlot=F)
  # # check all cgs covered at this sig threshold
  # gem.rnd.gxe.sin.res = read.table('./GEM_results/GEM_GEmodel_rnd_sindoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.rnd.gxe.sin.res$cpg))
  # gem.rnd.gxe.cos.res = read.table('./GEM_results/GEM_GEmodel_rnd_cosdoctheta.txt',header = T,sep = '\t')
  # length(unique(gem.rnd.gxe.cos.res$cpg))

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' create tables identifying top G and G x E effects for each CpG
#'
#' @return NULL
create_GEM_summary_results_tables = function(){
  
  setwd('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/')
  soc_model_sel = readRDS('./R_objects/soc_model_sel.RDS')
  ctrl_model_sel = readRDS('./R_objects/ctrl_model_sel.RDS')
  rnd_model_sel = readRDS('./R_objects/rnd_model_sel.RDS')
  
  # SoC-CpGs
  soc_res = data.frame(
    cpg = character(),
    which.E = character(),
    top.G=character(),
    top.G.pval=numeric(),
    top.G.fdr=numeric(),
    top.GxE=character(),
    top.GxE.pval=numeric(),
    top.GxE.fdr=numeric(),
    stringsAsFactors = F
  )
  
  
  # top.G
  res_G <- read.table("./GEM_results/GEM_Gmodel_SoC_p001.txt",head=T,stringsAsFactors = F) 
  res_G <- res_G[,c('cpg','snp','pvalue','FDR')]
  
  for (cpg in row.names(soc_model_sel)){ 
    G_snps <- res_G[which(res_G$cpg == cpg),]
    G_top_rownum <- which.min(G_snps$pvalue)
    soc_res[cpg,'cpg'] = cpg
    soc_res[cpg,'which.E'] = soc_model_sel[cpg,'which.E']
    soc_res[cpg,'top.G'] = G_snps[G_top_rownum,'snp']
    soc_res[cpg,'top.G.pval'] = G_snps[G_top_rownum,'pvalue']
    soc_res[cpg,'top.G.fdr'] = G_snps[G_top_rownum,'FDR']
  }
  rm(res_G)
  
  # top.GxE (combine results for CpGs with winning sin and cos terms)
  res_GE = rbind(
    read.table("./GEM_results/GEM_GEmodel_SoC_sindoctheta.txt",head=T,stringsAsFactors = F),
    read.table("./GEM_results/GEM_GEmodel_SoC_cosdoctheta.txt",head=T,stringsAsFactors = F)
  ) 
  res_GE <- res_GE[,c('cpg','snp','pvalue','FDR')]
  
  for (cpg in row.names(soc_model_sel)){ 
    GE_snps <- res_GE[which(res_GE$cpg == cpg),]
    GE_top_rownum <- which.min(GE_snps$pvalue) # note test p-val as FDRs computed separately for sin and cos
    soc_res[cpg,'top.GxE'] = GE_snps[GE_top_rownum,'snp']
    soc_res[cpg,'top.GxE.pval'] = GE_snps[GE_top_rownum,'pvalue']
    soc_res[cpg,'top.GxE.fdr'] = GE_snps[GE_top_rownum,'FDR']
  }
  rm(res_GE)
  
  # ctrl CpGs
  ctrl_res = data.frame(
    cpg = character(),
    which.E = character(),
    top.G=character(),
    top.G.pval=numeric(),
    top.G.fdr=numeric(),
    top.GxE=character(),
    top.GxE.pval=numeric(),
    top.GxE.fdr=numeric(),
    stringsAsFactors = F
  )
  
  # top.G
  res_G <- read.table("./GEM_results/GEM_Gmodel_ctrl_p001.txt",head=T,stringsAsFactors = F) 
  res_G <- res_G[,c('cpg','snp','pvalue','FDR')]
  
  for (cpg in row.names(ctrl_model_sel)){ 
    G_snps <- res_G[which(res_G$cpg == cpg),]
    G_top_rownum <- which.min(G_snps$pvalue)
    ctrl_res[cpg,'cpg'] = cpg
    ctrl_res[cpg,'which.E'] = ctrl_model_sel[cpg,'which.E']
    ctrl_res[cpg,'top.G'] = G_snps[G_top_rownum,'snp']
    ctrl_res[cpg,'top.G.pval'] = G_snps[G_top_rownum,'pvalue']
    ctrl_res[cpg,'top.G.fdr'] = G_snps[G_top_rownum,'FDR']
  }
  rm(res_G)
  
  # top.GxE (combine results for CpGs with winning sin and cos terms)
  res_GE = rbind(
    read.table("./GEM_results/GEM_GEmodel_ctrl_sindoctheta.txt",head=T,stringsAsFactors = F),
    read.table("./GEM_results/GEM_GEmodel_ctrl_cosdoctheta.txt",head=T,stringsAsFactors = F)
  ) 
  res_GE <- res_GE[,c('cpg','snp','pvalue','FDR')]
  
  # NB one ctrl SNP has no associated GxE SNP at recorded threshold (see GEM GxE above)
  for (cpg in row.names(ctrl_model_sel)){ 
    GE_snps <- res_GE[which(res_GE$cpg == cpg),]
    if (nrow(GE_snps)){
      GE_top_rownum <- which.min(GE_snps$pvalue)
      ctrl_res[cpg,'top.GxE'] = GE_snps[GE_top_rownum,'snp']
      ctrl_res[cpg,'top.GxE.pval'] = GE_snps[GE_top_rownum,'pvalue']
      ctrl_res[cpg,'top.GxE.fdr'] = GE_snps[GE_top_rownum,'FDR']
    }
  }
  rm(res_GE)

  
  # random CpGs
  rnd_res = data.frame(
    cpg = character(),
    which.E = character(),
    top.G=character(),
    top.G.pval=numeric(),
    top.G.fdr=numeric(),
    top.GxE=character(),
    top.GxE.pval=numeric(),
    top.GxE.fdr=numeric(),
    stringsAsFactors = F
  )
  
  # top.G
  res_G <- read.table("./GEM_results/GEM_Gmodel_rnd_p001.txt",head=T,stringsAsFactors = F) 
  res_G <- res_G[,c('cpg','snp','pvalue','FDR')]
  
  for (cpg in row.names(rnd_model_sel)){ 
    G_snps <- res_G[which(res_G$cpg == cpg),]
    G_top_rownum <- which.min(G_snps$pvalue)
    rnd_res[cpg,'cpg'] = cpg
    rnd_res[cpg,'which.E'] = rnd_model_sel[cpg,'which.E']
    rnd_res[cpg,'top.G'] = G_snps[G_top_rownum,'snp']
    rnd_res[cpg,'top.G.pval'] = G_snps[G_top_rownum,'pvalue']
    rnd_res[cpg,'top.G.fdr'] = G_snps[G_top_rownum,'FDR']
  }
  rm(res_G)
  
  # top.GxE (combine results for CpGs with winning sin and cos terms)
  res_GE = rbind(
    read.table("./GEM_results/GEM_GEmodel_rnd_sindoctheta.txt",head=T,stringsAsFactors = F),
    read.table("./GEM_results/GEM_GEmodel_rnd_cosdoctheta.txt",head=T,stringsAsFactors = F)
  ) 
  res_GE <- res_GE[,c('cpg','snp','pvalue','FDR')]
  
  for (cpg in row.names(rnd_model_sel)){ 
    GE_snps <- res_GE[which(res_GE$cpg == cpg),]
    GE_top_rownum <- which.min(GE_snps$pvalue)
    rnd_res[cpg,'top.GxE'] = GE_snps[GE_top_rownum,'snp']
    rnd_res[cpg,'top.GxE.pval'] = GE_snps[GE_top_rownum,'pvalue']
    rnd_res[cpg,'top.GxE.fdr'] = GE_snps[GE_top_rownum,'FDR']
  }
  rm(res_GE)
  
  # save all GEM results
  saveRDS(soc_res,'./R_objects/soc_topG_topE.RDS')
  saveRDS(rnd_res,'./R_objects/rnd_topG_topE.RDS')
  saveRDS(ctrl_res,'./R_objects/ctrl_topG_topE.RDS')
  

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Compare E, G and GxE models using 'winning' SNPs identified in GEM analysis
#'
#' @return NULL
#' @export
#'
#' @examples
compare_E_G_GxE_models = function(){
  
  setwd(dir = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/')
  
  # Load SNP and methylation data + covariates
  # SNPs
  gmb.snps = readRDS('./R_objects/GMB_SNPs_imputed_GEM.RDS')
  # covariates and predictors
  env_sin = read.table(file = './data/env_sin.txt',header = T,sep = '\t')
  env_cos = read.table(file = './data/env_cos.txt',header = T,sep = '\t')
  covs_only = read.table(file = './data/covs_only.txt',header = T,sep = '\t')
  covs_sin = read.table('./data/covs_sin.txt',header = T,sep = '\t')
  covs_cos = read.table('./data/covs_cos.txt',header = T,sep = '\t')
  # M-values
  soc.M = readRDS('./R_objects/soc_Mfil.RDS')
  ctrl.M = readRDS('./R_objects/ctrl_Mfil.RDS')
  rnd.M = readRDS('./R_objects/rnd_Mfil.RDS')
  # GEM results (winning SNPs for G and GxE models)
  gem.soc.res = readRDS('./R_objects/soc_topG_topE.RDS')
  gem.ctrl.res = readRDS('./R_objects/ctrl_topG_topE.RDS')
  gem.rnd.res = readRDS('./R_objects/rnd_topG_topE.RDS')
  
  # run regressions, choose best model, generate results table
  # design matrix for linear models incl covs and Fourier terms
  all.covs = rbind(covs_sin,env_cos) 
  
  fit.models <- function(gem.res = NULL,meth = NULL) {
    
    # gem.res is results object from GEM analysis:
    #    gem.soc.res - SoC CpGs
    #   gem.rand.res - RND controls
    #   gem.ctrl.res - HV controls
    # meth is meth object corresponding to gem.res
    
    message(paste('fitting E, G and GxE models for',deparse(substitute(gem.res))),'\n')
    
    for (cpg in gem.res$cpg){
      
      top.G = gem.res[cpg,'top.G']
      top.GxE = gem.res[cpg,'top.GxE']
      design = as.data.frame(
        t(
          rbind(
            meth[cpg,,drop=F],
            all.covs,
            gmb.snps[top.G,,drop=F],
            gmb.snps[top.GxE,,drop=F]
          )
        )
      ) # add M-values and G,GxE snps
      # remove special characters from snp names to make parseable for R forumula
      colnames(design)[which(colnames(design)=='sin(doc.theta)')] = 'sin.doc.theta'
      colnames(design)[which(colnames(design)=='cos(doc.theta)')] = 'cos.doc.theta'
      top.G.new = str_replace_all(top.G, "[^[:alnum:]]", "_")
      top.GxE.new = str_replace_all(top.GxE, "[^[:alnum:]]", "_")
      colnames(design)[which(colnames(design)==top.G)] = top.G.new
      colnames(design)[which(colnames(design)==top.GxE)] = top.GxE.new
      
      # fit NULL (cov-only) model
      lm.fit.null = lm(
        paste(
          cpg,'~',
          'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6'
        ), data = design
      )
      
      # fit G-model
      lm.fit.G = lm(
        paste(
          cpg,'~',
          'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +',
          top.G.new
        ), data = design
      )
      
      
      # fit E-model with most sig FT
      if (gem.res[cpg,'which.E'] == 'sin'){
        
        lm.fit.E = lm(
          paste(
            cpg,'~',
            'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + sin.doc.theta'
          ), data = design
        )
        
        lm.fit.GxE = lm(
          paste0(
            cpg,'~',
            'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + sin.doc.theta*',
            top.GxE.new
          ), data = design
        )
        
      } else if (gem.res[cpg,'which.E'] == 'cos'){
        
        lm.fit.E = lm(
          paste(
            cpg,'~',
            'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cos.doc.theta'
          ), data = design
        )
        
        if (!is.na(top.GxE.new)){ # one cpg has no signif GxE SNP at GEM threshold used for this analysis
          lm.fit.GxE = lm(
            paste0(
              cpg,'~',
              'SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + cos.doc.theta*',
              top.GxE.new
            ), data = design
          )
        } else {
          lm.fit.GxE = NA
        }
      }
      
      # adj R-squared
      gem.res[cpg,'null.adj.rsq'] = summary(lm.fit.null)$adj.r.squared
      gem.res[cpg,'E.adj.rsq'] = summary(lm.fit.E)$adj.r.squared
      gem.res[cpg,'G.adj.rsq'] = summary(lm.fit.G)$adj.r.squared
      if (!is.na(top.GxE.new)){ # one cpg has no signif GxE SNP at GEM threshold used for this analysis
        gem.res[cpg,'GxE.adj.rsq'] = summary(lm.fit.GxE)$adj.r.squared
      } else {
        gem.res[cpg,'GxE.adj.rsq'] = NA
      }
      
      # AIC
      gem.res[cpg,'null.AIC'] = AIC(lm.fit.null)
      gem.res[cpg,'E.AIC'] = AIC(lm.fit.E)
      gem.res[cpg,'G.AIC'] = AIC(lm.fit.G)
      if (!is.na(top.GxE.new)){ # one cpg has no signif GxE SNP at GEM threshold used for this analysis
        gem.res[cpg,'GxE.AIC'] = AIC(lm.fit.GxE)
      } else {
        gem.res[cpg,'GxE.AIC'] = NA
      }
      
    }
    
    # COMPARE MODELS
    
    gem.res$E.delta.adj.rsq = 
      as.numeric(gem.res$E.adj.rsq) - as.numeric(gem.res$null.adj.rsq)

    gem.res$G.delta.adj.rsq =
          as.numeric(gem.res$G.adj.rsq) - as.numeric(gem.res$null.adj.rsq)
    
    gem.res$GxE.delta.adj.rsq =
          as.numeric(gem.res$GxE.adj.rsq) - as.numeric(gem.res$null.adj.rsq)

    gem.res$winningR2 = 
      apply(
        gem.res[,c("E.adj.rsq","G.adj.rsq","GxE.adj.rsq")],1,
        function(x) {
          c("E.adj.rsq","G.adj.rsq","GxE.adj.rsq")[which.max(x)]
        }
      )
    
    gem.res$winningAIC = 
      apply(
        gem.res[,c("E.AIC","G.AIC","GxE.AIC")],1,
        function(x) {
          c("E.AIC","G.AIC","GxE.AIC")[which.min(x)]
        }
      )
    
    # summary stats
    message('** WINNING MODELS **\n')
    print(summary(as.factor(gem.res$winningR2)))
    print(summary(as.factor(gem.res$winningAIC)))
    cat('\n')
    
    return(gem.res)
    
  }
  
  # SoC-CpGs
  gem.soc.res = fit.models(gem.res = gem.soc.res, meth = soc.M)
  gem.ctrl.res = fit.models(gem.res = gem.ctrl.res, meth = ctrl.M)
  gem.rnd.res = fit.models(gem.res = gem.rnd.res, meth = rnd.M)
  
  # save results objects
  saveRDS(gem.soc.res,'./R_objects/gem.soc.res.RDS')
  saveRDS(gem.ctrl.res,'./R_objects/gem.ctrl.res.RDS')
  saveRDS(gem.rnd.res,'./R_objects/gem.rnd.res.RDS')

  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate location data for G1 and G2 SNPs
#'
#' @return
get_G1_G2_locations = function(){
  
  require(stringr)
  require(tidyverse)
  
  # load GSA manifest for missing SNP locaitons:
  gsa.manifest = read.csv(file = '/data/EMPHASIS_GSA/GSA-24v1-0_C1.csv',skip = 7, header = T)
  # load cg locations:
  cg.annot.df = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/R_objects/annotated.array.gr.RDS'))
  
  # SoC-CpGs ----
  gem.soc.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.RDS'))
  
  # Get CpG locations
  gem.soc.res$cpg.chr = left_join(gem.soc.res,cg.annot.df,by='cpg')$seqnames
  gem.soc.res$cpg.loc = left_join(gem.soc.res,cg.annot.df,by='cpg')$start
  
  # Get SNP locs  
  gem.soc.res$top.G.chr = 
    vapply(
      gem.soc.res$top.G, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.soc.res$top.G.loc = 
    vapply(
      gem.soc.res$top.G, FUN.VALUE = 1, function(x) as.numeric(unlist(strsplit(x,':'))[2])
    )
  gem.soc.res$top.GxE.chr = 
    vapply(
      gem.soc.res$top.GxE, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.soc.res$top.GxE.loc = 
    vapply(
      gem.soc.res$top.GxE, FUN.VALUE = 1, function(x) as.numeric(unlist(strsplit(x,':'))[2])
    )
  
  # get IDs corresponding to NA locations - these are all non-imputed loci since IMPUTE adds locs automatically
  snps.missing.locs.G = 
    gem.soc.res %>% 
    filter(is.na(top.G.loc)) %>%
    pull(top.G)
  
  snps.missing.locs.G.ids = str_replace(snps.missing.locs.G, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.G = filter(gsa.locs,Name %in% snps.missing.locs.G.ids)
  # check all locs found in GSA
  length(snps.missing.locs.G.ids)
  nrow(gsa.locs.for.snps.missing.locs.G)
  
  # add missing locs
  for (snp in snps.missing.locs.G.ids){
    gem.soc.res[grepl(pattern = snp, x = gem.soc.res$top.G),'top.G.loc'] = filter(gsa.locs.for.snps.missing.locs.G,Name==snp)$MapInfo
  }
  
  snps.missing.locs.GxE = 
    gem.soc.res %>% 
    filter(is.na(top.GxE.loc)) %>%
    pull(top.GxE)
  
  snps.missing.locs.GxE.ids = str_replace(snps.missing.locs.GxE, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.GxE = filter(gsa.locs,Name %in% snps.missing.locs.GxE.ids)
  # check all locs found in GSA
  length(snps.missing.locs.GxE.ids)
  nrow(gsa.locs.for.snps.missing.locs.GxE)
  
  # add missing locs
  for (snp in snps.missing.locs.GxE.ids){
    gem.soc.res[grepl(pattern = snp, x = gem.soc.res$top.GxE),'top.GxE.loc'] = filter(gsa.locs.for.snps.missing.locs.GxE,Name==snp)$MapInfo
  }
  
  sum(is.na(gem.soc.res)) # check no missing locs
  
  
  # ctrl CpGs ----
  gem.ctrl.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.ctrl.res.RDS'))
  
  # Get CpG locations
  gem.ctrl.res$cpg.chr = left_join(gem.ctrl.res,cg.annot.df,by='cpg')$seqnames
  gem.ctrl.res$cpg.loc = left_join(gem.ctrl.res,cg.annot.df,by='cpg')$start
  
  # Get SNP locs  
  gem.ctrl.res$top.G.chr = 
    vapply(
      gem.ctrl.res$top.G, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.ctrl.res$top.G.loc = 
    vapply(
      gem.ctrl.res$top.G, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,':'))[2]
    )
  gem.ctrl.res$top.GxE.chr = 
    vapply(
      gem.ctrl.res$top.GxE, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.ctrl.res$top.GxE.loc = 
    vapply(
      gem.ctrl.res$top.GxE, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,':'))[2]
    )

  # get IDs corresponding to NA locations - these are all non-imputed loci since IMPUTE adds locs automatically
  snps.missing.locs.G = 
    gem.ctrl.res %>% 
    filter(is.na(top.G.loc)) %>%
    pull(top.G)
  
  snps.missing.locs.G.ids = str_replace(snps.missing.locs.G, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.G = filter(gsa.locs,Name %in% snps.missing.locs.G.ids)
  # check all locs found in GSA
  length(unique(snps.missing.locs.G.ids))
  nrow(gsa.locs.for.snps.missing.locs.G)

  # add missing locs
  for (snp in snps.missing.locs.G.ids){
    gem.ctrl.res[grepl(pattern = snp, x = gem.ctrl.res$top.G),'top.G.loc'] = filter(gsa.locs.for.snps.missing.locs.G,Name==snp)$MapInfo
  }

  snps.missing.locs.GxE = 
    gem.ctrl.res %>% 
    filter(is.na(top.GxE.loc)) %>%
    pull(top.GxE)
  
  snps.missing.locs.GxE.ids = str_replace(snps.missing.locs.GxE, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.GxE = filter(gsa.locs,Name %in% snps.missing.locs.GxE.ids)
  # check all locs found in GSA
  length(unique(snps.missing.locs.GxE.ids))
  nrow(gsa.locs.for.snps.missing.locs.GxE) # note one cg has no GxE snp
  
  # add missing locs
  for (snp in snps.missing.locs.GxE.ids){
    if (!is.na(snp)){
      gem.ctrl.res[grepl(pattern = snp, x = gem.ctrl.res$top.GxE),'top.GxE.loc'] = filter(gsa.locs.for.snps.missing.locs.GxE,Name==snp)$MapInfo
    }
  }
  
  sum(is.na(gem.ctrl.res)) # check no missing locs

  
  # random CpGs ----
  gem.rand.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.rnd.res.RDS'))
  
  # Get CpG locations
  gem.rand.res$cpg.chr = left_join(gem.rand.res,cg.annot.df,by='cpg')$seqnames
  gem.rand.res$cpg.loc = left_join(gem.rand.res,cg.annot.df,by='cpg')$start
  
  # Get SNP locs  
  gem.rand.res$top.G.chr = 
    vapply(
      gem.rand.res$top.G, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.rand.res$top.G.loc = 
    vapply(
      gem.rand.res$top.G, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,':'))[2]
    )
  gem.rand.res$top.GxE.chr = 
    vapply(
      gem.rand.res$top.GxE, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,'_'))[1]
    )
  gem.rand.res$top.GxE.loc = 
    vapply(
      gem.rand.res$top.GxE, FUN.VALUE = 'chrx', function(x) unlist(strsplit(x,':'))[2]
    )
  
  # get IDs corresponding to NA locations - these are all non-imputed loci since IMPUTE adds locs automatically
  snps.missing.locs.G = 
    gem.rand.res %>% 
    filter(is.na(top.G.loc)) %>%
    pull(top.G)
  
  snps.missing.locs.G.ids = str_replace(snps.missing.locs.G, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.G = filter(gsa.locs,Name %in% snps.missing.locs.G.ids)
  # check all locs found in GSA
  length(snps.missing.locs.G.ids)
  nrow(gsa.locs.for.snps.missing.locs.G)
  
  # add missing locs
  # note chr4_rs5009180 is not in manifest - not sure why; add snp location info from dbSNP
  for (snp in snps.missing.locs.G.ids[snps.missing.locs.G.ids!='rs5009180']){
    gem.rand.res[grepl(pattern = snp, x = gem.rand.res$top.G),'top.G.loc'] = filter(gsa.locs.for.snps.missing.locs.G,Name==snp)$MapInfo
  }
  
  gem.rand.res[gem.rand.res$top.G=='chr4_rs5009180','top.G.loc'] = 86645644
    
  snps.missing.locs.GxE = 
    gem.rand.res %>% 
    filter(is.na(top.GxE.loc)) %>%
    pull(top.GxE)
  
  snps.missing.locs.GxE.ids = str_replace(snps.missing.locs.GxE, "chr[0-9]+_", "")
  # get these from GSA manifest
  gsa.locs = select(gsa.manifest,Name,Chr,MapInfo)
  gsa.locs$Name = as.character(gsa.locs$Name)
  gsa.locs.for.snps.missing.locs.GxE = filter(gsa.locs,Name %in% snps.missing.locs.GxE.ids)
  # check all locs found in GSA
  length(snps.missing.locs.GxE.ids)
  nrow(gsa.locs.for.snps.missing.locs.GxE)
  
  # add missing locs
  for (snp in snps.missing.locs.GxE.ids){
    gem.rand.res[grepl(pattern = snp, x = gem.rand.res$top.GxE),'top.GxE.loc'] = filter(gsa.locs.for.snps.missing.locs.GxE,Name==snp)$MapInfo
  }
  
  sum(is.na(gem.rand.res)) # check no missing locs


  
  # save G1 / G2 SNP annotations
  saveRDS(gem.soc.res,'~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.incl.locs.RDS')
  saveRDS(gem.ctrl.res,'~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.ctrl.res.incl.locs.RDS')
  saveRDS(gem.rand.res,'~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.rand.res.incl.locs.RDS')
  
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate annotation tables for G1 and G2 SNPs
#'
#' @return NULL
generate_G1_G2_annotation_tables = function(){
  
  require(GenomicRanges)
  require(tidyverse)
  
  gem.soc.res = read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.incl.locs.RDS')
  gem.ctrl.res = read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.ctrl.res.incl.locs.RDS')
  gem.rand.res = read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.rand.res.incl.locs.RDS')
  
  # generate G1 / G2 annotation tables
  g1.g2.table = gem.soc.res %>%
    select(
      cpg,cpg.chr,cpg.loc,
      winning.G1=top.G,winning.G1.chr=top.G.chr,winning.G1.loc=top.G.loc,winning.G1.pval=top.G.pval,winning.G1.fdr=top.G.fdr,
      winning.G2=top.GxE,winning.G2.chr=top.GxE.chr,winning.G2.loc=top.GxE.loc,winning.G2.pval=top.GxE.pval,winning.G2.fdr=top.GxE.fdr
    )
  # order by chr, location
  g1.g2.table = g1.g2.table %>%
    mutate_at(c('cpg.chr','winning.G1.chr','winning.G2.chr'),str_remove,'chr') %>%
    arrange(as.numeric(cpg.chr),as.numeric(cpg.loc))
  
  # clean up SNP IDs
  g1.snp.id = str_remove_all(g1.g2.table$winning.G1,'chr\\d++_')
  g1.snp.id = str_replace(g1.snp.id,"(rs\\d+).*",replacement = "\\1")
  g1.g2.table$winning.G1 = g1.snp.id
  g2.snp.id = str_remove_all(g1.g2.table$winning.G2,'chr\\d++_')
  g2.snp.id = str_replace(g2.snp.id,"(rs\\d+).*",replacement = "\\1")
  g1.g2.table$winning.G2 = g2.snp.id
  
  write.csv(
    g1.g2.table,
    file = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/results_tables/g1.g2.soc.cpg.results.csv',
  )
  
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Investigate G1 and G2 SNP locations
#' 
#' cpg-snp distances; cis-trans; clustering etc
#'
#' @param snp.clust.dist integer distance over which to check for G1/G2 SNP clusters
#'
#' @return
characterise_G1_G2_SNPs = function(snp.clust.dist=NULL){
  
  require(tidyverse)
  
  gem.soc.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.incl.locs.RDS'))
  gem.ctrl.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.ctrl.res.incl.locs.RDS'))
  gem.rand.res = as_tibble(read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.rand.res.incl.locs.RDS'))
  

  cat('\n***',length(unique(gem.soc.res$top.G)),
      'unique G1 SNPs associated with',nrow(gem.soc.res),'SoC-CpGs identified ***\n\n')
  cat('\n***',length(unique(gem.soc.res$top.GxE)),
      'unique G2 SNPs associated with',nrow(gem.soc.res),'SoC-CpGs identified ***\n\n')

  # cis-trans ----
  # check cis/trans distributions
  message('poportions of G1/G2 SNPs on same chromosome (defined as cis)\n')
  message('SoC-CpGs:\n')
  cat(
    sum(gem.soc.res$cpg.chr==gem.soc.res$top.G.chr)/nrow(gem.soc.res),
    'proprtion G1 in cis\n'
  )
  cat(
    sum(gem.soc.res$cpg.chr==gem.soc.res$top.GxE.chr)/nrow(gem.soc.res),
    'proportion G2 in cis\n\n'
  )
  
  message('ctrl CpGs:\n')
  cat(
    sum(gem.ctrl.res$cpg.chr==gem.ctrl.res$top.G.chr)/nrow(gem.ctrl.res),
    'proprtion G1 in cis\n'
  )
  cat(
    sum(gem.ctrl.res$cpg.chr==gem.ctrl.res$top.GxE.chr)/nrow(gem.ctrl.res),
    'proportion G2 in cis\n\n'
  )
  
  message('rand CpGs:\n')
  cat(
    sum(gem.rand.res$cpg.chr==gem.rand.res$top.G.chr)/nrow(gem.rand.res),
    'proprtion G1 in cis\n'
  )
  cat(
    sum(gem.rand.res$cpg.chr==gem.rand.res$top.GxE.chr)/nrow(gem.rand.res),
    'proportion G2 in cis\n\n'
  )

  
  message('ALL FOLLOWING REFER TO SOC-CPGS ONLY\n')
  
  g1.g2.table = read_csv(
    '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/results_tables/g1.g2.soc.cpg.results.csv'
  )
  
  # same G1 and G2 ----
  # check for same G1 and G2
  if (nrow(filter(g1.g2.table,winning.G1==winning.G2))==0){
    cat('\n*** No G1 and G2 SNPs match***\n\n')
  }

  # check for G1 or G2 driving more than one association
  # G1 influencing multiple CpGs ----
  g1.duplicates = table(g1.g2.table$winning.G1)[table(g1.g2.table$winning.G1)>1]
  cat('***',length(g1.duplicates),'mQTLs have more than 1 associated CpG:\n')
  print(
    g1.g2.table %>%
      filter(winning.G1 %in% names(g1.duplicates)) %>%
      arrange(winning.G1) %>%
      select(cpg,cpg.chr,cpg.loc,winning.G1,winning.G1.chr,winning.G1.loc,winning.G1.pval,winning.G1.fdr)
  )
  cat('\n')
  
  # G2 influencing multiple CpGs ----
  g2.duplicates = table(g1.g2.table$winning.G2)[table(g1.g2.table$winning.G2)>1]
  cat('***',length(g2.duplicates),'GxE SNPs have more than 1 associated CpG. Number of associated CpGs as follows:\n')
  print(
    g1.g2.table %>%
      filter(winning.G2 %in% names(g2.duplicates)) %>%
      arrange(winning.G2) %>%
      select(cpg,cpg.chr,cpg.loc,winning.G2,winning.G2.chr,winning.G2.loc,winning.G2.pval,winning.G2.fdr)
  )
  cat('\n')
  
  # check for SNP clusters driving multiple replicated signals
  # G1 SNP clusters ----
  cat('\n\nCheck G1 clusters within',snp.clust.dist,'bp\n')
  g1.g2.table$winning.G1.loc.end = g1.g2.table$winning.G1.loc + 1
  g1.gr =
    makeGRangesFromDataFrame(filter(g1.g2.table,!is.na(winning.G1.chr)),
                             keep.extra.columns=TRUE,
                             ignore.strand=TRUE,
                             seqinfo = levels(g1.g2.table$winning.G1.chr)[2:length(levels(g1.g2.table$winning.G1.chr))],
                             seqnames.field="winning.G1.chr",
                             start.field="winning.G1.loc",
                             end.field="winning.G1.loc.end"
    )
  
  g1.gr.clust = GenomicRanges::reduce(g1.gr,min.gapwidth=snp.clust.dist)
  g1.gr.clust$g1.width = width(g1.gr.clust)
  g1.gr.clust.df = data.frame(g1.gr.clust)
  
  cat(length(unique(g1.gr$winning.G1)),'G1 SNPs reduce to',length(g1.gr.clust),'clusters within',snp.clust.dist,'bp\n')
  cat('CpGs overlapping SNP clusters using an inter-snp distance of',snp.clust.dist,'bp:\n')
  
  # print SNPs and CpGs in identified clusters
  cg.clust.all.df = 
    data.frame(
      cluster=numeric(),
      cpg=character(), 
      cpg.chr=character(),
      cpg.loc=numeric(),
      winning.G1=character(),
      winning.G1.chr=character(),
      winning.G1.loc=numeric()
    )
  for (row in seq(length(g1.gr.clust[g1.gr.clust$g1.width>2]))){
    clust = g1.gr.clust[g1.gr.clust$g1.width>2][row]
    cg.clust = g1.gr[subjectHits(findOverlaps(clust,g1.gr))]
    cg.clust$winning.G1.chr = seqnames(cg.clust)
    cg.clust$winning.G1.loc = start(cg.clust)
    clust.df = as.data.frame(cg.clust) %>%
      select(
        cpg,cpg.chr,cpg.loc,winning.G1,winning.G1.chr,winning.G1.loc
      )
    cg.clust.all.df = rbind(cg.clust.all.df,cbind(cluster=row,clust.df))
  }
  print(
    cg.clust.all.df %>%
      arrange(cluster,winning.G1.loc),
    row.names = F
  )

  # G2 SNP clusters ----
  # print SNPs and CpGs in identified clusters
  cat('\n\nCheck G1 clusters within',snp.clust.dist,'bp\n')
  # G2xE
  g1.g2.table$winning.G2.loc.end = g1.g2.table$winning.G2.loc + 1
  g2.gr =
    makeGRangesFromDataFrame(filter(g1.g2.table,!is.na(winning.G2.chr)),
                             keep.extra.columns=TRUE,
                             ignore.strand=TRUE,
                             seqinfo = levels(g1.g2.table$winning.G2.chr)[2:length(levels(g1.g2.table$winning.G2.chr))],
                             seqnames.field="winning.G2.chr",
                             start.field="winning.G2.loc",
                             end.field="winning.G2.loc.end"
    )
  g2.gr.clust = GenomicRanges::reduce(g2.gr,min.gapwidth=snp.clust.dist)
  g2.gr.clust$g2.width = width(g2.gr.clust)
  g2.gr.clust.df = data.frame(g2.gr.clust)
  cat(length(unique(g2.gr$winning.G2)),'G2 SNPs reduce to',length(g2.gr.clust),'clusters within',snp.clust.dist,'bp\n')
  cat('CpGs overlapping SNP clusters using an inter-snp distance of',snp.clust.dist,'bp:\n')
  
  # print SNPs and CpGs in identified clusters
  cg.clust.all.df = 
    data.frame(
      cluster=numeric(),
      cpg=character(), 
      cpg.chr=character(),
      cpg.loc=numeric(),
      winning.G2=character(),
      winning.G2.chr=character(),
      winning.G2.loc=numeric()
    )
  for (row in seq(length(g2.gr.clust[g2.gr.clust$g2.width>2]))){
    clust = g2.gr.clust[g2.gr.clust$g2.width>2][row]
    cg.clust = g2.gr[subjectHits(findOverlaps(clust,g2.gr))]
    cg.clust$winning.G2.chr = seqnames(cg.clust)
    cg.clust$winning.G2.loc = start(cg.clust)
    clust.df = as.data.frame(cg.clust) %>%
      select(
        cpg,cpg.chr,cpg.loc,winning.G2,winning.G2.chr,winning.G2.loc
      )
    cg.clust.all.df = rbind(cg.clust.all.df,cbind(cluster=row,clust.df))
  }

  print(
    cg.clust.all.df %>%
      arrange(cluster,winning.G2.loc),
    row.names = F
  )
  
    
}




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Test association between all G1 and G2 SNPs and SoC
#'
#' @return NULL
G1_G2_SoC_association_tests = function(){
  
  require(SNPassoc)
  require(tidyverse)
  
  # get G1/G2 SNPs
  replicated.all.locs = read_rds('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.incl.locs.RDS')
  G1.snps = unique(replicated.all.locs$top.G)
  G2.snps = unique(replicated.all.locs$top.GxE)
  snps.of.interest = unique(union(G1.snps,G2.snps))
  
  # get genotypes
  all.snps = as_tibble(
    read_rds( '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/GMB_SNPs_imputed_GEM.RDS'),
    rownames = 'snp'
  )
  all.snps = all.snps %>% filter(snp %in% snps.of.interest)
  # convert to SNPassoc format
  snp.names = all.snps$snp
  all.snps$snp = NULL
  all.snps = t(all.snps)
  colnames(all.snps) = snp.names
  all.snps = as_tibble(all.snps,rownames = 'ID')
  all.snps[all.snps==0] = 'A/A'; all.snps[all.snps==1] = 'A/T'; all.snps[all.snps==2] = 'T/T'
  all.snps = mutate_all(all.snps,as.factor)
  g1.snps.geno = all.snps[,c('ID',G1.snps)]
  g2.snps.geno = all.snps[,c('ID',G2.snps)]
  cat('genotype data for',ncol(g1.snps.geno),'G1 and',ncol(g2.snps.geno),'G2 SNPs loaded\n')

  # get SoC data - use dichotomised 'Moore' SoC
  emph.pdata = readRDS('~/projects/soc_enid_emphasis/cohort_data/EMPHASIS_pdata.RDS')
  emph.pdata = emph.pdata %>%
    select(ID=Subject_ID,MooreSoC)
  
  # match genotypes to indivdiuals with methylation data (4 more genotyped subjects than methylation)
  g1.snps.geno = g1.snps.geno[match(emph.pdata$ID,g1.snps.geno$ID),]
  g1.data = as_tibble(left_join(emph.pdata,g1.snps.geno))
  g1.data = setupSNP(g1.data,colSNPs = seq(3,ncol(g1.data)),sep = '/')
  g2.snps.geno = g2.snps.geno[match(emph.pdata$ID,g2.snps.geno$ID),]
  g2.data = as_tibble(left_join(emph.pdata,g2.snps.geno))
  g2.data = setupSNP(g2.data,colSNPs = seq(3,ncol(g2.data)),sep = '/')
  
  # association analysis ----
  g1.res = WGassociation(MooreSoC,data = g1.data,quantitative = FALSE)
  g2.res = WGassociation(MooreSoC,data = g2.data,quantitative = FALSE)
  
  # add HWE
  g1.hwe = tableHWE(g1.data)
  g2.hwe = tableHWE(g2.data)

  g1.res.tab = cbind(snp.type='G1',g1.res,hwe.p=g1.hwe[row.names(g1.res),])
  g2.res.tab = cbind(snp.type='G2',g2.res,hwe.p=g2.hwe[row.names(g2.res),])
  g1.res.tab$comments = NULL
  g2.res.tab$comments = NULL
  res.tab.all = rbind(g1.res.tab,g2.res.tab)
  
  # clean up SNP IDs
  snp.id = str_remove_all(row.names(res.tab.all),'chr\\d++_')
  snp.id = str_replace(snp.id,"(rs\\d+).*",replacement = "\\1")
  row.names(res.tab.all) = snp.id
  
  # plot p-value distributions
  g1.pvals = res.tab.all %>%
    filter(snp.type=='G1') %>%
    select(-(c('hwe.p','snp.type')))
  hist(unlist(g1.pvals))
  g2.pvals = res.tab.all %>%
    filter(snp.type=='G2') %>%
    select(-(c('hwe.p','snp.type')))
  hist(unlist(g2.pvals))
  
  # print min p-vals for each model
  cat('G1 min p-vals:\n',min(g1.pvals),'\n')
  cat('G2 min p-vals:\n',min(g2.pvals),'\n')

  # QQ plot
  plot.dir = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/plots/'
  plt.width = 8 # plot width in cms
  plt.height = 6 # plot height in cms
  plt.dpi = 600
  source('~/src/utils/misc/qqPlot.R')
  gg_qqplot(unlist(g1.pvals),base.size = 10)
  ggsave(filename = paste0(plot.dir,"G1_qqPlot.pdf"),
         dpi = plt.dpi,width = plt.width,height = plt.height,units = 'cm'
  )
  gg_qqplot(unlist(g2.pvals))
  ggsave(filename = paste0(plot.dir,"G2_qqPlot.pdf"),
         dpi = plt.dpi,width = plt.width,height = plt.height,units = 'cm'
  )

  # save
  write.csv(res.tab.all,file = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/results_tables/g1_g2_soc_assoc_res.csv')
  
  
}

















