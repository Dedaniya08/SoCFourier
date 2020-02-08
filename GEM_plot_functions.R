#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# produce plots from GEM analysis
# functions called from SoCFourier_GEM_analysis.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Define plot parameters, size, DPI etc
#'
#' @param plot.dir directory to save plots
#'
#' @return list of plot parameters
get_figure_params = function(plot.dir=NULL,width=NULL,height=NULL,dpi=600){
  
  dir = plot.dir
  width = 8 # plot width in cms
  height = 6 # plot height in cms
  dpi = 600
  return(
    list(dir=dir,width=width,height=height,dpi=dpi)
  )
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot % of winning models that are E, G, GxE for all cpg sets
#'
#' @param soc.res DF SoC-CpG GEM modelling results
#' @param ctrl.res DF matched control GEM modelling results
#' @param rnd.res DF random control GEM modelling results
#' @param plot.params list of plot dimensions etc
#'
#' @return NULL
winning_model_pie_charts = function(soc.res=NULL,ctrl.res=NULL,rnd.res=NULL,plot.params = NULL){
  
  require(tidyverse)
  
  # plot colours
  colpal1 <- c("E" = "#FBB113","G1" = "#501C34","G2xE" = "#FD6621")
  colpal2 <- c("0" = "#014A33","1" = "#B33685","2" = "#42BFFF")
  
  percent_win_model = 
    as.data.frame(
      rbind(
        c(win.mod = 'SoC-CpG',prop.table(table(soc.res$winningAIC))),
        c(win.mod = 'matched\ncontrol',prop.table(table(ctrl.res$winningAIC))),
        c(win.mod = 'random\ncontrol',prop.table(table(rnd.res$winningAIC)))
      )
    )
  
  percent_win_model = 
    percent_win_model %>% 
      select(win.mod,G1=G.AIC,G2xE=GxE.AIC) %>%
      gather(G1,G2xE,key = 'model',value='percent') %>%
      select(cpgset = win.mod,model,percent)
  
  percent_win_model$cpgset = fct_relevel(factor(percent_win_model$cpgset,levels=c('SoC-CpG','matched\ncontrol','random\ncontrol')))
  
  
  print(
    ggplot(percent_win_model, aes(x="",y=as.numeric(percent),fill=model)) +
      geom_bar(width = 1, stat = "identity") + 
      coord_polar("y", start=0) +
      facet_wrap(cpgset ~ .) +  
      scale_fill_manual("winning model", values = colpal1) +
      xlab("") +
      ylab("") +
      theme_minimal() +
      theme(axis.text.x=element_blank())
  )
  
  ggsave(filename = paste0(plot.params$dir,"G1vG2E_pie.pdf"),
         dpi = plot.params$dpi,width = 1.2*plot.params$width,height = 1.2*plot.params$height,units = 'cm'
  )
  
  print(percent_win_model)
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Barplots of delta adjusted R2 from null explained by each model for all cpg sets
#'
#' @param soc.res DF SoC-CpG GEM modelling results
#' @param ctrl.res DF matched control GEM modelling results
#' @param rnd.res DF random control GEM modelling results
#' @param plot.params list of plot dimensions etc
#'
#' @return NULL
delta_adj_r2_barplots = function(soc.res=NULL,ctrl.res=NULL,rnd.res=NULL,plot.params = NULL){

  require(tidyverse)
  require(bootstrap)

  # plot colours
  colpal1 <- c("E" = "#FBB113","G1" = "#501C34","G2xE" = "#FD6621")
  colpal2 <- c("0" = "#014A33","1" = "#B33685","2" = "#42BFFF")
  
  # get adj R2 values
  ALL_lm_res_rsq = 
    as_tibble(
      rbind(
        data.frame(win.mod = 'SoC-CpG',soc.res[,c("null.adj.rsq","E.adj.rsq","G.adj.rsq","GxE.adj.rsq")]),
        data.frame(win.mod = 'matched\ncontrol',ctrl.res[,c("null.adj.rsq","E.adj.rsq","G.adj.rsq","GxE.adj.rsq")]),
        data.frame(win.mod = 'random\ncontrol',rnd.res[,c("null.adj.rsq","E.adj.rsq","G.adj.rsq","GxE.adj.rsq")])
      )
    )
  
  # compute delta adj R2 values
  delta.adj.r2 =
    ALL_lm_res_rsq %>%
      mutate(
        E = E.adj.rsq-null.adj.rsq,
        G1 = G.adj.rsq-null.adj.rsq,
        G2xE = GxE.adj.rsq-null.adj.rsq
      ) %>%
      select(win.mod,E,G1,G2xE)
  
  # generate bootstraps for median delta r-squared
  ci.95 <- function(x){quantile(x, c(.05,.95))}
  
  bootstrap.median = function(test.model=NULL,locus=NULL){
    adj.rsq = delta.adj.r2 %>% filter(win.mod==locus) %>% pull(test.model)
    median.val = median(adj.rsq,na.rm = T)
    boot.ci = bootstrap::bootstrap(x = na.exclude(adj.rsq),theta = median,nboot = 1000,func = ci.95)$func.thetastar
    # (exclude single missing value in ctrl.HV GxE set)
    # cat('locus:',locus,'; test.model:',test.model,'; median:',median.val,'; boot.ci:',boot.ci,'\n')
    return(list(
              locus = locus,
              test.model = test.model,
              median=median.val,
              ci.lower = as.vector(boot.ci[1]),
              ci.upper = as.vector(boot.ci[2])
          )
    )
  }
  
  res.df = data.frame(
    locus = character(),
    test.model = character(),
    median = numeric(),
    ci.lower = numeric(),
    ci.upper = numeric()
  )
  
  for (locus in levels(delta.adj.r2$win.mod)){
    for (test.model in c('E','G1','G2xE')){
      res.df = rbind(res.df,data.frame(bootstrap.median(test.model = test.model,locus = locus)))
    }
  }
  
  print(res.df)
  
  # plot
  print(
    ggplot(data = res.df, aes(x = locus, y = 100*median, fill = test.model)) +
    geom_col(position = position_dodge()) +
    geom_errorbar(aes(ymin = 100*ci.lower, ymax = 100*ci.upper), width=0.2, position = position_dodge(0.9),colour='black') +
    labs(
      fill='',
      x='',
      y=expression("median "*Delta~adjR^2)
    ) +
    scale_fill_manual("model",values = colpal1) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(colour = "black",size = 10),
      axis.text.y = element_text(colour = "black")
    )
  )
  
  ggsave(filename = paste0(plot.params$dir,"GxE_delta_adj_Rsquared.pdf"),
         dpi = plot.params$dpi,width = 1.25*plot.params$width,height = plot.params$height,units = 'cm'
  )

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot SoC Fourier curves stratified by genotype
#'
#' @param soc.res DF SoC-CpG GEM modelling results
#'
#' @return NULL
plot_soc_fourier_stratified_by_genotype = function(soc.res=NULL,plot.params=NULL,log.to.file=F,sink.file=NULL){
  
  require(tidyverse)
  require(gridExtra)
  require(cowplot)

  setwd('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/')  
  
  if (log.to.file) sink(sink.file)
  
  # load genotypes
  if (!exists('gmb.snps')) gmb.snps  = readRDS('./R_objects/GMB_SNPs_imputed_GEM.RDS')

  # M-values
  if (!exists('soc_Mfil')) soc_Mfil = readRDS('./R_objects/soc_Mfil.RDS')
  
  # SoC Fourier modelling results
  soc.res = readRDS('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/R_objects/gem.soc.res.RDS')
  
  # covariates and predictors
  env_sin = read.table(file = './data/env_sin.txt',header = T,sep = '\t')
  env_cos = read.table(file = './data/env_cos.txt',header = T,sep = '\t')
  covs_only = read.table(file = './data/covs_only.txt',header = T,sep = '\t')
  covs_sin = read.table('./data/covs_sin.txt',header = T,sep = '\t')
  covs_cos = read.table('./data/covs_cos.txt',header = T,sep = '\t')
  # create design matrix
  soc_min_cov = cbind(t(covs_sin),t(env_cos))
  # add doc.theta
  pdata = readRDS('../../cohort_data/EMPHASIS_pdata.RDS')
  soc_min_cov = 
    cbind(
      soc_min_cov,
      doc.theta=pdata[match(row.names(soc_min_cov),pdata$Subject_ID),'doc.theta']
    )
  
  colnames(soc_min_cov)[10:11] = c('sin(doc.theta)','cos(doc.theta)')
  
  # 1. SoC-CpGs: winning GxE
  
  soc.win.gxe = 
    soc.res %>%
    filter(winningAIC=='GxE.AIC') %>%
    select(cpg,which.E,top.GxE,GxE.delta.adj.rsq)
  
  # pick cg-snp pairs with highest delta adj r2 for gxe model
  soc.win.gxe.top = arrange(soc.win.gxe,desc(GxE.delta.adj.rsq))[1:20,]
  soc.gxe.cgs = soc.win.gxe.top$cpg
  
  for (cg in soc.gxe.cgs){

    snp = as.character(filter(soc.res,cpg==cg)$top.GxE)
    combined.stratified.socFourier.plot(
      cg = cg,snp = snp,g.model = 'G2xE', title.add = 'repl_GxE_winning',plot.params = plot.params
    )
    
  }
  
  # 2. SoC-CpGs: winning G (mQTL)
  soc.win.g = 
    soc.res %>%
    filter(winningAIC=='G.AIC') %>%
    select(cpg,which.E,top.G,G.delta.adj.rsq)
  
  # pick cg-snp pairs with highets delta adj r2 for gxe model; pick only one of the 3 cgs mapping to IGF1R
  soc.win.g.top = arrange(soc.win.g,desc(G.delta.adj.rsq))[1:10,]
  soc.g.cgs = soc.win.g.top$cpg
  
  for (cg in soc.g.cgs){
    
    snp = as.character(filter(soc.res,cpg==cg)$top.G)
    combined.stratified.socFourier.plot(
      cg = cg,snp = snp,g.model = 'G1.E', title.add = 'repl_G_winning',plot.params = plot.params
    )
    
  }
  
  if (log.to.file) sink()  
  
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate plots visualising SoC Fourier curves stratified by genotype
#' 
#' called from plot_soc_fourire_stratified_by_genotype()
#'
#' @param cg 
#' @param snp 
#' @param g.model 
#' @param title.add 
#' @param dominant.model 
#'
#' @return NULL
combined.stratified.socFourier.plot = function(cg=NULL,snp=NULL,g.model=NULL,title.add = NULL,plot.params=NULL){

  plot.title=paste0(title.add,'_',cg)

  data.df.G = get_MxG_data(cg = cg, snp = snp, model = g.model)
  colnames(data.df.G)[which(colnames(data.df.G)==snp)] = 'snp'
  data.df.E = get_MxG_data(cg = cg, snp = snp, model = 'E')
  data.df.E$snp = 'E'
  data.df = rbind(data.df.G,data.df.E)
  
  # simplify SNP ID
  snp.id = str_remove(snp,'chr\\d++_')
  snp.id = str_replace(snp.id,"(rs\\d+).*",replacement = "\\1")

  plt.1 = 
    ggplot(data.df, aes(x=doc.theta,y=100 * betahat, color=snp)) + 
    geom_line(aes(linetype=snp)) +
    scale_x_continuous(name = 'month of conception',
                       breaks = seq(0+pi/12,2*pi-pi/12,length.out = 12),
                       labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
    {
      if(!any(data.df.G$snp==2)) {
        scale_color_manual(
          name=paste0(snp.id,'\nallele count'),
          labels=c(
            paste0('0 (n=',sum(data.df$snp==0),')'),
            paste0('1 (n=',sum(data.df$snp==1),')'),
            'E'
          ),
          values = c(colpal2[1:2],'E'='red')
        )
      } else {
        scale_color_manual(
          name=paste0(cg,'\n',snp.id),
          labels=c(
            paste0('AA (n=',sum(data.df$snp==0),')'),
            paste0('Aa (n=',sum(data.df$snp==1),')'),
            paste0('aa (n=',sum(data.df$snp==2),')'),
            'E'
          ),
          values = c(colpal2,'E'='red')
        )
      } 
    } +
    scale_linetype_manual(
      name=snp.id,
      values = c(rep("dotdash", 3), "solid"),
      guide=FALSE
    ) +
    guides(colour = guide_legend(override.aes = list(linetype = c(rep("dotdash", 3), "solid")))) +
    labs(
      x = "date of conception",
      y = "% methylation\n(mean-centred)"
      # title = 'full model'
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(
        colour = c(rep(x = '#FF9933',6),rep('#33CC66',6))
      )
    )

  # scatter plot with dichotomise SoC using meth beta-values adjusted for
  # covariates only - i.e. without fourier terms from SoC model
  data.df.cov.adj = get_MxG_data(cg = cg, snp = snp, model = 'G.cov.adj')
  colnames(data.df.cov.adj)[which(colnames(data.df.cov.adj)==snp)] = 'snp'
  data.df.cov.adj[data.df.cov.adj$doc.theta<=pi,'SoC'] = 'dry'
  data.df.cov.adj[data.df.cov.adj$doc.theta>pi,'SoC'] = 'rainy'

  # irrespective of SNP allele count
  plt.3 = 
  ggplot(data.df.cov.adj, aes(x=SoC,y=100 * betahat, colour=SoC)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.25),size=0.5) +
    stat_summary(
      data = subset(data.df.cov.adj,SoC=='rainy'),colour='black',
      fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.25, size=0.25
                 # position = position_nudge(x = 0.18)
    ) +
    stat_summary(
      data = subset(data.df.cov.adj,SoC=='dry'),colour='black',
      fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.25, size=0.25
                 # position = position_nudge(x = -0.18)
    ) +
    scale_color_manual(
      name='SoC',
      values = c('dry'='#FF9933','rainy'='#33CC66')
    ) +
    labs(
      y='% methylation'
      # title = 'un-modelled adjusted beta values'
    ) +
    # coord_equal(ratio = 2) +
    theme_bw(base_size = 10) + 
    theme(
      legend.position = 'none',
      axis.title.x=element_blank(),
      axis.text.x = element_text(size=11,colour='black')
      # aspect.ratio = 2
    )
  
  # stratified by allelecount
  plt.4 = 
  ggplot(data.df.cov.adj, aes(x=snp,y=100 * betahat, colour=SoC)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.25),size=0.5) +
    stat_summary(
      data = subset(data.df.cov.adj,SoC=='rainy'),colour='black',
      fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.25, size=0.25,
                 position = position_nudge(x = 0.18)
    ) +
    stat_summary(
      data = subset(data.df.cov.adj,SoC=='dry'),colour='black',
      fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                 geom = "crossbar", width = 0.25, size=0.25,
                 position = position_nudge(x = -0.18)
    ) +
    scale_color_manual(
      name='SoC',
      values = c('dry'='#FF9933','rainy'='#33CC66')
    ) +
    labs(
      x='allele count'
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.y = element_blank()
    )
  

  plt = plot_grid(plt.1,plt.3,plt.4, align = "h", nrow = 1, rel_widths = c(1/2, 3/16, 1/3))
  print(plt)

  ggsave(filename = paste0(plot.params$dir,'stratified_fourier_plots/',plot.title,'.pdf'),
         dpi = plot.params$dpi,width = 20,height = 5,units = 'cm'
  )
  
    
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate methylation x genotype data for stratified SoCFourier plots
#' 
#' called from combined.stratified.socFourier.plot()
#'
#' @param cg 
#' @param snp 
#' @param model 
#'
#' @return DF with fitted values stratified by genotype
get_MxG_data = function(cg=NULL,snp=NULL,model=NULL){
  
  source('~/src/soc_multi_cohort_analysis/disc_repl/SoCFourier_modelling_functions.R') # m2beta
  
  cat(
    '****\n',
    'ANALYSIS FOR: ',cg,'and',snp,'\n'
  )
  
  cg_cpg <- as.data.frame(t(soc_Mfil[cg,,drop=F]))
  cg_cov <- as.data.frame(
    cbind(
      soc_min_cov,
      cg_cpg[row.names(soc_min_cov),cg]
    )
  )
  colnames(cg_cov)[ncol(cg_cov)] = cg
  
  if (grepl(pattern = 'G',model)){
    
    cg_cov = cbind(cg_cov,t(gmb.snps[snp,row.names(cg_cov)]))
    cg_cov[,snp] <- as.numeric(cg_cov[,snp])
    
    # remove special characters from snp names to make parseable for R forumula
    snp.new = str_replace_all(snp, "[^[:alnum:]]", "_")
    colnames(cg_cov)[which(colnames(cg_cov)==snp)] = snp.new
    snp.orig = snp
    snp = snp.new    
    
    row.names(cg_cov) = NULL
    
  }  
  
  # run regression models and obtain fitted methylation values
  if (model=='E'){
    
    # generate fitted M-values adjusted for covariates
    lm.formula = as.formula(
      paste0(
        cg,'~','SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6',
        '+ `sin(doc.theta)` + `cos(doc.theta)`'
      )
    )
    lm.cg <- lm(lm.formula,cg_cov)
    print(summary(lm.cg))
    
    # get fitted M and beta
    cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`")] %*% 
      t(model.matrix(lm.cg)[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`")])
    
  } else if (model=='G.cov.adj'){
    
    # generate M-value residuals adjusted for covariates only
    lm.formula = as.formula(
      paste0(
        cg,'~','SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6'
      )
    )
    lm.cg <- lm(lm.formula,cg_cov)
    print(summary(lm.cg))
    
    # get fitted M and beta
    # unadjusted
    # cg_Mfit = t(cg_cpg)
    # row.names(cg_Mfit) = NULL
    # cov adjusted
    cg_Mfit <- t(as.data.frame(coef(lm.cg)[1] + lm.cg$residuals))
    row.names(cg_Mfit) = NULL
    
  } else if (model=='G1.E'){
    
    # model effect of G1 (mQTL) and E only
    lm.formula = as.formula(
      paste0(
        cg,'~','SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6',
        '+ `sin(doc.theta)` + `cos(doc.theta)` + `', snp,'`'
      )
    )
    lm.cg <- lm(lm.formula,cg_cov)
    print(summary(lm.cg))
    
    # get fitted M and beta
    if (grepl(pattern = '-',snp)){
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                   paste0('`',snp,'`'))] %*%
        t(model.matrix(lm.cg)[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                 paste0('`',snp,'`'))])
    } else {
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                   snp)] %*% 
        t(model.matrix(lm.cg)[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                 snp)])
    }
    
  } else if (model=='G1'){
    
    # model effect of G1 (mQTL) only
    lm.formula = as.formula(
      paste0(
        cg,'~','SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6',
        '+ `', snp,'`'
      )
    )
    lm.cg <- lm(lm.formula,cg_cov)
    print(summary(lm.cg))
    
    # get fitted M and beta
    if (grepl(pattern = '-',snp)){
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)", paste0('`',snp,'`'))] %*%
        t(model.matrix(lm.cg)[,c("(Intercept)",paste0('`',snp,'`'))])
    } else {
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)",snp)] %*% 
        t(model.matrix(lm.cg)[,c("(Intercept)",snp)])
    }
    
  } else if (grepl('x',model)){ # interaction model
    
    lm.formula = as.formula(
      paste0(
        cg,'~','SexM + Age + MasterGroupNo2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6',
        '+ `sin(doc.theta)` + `cos(doc.theta)` + `', snp,'`+ `sin(doc.theta)`:`',snp,'`',
        '+ `cos(doc.theta)`:`',snp,'`'
      )
    )
    
    lm.cg <- lm(lm.formula,cg_cov)
    print(summary(lm.cg))
    
    # get fitted M and beta (explained by season + genotype + their interaction only) 
    if (grepl(pattern = '-',snp)){
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                   paste0('`',snp,'`'),paste0("`sin(doc.theta)`:",'`',snp,'`'),paste0("`cos(doc.theta)`:`",snp,'`'))] %*%
        t(model.matrix(lm.cg)[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                 paste0('`',snp,'`'),paste0("`sin(doc.theta)`:",'`',snp,'`'),paste0("`cos(doc.theta)`:`",snp,'`'))])
    } else {
      cg_Mfit <- t(coef(lm.cg))[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                   snp,paste0("`sin(doc.theta)`:",snp),paste0("`cos(doc.theta)`:",snp))] %*% 
        t(model.matrix(lm.cg)[,c("(Intercept)","`sin(doc.theta)`","`cos(doc.theta)`",
                                 snp,paste0("`sin(doc.theta)`:",snp),paste0("`cos(doc.theta)`:",snp))])
    }
    
    
  } else {
    
    cat('model not supported!\n')
    
  }
  
  cg_Bfit <- as.data.frame(t(m2beta(cg_Mfit)))
  # mean-centre (modelled betas only)
  if (model!='G.cov.adj') cg_Bfit$V1 = scale(cg_Bfit$V1,center = T,scale = F)
  
  if (grepl(pattern = 'G',model)){
    cg_Bfit <- left_join(rownames_to_column(cg_cov[,snp,drop=F]),
                         rownames_to_column(cg_Bfit), by="rowname")
    rownames(cg_Bfit) <- NULL
    cg_Bfit <- column_to_rownames(cg_Bfit,var="rowname")
  }
  cg_Bfit <- cg_Bfit %>% na.omit()
  if (grepl(pattern = 'G',model)){
    colnames(cg_Bfit)[2] <- "Beta"
  } else {
    colnames(cg_Bfit)[1] <- "Beta"
  }
  head(cg_Bfit)
  
  cg_df <- 
    left_join(rownames_to_column(as.data.frame(cg_cov[,"doc.theta",drop=F])),
              rownames_to_column(as.data.frame(t(cg_Mfit))), by="rowname")
  colnames(cg_df)[3] <- "Mhat"
  cg_df <- left_join(cg_df,rownames_to_column(cg_Bfit),by="rowname") 
  if (grepl(pattern = 'G',model)){
    colnames(cg_df)[5] <- "betahat"
    cg_df[,snp] <- as.factor(cg_df[,snp])
  } else {
    colnames(cg_df)[4] <- "betahat"
  }
  cg_df <- cg_df[complete.cases(cg_df),]
  if (grepl(pattern = 'G',model)){
    colnames(cg_df)[which(colnames(cg_df)==snp)] = snp.orig
  }
  
  return(cg_df)
  
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot G1 and G2 SNP locations as karyograms
#'
#' @return NULL
chromosome_graph_plot = function(plot_file.gr=NULL,file_name_1=NULL,colorset_code=NULL){
 
  require(magick)
  require(RIdeogram)
  data(human_karyotype, package="RIdeogram")
  data(gene_density, package="RIdeogram")
  human_karyotype <- human_karyotype[1:22,]
  gene_density$Value = "0"

  chr_plot.gr <- plot_file.gr 
  chr_plot_replicate_true <- as.data.frame(chr_plot.gr[which(chr_plot.gr$replicated==TRUE,)])
  chr_1<-data.frame(chr_plot_replicate_true$seqnames,chr_plot_replicate_true$start-1600000,chr_plot_replicate_true$start-1400000)
  chr_1$chr_plot_replicate_true.seqnames<-gsub("chr","",chr_1$chr_plot_replicate_true.seqnames)
  chr_1$Value = "200"
  chr_1 <- setNames(chr_1,  c("Chr", "Start", "End","Value"))

  chr_loc_marked_for_plot <-rbind(gene_density,chr_1,stringsAsFactors=TRUE)
  chr_loc_marked_for_plot<-chr_loc_marked_for_plot[which(chr_loc_marked_for_plot$Chr != "X"), ]
  chr_loc_marked_for_plot<-chr_loc_marked_for_plot[which(chr_loc_marked_for_plot$Chr != "Y"), ]
  chr_loc_marked_for_plot<-chr_loc_marked_for_plot[order(chr_loc_marked_for_plot$Chr,chr_loc_marked_for_plot$Start), ]
  rownames(chr_loc_marked_for_plot) <- 1 : length(rownames(chr_loc_marked_for_plot))
  chr_loc_marked_for_plot$Chr<-as.numeric(chr_loc_marked_for_plot$Chr)
  chr_loc_marked_for_plot$Value<-as.numeric(chr_loc_marked_for_plot$Value)

  ideogram(karyotype = human_karyotype,label = NULL, label_type = NULL, overlaid = chr_loc_marked_for_plot, colorset1 = colorset_code,width = 170,Lx=160,Ly=35)
  convertSVG("chromosome.svg", file = '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/chromosome_raw_graph.png', device = "png")

  chr_image_test <- image_read('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/chromosome_raw_graph.png'))
  chr_image_test_1 <- image_trim(chr_image_test,fuzz=0)
  chr_image_test_2 <- image_border(image_annotate(chr_image_test_1, "CONFIDENTIAL", size = 100, color = "white", boxcolor = "white",gravity="northeast",location="+40+90"),color = "white","20X20") 
  image_write(chr_image_test_2, path = paste0('~/projects/soc_enid_emphasis/analysis/emphasis_GEM/',file_name_1,'.png'), format = "png")
  file.remove(paste0(~/projects/soc_enid_emphasis/analysis/emphasis_GEM/'chromosome.svg'))
  file.remove(paste0(~/projects/soc_enid_emphasis/analysis/emphasis_GEM/'chromosome_raw_graph.png'))
  rm(chr_plot.gr,chr_plot_replicate_true,chr_1,chr_loc_marked_for_plot,chr_image_test,chr_image_test_1,chr_image_test_2)
  #' please check if file labelled as G1 or G2 is only created in ~/projects/soc_enid_emphasis/analysis/emphasis_GEM/
}

plot_G1_G2_locations_karyogram = function(){
  
  require(GenomicRanges)

  g1.g2.table = read_csv(
    '~/projects/soc_enid_emphasis/analysis/emphasis_GEM/results_tables/g1.g2.soc.cpg.results.csv'
  )
  
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
  
  chromosome_graph_plot(plot_file.gr=g1.gr,file_name_1=G1.loc,colorset_code= c( "#ededed", "#e5432d"))
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
  
   chromosome_graph_plot(plot_file.gr=g2.gr,file_name_1=G2.loc,colorset_code= c( "#ededed", "#3CB371")) 
}









