#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Functions for plotting characteristics of loci identified in
# the SoC discovery-relication analysis.
# These are called from SoCFourier_main_analysis.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



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

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Summarise numbers of features overlapping discovery/replicated CpGs
#'
#' @param annot.gr 
#'
#' @return NULL
summarise_overlapping_features = function(annot.gr=NULL){
  
  cat('Features overlapping',sum(annot.gr$discovery),'discovery cgs:\n')
  cat('ME:',sum(annot.gr$discovery & annot.gr$ME),'\n')
  cat('ME100bp:',sum(annot.gr$discovery & annot.gr$ME.100bp),'\n')
  cat('oo.gDMR:',sum(annot.gr$discovery & annot.gr$oo.gDMR),'\n')
  cat('oo.blast.gDMR:',sum(annot.gr$discovery & annot.gr$oo.blast.gDMR),'\n')
  cat('oo.blast.plac.gDMR:',sum(annot.gr$discovery & annot.gr$oo.blast.plac.gDMR),'\n')
  cat('sp.gDMR:',sum(annot.gr$discovery & annot.gr$sp.gDMR),'\n')
  cat('sp.blast.gDMR:',sum(annot.gr$discovery & annot.gr$sp.blast.gDMR),'\n')
  cat('\n\n')
  cat('Features overlapping',sum(annot.gr$replicated),'replicated cgs:\n')
  cat('ME:',sum(annot.gr$replicated & annot.gr$ME),'\n')
  cat('ME100bp:',sum(annot.gr$replicated & annot.gr$ME.100bp),'\n')
  cat('oo.gDMR:',sum(annot.gr$replicated & annot.gr$oo.gDMR),'\n')
  cat('oo.blast.gDMR:',sum(annot.gr$replicated & annot.gr$oo.blast.gDMR),'\n')
  cat('oo.blast.plac.gDMR:',sum(annot.gr$replicated & annot.gr$oo.blast.plac.gDMR),'\n')
  cat('sp.gDMR:',sum(annot.gr$replicated & annot.gr$sp.gDMR),'\n')
  cat('sp.blast.gDMR:',sum(annot.gr$replicated & annot.gr$sp.blast.gDMR),'\n')
  cat('\n\n')

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Combined line plot of Fourier regression curve for each CpG
#'
#' MEs and non-MEs are distinguished by line colour 
#'
#' @param annot.gr 
#' @param fitted.vals 
#' @param plot.legend 
#'
#' @return ggplot object
plot_SoC_Fourier_curves_discovery = function(
  annot.gr=NULL,fitted.vals=NULL,plot.legend=FALSE,plot.save=FALSE,plot.filename=NULL,plot.params=NULL){
  
  require(tidyverse)

  fitted.vals = fitted.vals[annot.gr[annot.gr$discovery]$cpg,]
  fitted.vals$locus = 'other'
  fitted.vals[row.names(fitted.vals) %in% annot.gr[annot.gr$ME]$cpg,'locus'] = 'ME'
  fitted.vals$cpg = row.names(fitted.vals)
  fitted.vals.m = gather(fitted.vals,key = 'doc.theta',value = 'beta.change',-cpg,-locus,convert = T)
  fitted.vals.m$doc.theta = as.numeric(fitted.vals.m$doc.theta)
  fitted.vals.m$beta.change = as.numeric(fitted.vals.m$beta.change)

  plot.alpha = 0.5
  plt = ggplot(data = fitted.vals.m, aes(x=doc.theta,y=100*beta.change,fill=cpg)) +
    geom_line(data = filter(fitted.vals.m,locus=='other'),aes(colour='other'),alpha=plot.alpha,size=0.3) +
    geom_line(data = filter(fitted.vals.m,locus=='ME'),aes(colour='ME'),alpha=plot.alpha,size=0.3) +
    scale_colour_manual("locus",
      breaks=c('ME','other'),
      values=c('red','royalblue2')) +
    scale_x_continuous(name = 'month of conception',
      breaks = seq(0+pi/12,2*pi-pi/12,length.out = 12),
      labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
    coord_cartesian(ylim = c(-8,8)) +
    theme_classic(base_size = 12) +
    labs(y='seasonal change (meth %)')

  if (!plot.legend) plt = plt + theme(legend.position="none")
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,plot.filename),
      dpi = plot.params$dpi,width = plot.params$width,height = plot.params$height,units = 'cm'
    )
  }

  return(plt)

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Scatter plot of DoC at methylation maximum across cohorts
#'
#' @param plot.save 
#' @param plot.params 
#'
#' @return NULL
plot_cross_cohort_doc_max = function(plot.save=FALSE,plot.params=NULL){
  
  require(gridExtra)
  require(cowplot)
  
  plot.max = as.data.frame(mcols(annot.gr)) %>%
    filter(discovery | ctrl) %>% 
    select(cpg,discovery,replicated,ctrl,enid.max.doc.theta,emph.max.doc.theta) %>%
    gather(enid.max.doc.theta,emph.max.doc.theta,
           key = 'max.min',
           value = 'doc.theta'
    )
  plot.max[grepl('enid',plot.max$max.min),'cohort'] = 'ENID'
  plot.max[grepl('emph',plot.max$max.min),'cohort'] = 'EMPHASIS'
  
  scatter_doc_max_across_cohorts = function(df=NULL,plt.title=NULL){
    
    
    plt = ggplot(data=comp.effects.df,aes(x=ENID,y=EMPHASIS)) +
      geom_point(size=0.4) + 
      scale_x_continuous(name ='methylation maximum\nDISCOVERY',
                         breaks = seq(0+pi/12,2*pi-pi/12,length.out = 12), 
                         labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
      scale_y_continuous(name = 'methylation maximum\nREPLICATION',
                         breaks = seq(0+pi/12,2*pi-pi/12,length.out = 12), 
                         labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
      theme_classic(base_size = 11) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
      coord_fixed(ratio = 1,xlim = c(0,2*pi),ylim = c(0,2*pi)) +
      annotate('label',x = 6.5,y = 0,label=plt.title,hjust=1,size=2)
    
    return(plt)
    
  }
  
  # SoC-CpGs
  comp.effects.df = plot.max %>%
    filter(replicated) %>%
    select(cpg,doc.theta,cohort) %>%
    spread(cohort,doc.theta)
  # calculate Spearman correlation of DoC max across cohorts
  corr = cor.test(comp.effects.df$ENID,comp.effects.df$EMPHASIS,method = 'spearman')
  plt.title = paste0('SoC-CpGs: rho=',round(corr$estimate,2),'; p=',formatC(corr$p.value,format='e',digits=2))
  plt1 = scatter_doc_max_across_cohorts(df = comp.effects.df,plt.title = plt.title)
  
  # discovery CpGs
  comp.effects.df = plot.max %>%
    filter(discovery) %>%
    select(cpg,doc.theta,cohort) %>%
    spread(cohort,doc.theta)
  # calculate Spearman correlation of DoC max across cohorts
  corr = cor.test(comp.effects.df$ENID,comp.effects.df$EMPHASIS,method = 'spearman')
  plt.title = paste0('Discovery CpGs: rho=',round(corr$estimate,2),'; p=',formatC(corr$p.value,format='e',digits=2))
  plt2 = scatter_doc_max_across_cohorts(df = comp.effects.df,plt.title = plt.title)
  
  # matched controls
  comp.effects.df = plot.max %>%
    filter(ctrl) %>%
    select(cpg,doc.theta,cohort) %>%
    spread(cohort,doc.theta)
  # calculate Spearman correlation of DoC max across cohorts
  corr = cor.test(comp.effects.df$ENID,comp.effects.df$EMPHASIS,method = 'spearman')
  plt.title = paste0('Matched controls: rho=',round(corr$estimate,2),'; p=',round(corr$p.value,2))
  plt3 = scatter_doc_max_across_cohorts(df = comp.effects.df,plt.title = plt.title)
  
  # plot soc-cpgs and matched controls
  plt.soc.vs.matched = 
    grid.arrange(
      plt1,plt3, nrow = 1,
      bottom='methylation maximum\nDISCOVERY',
      left='methylation maximum\nREPLICATION'
    )
  ggdraw(plt.soc.vs.matched)
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,'enid_vs_emph_max_doc_soc_vs_matched.pdf'),
      dpi = plot.params$dpi,width = 1.5*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
  
  # plot discovery cpgs
  print(plt2)
  if(plot.save){
    ggsave(plt2,
           filename = paste0(plot.params$dir,'enid_vs_emph_max_doc_discovery.pdf'),
           dpi = plot.params$dpi,width = 1.5*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot date of conception at methylation maxima/minima for various loci
#'
#' @param annot.gr 
#' @param plot.max.doc.only logical plot doc at meth max only?
#' @param plot.save save plot to file?
#' @param plot.params plot size params
#'
#' @return NULL
plot_doc_maxmin = function(annot.gr=NULL,plot.max.doc.only=TRUE,plot.save=FALSE,plot.params=NULL){
  
  plot.maxmin = as.data.frame(mcols(annot.gr)) %>%
    filter(replicated | ctrl | rand | (ME & discovery)) %>% # replicated, controls and MEs in discovery set
    select(cpg,ME,discovery,replicated,rand,ctrl,enid.max.doc.theta,enid.min.doc.theta,emph.max.doc.theta,emph.min.doc.theta) %>%
    gather(enid.max.doc.theta,enid.min.doc.theta,emph.max.doc.theta,emph.min.doc.theta,
        key = 'max.min',
        value = 'doc.theta'
      )
  plot.maxmin[grepl('enid',plot.maxmin$max.min),'cohort'] = 'ENID'
  plot.maxmin[grepl('emph',plot.maxmin$max.min),'cohort'] = 'EMPHASIS'
  plot.maxmin[plot.maxmin$discovery & plot.maxmin$ME,'locus'] = 'ME\nnot replicated'
  plot.maxmin[plot.maxmin$replicated & plot.maxmin$ME,'locus'] = 'SoC-CpGs\nMEs'
  plot.maxmin[plot.maxmin$replicated & !plot.maxmin$ME,'locus'] = 'SoC-CpGs\nnon-MEs'
  plot.maxmin[plot.maxmin$ctrl,'locus'] = 'matched\ncontrols'
  plot.maxmin[plot.maxmin$rand,'locus'] = 'random\ncontrols'
  plot.maxmin$max.min[grepl(pattern = 'min',x = plot.maxmin$max.min)] = 'methylation miniumum'
  plot.maxmin$max.min[grepl(pattern = 'max',x = plot.maxmin$max.min)] = 'methylation maximum'
  
  if (plot.max.doc.only) plot.maxmin = filter(plot.maxmin,max.min=='methylation maximum') # plot maxima only
  # table(plot.maxmin$locus)/2
  
  # specify plot order
  plot.maxmin$locus = 
    recode_factor(plot.maxmin$locus, 
        'SoC-CpGs\nMEs' = 'SoC-CpGs\nMEs',
        'SoC-CpGs\nnon-MEs' = 'SoC-CpGs\nnon-MEs',
        'ME\nnot replicated' = 'ME\nnot replicated',
        'matched\ncontrols' = 'matched\ncontrols',
        'random\ncontrols' = 'random\ncontrols'
  )
  plot.maxmin$cohort = recode_factor(plot.maxmin$cohort,
        'ENID' = 'ENID',
        'EMPHASIS' = 'EMPHASIS'
  )

  if (plot.max.doc.only) {

    plt = ggplot(
      data = filter(plot.maxmin, locus != 'other SoC\nnot replicated'),
      aes(x = locus,y = doc.theta,fill = locus,colour = locus)
    ) +
      geom_jitter(position = position_jitterdodge(dodge.width = 0.9),alpha = 1,size = 0.5) +
      geom_violin(colour = 'black', alpha = 0) +
      facet_wrap(. ~ cohort, strip.position = 'top') +
      scale_y_continuous(
        name = 'month of conception',
        breaks = seq(0 + pi / 12, 2 * pi - pi / 12, length.out = 12),
        labels = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')
      ) +
      scale_colour_manual("locus",
        values = c('red','royalblue2','pink','lightblue','lightgrey','darkgrey')) +
      theme_classic(base_size = 12) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(
          colour = 'black',
          size = 10,
          angle = 45,
          hjust = 1
        ),
        axis.title.x = element_blank(),
        panel.grid.major.x =  element_blank(),
        panel.grid.major.y =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        legend.position = "none"
      )
    
    print(plt)
    
    if(plot.save){
      ggsave(
        filename = paste0(plot.params$dir,'doc_max_ENID_EMPHASIS.pdf'),
        dpi = plot.params$dpi,width = 2*plot.params$width,height = 1.25*plot.params$height,units = 'cm'
      )
      
    }

  } else {

    plt = ggplot(data = filter(plot.maxmin,locus != 'other SoC\nnot replicated' & cohort == 'ENID'),
      aes(x = locus,y = doc.theta,fill = locus,colour = locus)) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9),alpha = 1,size = 0.5) +
    geom_violin(colour = 'black', alpha = 0) +
    facet_wrap(. ~ max.min, strip.position = 'top') +
    scale_y_continuous(
      name = 'month of conception',
      breaks = seq(0 + pi / 12, 2 * pi - pi / 12, length.out = 12),
      labels = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')) +
    scale_colour_manual("locus",
      values = c('red','royalblue2','pink','lightblue','lightgrey','darkgrey')) +
    labs(title = 'DISCOVERY (ENID)') +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(
        colour = 'black',
        size = 10,
        angle = 45,
        hjust = 1
      ),
      axis.title.x = element_blank(),
      panel.grid.major.x =  element_blank(),
      panel.grid.major.y =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      legend.position = "none"
    )
    print(plt)

    if(plot.save){
      ggsave(
        filename = paste0(plot.params$dir,'doc_maxmin_ENID.pdf'),
        dpi = plot.params$dpi,width = 2*plot.params$width,height = 1.25*plot.params$height,units = 'cm'
      )
      
    }

    plt = ggplot(
      data = filter(plot.maxmin,locus != 'other SoC\nnot replicated' & cohort == 'EMPHASIS'),
      aes(x = locus,y = doc.theta,fill = locus,colour = locus)
    ) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.9),alpha = 1,size = 0.5) +
    geom_violin(colour = 'black', alpha = 0) +
    facet_wrap(. ~ max.min, strip.position = 'top') +
    scale_y_continuous(
      name = 'month of conception',
      breaks = seq(0 + pi / 12, 2 * pi - pi / 12, length.out = 12),
      labels = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D')
    ) +
    scale_colour_manual("locus",
      values = c('red','royalblue2','pink','lightblue','lightgrey','darkgrey')) +
    labs(title = 'REPLICATION (EMPHASIS)') +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.text.x = element_text(
        colour = 'black',
        size = 10,
        angle = 45,
        hjust = 1
      ),
      axis.title.x = element_blank(),
      panel.grid.major.x =  element_blank(),
      panel.grid.major.y =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      legend.position = "none"
    )
    print(plt)
    
    if(plot.save){
      ggsave(
        filename = paste0(plot.params$dir,'doc_maxmin_EMPHASIS.pdf'),
        dpi = plot.params$dpi,width = 2*plot.params$width,height = 1.25*plot.params$height,units = 'cm'
      )
      
    }
  }
  

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot amplitude (distance between meth max and min) for both cohorts
#'
#' @param annot.gr 
#' @param plot.save save plot to file?
#' @param plot.params plot dimensions
#'
#' @return NULL
plot_amplitude = function(annot.gr=NULL,plot.save=FALSE,plot.params=NULL){

  plot.amplitude = as.data.frame(mcols(annot.gr)) %>%
    filter(replicated | ctrl | rand | (ME & discovery)) %>% # replicated, controls and MEs in discovery set
    select(cpg,ME,discovery,replicated,rand,ctrl,enid.amplitude,emph.amplitude,replicated) %>%
    gather(enid.amplitude,emph.amplitude,
           key = 'cohort',
           value = 'amplitude'
    )
  plot.amplitude$cpg = as.character(plot.amplitude$cpg)
  plot.amplitude[grepl('enid',plot.amplitude$cohort),'cohort'] = 'ENID'
  plot.amplitude[grepl('emph',plot.amplitude$cohort),'cohort'] = 'EMPHASIS'
  plot.amplitude[plot.amplitude$discovery & plot.amplitude$ME,'locus'] = 'ME\nnot replicated'
  plot.amplitude[plot.amplitude$replicated & plot.amplitude$ME,'locus'] = 'SoC-CpGs\nMEs'
  plot.amplitude[plot.amplitude$replicated & !plot.amplitude$ME,'locus'] = 'SoC-CpGs\nnon-MEs'
  plot.amplitude[plot.amplitude$ctrl,'locus'] = 'matched\ncontrols'
  plot.amplitude[plot.amplitude$rand,'locus'] = 'random\ncontrols'
  #  table(plot.amplitude$locus)/2
  
  # specify plot order
  plot.amplitude$locus = 
    recode_factor(plot.amplitude$locus, 
        'ME\nreplicated' = 'SoC-CpGs\nMEs',
        'SoC\nreplicated' = 'SoC-CpGs\nnon-MEs',
        'ME-SoC' = 'ME\nnot replicated',
        'matched\ncontrols' = 'matched\ncontrols',
        'random\ncontrols' = 'random\ncontrols'
    )
  plot.amplitude$cohort = recode_factor(plot.amplitude$cohort,
        'ENID' = 'ENID',
        'EMPHASIS' = 'EMPHASIS'
  )
  
  plt = ggplot(data = plot.amplitude, aes(x=locus, y=amplitude*100,colour=locus)) +
    geom_jitter(width=0.1,alpha=1,size=0.5) +
    geom_boxplot(width=0.2,colour='black',alpha=0, outlier.shape = NA) +
    facet_wrap(.~cohort,strip.position = 'top') +
    labs(
      x='',
      y='seasonal amplitude (%)'
    ) +
    scale_colour_manual("locus",
      values=c('red','royalblue2','pink','lightblue','lightgrey','darkgrey')) +
    coord_cartesian(ylim = c(0,17.5)) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside",
      axis.text.x = element_text(colour = 'black',size = 10,angle = 45,hjust = 1),
      legend.position="none"
    )
  
  print(plt)

  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,'amplitude.pdf'),
      dpi = plot.params$dpi,width = 2*plot.params$width,height = 1.35*plot.params$height,units = 'cm'
    )
  }

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot DoC at meth max vs amplitude
#' 
#' This is done for discovery (ENID) and replication (EMPHASIS) cohorts
#' replicated loci distinguished by colour
#'
#' @param annot.gr 
#' @param plot.save logical save plot to file?
#' @param plot.params list of params specifying plot dir and dimensions etc
#'
#' @return
plot_doc_vs_amplitude = function(annot.gr=NULL,plot.save=FALSE,plot.params=NULL){

  plot.doc.max.amp.enid = as.data.frame(mcols(annot.gr)) %>%
    filter(discovery) %>% # discovery and replicated only
    select(cpg,replicated,doc.max=enid.max.doc.theta,amplitude=enid.amplitude)
  plot.doc.max.amp.enid$cohort = 'ENID'
  plot.doc.max.amp.emph = as.data.frame(mcols(annot.gr)) %>%
    filter(discovery) %>% # discovery and replicated only
    select(cpg,replicated,doc.max=emph.max.doc.theta,amplitude=emph.amplitude)
  plot.doc.max.amp.emph$cohort = 'EMPHASIS'
  plot.doc.max.amp.all = rbind(plot.doc.max.amp.enid,plot.doc.max.amp.emph)  
  # specify plot order
  plot.doc.max.amp.all$cohort = recode_factor(plot.doc.max.amp.all$cohort,
    'ENID' = 'ENID',
    'EMPHASIS' = 'EMPHASIS'
  )
  
  plt = ggplot(plot.doc.max.amp.all,aes(x=100*amplitude,y=doc.max,colour=replicated)) +
    geom_point(data = filter(plot.doc.max.amp.all,replicated==FALSE),size=1) + 
    geom_point(data = filter(plot.doc.max.amp.all,replicated==TRUE),size=1) + 
    facet_wrap(~cohort) +
    theme_classic() +
    scale_y_continuous(name = 'month of conception', 
                       breaks = seq(0+pi/12,2*pi-pi/12,length.out = 12), 
                       labels = c('J','F','M','A','M','J','J','A','S','O','N','D')) +
    scale_colour_manual("replicated",
      breaks=c('FALSE','TRUE'),
      values=c('lightblue','royalblue2'),
      limits=c('FALSE','TRUE')
    ) +
    labs(
      x='seasonal amplitude',
      y='month of conception at\nmethylation maximum'
    ) +
    theme(
      strip.background = element_blank(), 
      strip.placement = "outside"
    )
  print(plt)
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,'doc.max.vs.amplitude.pdf'),
      dpi = plot.params$dpi,width = 2*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot mean or median methylation Beta distribution for SoC-CpGs vs 
#' array background
#' 
#' Beta values are from both cohorts combined. SoC-CpGs are 
#' stratified by ME status
#'
#' @param annot.gr 
#' @param summary.stat string which summary stat to plot: 'mean'/'median'/'var'
#' @param plot.save logical save plot to file?
#' @param plot.params list of params specifying plot dir and dimensions etc
#'
#' @return
plot_beta_summary_stat_vs_bgd = function(annot.gr=NULL,summary.stat=NULL,plot.save=FALSE,plot.params=NULL){
 
  plot.beta.stat = as.data.frame(mcols(annot.gr)) %>%
    select(replicated,ME,
           !!paste0('enid.',summary.stat,'.beta'),
           !!paste0('emph.',summary.stat,'.beta')
    ) %>%
    gather(!!paste0('enid.',summary.stat,'.beta'),!!paste0('emph.',summary.stat,'.beta'),
           key = 'cohort',
           value = 'stat'
    )
  
  plot.beta.stat$locus = 'array\nbackground'
  plot.beta.stat[plot.beta.stat$replicated,'locus'] = 'SoC-CpGs\nnon-MEs'
  plot.beta.stat[plot.beta.stat$replicated & plot.beta.stat$ME ,'locus'] = 'SoC-CpGs\nMEs'

  # reduce number of array background points to aid plot rendering as pdf (2*array size for both cohorts)
  plot.beta.stat.bgd = filter(plot.beta.stat,locus=='array\nbackground')
  plot.beta.stat.bgd.sample = 
    plot.beta.stat.bgd[sample(seq(nrow(plot.beta.stat.bgd)),size = nrow(plot.beta.stat.bgd)/16,replace = F),]
  plot.beta.stat = rbind(
    filter(plot.beta.stat,locus!='array\nbackground'),
    plot.beta.stat.bgd.sample
  )

  # specify plot order
  plot.beta.stat$locus = recode_factor(plot.beta.stat$locus, 
      'SoC-CpGs\nMEs' = 'SoC-CpGs\nMEs',
      'SoC-CpGs\nnon-MEs' = 'SoC-CpGs\nnon-MEs',
      'array\nbackground' = 'array\nbackground'
  )

  plt = ggplot(data=plot.beta.stat, aes(x=locus,y=100*stat,colour=locus,alpha=locus,size=locus)) +
    geom_jitter(width=0.1) +
    geom_violin(colour='black',alpha=0,size=0.5) +
    theme_bw() +
    scale_colour_manual(values=c('red','royalblue2','black')) +
    scale_alpha_manual(values=c(0.5,0.5,0.01)) +
    scale_size_manual(values=c(0.5,0.5,0.01)) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle=30,hjust = 1,colour='black',size='9'),
      legend.position="none",
      strip.background = element_rect(fill = 'white')
    ) +
    labs(x='',y=paste0(summary.stat,'\nmethylation (%)'))
  print(plt)

  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,summary.stat,'.betas_vs_bgd.pdf'),
      dpi = plot.params$dpi,width = plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
    
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot mean or median methylation Beta distribution for SoC-CpGs vs 
#' matched controls
#' 
#' Beta values are from both cohorts combined. 
#'
#' @param annot.gr 
#' @param summary.stat string which summary stat to plot: 'mean'/'median'/'var'
#' @param plot.save logical save plot to file?
#' @param plot.params list of params specifying plot dir and dimensions etc
#'
#' @return
plot_beta_summary_stat_vs_matched = function(annot.gr=NULL,summary.stat=NULL,plot.save=FALSE,plot.params=NULL){
  
  plot.beta.stat = as.data.frame(mcols(annot.gr)) %>%
    filter(replicated | ctrl) %>%
    select(replicated,ME,ctrl,
           !!paste0('enid.',summary.stat,'.beta'),
           !!paste0('emph.',summary.stat,'.beta')
    ) %>%
    gather(!!paste0('enid.',summary.stat,'.beta'),!!paste0('emph.',summary.stat,'.beta'),
           key = 'cohort',
           value = 'stat'
    )
  
  plot.beta.stat[grepl('enid',plot.beta.stat$cohort),'cohort'] = 'ENID'
  plot.beta.stat[grepl('emph',plot.beta.stat$cohort),'cohort'] = 'EMPHASIS'
  plot.beta.stat[plot.beta.stat$replicated,'locus'] = 'SoC-CpGs'
  plot.beta.stat[plot.beta.stat$ctrl,'locus'] = 'matched\ncontrols'
  
  # specify plot order
  plot.beta.stat$locus = recode_factor(plot.beta.stat$locus, 
   'SoC-CpGs' = 'SoC-CpGs',
   'matched\ncontrols' = 'matched\ncontrols'
  )
  
  plt = ggplot(data=plot.beta.stat, aes(x=locus,y=100*stat,colour=locus)) +
    geom_jitter(width=0.1,alpha=0.5, size=0.5) +
    geom_violin(colour='black',alpha=0) +
    theme_bw() +
    scale_colour_manual(values=c('red','royalblue2','lightgrey','darkgrey')) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle=30,hjust = 1,colour='black',size='9'),
      legend.position="none",
      strip.background = element_rect(fill = 'white')
    ) +
    labs(x='',y=paste0(summary.stat,'\nmethylation (%)'))
  
  print(plt)
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,summary.stat,'.betas_vs_matched.pdf'),
      dpi = plot.params$dpi,width = 2*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot Beta distributions of KS-matched controls
#'
#' Each plot compares Beta distribution of query CpG with 
#' all matched target CpGs. A total of n.query.cgs.to.plot
#' query CpGs and corresonding matches are plotted
#'
#' @param matches List of KS matches mapping query to targets
#' @param n.query.cgs.to.plot Number of query CpGs to plot
#' @param plot.save logical
#' @param betas DF methylation betas
#'
#' @return NULL
plot_ks_matches = function(matches=NULL,n.query.cgs.to.plot=1,betas=NULL,plot.save=FALSE,plot.params=NULL){
  
  require(tidyverse)
  
  # plot some matches
  for (query.cg in names(matches)[1:n.query.cgs.to.plot]){
    
    betas.plot = as.data.frame(betas[c(query.cg,matches[[query.cg]]),])
    betas.plot$cg = row.names(betas.plot)
    betas.plot = betas.plot %>% 
      gather(key='id',value = 'beta',-cg) %>%
      mutate(locus = ifelse(cg==query.cg, 'query','ks.match'))
    betas.plot$locus = recode_factor(betas.plot$locus,
       'query' = 'query',
       'ks.match' = 'ks.match'
    )
    
    print(
      ggplot() + 
        geom_density(data = betas.plot, aes(x=beta,group=cg,colour=locus,size=locus),show.legend = FALSE) + 
        stat_density(data = betas.plot, aes(x=beta,group=cg,colour=locus,size=locus),geom="line",position="identity") + 
        scale_color_manual(values = c('red','black')) +
        scale_size_manual(values = c(1,0.2)) +
        labs(title = query.cg) +
        coord_cartesian(xlim = c(0,1)) +
        theme_bw()
    )
    
    if(plot.save){
      ggsave(
        filename = paste0(plot.params$dir,'ks_match_beta_distributions/KS-matches_',query.cg,'.pdf'),
        dpi = plot.params$dpi,width = 1.25*plot.params$width,height = 1.1*plot.params$height,units = 'cm'
      )
    }
    
  }
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot SoC-CpG locations in karyogramm
#'
#' @param annot.gr 
#' @param plot.save 
#' @param plot.params 
#'
#' @return NULL
plot_cpg_locations_karyo = function(annot.gr=NULL){

  require(magick)
  require(GenomicRanges)
  require(RIdeogram)
  data(human_karyotype, package="RIdeogram")
  data(gene_density, package="RIdeogram")
  human_karyotype <- human_karyotype[1:22,]
  gene_density$Value = "0"
  
  chr_plot.gr <- annot.gr  #'please check weather the file is correct ( in mail the file name was annotated.array.gr.rds)
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
  
  ideogram(karyotype = human_karyotype,label = NULL, label_type = NULL, overlaid = chr_loc_marked_for_plot, colorset1 = c( "#ededed", "#e5432d"),width = 170,Lx=160,Ly=35)
  convertSVG("chromosome.svg", file = paste0(plot.dir,'chromosome_raw_graph'), device = "png")
    
  chr_image_test <- image_read(paste0(plot.dir,'chromosome_raw_graph.png'))
  chr_image_test_1 <- image_trim(chr_image_test,fuzz=0)
  chr_image_test_2 <- image_border(image_annotate(chr_image_test_1, "CONFIDENTIAL", size = 100, color = "white", boxcolor = "white",gravity="northeast",location="+40+90"),color = "white","20X20") 
  image_write(chr_image_test_2, path = paste0(plot.dir,'replicated.cg.karyogram.png'), format = "png")
  file.remove(paste0(plot.dir,'chromosome_raw_graph.svg'))
  file.remove(paste0(plot.dir,'chromosome_raw_graph.png'))
  rm(chr_plot.gr,chr_plot_replicate_true,chr_1,chr_loc_marked_for_plot,chr_image_test,chr_image_test_1,chr_image_test_2)
  #' please check if file labelled as replicated.cg.karyogram.png is only created in plot.dr
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Plot genomic context wrt CG Islands and gene model
#' 
#' SoC-CpGs compared to controls and array background.
#' Plots include bootstrap CIs
#'
#' @param annot.gr 
#' @param plot.save 
#' @param plot.params 
#' @param verbose logical verbose output?
#'
#' @return NULL
plot_genomic_context = function(annot.gr=NULL,plot.save=FALSE,plot.params=NULL,verbose=F){
  
  require(tidyverse)
  require(GenomicRanges)

  # PARAMS
  n.bootstrap = 200
  cluster = FALSE # collapse clusters?
  cluster.max.gap = 1000 # min inter-cg distance for clustering
  
  # annotations from 450k manifest
  annot.df = as.data.frame(annot.gr) %>%
    select(replicated,ctrl,gene_annot=UCSC_RefGene_Group,CGI_annot=Relation_to_UCSC_CpG_Island,enhancer=Enhancer)
  
  # define inter CGI cgs and non-enhancers
  annot.df[is.na(annot.df$CGI_annot),'CGI_annot'] = "NA"
  annot.df[annot.df$CGI_annot=="",'CGI_annot'] = 'inter CGI'
  annot.df[is.na(annot.df$enhancer),'enhancer'] = FALSE
  
  # gene annotations
  annot.df$body = grepl('Body',annot.df$gene_annot)
  annot.df$TSS200 = grepl('TSS200',annot.df$gene_annot)
  annot.df$TSS1500 = grepl('TSS1500',annot.df$gene_annot)
  annot.df$UTR5 = grepl('5\'UTR',annot.df$gene_annot)
  annot.df$UTR3 = grepl('3\'UTR',annot.df$gene_annot)
  annot.df$first.exon = grepl('1stExon',annot.df$gene_annot)
  annot.df$gene_annot = NULL
  
  annot.df$array.background = T
  
  # CGI ANNOTATIONS
  source('~/src/utils/misc/bootstrap.proportion.CI.R')

  # convert to proportions with bootstrap CI
  plot.cgi = data.frame(
    locus = character(),
    cgi.type = character(),
    freq = numeric(),
    lower.ci = numeric(),
    upper.ci = numeric()
  )
  
  if (verbose) cat('CGI annotations\n')
  
  for (test.set in c('replicated','ctrl','array.background')){
    
    if (verbose) cat('bootstrapping',test.set,'...\n')
    cgi.test.set = annot.df[annot.df[[test.set]]==T,'CGI_annot']
    
    for (cgi.type in unique(cgi.test.set)){
      plot.cgi = rbind(plot.cgi,
         data.frame(
           locus=test.set,
           cgi.type=cgi.type,
           freq=sum(cgi.test.set==cgi.type)/length(cgi.test.set),
           lower.ci=alpha.quantiles.lower(bootstrap.prop(cgi.test.set==cgi.type,n.bootstrap),0.05),
           upper.ci=alpha.quantiles.upper(bootstrap.prop(cgi.test.set==cgi.type,n.bootstrap),0.05)
         )
      )
    }
  }
  cat('done\n')
  
  cat('CGI stats:\n')
  print(plot.cgi)
  
  plot.cgi$locus = dplyr::recode_factor(plot.cgi$locus,
      'replicated' = 'SoC-CpGs',
      'ctrl' = 'matched\ncontrols',
      'array.background' = 'array\nbackground'
  )
  plot.cgi$cgi.type = recode_factor(plot.cgi$cgi.type,
      'Island' = 'Island',
      'N_Shelf' = 'N Shelf',
      'N_Shore' = 'N Shore',
      'S_Shelf' = 'S Shelf',
      'S_Shore' = 'S Shore',
      'inter CGI' = 'open sea'
  )
  
  print(
    ggplot(filter(plot.cgi,cgi.type!='NA'),aes(x = cgi.type, y=as.numeric(freq), fill = locus)) + 
      geom_bar(position = "dodge",stat = "identity") +
      geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.2,position=position_dodge(.9)) +
      scale_fill_manual("locus",values=c('royalblue2','lightgrey','darkgrey')) +
      theme_classic(base_size = 10) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1,colour='black',size='11'),
        legend.title = element_blank()
      ) +
      labs(y='proportion')
  )
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,'CGI_enrichment.pdf'),
      dpi = plot.params$dpi,width = 1.25*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }
  
  
  # GENE ANNOTATIONS
  plot.gene = data.frame(
    locus = character(),
    annot.type = character(),
    freq = numeric(),
    lower.ci = numeric(),
    upper.ci = numeric()
  )
  
  
  if (verbose) cat('\nGene position annotations\n')
  
  for (test.set in c('replicated','ctrl','array.background')){
    
    if (verbose) cat('bootstrapping',test.set,'...\n')
  
    
    for (gene.type in c('enhancer','TSS1500','TSS200','UTR5','body','UTR3')){
      gene.test.set = as.logical(annot.df[annot.df[[test.set]]==T,gene.type])
      plot.gene = rbind(plot.gene,
         data.frame(
           locus=test.set,
           gene.type=gene.type,
           freq=sum(gene.test.set)/length(gene.test.set),
           lower.ci=alpha.quantiles.lower(bootstrap.prop(gene.test.set,n.bootstrap),0.05),
           upper.ci=alpha.quantiles.upper(bootstrap.prop(gene.test.set,n.bootstrap),0.05)
        )
      )
    }
  }
  
  print(plot.gene)

  plot.gene$locus = dplyr::recode_factor(plot.gene$locus,
     'replicated' = 'SoC-CpGs',
     'ctrl' = 'matched\ncontrols',
     'array.background' = 'array\nbackground'
  )
  
  print(
    ggplot(plot.gene,aes(x = gene.type, y=as.numeric(freq), fill = locus)) + 
      geom_bar(position = "dodge",stat = "identity") +
      geom_errorbar(aes(ymin=lower.ci, ymax=upper.ci), width=.2,position=position_dodge(.9)) +
      scale_fill_manual("locus",values=c('royalblue2','lightgrey','darkgrey')) +
      theme_classic(base_size = 10) +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1,colour='black',size='11'),
        legend.title = element_blank()
      ) +
      labs(y='proportion')
  )
  
  if(plot.save){
    ggsave(
      filename = paste0(plot.params$dir,'gene_position_enrichment.pdf'),
      dpi = plot.params$dpi,width = 1.25*plot.params$width,height = plot.params$height,units = 'cm'
    )
  }

}





#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#' Plot distributions of cpg pairwise correlations in 
#' discovery and replication cohorts
#' 
#' Distributions are plotted for discovery, SoC-CpGs and controls
#' Option to de=cluster GRs
#'
#' @param annot.gr GR
#' @param disc.fitted.vals DF fitted values for discovery cohort
#' @param repl.fitted.vals DF fitted values for replication cohort
#' @param plot.save 
#' @param plot.params 
#' @param cluster.adj logical de-cluster cgs first?
#' @param cluster.size integer bp distance to define cluster
#'
#' @return NULL
plot_pairwise_cpg_correlations = function(
  annot.gr=NULL,disc.fitted.vals=NULL,repl.fitted.vals=NULL,cluster.adj=F,cluster.size=NULL,
  plot.save=FALSE,plot.params=NULL
){
  
  require(tidyverse)
  
  discovery.gr = annot.gr[annot.gr$discovery]
  replicated.gr = annot.gr[annot.gr$replicated]
  ctrl.gr = annot.gr[annot.gr$ctrl]
  rnd.gr = annot.gr[annot.gr$rand]
  
  if (cluster.adj){
    source('~/src/soc_enid_emphasis/SoCFourier_stats_functions.R')
    discovery.gr = de_cluster_gr(test.gr = discovery.gr,clustergap = cluster.size,verbose = F)
    replicated.gr = de_cluster_gr(test.gr = replicated.gr,clustergap = cluster.size,verbose = F)
    rnd.gr = de_cluster_gr(test.gr = rnd.gr,clustergap = cluster.size,verbose = F)
    ctrl.gr = de_cluster_gr(test.gr = ctrl.gr,clustergap = cluster.size,verbose = F)
  }
  
  # DISCOVERY (ENID) SPEARMAN CORRS
  corrs.disc = list()
  cors = cor(t(disc.fitted.vals[discovery.gr$cpg,]),method = 'spearman')
  corrs.disc[['discovery']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(disc.fitted.vals[replicated.gr$cpg,]),method = 'spearman')
  corrs.disc[['replicated']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(disc.fitted.vals[ctrl.gr$cpg,]),method = 'spearman')
  corrs.disc[['ctrl']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(disc.fitted.vals[rnd.gr$cpg,]),method = 'spearman')
  corrs.disc[['rand']] = cors[lower.tri(cors,diag = F)]
  
  corrs.disc = rbind(
    data.frame(locus='discovery',corr=corrs.disc[['discovery']]),
    data.frame(locus='replicated',corr=corrs.disc[['replicated']]),
    data.frame(locus='ctrl',corr=corrs.disc[['ctrl']]),
    data.frame(locus='rand',corr=corrs.disc[['rand']])
  )
  
  corrs.disc$cohort = 'DISCOVERY'
  
  # REPLICATION (EMPHASIS) SPEARMAN CORRS
  corrs.repl = list()
  cors = cor(t(repl.fitted.vals[discovery.gr$cpg,]),method = 'spearman')
  corrs.repl[['discovery']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(repl.fitted.vals[replicated.gr$cpg,]),method = 'spearman')
  corrs.repl[['replicated']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(repl.fitted.vals[ctrl.gr$cpg,]),method = 'spearman')
  corrs.repl[['ctrl']] = cors[lower.tri(cors,diag = F)]
  cors = cor(t(repl.fitted.vals[rnd.gr$cpg,]),method = 'spearman')
  corrs.repl[['rand']] = cors[lower.tri(cors,diag = F)]
  
  corrs.repl = rbind(
    data.frame(locus='discovery',corr=corrs.repl[['discovery']]),
    data.frame(locus='replicated',corr=corrs.repl[['replicated']]),
    data.frame(locus='ctrl',corr=corrs.repl[['ctrl']]),
    data.frame(locus='rand',corr=corrs.repl[['rand']])
  )
  
  corrs.repl$cohort = 'REPLICATION'

  corrs.all = rbind(corrs.disc,corrs.repl)
  
  # PLOT
  corrs.all$cohort = recode_factor(corrs.all$cohort,
       'DISCOVERY' = 'DISCOVERY (ENID)',
       'REPLICATION' = 'REPLICATION (EMPHASIS)'
  )
  corrs.all$locus = recode_factor(corrs.all$locus,
        'discovery' = 'discovery\nCpGs',
        'replicated' = 'SoC\nCpGs',
        'ctrl' = 'matched\ncontrols',
        'rand' = 'random\ncontrols'
  )
  
  print(
    ggplot(data = corrs.all, aes(x=locus,y=corr,fill=locus)) +
      geom_boxplot(outlier.size=0.5) +
      geom_hline(aes(yintercept = 0),linetype='dashed',colour='grey') +
      facet_wrap(.~cohort) +
      scale_fill_manual("locus",values=c('lightblue','royalblue2','lightgrey','darkgrey')) +
      labs(
        x='',
        y='Spearman\npairwise correlation'
      ) +
      theme_bw() +
      theme(
        legend.position = 'none',
        strip.background = element_blank(), 
        strip.placement = "outside"
      )
  )
  
  if(plot.save){
    if (cluster.adj){
      ggsave(
        filename = paste0(plot.params$dir,'pairwise_correlation_boxplot_CLUSTER_ADJUSTED.pdf'),
        dpi = plot.params$dpi,width = 2*plot.params$width,height = plot.params$height,units = 'cm'
      )
    } else {
      ggsave(
        filename = paste0(plot.params$dir,'pairwise_correlation_boxplot.pdf'),
        dpi = plot.params$dpi,width = 2*plot.params$width,height = plot.params$height,units = 'cm'
      )
    }
  }
  
}








