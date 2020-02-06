#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# Helper functions use to perform Fourier regression, likelihood ratio tests,
# generate fitted values, and produce model fit summary stats
# Functions are called from Fourier_lrt_model_fit.R
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Compute LRT pvals from Fourier regression model fit
#'
#' LRT computed by comparing Fourier model with 1 pair of Fourier terms
#' against baseline model with adjustment covariates only
#'
#' @param meth df methylation values (cgs x samples)
#' @param design df model design matrix
#'
#' @return df with cg.id, LRT pval, number of Fourier terms modelled and FDR qval
#'
generate_fourier_lrt_pvals = function(meth=NULL,design=NULL){

  require(lmtest)
  require(lubridate)
  require(doParallel) # recommended for multi-thread execution of array-wide LRTs
  require(tidyverse)

  cat('\n** fitting model...\n')
  cat('design matrix covariates:\n',colnames(design),'\n')

  # covert date of conception to radians (doc.theta) in [0,2*pi]
  day.of.yr = yday(design$doc)
  design$doc.theta = day.of.yr/366*2*pi
  design$sin.doc.theta = sin(design$doc.theta)
  design$cos.doc.theta = cos(design$doc.theta)
  design$sin.2doc.theta = sin(2*design$doc.theta)
  design$cos.2doc.theta = cos(2*design$doc.theta)
  design$doc = NULL
  design$doc.theta = NULL

  # create nested models for fourier fit
  design.base = design %>% # baseline model (covariates only)
    dplyr::select(
      -c(sin.doc.theta,cos.doc.theta,sin.2doc.theta,cos.2doc.theta)
    )
  design.fourier = design %>% # covariates + 1pr FTs
    dplyr::select(
      -c(sin.2doc.theta,cos.2doc.theta)
    )

  cat('performing likelihood ratio tests to identify
      models with significant seasonality...\n')

  lrt.pval = function(h0,h1){
    test = lrtest(h0,h1)
    return(test$`Pr(>Chisq)`[2])
  }

  cl <- makeCluster(16)
  registerDoParallel(cl)

  lrt.results =

    foreach (cg = row.names(meth), .combine='c', .packages = 'lmtest') %dopar%

    {

      # compare best Fourier model against baseline
      fit.base = lm(meth[cg,]~.,data=design.base)
      fit.fourier = lm(meth[cg,]~.,data=design.fourier)
      lrt.p = lrt.pval(fit.base,fit.fourier)
      c(cg,lrt.p)

    }
  stopCluster(cl)

  cat('lrt completed\n\n')

  lrt.results = as.data.frame(matrix(lrt.results,ncol = 2,byrow = T),stringsAsFactors = F)
  colnames(lrt.results) = c('cg','lrt.pval')
  row.names(lrt.results) = lrt.results$cg
  lrt.results$lrt.pval = as.numeric(lrt.results$lrt.pval)
  lrt.results$lrt.fdr = p.adjust(lrt.results$lrt.pval,method = 'fdr')

  return(lrt.results)
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Transform methylation M-values to Beta values
#'
#' @param M df or matrix of M-values
#'
#' @return df or matrix of Beta values
m2beta = function (M)
{
  beta <- 2^M/(2^M + 1)
  return(beta)
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Generate fitted values from Fourier regression model over [0,2*pi]
#'
#' @param model.fit Fourier model fit object of class 'lm'
#' @param design df model design matrix
#'
#' @return df of fitted values; rows = cpgs; cols = theta in radians
#'
get_fourier_fitted_vals <- function(model.fit = NULL, design = NULL) {

  # generate fitted M-values with intercept
  # do this in the interveal [0,2pi] for plotting, calculating amplitudes etc
  doc.theta.vals = seq(0,2*pi,0.02)
  design.fitted.vals = data.frame(
    intercept = rep(1,length(doc.theta.vals)),
    sin.doc.theta = sin(doc.theta.vals),
    cos.doc.theta = cos(doc.theta.vals)
  )

  coefs = t(coef(model.fit))
  fitted.vals =
    coefs[,c("`(Intercept)`","`sin(doc.theta)`","`cos(doc.theta)`")] %*%
    t(design.fitted.vals[,c('intercept','sin.doc.theta','cos.doc.theta')])

  # convert to Betas and centre
  fitted.vals = m2beta(fitted.vals)
  fitted.vals=as.data.frame(fitted.vals-rowMeans(fitted.vals))
  colnames(fitted.vals) = doc.theta.vals

  return(fitted.vals)

}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Get summary stats for fitted Fourier regression model
#'
#' @param fitted.vals
#' @param lrt.results
#' @param discovery.cgs
#' @param cn.var.cgs
#' @param cn.rand.cgs
#'
#' @return summary stats
#'
generate_model_fit_stats <- function(
  fitted.vals=NULL, lrt.results=NULL, discovery=NULL, ctrl.cgs=NULL, rand.cgs=NULL, soc=NULL)
{
  fitted.vals.max = apply(fitted.vals,1,function(x) max(as.numeric(x[1:(length(x)-2)])))
  fitted.vals.max.id = apply(fitted.vals,1,function(x) which.max(as.numeric(x[1:(length(x)-2)])))

  fitted.vals.min = apply(fitted.vals,1,function(x) min(as.numeric(x[1:(length(x)-2)])))
  fitted.vals.min.id = apply(fitted.vals,1,function(x) which.min(as.numeric(x[1:(length(x)-2)])))

  model.stats = data.frame(
    cg = row.names(fitted.vals),
    max.doc.theta =
      colnames(fitted.vals)[fitted.vals.max.id],
    min.doc.theta =
      colnames(fitted.vals)[fitted.vals.min.id],
    amplitude = abs(fitted.vals.max - fitted.vals.min),
    stringsAsFactors = F
  )

  # add LRT pvals
  model.stats$lrt.pval = lrt.results[row.names(model.stats),'lrt.pval']
  model.stats$lrt.fdr = lrt.results[row.names(model.stats),'lrt.fdr']

  model.stats$locus = NA
  model.stats[model.stats$cg %in% discovery.cgs,'locus'] = 'discovery'
  model.stats[model.stats$cg %in% ctrl.cgs,'locus'] = 'ctrl'
  model.stats[model.stats$cg %in% rand.cgs,'locus'] = 'rand'

  model.stats$replicated = FALSE
  model.stats[model.stats$cg %in% soc,'replicated'] = TRUE

  return(model.stats)
}
