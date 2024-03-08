# Rupert's function
# See original code: https://github.com/uclahs-cds/project-ProstateCancer-m6A/blob/rhughwhite/somatic_driver_events/functions/RADAR_diffIP_parallel_handle_NA.R
diffIP_parallel_NA <- function(object, exclude = NULL, maxPsi = 100, fdrBy = 'fdr', thread = 8) {
  allY <- object@ip_adjExpr_filtered
  psi <- 10 # start point
  
  ## convert predictor variable if it is not numeric
  if (is.numeric(variable(object)[,1]) ) {
    X <- variable(object)[,1]
  } else {
    tmp <-  as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] )
    names(tmp) <- as.character(variable(object)[,1])
    cat('The predictor variable has been converted:\n')
    print(tmp)
    X <-  as.integer(as.character(variable(object)[,1]) == unique(as.character(variable(object)[,1]) )[2] ) # convert categorical variable into numerical variable.
  }
  
  #if( ncol(variable(object)) == 1 ){
  
  cat('running PoissonGamma test at single beta mode\n')
  
  start_time <- Sys.time()  # track run time
  ## register cluster for hyperthread computing
  doParallel::registerDoParallel(cores=thread)
  cat(paste('Hyper-thread registered:',getDoParRegistered(),'\n'))
  cat(paste('Using',getDoParWorkers(),'thread(s) to run PoissonGamma test...\n'))
  
  error.id <- NULL
  all.est <- foreach(kk = 1:nrow(allY), .combine = rbind, .errorhandling = 'remove') %dopar% {
    Y <- unlist(allY[kk, ])
    non.missing.data <- which(!is.na(Y));
    if (min(table(X[non.missing.data])) < 3) {
      est <- setNames(
        rep(NA, 8),
        c('beta', 'p_value', 'psi', 'alt_psi', 'mu2', 'alt_mu2', 'likelihood', 'alt_likelihood')
      );
      return(est)
      next
    }
    model1 <- glm(Y[non.missing.data] ~ X[non.missing.data], family = poisson(link = 'log'))
    coef <- model1$coefficients
    mu2 <- coef[1]
    beta <- coef[2]
    est <- try(unlist(PoissionGamma(Y[non.missing.data], X[non.missing.data], beta, psi, mu2, gamma = 0.75, steps = 50, down = 0.1,psi_cutoff = maxPsi)))
    if(class(est) == 'try-error'){
      error.id <- c(error.id, kk)
    }
    return(est)
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste('Time used to run PoissonGamma test:',difftime(end_time, start_time, units = 'mins'),'mins... \n'))
  
  all.id <- which(! 1:nrow(allY) %in% error.id )
  rownames(all.est) <- rownames(allY)[all.id]
  
  if (fdrBy == 'qvalue') {
    fdr <- qvalue::qvalue(all.est[,'p_value'] )$qvalue
    object@fdr.method = 'qvalue'
  } else {
    fdr <- p.adjust(all.est[,'p_value'],method = fdrBy )
    object@fdr.method = 'Benjamini & Hochberg'
  }
  
  object@test.est <- cbind(all.est, fdr)
  object@test.method <- 'PoissonGamma test (RADAR)'
  cat('\n')
  return(object)
}