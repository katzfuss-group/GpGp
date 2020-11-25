#' Fisher scoring algorithm
#' 
#' @param likfun likelihood function, returns likelihood, gradient, and hessian
#' @param start_parms starting values of parameters
#' @param link link function for parameters (used for printing)
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on step dot grad
#' @param max_iter maximum number of Fisher scoring iterations
fisher_scoring_modified <- function( likfun, start_parms, 
                            silent = FALSE, convtol = 1e-4, max_iter = 40 ){
  
  # evaluate function at initial values
  parms <- start_parms
  likobj <- likfun(parms)
  
  # test likelihood object    
  if( !test_likelihood_object(likobj) ){
    parms <- 0.1*parms
    likobj <- likfun(parms)
  }
  
  # assign loglik, grad, and info
  loglik <- likobj$loglik        
  grad <- likobj$grad
  info <- likobj$info
  
  # add a small amount of regularization
  diag(info) <- diag(info) + 0.1*min(diag(info))
  
  # print some stuff out
  if(!silent){
    cat(paste0("Iter ",0,": \n"))
    cat("pars = ",  paste0(round(parms ,4)), "  \n" )
    cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
    cat("grad = ")
    cat(as.character(round(-grad,3)))
    cat("\n\n")
  }
  
  for(j in 1:max_iter){
    
    likobj0 <- likobj
    
    # if condition number of info matrix large, 
    # then gradient descent
    tol <- 1e-10
    if (condition_number(info) > 1 / tol) {
      if (!silent) cat("Cond # of info matrix > 1/tol \n")
      #info <- 1.0*max(likobj0$info)*diag(nrow(likobj0$info))
      # regularize
      diag(info) <- diag(info) + tol*max(diag(info))
    }
    
    # calculate fisher step 
    step <- -solve(info, grad)
    
    # if step size large, then make it smaller
    if (mean(step^2) > 1) {
      if(!silent) cat("##\n")
      step <- step/sqrt(mean(step^2))
    }
    
    # take step and calculate loglik, grad, and info
    newparms <- parms + step
    likobj <- likfun(newlogparms)
    
    # check for Inf, NA, or NaN
    cnt <- 1
    while (!test_likelihood_object(likobj)) {
      if (!silent) cat("inf or na or nan in likobj\n")
      step <- 0.5 * step
      newparms <- parms + step
      likobj <- likfun(newlogparms)
      if (cnt == 10) { stop("could not avoid inf, na or nan\n") }
    }
    
    # Check the wolfe conditions
    # if not satisfied, shrink fisher step
    # after some iterations, switch to gradient
    cnt <- 1
    no_decrease <- FALSE
    both <- FALSE
    mult <- 1.0
    stepgrad <- c(crossprod(step,grad))
    
    # redefine parms, loglik, grad, info
    parms <- parms + step
    loglik <- likobj$loglik        
    grad <- likobj$grad
    info <- likobj$info
    
    # print some stuff out
    if(!silent){
      cat(paste0("Iter ",j,": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("loglik = ", round(-loglik,6),         "  \n"))
      cat("grad = ")
      cat(as.character(round(grad,4)),"\n")
      cat("step dot grad = ",stepgrad,"\n")
      cat("\n")
    }
    
    # if gradient is small, then break and return the results        
    if( abs(stepgrad) < convtol || no_decrease ){
      break
    }
  }
  
  
  ret <- list(
    covparms = parms, 
    loglik = loglik,
    grad = likobj$grad,
    info = likobj$info
  )
  return(ret)
}






