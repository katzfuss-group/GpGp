#' Coordinate descent using the 2nd order approximation of the log-likelihood 
#' 
#' @param likfun likelihood function, returns log-likelihood, gradient, and FIM
#' @param start_parms starting values of parameters
#' @param lambda L1 penalty parameter
#' @param pen_idx the indices of parameters that will be penalized by the L1
#' @param epsl step size when coordinate descent does not reduce the obj func
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on the norm of grad
#' @param convtol2 convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of 2nd order approximations
#' @param max_iter2 maximum number of epochs in coordinate descent
coord_grad_L1_opt <- function(likfun, start_parms, lambda, pen_idx, epsl, 
                              silent = FALSE, convtol = 1e-4, convtol2 = 1e-4, 
                              max_iter = 40, max_iter2 = 40)
{
  parms <- start_parms
  if(lambda < 0 || epsl < 0)
    stop("lambda and epsl should both be greater than zero\n")
  for(i in 1 : max_iter)
  {
    likobj <- likfun(parms)
    
    # check for Inf, NA, or NaN
    if( !test_likelihood_object(likobj) ){
      stop("inf or na or nan in likobj\n")
    }
    
    obj <- -likobj$loglik + lambda * sum(abs(parms[pen_idx]))
    grad <- likobj$grad
    info <- likobj$info
    
    if(sqrt(sum(grad^2)) < convtol)
      break;
    
    # coordinate descent
    parms_new <- parms
    for(k in 1 : max_iter2)
    {
      parms_new_cp <- parms_new
      for(j in 1 : length(parms_new_cp))
      {
        
      }
      if(sqrt(sum((parms_new - parms_new_cp)^2)) < convtol2)
        break
    }
    
    
  }
}