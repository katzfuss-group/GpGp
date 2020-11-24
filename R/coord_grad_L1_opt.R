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
  compute_likobj = T
  for(i in 1 : max_iter)
  {
    if(compute_likobj)
      likobj <- likfun(parms)
    
    # check for Inf, NA, or NaN
    if( !test_likelihood_object(likobj) ){
      stop("inf or na or nan in likobj\n")
    }
    
    obj <- -likobj$loglik + lambda * sum(abs(parms[pen_idx]))
    grad <- likobj$grad
    info <- likobj$info
    
    if(!silent)
    {
        cat(paste0("Iter ", i, ": \n"))
        cat("pars = ",  paste0(round(parms, 4)), "  \n" )
        cat(paste0("obj = ", round(obj, 6), "  \n"))
        cat("grad = ")
        cat(as.character(round(grad, 3)))
        cat("\n")
    }
    
    gradCp <- grad
    gradCp[pen_idx] <- gradCp[pen_idx] - lambda
    if(sqrt(sum(gradCp^2)) < convtol)
      break
    
    # coordinate descent
    b <- - grad
    b[pen_idx] <- b[pen_idx] + lambda
    b <- b - as.vector(info %*% parms)
    coord_des_obj <- coord_quad_posi_domain(info, b, parms, silent, 
                                            convtol2, max_iter2, 1e6)
    # check if obj func decreases
    if(coord_des_obj$code < 2) # parms_new is valid
    {
        likobj_new <- likfun(coord_des_obj$parms)
        obj_new <- -likobj_new$loglik + lambda * sum(abs(coord_des_obj$parms[pen_idx]))
        if(obj_new < obj)
            grad_des <- F
        else
            grad_des <- T
    }
    else
        grad_des <- T
    
    if(grad_des)
    {
        if(!silent)
        {
            cat("Gradient descent is used\n")
        }
        parms <- parms - grad * epsl
        parms[parms < 0] <- 0
        compute_likobj <- T
    }else
    {
        parms <- coord_des_obj$parms
        likobj <- likobj_new
        compute_likobj <- F
    }
    if(!silent)
        cat("\n")
  }
  likobj <- likfun(parms)
  return(list(covparms = parms))
}

#' Coordinate descent for a quadratic function in the positive domain 
#' 
#' @param A 2nd order coefficient matrix (\frac{1}{2} x^\top A x)
#' @param b 1st order coefficient vector
#' @param start_parms starting values of parameters
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of epochs in coordinate descent
#' @param max_parm maximum parameter value
#' 
#' Return a list of two
#' 
#' @return code 0 if convtol is reached, 1 if max number of epochs reached, 2 parms become invalid
#' @return parms new parameter values
coord_quad_posi_domain <- function(A, b, start_parms, silent, convtol, max_iter, max_parm)
{
    parms_new <- start_parms
    for(k in 1 : max_iter)
    {
        parms_new_cp <- parms_new
        for(j in 1 : length(parms_new))
        {
            chg <- - sum(A[j, -j] * parms_new[-j])
            parms_new[j] <- max((- b[j] + chg) / A[j, j], 0)
            if(parms_new[j] > max_parm || is.na(parms_new[j]) || is.infinite(parms_new[j]))
                return(list(code = 2, parms = parms_new))
        }
        if(sqrt(sum((parms_new - parms_new_cp)^2)) < convtol)
            return(list(code = 0, parms = parms_new))
    }
    return(list(code = 1, parms = parms_new))
}

test_coord_quad_posi_domain <- function()
{
    set.seed(120)
    n <- 10
    tmpM <- matrix(runif(n * n), n, n)
    covM <- tmpM %*% t(tmpM)
    diagV <- diag(covM)
    diagV <- 1 / sqrt(diagV)
    covM <- diag(diagV) %*% covM %*% diag(diagV)
    
    b <- runif(n, min = -100, max = 100)
    
    obj_func <- function(x, A, b)
    {
        0.5 * t(x) %*% A %*% x + sum(x * b)
    }
    d_obj_func <- function(x, A, b)
    {
        as.vector(A %*% x) + b
    }
    
    coord_obj1 <- coord_quad_posi_domain(covM, b, rep(0.1, n), 2, 1e-4, 40, 1e6)
    coord_obj2 <- coord_quad_posi_domain(covM, b, rep(0.1, n), 2, 1e-4, 4, 1e6)
    coord_obj3 <- coord_quad_posi_domain(covM, b, rep(0.1, n), 2, 1e-4, 40, 10)
    
    optim_obj <- optim(rep(0.1, n), obj_func, d_obj_func, A = covM, b = b, method = "L-BFGS-B",
                       lower = rep(0, n))
}

# -(x^4 + y^4 + z^4 - 8*x*y*z - 3*x*y - 4*y*z - 5*x*z - x^3 - 2*y^3 - 3*z^3)
# -(4*x^3 - 8*y*z - 3*y - 5*z - 3*x^2, 4*y^3 - 8*x*z - 3*x - 4*z - 6*y^2, 4*z^3 - 8*x*y - 4*y - 5*x - 9*z^2)
# -(12*x^2 - 6*x, -8*z - 3, -8*y - 5)
# -(-8*z - 3, 12*y^2 - 12*y, -8*x - 4)
# -(-8*y - 5, -8*x - 4, 12*z^2 - 18*z)
test_coord_grad_L1_opt <- function()
{
    likfun <- function(x)
    {
        loglik <- -(x[1]^4 + x[2]^4 + x[3]^4 - 8*x[1]*x[2]*x[3] - 3*x[1]*x[2] - 
              4*x[2]*x[3] - 5*x[1]*x[3] - x[1]^3 - 2*x[2]^3 - 3*x[3]^3)
        grad <- -c(4*x[1]^3 - 8*x[2]*x[3] - 3*x[2] - 5*x[3] - 3*x[1]^2, 
                   4*x[2]^3 - 8*x[1]*x[3] - 3*x[1] - 4*x[3] - 6*x[2]^2, 
                   4*x[3]^3 - 8*x[1]*x[2] - 4*x[2] - 5*x[1] - 9*x[3]^2)
        info <- matrix(c(12*x[1]^2 - 6*x[1], -8*x[3] - 3, -8*x[2] - 5, 
                         -8*x[3] - 3, 12*x[2]^2 - 12*x[2], -8*x[1] - 4, 
                         -8*x[2] - 5, -8*x[1] - 4, 12*x[3]^2 - 18*x[3]), 
                       3, 3, byrow = T)
        return(list(loglik = loglik, grad = grad, info = info))
    }
    
    func <- function(x, lambda)
    {
        x[1]^4 + x[2]^4 + x[3]^4 - 8*x[1]*x[2]*x[3] - 3*x[1]*x[2] - 
            4*x[2]*x[3] - 5*x[1]*x[3] - x[1]^3 - 2*x[2]^3 - 3*x[3]^3 + 
            lambda * sum(x)
    }
    dfunc <- function(x, lambda)
    {
        c(4*x[1]^3 - 8*x[2]*x[3] - 3*x[2] - 5*x[3] - 3*x[1]^2, 
           4*x[2]^3 - 8*x[1]*x[3] - 3*x[1] - 4*x[3] - 6*x[2]^2, 
           4*x[3]^3 - 8*x[1]*x[2] - 4*x[2] - 5*x[1] - 9*x[3]^2) + lambda
    }
    
    parms_init <- runif(3)
    lambda <- 1
    pen_idx <- c(1 : 3)
    epsl <- 1e-3
    silent <- F
    convtol <- 1e-4
    convtol2 <- 1e-4
    max_iter <- 40
    max_iter2 <- 40
    coord_obj <- coord_grad_L1_opt(likfun, parms_init, lambda, pen_idx, 
                      epsl, silent, convtol, convtol2, max_iter, max_iter2)
    
    optim_obj <- optim(parms_init, func, dfunc, lambda = lambda, method = "L-BFGS-B",
                       lower = rep(0, n))
}

# -(x^2 + y^2 + z^2 - x*y - 0.3*y*z - 0.5*z*x - 3*x - 4*y - 5*z)
# -(2*x - y - 0.5*z - 3, 2*y - x - 0.3*z - 4, 2*z - 0.3*y - 0.5*x - 5)
# -(2, -1, -0.5)
# -(-1, 2, -0.3)
# -(-0.5, -0.3, 2)
test_coord_grad_L1_opt_2nd <- function()
{
    likfun <- function(x)
    {
        loglik <- -(x[1]^2 + x[2]^2 + x[3]^2 - x[1]*x[2] - 0.3*x[2]*x[3] - 0.5*x[3]*x[1] - 3*x[1] - 4*x[2] - 5*x[3])
        grad <- -c(2*x[1] - x[2] - 0.5*x[3] - 3, 2*x[2] - x[1] - 0.3*x[3] - 4, 2*x[3] - 0.3*x[2] - 0.5*x[1] - 5)
        info <- matrix(c(2, -1, -0.5, 
                         -1, 2, -0.3, 
                         -0.5, -0.3, 2), 
                       3, 3, byrow = T)
        return(list(loglik = loglik, grad = grad, info = info))
    }
    
    func <- function(x, lambda)
    {
        x[1]^2 + x[2]^2 + x[3]^2 - x[1]*x[2] - 0.3*x[2]*x[3] - 0.5*x[3]*x[1] - 3*x[1] - 4*x[2] - 5*x[3] + 
            lambda * sum(x)
    }
    dfunc <- function(x, lambda)
    {
        c(2*x[1] - x[2] - 0.5*x[3] - 3, 2*x[2] - x[1] - 0.3*x[3] - 4, 2*x[3] - 0.3*x[2] - 0.5*x[1] - 5) + lambda
    }
    
    parms_init <- runif(3)
    lambda <- 1
    pen_idx <- c(1 : 3)
    epsl <- 1e-3
    silent <- F
    convtol <- 1e-4
    convtol2 <- 1e-4
    max_iter <- 40
    max_iter2 <- 40
    coord_obj <- coord_grad_L1_opt(likfun, parms_init, lambda, pen_idx, 
                                   epsl, silent, convtol, convtol2, max_iter, max_iter2)
    
    optim_obj <- optim(parms_init, func, dfunc, lambda = lambda, method = "L-BFGS-B",
                       lower = rep(0, n))
}





















