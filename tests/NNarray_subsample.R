library(GpGp)

n1 <- 20
n2 <- 20
n <- n1*n2
locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
covparms <- c(2, 0.2, 0.75, 0)
y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
ord <- order_maxmin(locs)
NNarray <- find_ordered_nn(locs, 20)
loglik <- vecchia_meanzero_loglik_grad_info( covparms, "matern_isotropic",
    y, locs, NNarray )



loglik1 <- vecchia_meanzero_loglik_grad_info( covparms, "matern_isotropic",
                                              y, locs, NNarray[1 : (n / 2), ])
loglik2 <- vecchia_meanzero_loglik_grad_info( covparms, "matern_isotropic",
                                              y, locs, NNarray[(n / 2 + 1) : n, ])

loglik$loglik - loglik1$loglik - loglik2$loglik
loglik$grad - loglik1$grad - loglik2$grad
loglik$info - loglik1$info - loglik2$info