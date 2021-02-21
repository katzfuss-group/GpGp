library(GpGp)
n <- 1e3
d <- 1e2

set.seed(123)
parmsRange <- runif(d + 2)
parmsRel <- parmsRange
parmsRel[2 : (d + 1)] <- 1 / parmsRel[2 : (d + 1)]
parmsSqRel <- parmsRel
parmsSqRel[2 : (d + 1)] <- parmsSqRel[2 : (d + 1)]^2
locs <- lhs::randomLHS(n, d)

covMRange <- GpGp::matern25_scaledim(parmsRange, locs)
covMRel <- GpGp::matern25_scaledim_relevance(parmsRel, locs)
covMSqRel <- GpGp::matern25_scaledim_sqrelevance(parmsSqRel, locs)
cat("Difference Norm", norm(covMRange - covMRel), 
    norm(covMRange - covMSqRel), "\n")

dcovMRange <- GpGp::d_matern25_scaledim(parmsRange, locs)
dcovMRel <- GpGp::d_matern25_scaledim_relevance(parmsRel, locs)
dcovMSqRel <- GpGp::d_matern25_scaledim_sqrelevance(parmsSqRel, locs)
for(j in 2 : (d + 1))
{
  cat("Difference norm", 
      norm(dcovMRange[, , j] + dcovMRel[, , j] * parmsRel[j]^2), 
      norm(dcovMRange[, , j] + 2 * dcovMSqRel[, , j] * parmsRel[j]^3),
      "\n")
}



















