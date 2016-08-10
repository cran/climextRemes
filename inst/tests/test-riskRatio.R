require(testthat)

# POT and binom for yearly data

context("Testing risk ratio calculations.")

set.seed(0)
n <- 500
y1 <- rnorm(n, 0.25, 1)
y0 <- rnorm(n, 0, 1)

cutoff <- 2

yC <- c(sum(y1 > cutoff), sum(y0 > cutoff))

thr <- 1.4
wh1 <- y1 > thr
wh0 <- y0 > thr
y1e <- y1[wh1]
y0e <- y0[wh0]
bi1 <- (1:n)[wh1]
bi2 <- (1:n)[wh0]

outb <- calc_riskRatio_binom(yC, rep(n, 2), bootSE = TRUE, lrtCI = TRUE)
outp <- calc_riskRatio_pot(cutoff, y1e, y0e, threshold1 = thr, nBlocks1 = n,
                           nBlocks2 = n, bootSE = FALSE, lrtCI = TRUE)

test_that(paste0("binom RR uncertainty methods"), {
    expect_lt(max(abs(outb$ci_riskRatio - outb$ci_riskRatio_lrt)), (.25))
})
test_that(paste0("binom RR uncertainty methods"), {
    expect_lt(max(abs(outb$ci_riskRatio - outb$ci_riskRatio_boot)), (.5))
})
test_that(paste0("POT stationary RR uncertainty methods"), {
    expect_lt(max(abs(outp$ci_riskRatio - outp$ci_riskRatio_lrt)), (.25))
})
# boot getting -Inf for some logReturnProbs so get NaN for SE

# try repeated sampling

m <- 500
set.seed(0)

rrp <- se_rrp <- rep(NA, m)
rrb <- se_rrb <- se_rrb_b <- rep(NA, m)

for( i in 1:m ) {
    set.seed(i)
    n <- 300
    y1 <- rnorm(n, 0.25, 1)
    y0 <- rnorm(n, 0, 1)
    
    cutoff <- 2
    
    yC <- c(sum(y1 > cutoff), sum(y0 > cutoff))
    
    thr <- 1.4
    wh1 <- y1 > thr
    wh0 <- y0 > thr
    y1e <- y1[wh1]
    y0e <- y0[wh0]
    bi1 <- (1:n)[wh1]
    bi2 <- (1:n)[wh0]
    
    outb <- calc_riskRatio_binom(yC, rep(n, 2), bootSE = TRUE, lrtCI = FALSE)
    outp <- calc_riskRatio_pot(cutoff, y1e, y0e, threshold1 = thr, nBlocks1 = n,
                               nBlocks2 = n, bootSE = FALSE, lrtCI = FALSE)
    rrb[i] <- outb$logRiskRatio
    se_rrb[i] <- outb$se_logRiskRatio
    se_rrb_b[i] <- outb$se_logRiskRatio_boot
    rrp[i] <- outp$logRiskRatio
    se_rrp[i] <- outp$se_logRiskRatio
}

wh <- rrb < Inf & rrb > 0
sd <- sd(rrb[wh], na.rm = TRUE)
se <- mean(se_rrb[wh], na.rm = TRUE)
test_that(paste0("test delta se(logRR) for binom"), {
    expect_lt(max(abs(sd - se)), (.08)) 
})
se <- mean(se_rrb_b[wh], na.rm = TRUE)
test_that(paste0("test boot se(logRR) for binom"), {
    expect_lt(max(abs(sd - se)), (.04)) 
})
wh <- !is.na(rrp)
sd <- sd(rrp[wh], na.rm = TRUE)
se <- mean(se_rrp[wh], na.rm = TRUE)
test_that(paste0("test delta se(logRR) for POT"), {
    expect_lt(max(abs(sd - se)), (.02)) 
})


set.seed(1)
cutoff <- 3

nT <- 300
nObs <- 100
mn <- rep(c(-0.25,.25), each = nT*nObs/2)
y <- matrix(rnorm(nT*nObs, mn),  nrow = nObs)

grps <- rep(c(1,2), each = nT*nObs/2)
yrgrps <- rep(c(1,2), each = nT/2)
yMax <- apply(y, 2, max)

blockIndex <- matrix(rep(1:nT, nObs), nrow = nObs, byrow = TRUE)
thr <- quantile(y, 0.97)
wh <- y > thr
yExc <- y[wh]
grps <- grps[wh]
blockIndexObs <- blockIndex[wh]
index <- (1:(nT*nObs))[wh]

p1 = mean(yMax[yrgrps==1] > cutoff)
p2 = mean(yMax[yrgrps==2] > cutoff)
rpd = log(p1)-log(p2)
se1 = sqrt(p1*(1-p1)/(nT/2))
se2 = sqrt(p2*(1-p2)/(nT/2))

yC <- c(sum(y[grps == 1] > cutoff), sum(y[grps==2] > cutoff))

rp = 20
rv1 = quantile(yMax[yrgrps==1], 1-1/rp)
rv2 = quantile(yMax[yrgrps==2], 1-1/rp)
rvd = rv1 - rv2

w <- rnorm(nT/2)
z <- rnorm(nT/2)

wNew = c(.001, -.001)
zNew = c(.001, -.001)

weights1 = runif(nT/2, 0.5, 1.5)
weights2 = runif(nT/2, 0.5, 1.5)
pm1 = runif(nT/2, 0, 0.3)
pm2 = runif(nT/2, 0, 0.3)

outg <- calc_riskRatio_gev(yMax[yrgrps == 2], yMax[yrgrps == 1], x1 = data.frame(w = w),
                   x2 = data.frame(w= w, z = z),
                   locationFun1 = ~ w, scaleFun1 = ~w, shapeFun1 = ~1,
                   locationFun2 = ~ w+z, scaleFun2 = ~w+z, shapeFun2 = ~w+z,
                   xNew1 = data.frame(w = wNew), xNew2 = data.frame(w = wNew, z = zNew),
                   returnValue = cutoff, bootSE = TRUE,
                   lrtCI = TRUE, optimArgs = list(method = 'BFGS'))
test_that(paste0("GEV nonstationary RR methods"), {
    expect_lt(max(abs(outg$ci_riskRatio[1,] - outg$ci_riskRatio_lrt[1, ])), (.4))
})

outg <- calc_riskRatio_gev(yMax[yrgrps == 2], yMax[yrgrps == 1], x1 = NULL, x2 = data.frame(z = z),
                   locationFun1 = ~ 1, scaleFun1 = ~1, shapeFun1 = ~1,
                   locationFun2 = ~ z, scaleFun2 = ~z, shapeFun2 = ~z,
                   xNew1 = NULL, xNew2 = data.frame(z = zNew[1]),
                   returnValue = cutoff, bootSE = TRUE,
                   lrtCI = TRUE)

test_that(paste0("GEV nonstationary RR methods, one stationary"), {
    expect_lt(max(abs(outg$ci_riskRatio - outg$ci_riskRatio_lrt)), (.5))
})

outp <- calc_riskRatio_pot(yExc[grps == 1], yExc[grps == 2],
                   x1 = data.frame(w = w), x2 = data.frame(w= w, z = z), nBlocks1 = nT/2, nBlocks2 = nT/2,
                   blockIndex1 = blockIndexObs[grps == 1], blockIndex2 = blockIndexObs[grps==2]-(nT/2),
                   threshold1 = thr, 
                   locationFun1 = ~ w, scaleFun1 = ~w, shapeFun1 = ~1,
                   locationFun2 = ~ w+z, scaleFun2 = ~w+z, shapeFun2 = ~w+z,
                   xNew1 = data.frame(w = wNew), xNew2 = data.frame(w = wNew, z = zNew),
                   returnValue = cutoff, bootSE = TRUE,
                   lrtCI = TRUE, optimArgs = list(method = 'BFGS'))
test_that(paste0("POT nonstationary RR methods"), {
    expect_lt(max(abs(outp$ci_riskRatio[1,] - outp$ci_riskRatio_lrt[1, ])), (.01))
})
test_that(paste0("POT nonstationary RR methods"), {
    expect_lt(max(abs(outp$ci_riskRatio[1,] - outp$ci_riskRatio_boot[1, ])), (.03))
})

# bootstrap has -Inf for group 1 so not doing
outp <- calc_riskRatio_pot(yExc[grps == 1], yExc[grps == 2],
                   x1 = NULL, x2 = data.frame(w= w, z = z), nBlocks1 = nT/2, nBlocks2 = nT/2,
                   blockIndex1 = blockIndexObs[grps == 1], blockIndex2 = blockIndexObs[grps==2]-(nT/2),
                   threshold1 = thr, 
                   locationFun1 = ~ 1, scaleFun1 = ~1, shapeFun1 = ~1,
                   locationFun2 = ~ w+z, scaleFun2 = ~w+z, shapeFun2 = ~w+z,
                   xNew1 = NULL, xNew2 = data.frame(w = wNew[1], z = zNew[1]),
                   returnValue = cutoff, bootSE = TRUE,
                   lrtCI = TRUE, optimArgs = list(method = 'BFGS'))
test_that(paste0("POT nonstationary RR methods, one stationary"), {
    expect_lt(max(abs(outp$ci_riskRatio - outp$ci_riskRatio_lrt)), (.01))
})
test_that(paste0("POT nonstationary RR methods, one stationary"), {
    expect_lt(max(abs(outp$ci_riskRatio - outp$ci_riskRatio_boot)), (.01))
})


# try with pm, weights
outp <- calc_riskRatio_pot(yExc[grps == 1], yExc[grps == 2],
                   x1 = NULL, x2 = data.frame(w= w, z = z), nBlocks1 = nT/2, nBlocks2 = nT/2,
                   blockIndex1 = blockIndexObs[grps == 1], blockIndex2 = blockIndexObs[grps==2]-(nT/2),
                   threshold1 = thr, weights1 = weights1, weights2 = weights2,
                   proportionMissing1 = pm1, proportionMissing2 = pm2,
                   locationFun1 = ~ 1, scaleFun1 = ~1, shapeFun1 = ~1,
                   locationFun2 = ~ w+z, scaleFun2 = ~w+z, shapeFun2 = ~w+z,
                   xNew1 = NULL, xNew2 = data.frame(w = wNew[1], z = zNew[1]),
                   returnValue = cutoff, bootSE = TRUE,
                   lrtCI = TRUE, optimArgs = list(method = 'BFGS'))
test_that(paste0("POT nonstationary RR methods, one stationary"), {
    expect_lt(max(abs(outp$ci_riskRatio - outp$ci_riskRatio_lrt)), (.01))
})

