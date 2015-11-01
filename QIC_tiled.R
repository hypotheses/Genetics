    
    fname<-"gaussian"
    geemod<-gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="exchangeable")
    W <- geemod$naive.variance
#   if (fname == "binomial")
#     W <- summary(glm(formula, family = quasibinomial, data = data))$cov.scaled
    N <- geemod$nobs

    ## compute QIC:
    Y <- geemod$y
    MU <- geemod$fitted.values
    Qlik <- switch(fname,
                   "gaussian" = -sum((Y - MU)^2)/2,
                   "binomial" = sum(Y*log(MU/(1 - MU)) + log(1 - MU)),
                   "poisson" = sum(Y*log(MU) - MU),
                   "Gamma" = sum(Y/MU + log(MU)),
                   "inverse.gaussian" = sum(-Y/(2*MU^2) + 1/MU))
    Ai <- gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="independence")$naive.variance
    QIC <- -2*Qlik + 2*sum(diag(solve(Ai) %*% W))

geenull<-gee(breaks ~ 1, id=wool, data=warpbreaks, corstr="independence")
    W.null <- geenull$naive.variance
    Y <- geenull$y
    MU <- geenull$fitted.values
    Qlik.null <- switch(fname,
                   "gaussian" = -sum((Y - MU)^2)/2,
                   "binomial" = sum(Y*log(MU/(1 - MU)) + log(1 - MU)),
                   "poisson" = sum(Y*log(MU) - MU),
                   "Gamma" = sum(Y/MU + log(MU)),
                   "inverse.gaussian" = sum(-Y/(2*MU^2) + 1/MU))
    Ai.null <- gee(breaks ~ 1, id=wool, data=warpbreaks, corstr="independence")$naive.variance
    QIC.null <- -2*Qlik.null + 2*sum(diag(solve(Ai.null) %*% W.null))



stat<-2*(QIC-QIC.null)
1-2*pchisq(abs(stat),1,lower.tail=F)


  