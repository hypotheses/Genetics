#mlinfos <- system("cat ~/newhome/aa_imputed_genotypes/mlinfofile.list",intern=T )
#mldoses <- system("cat ~/newhome/aa_imputed_genotypes/mldosefile.list",intern=T)
## smallest file = 3_31 ## 61
#i=61
#for i in `seq 1 1`
#   do
#  echo "mlinfo= ${mlinfos[$i]}"
#  echo "mldose= ${mldoses[$i]}"
## mach.lme.batch(
##         phenfile="/research/2nd_imp_data/users/lyanek/crpphen.csv";
## 	genfile="/research/2nd_imp_data/users/bhoom/aa_imputed_genotypes/c3/machc3_31.out.mldose.GeneSTAR.gz";
## 	pedfile="/research/2nd_imp_data/users/lyanek/crpped.csv";
## 	phen="CRPRES";
## 	kinmat="/research/2nd_imp_data/users/lyanek/gecrp.kinship.Rdata";
##         mapfile="/research/2nd_imp_data/users/bhoom/aa_imputed_genotypes/c3/machc3_31.out.mlinfo.gz";
## 	model="a";
## 	outfile="/tmp/crplme-black-63.csv";
## 	col.names=T;
## 	sep.ped=",";
## 	sep.phe=",";
## 	sep.gen=" "
## x <- read.in.data(phenfile,genfile,pedfile,mapfile,sep.ped,sep.phe,sep.gen)
## y <- data.frame(SNP=x$snps,N=numeric(length(x$snps)))
print("my.lme.batch.R revision 11/18/10 -- fixed read.in.data")
read.in.data <- function(phenfile, genfile, pedfile, mapfile,sep.ped=",", sep.phe=",", sep.gen=" ") {
        print("Reading in Data")
        ped.dat <- read.table(genfile, header = FALSE, sep = sep.gen, as.is=T)
#        ped.dat[,1] <- matrix(unlist(strsplit(ped.dat[,1],"->")),byrow=T,ncol=2)[,2]
#        ped.dat[,2] <- NULL
        mlinfo <- read.table(mapfile,header=T,as.is=T,)
#        mlinfo[mlinfo[,4]>mlinfo[,5],"MA"] <- mlinfo[(mlinfo[,4]>mlinfo[,5]),3]
#        mlinfo[mlinfo[,4]<=mlinfo[,5],"MA"] <- mlinfo[(mlinfo[,4]<=mlinfo[,5]),2]
        mlinfo[mlinfo[,2]==1,2] <- "A"
        mlinfo[mlinfo[,2]==2,2] <- "C"
        mlinfo[mlinfo[,2]==3,2] <- "G"
        mlinfo[mlinfo[,2]==4,2] <- "T"
        mlinfo[mlinfo[,3]==1,3] <- "A"
        mlinfo[mlinfo[,3]==2,3] <- "C"
        mlinfo[mlinfo[,3]==3,3] <- "G"
        mlinfo[mlinfo[,3]==4,3] <- "T"
        snp.names <- paste(gsub("X","",mlinfo[,1]),mlinfo[,2],sep="_")
        names(ped.dat) <- c("id",snp.names)
        pedigree <- read.table(pedfile, header = TRUE, sep = sep.ped)
        gntp.all <- merge(pedigree, ped.dat, by = "id")
        phen.dat = read.table(phenfile, header = TRUE, sep = sep.phe)
        phen.name = colnames(phen.dat)[-1]
        n.snp = length(names(gntp.all))
        if (length(grep("^sex$", colnames(phen.dat))) == 0) {
            phensnp.dat <- merge(gntp.all, phen.dat, by = c("id"))
        }
        else {
            phensnp.dat <- merge(gntp.all, phen.dat, by = c("id", 
                "sex"))
        }
        print("Done reading in data")
        return(list(data = phensnp.dat, snps = snp.names, phen.name = phen.name))
    }
## mincount implemented here the count of all genotypes coded as 1 and 2 (where 2 is rare homozygous genotype)
mach.lme.batch <- function (phenfile, genfile, pedfile, phen,mapfile, kinmat, model = "a",
                          covars = NULL, outfile, col.names = T, sep.ped = ",", sep.phe = ",",
                          sep.gen = " ",mincount=10)
{
    if (!model %in% c("g", "r", "a", "d")) 
        stop("please specify model as \"a\",\"g\",\"r\" or \"d\" only")
    trykin <- try(load(kinmat))
    if (inherits(trykin, "try-error")) 
        stop(paste("kinship matrix does not exist at ", kinmat))
    cor.snp <- function(y, x) {
        if (!is.numeric(y)) 
            y <- as.numeric(as.factor(y))
        return(sd(y) == 0 || abs(cor(y, x, use = "complete")) > 
            0.99999999)
    }
    assign("phen", phen, env = .GlobalEnv, inherits = T)
    phensnp.dat <- read.in.data(phenfile, genfile, pedfile,mapfile,sep.ped,sep.phe,sep.gen)
    snplist <- phensnp.dat$snps
    if (is.null(covars)) 
        phenlist <- phensnp.dat$phen.name
    else if (!is.null(covars) & sum(phensnp.dat$phen.name %in% 
        covars) == length(covars)) 
        phenlist <- phensnp.dat$phen.name[!phensnp.dat$phen.name %in% 
            covars]
    else stop("some covariates are not available")
    test.dat <- phensnp.dat$data
    assign("test.dat", test.dat, env = .GlobalEnv, inherits = T)
    if (!is.null(covars) & sum(snplist %in% covars) >= 1) {
        names(test.dat)[which(names(test.dat) == paste(snplist[snplist %in% 
            covars], ".x", sep = ""))] <- snplist[snplist %in% 
            covars]
        covars[covars %in% snplist] <- paste(covars[covars %in% 
            snplist], ".y", sep = "")
    }
    if (!is.null(covars)) {
        covars.dat <- na.omit(test.dat[, covars])
        single.cov <- F
        if (length(covars) == 1) 
            single.cov <- var(covars.dat) == 0
        else {
          ## <note> Check if variance of any covariates is "ZERO"
            single.cov <- any(apply(covars.dat, 2, var) == 0)
            if (single.cov) 
                stop(paste("Single category in covariates!"))
            ## <note> Check if correlation between each covariate is 0.99999999 or more
            for (i in covars) {
                cov1 <- covars.dat[, i]
                if (!is.numeric(cov1)) 
                  cov1 <- as.numeric(as.factor(cov1))
                for (j in covars[covars != i]) {
                  cov2 <- covars.dat[, j]
                  if (!is.numeric(cov2)) 
                    cov2 <- as.numeric(as.factor(cov2))
                  if (abs(cor(cov1, cov2)) > 0.99999999) 
                    stop(paste("Highly correlated covariates ", 
                      i, " and ", j, "!!", sep = ""))
                }
            }
        }
    }
    idlab <- "id"
    result <- NULL
    for (i in snplist) {
        assign("i", i, env = .GlobalEnv, inherits = T)
        if (is.null(covars)) 
            test2.dat <- na.omit(test.dat[, c(i, phen, idlab)])
        else {
            test2.dat <- na.omit(test.dat[, c(i, phen, idlab, 
                covars)])
            x.covar <- as.matrix(test2.dat[, covars])
            assign("x.covar", x.covar, env = .GlobalEnv, inherits = T)
        }
        id <- test2.dat[, idlab]
        assign("test2.dat", test2.dat, env = .GlobalEnv, inherits = T)
        assign("id", id, env = .GlobalEnv, inherits = T)
        if (is.null(covars)) 
            v.cov <- sum(try(lmekin(test2.dat[, phen] ~ 1, random = ~1 | 
                id, varlist = kmat, na.action = na.omit))$theta)
        else {
            lme.cov.out <- try(lmekin(test2.dat[, phen] ~ x.covar, 
                random = ~1 | id, varlist = kmat, na.action = na.omit))
            v.cov <- sum(lme.cov.out$theta)
        }
        ## <note> modify count ---> round it to fit dosage data
        count <- table(round(test2.dat[, i]))
        gntps <- names(count)
        count1 <- rep(0, 3)
        count1[as.numeric(gntps) + 1] <- round(count)
        if (!is.null(covars) & length(count) == 1) 
            colinear <- F
        else if (!is.null(covars) & length(covars) > 1 & length(count) > 
            1)
          ## <note> test for colinear with SNP i
            colinear <- apply(x.covar, 2, cor.snp, x = test2.dat[, 
                i])
        else if (!is.null(covars) & length(covars) == 1 & length(count) > 
            1) 
            colinear <- cor.snp(x.covar, test2.dat[, i])
        else if (is.null(covars)) 
            colinear <- F
        ## <note> summary if there is any colinearity in the data
        if (sum(colinear) > 0) {
            if (model %in% c("a", "r", "d")) 
                result <- rbind(result, c(phen, i, count1, rep(NA, 
                  7)))
            else result <- rbind(result, c(phen, i, count1, rep(NA,11)))
        }
        else {
            if (sort(count1)[1] + sort(count1)[2] < mincount) {
                if (model %in% c("a", "r", "d")) 
                  result <- rbind(result, c(phen, i, count1, 
                    rep(NA, 7)))
                else result <- rbind(result, c(phen, i, count1, 
                  rep(NA, 11)))
            }
            else {
                if (model == "a") {
                  mod.lab <- "additive"
                  if (is.null(covars)) 
                    lme.out <- try(lmekin(test2.dat[, phen] ~ 
                      test2.dat[, i], random = ~1 | id, varlist = kmat, 
                      na.action = na.omit))
                  else lme.out <- try(lmekin(test2.dat[, phen] ~ 
                    test2.dat[, i] + x.covar, random = ~1 | id, 
                    varlist = kmat, na.action = na.omit))
                  chisq <- lme.out$ctable[2, 1]^2/lme.out$var[2, 
                    2]
                  tmp <- c(max(v.cov - sum(lme.out$theta), 0)/var(test2.dat[, 
                    phen]), lme.out$ctable[2, 1], sqrt(lme.out$var[2, 
                    2]), chisq, 1, mod.lab, pchisq(chisq, 1, 
                    lower.tail = F))
                }
                else if (model == "g") {
                  if (min(count1) < 10) {
                    snp.i <- test2.dat[, i]
                    snp.i[snp.i == 2] <- 1
                    assign("snp.i", snp.i, env = .GlobalEnv, 
                      inherits = T)
                    if (is.null(covars)) 
                      lme.out <- try(lmekin(test2.dat[, phen] ~ 
                        snp.i, random = ~1 | id, varlist = kmat, 
                        na.action = na.omit))
                    else lme.out <- try(lmekin(test2.dat[, phen] ~ 
                      snp.i + x.covar, random = ~1 | id, varlist = kmat, 
                      na.action = na.omit))
                    chisq <- lme.out$ctable[2, 1]^2/lme.out$var[2, 
                      2]
                    tmp <- c(max(v.cov - sum(lme.out$theta), 
                      0)/var(test2.dat[, phen]), lme.out$ctable[2, 
                      1], NA, NA, sqrt(lme.out$var[2, 2]), NA, 
                      NA, chisq, 1, "dominant", pchisq(chisq, 
                        1, lower.tail = F))
                  }
                  else {
                    if (is.null(covars)) 
                      lme.out <- try(lmekin(test2.dat[, phen] ~ 
                        factor(test2.dat[, i]), random = ~1 | 
                        id, varlist = kmat, na.action = na.omit))
                    else lme.out <- try(lmekin(test2.dat[, phen] ~ 
                      factor(test2.dat[, i]) + x.covar, random = ~1 | 
                      id, varlist = kmat, na.action = na.omit))
                    chisq <- lme.out$ctable[2:3, 1] %*% solve(lme.out$var[2:3, 
                      2:3]) %*% lme.out$ctable[2:3, 1]
                    tmp <- c(max(v.cov - sum(lme.out$theta), 
                      0)/var(test2.dat[, phen]), lme.out$ctable[2:3, 
                      1], lme.out$ctable[3, 1] - lme.out$ctable[2, 
                      1], sqrt(diag(lme.out$var)[2:3]), sqrt(diag(lme.out$var)[2] + 
                      diag(lme.out$var)[3] - 2 * lme.out$var[2, 
                      3]), chisq, 2, "general", pchisq(chisq, 
                      2, lower.tail = F))
                  }
                }
                else if (model == "d") {
                  snp.i <- test2.dat[, i]
                  snp.i[snp.i == 2] <- 1
                  mod.lab <- "dominant"
                  assign("snp.i", snp.i, env = .GlobalEnv, inherits = T)
                  if (is.null(covars)) 
                    lme.out <- try(lmekin(test2.dat[, phen] ~ 
                      snp.i, random = ~1 | id, varlist = kmat, 
                      na.action = na.omit))
                  else lme.out <- try(lmekin(test2.dat[, phen] ~ 
                    snp.i + x.covar, random = ~1 | id, varlist = kmat, 
                    na.action = na.omit))
                  chisq <- lme.out$ctable[2, 1]^2/lme.out$var[2, 
                    2]
                  tmp <- c(max(v.cov - sum(lme.out$theta), 0)/var(test2.dat[, 
                    phen]), lme.out$ctable[2, 1], sqrt(lme.out$var[2, 
                    2]), chisq, 1, mod.lab, pchisq(chisq, 1, 
                    lower.tail = F))
                }
                else if (model == "r") {
                  if (count1[3] < 10) 
                    tmp <- rep(NA, 7)
                  else {
                    snp.i <- test2.dat[, i]
                    snp.i[snp.i == 1] <- 0
                    snp.i[snp.i == 2] <- 1
                    mod.lab <- "recessive"
                    assign("snp.i", snp.i, env = .GlobalEnv, 
                      inherits = T)
                    if (is.null(covars)) 
                      lme.out <- try(lmekin(test2.dat[, phen] ~ 
                        snp.i, random = ~1 | id, varlist = kmat, 
                        na.action = na.omit))
                    else lme.out <- try(lmekin(test2.dat[, phen] ~ 
                      snp.i + x.covar, random = ~1 | id, varlist = kmat, 
                      na.action = na.omit))
                    chisq <- lme.out$ctable[2, 1]^2/lme.out$var[2, 
                      2]
                    tmp <- c(max(v.cov - sum(lme.out$theta), 
                      0)/var(test2.dat[, phen]), lme.out$ctable[2, 
                      1], sqrt(lme.out$var[2, 2]), chisq, 1, 
                      mod.lab, pchisq(chisq, 1, lower.tail = F))
                  }
                }
                if (class(tmp) == "try-error") {
                  if (model %in% c("a", "r", "d")) 
                    result <- rbind(result, c(phen, i, count1, 
                      rep(NA, 7)))
                  else result <- rbind(result, c(phen, i, count1, 
                    rep(NA, 11)))
                }
                else result <- rbind(result, c(phen, i, count1, 
                  tmp))
            }
        }
    }
    if (model %in% c("a", "d", "r")) {
        colnames(result) <- c("phen", "snp", "n0", "n1", "n2", 
            "h2q", "beta", "se", "chisq", "df", "model", "pval")
    }
    else colnames(result) <- c("phen", "snp", "n0", "n1", "n2", 
        "h2q", "beta10", "beta20", "beta21", "se10", "se20", 
        "se21", "chisq", "df", "model", "pval")
    write.table(result, outfile, quote = F, row.names = F, col.names = T, 
        sep = ",", na = "", append = T)
}
