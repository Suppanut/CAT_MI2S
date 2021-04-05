# Data generation
genData_check_write <- function(maindir, nCat, thDist, N, repno, seed = NULL, writedat = TRUE, missflag = "9", 
                                sourcedir = NULL, genMethod = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  set.seed(NULL)
  proper <- FALSE
  i <- 0    
  while(isFALSE(proper)) { # only generate data with proper number of categories in each variable
    if(!is.null(seed)) seed <- seed + i
    i <- i + 1
    
    genData <- get(paste0("genData_Method", genMethod))
    missvar <- get(paste0("missvar", genMethod))

    dat_all <- genData(nCat, thDist, N, repno, seed)
    
    filename_comp <- genPath(maindir, nCat, thDist, N, repno, 0, createdir = TRUE)
    filename_miss20 <- genPath(maindir, nCat, thDist, N, repno, 20, createdir = TRUE)
    filename_miss40 <- genPath(maindir, nCat, thDist, N, repno, 40, createdir = TRUE)
    
    check_comp <- checkDat(dat_all$comp, nCat, comp = TRUE)
    check_miss20 <- checkDat(dat_all$miss20, nCat, comp = FALSE, missvar = missvar)
    check_miss40 <- checkDat(dat_all$miss40, nCat, comp = FALSE, missvar = missvar)
    
    proper <- all(c(check_comp$proper, check_miss20$proper, check_miss40$proper))
    if (isFALSE(proper)) cat("Regenerate datasets:", basename(filename_comp), "\n")
    if (i > 10) stop("too many iterations")
  }
  
  if (isTRUE(writedat)) {
    write.table(dat_all$comp, filename_comp, sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(dat_all$miss20, filename_miss20, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
    write.table(dat_all$miss40, filename_miss40, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
  } 

  output <- list(comp = filename_comp, 
                 miss20 = filename_miss20, 
                 miss40 = filename_miss40, 
                 seed = dat_all$seed, 
                 proper = proper, 
                 missrate20 = check_miss20$missrate, 
                 missrate40 = check_miss40$missrate, 
                 iteration = i)
  output
}

missvar1 <- c(0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,1,1,1) # missing in last 6 (Y) out of 18 items
missvar2 <- c(0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,1,1,1) # missing in last 6 (Y) out of 18 items
missvar3 <- c(0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,1,1,1) # missing in last 6 (Y) out of 18 items
missvar4 <- c(0,0,0,0,0,0, 1,1,1,1,1,1, 1,1,1,1,1,1) # missing in last 6 (Y) out of 18 items
missvar5 <- c(0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1) # missing in last 3 of each factor

checkDat <- function(dat = NULL, nCat = NULL, comp = NULL, missvar = NULL) {
  check_nCat <- sapply(lapply(dat, unique), length)[-1] # exclude repno column
  if (isTRUE(comp)) {
    proper <- all((check_nCat - nCat) == 0)
  } else {
    # missvar <- c(0,0,0,0,0,0, 0,0,0,0,0,0, 1,1,1,1,1,1) # missing in last 6 (Y) out of 18 items
    # missvar <- c(0,0,0,0,0,0, 1,1,1,1,1,1, 1,1,1,1,1,1) # missing MY
    # missvar <- c(0,0,0,1,1,1, 0,0,0,1,1,1, 0,0,0,1,1,1) # missing in last 3 of each factor
    proper <- all((check_nCat - rep(nCat, 18) - missvar) == 0)
  }
  missrate <- sapply(lapply(dat, is.na), mean)[-1]
  list(proper = proper, missrate = missrate)
}

genPath <- function(maindir, nCat, thDist, N, repno, propMiss, MI = FALSE, check = FALSE, createdir = FALSE) {
  if (is.numeric(nCat)) nCat <- paste0("C", nCat)
  if (is.numeric(N)) N <- paste0("N", N)
  between <- c(nCat, thDist, N)

  if (is.numeric(propMiss)) propMiss <- paste0("miss", propMiss)

  folname <- file.path(maindir, paste(between, collapse = "_"))
  filename <- paste(c(between, propMiss, repno), collapse = "_")

  if (isFALSE(MI)) { 
    folname <- file.path(folname, propMiss)
    filename <- paste0(filename, ".csv")
  } else if (isTRUE(MI)) {
    folname <- file.path(folname, sub("miss" , "MI", propMiss), paste0("rep", repno))
    filename <- paste0(filename, "_imp", repno, ".dat")
  }

  abspath <- file.path(folname, filename)
  
  if (isTRUE(createdir) & !dir.exists(folname)) dir.create(folname, recursive=TRUE)

  if (isTRUE(check)) {
    cat("dir.exists(folname):", dir.exists(folname), "\n")
    cat("file.exists(abspath):", file.exists(abspath), "\n")
  }

  # list (abspath = abspath, folname = folname, filename = filename)
  abspath
}

# Catch() is used to save errors and warnings
# Modified from:
# https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
Catch <- function( ... ) {
  warn <- err <- NULL
  # withCallingHandlers() is used to handle warnings
  # tryCatch() is used to handle errors
  res <- withCallingHandlers(tryCatch( ... , 
                                       error=function(e) {
                                         err <<- conditionMessage(e)
                                         NULL 
                                       }), 
                             warning=function(w) {
                               warn <<- append(warn, conditionMessage(w))
                               invokeRestart("muffleWarning") # the warning is not printed
                             })
  list(res=res, err=err, warn=warn)
}



anaComp <- function(maindir, nCat, thDist, N, repno, propMiss = 0, anaModel, est, sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  cat_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))
  
  model <- get(anaModel) # CM1, ICM1, BM

  # Import data
  filename <- genPath(maindir, nCat, thDist, N, repno, propMiss)
  data <- data.frame(read.table(filename, sep = ",", header = FALSE, na.strings="9"))
  colnames(data) <- c("repno", cat_items)
  
  # Analysis of complete data
  fit <- Catch(sem(model            = model,
                   data             = data,
                   ordered          = cat_items,
                   std.lv           = TRUE,
                   parameterization = "delta",
                   estimator        = est,
                   missing          = "pairwise"))
  
  err_and_warn <- fit[2:3]
  
  fit <- fit[[1]]
  param <- try(lavaan::coef(fit), silent = TRUE)
  fitstat <- try(fitMeasures(fit), silent = TRUE)
  teststat <- try(calculateT(fit), silent = TRUE)
  nobs <- try(inspect(fit, "nobs"))

  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, 
             anaModel = anaModel, est = est)

  RES <- list(conds     = conds,
              param     = param, 
              fitstat   = fitstat,
              teststat  = teststat,
              nobs      = nobs,
              err       = err_and_warn$err,
              warn      = err_and_warn$warn)
  RES
}

# Model CM1: 3-factor CFA for X, M, Y with 6 items per factor
CM1 <- '
  X =~ NA*X1 + X2 + X3 + X4 + X5 + X6
  M =~ NA*M1 + M2 + M3 + M4 + M5 + M6
  Y =~ NA*Y1 + Y2 + Y3 + Y4 + Y5 + Y6
  X ~~ 1*X + M + Y
  M ~~ 1*M + Y
  Y ~~ 1*Y'

# Model ICM1: X --> M --> Y
ICM1 <- '
  X =~ NA*X1 + X2 + X3 + X4 + X5 + X6
  M =~ NA*M1 + M2 + M3 + M4 + M5 + M6
  Y =~ NA*Y1 + Y2 + Y3 + Y4 + Y5 + Y6
  Y ~ M
  M ~ X
  Y ~ 0*X'

# Model baseline model: thresholds only
BM <- '
  X1~~0*X1; X2~~0*X2; X3~~0*X3; X4~~0*X4; X5~~0*X5; X6~~0*X6 
  M1~~0*M1; M2~~0*M2; M3~~0*M3; M4~~0*M4; M5~~0*M5; M6~~0*M6 
  Y1~~0*Y1; Y2~~0*Y2; Y3~~0*Y3; Y4~~0*Y4; Y5~~0*Y5; Y6~~0*Y6'


# Model CM2: 3-factor CFA for X, M, Y with 5 items per factor
CM2 <- '
  X =~ NA*X1 + X2 + X3 + X4 + X5
  M =~ NA*M1 + M2 + M3 + M4 + M5
  Y =~ NA*Y1 + Y2 + Y3 + Y4 + Y5
  X ~~ 1*X + M + Y
  M ~~ 1*M + Y
  Y ~~ 1*Y'

# Model CM3: 3-factor CFA for X, M, Y with 3 items per factor
CM3 <- '
  X =~ NA*X1 + X2 + X3 
  M =~ NA*M1 + M2 + M3 
  Y =~ NA*Y1 + Y2 + Y3 
  X ~~ 1*X + M + Y
  M ~~ 1*M + Y
  Y ~~ 1*Y'

calculateT <- function(fit) {
  # Browne (1984)
  # Or see eq. 3.1-3.2 in Chun, Browne, & Shapiro (2018)
  N <- lavaan::inspect(fit, "nobs")
  df <- fit@Fit@test$standard$df
  G <- lavaan::lavTech(fit, "gamma")[[1]]
  D <- lavaan:::computeDelta(lavmodel = fit@Model)[[1]]
  e <- matrix(c(residuals(fit)$th, lav_matrix_vech(residuals(fit)$cov, diag = FALSE)), ncol = 1)
  G.inv <- try(solve(G), silent = TRUE)
  if ("try-error" %in% class(G.inv)) G.inv <- MASS::ginv(G)
  U <- try(G.inv - (G.inv %*% D %*% solve(t(D) %*% G.inv %*% D) %*% t(D) %*% G.inv), silent = TRUE)
  if ("try-error" %in% class(U)) U <- G.inv - (G.inv %*% D %*% MASS::ginv(t(D) %*% G.inv %*% D) %*% t(D) %*% G.inv)
  # U <- G.inv - (G.inv %*% D %*% MASS::ginv(t(D) %*% G.inv %*% D) %*% t(D) %*% G.inv)
  TB <- as.numeric(N * (t(e) %*% U %*% e))
  TYB <- TB/(1+N*TB/(N-1)^2) # (see Yuan & Bentler, 1998, eq. 7)

  # See Savalei & Rhemtulla (2013) & SB1994 & AM2010
  T <- fit@Fit@test$standard$stat
  df <- fit@Fit@test$standard$df
  UG <- lavInspect(fit, "UGamma")
  UG.tr <- sum(diag(UG))
  # SB1994
  c <- df/UG.tr
  TM <- c * T
  # AM2010
  a <- sqrt(df/sum(diag(UG %*% UG)))
  b <- df - a * UG.tr
  TMV <- a*T+b
  output <- c("T" = T, "df" = df, 
            "TM" = TM, "pvalue.TM" = 1-pchisq(TM,df),
            "TMV" = TMV, "pvalue.TMV" = 1-pchisq(TMV,df),
            "TB" = TB, "pvalue.TB" = 1-pchisq(TB,df),
            "TYB" = TYB, "pvalue.TYB" = 1-pchisq(TYB,df))
  return(output)
}

rmsea <- function(T, df, N) {
  rmsea <- sqrt(max(c(T-df, 0))/(df*(N-1)))
  c(rmsea = rmsea)
}

cfitli <- function(T, df, TI, dfI) {
  cfi <- 1 - max(c(T - df, 0))/max(c(TI - dfI, 0))
  tli <- 1 - ((T - df) * dfI) / ((TI - dfI) * df)
  c(cfi = cfi, tli = tli)
}

### 1
genData_Method1 <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
  if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
  set.seed(seedused)
  # Conditions
  nCat <- nCat # Number of categories: C2, C5
  thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
  N <- N # Sample size of the observed data: N = 250, 500, 1000

  # Factor loadings and factor correlations
  factorLoading <- 0.8
  factorCorr <- 0.4

  # Generate observed data on the indicators 

  # Population loading matrix
  nItemPerFactor <- 6

  lambda <- matrix(0, nrow = nItemPerFactor*3, ncol=3)
  lambda[1:nItemPerFactor,1] <- factorLoading
  lambda[(nItemPerFactor+1):(nItemPerFactor+nItemPerFactor),2] <- factorLoading
  lambda[(nItemPerFactor*2+1):(nItemPerFactor*2+nItemPerFactor),3] <- factorLoading

  # Population factor covariance matrix
  psi <- matrix(c(1,factorCorr,factorCorr,
                  factorCorr,1,factorCorr,
                  factorCorr,factorCorr,1),
                nrow=3,ncol=3,byrow = TRUE)

  # Population covariance matrix of the latent responses of the ordinal items
  Sigma <- lambda %*% psi %*% t(lambda)
  diag(Sigma) <- 1
    
  # Generate data on the latent responses
  contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*nItemPerFactor), Sigma = Sigma)
  contDat <- as.data.frame(contDat)
  colnames(contDat) <- c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor), paste0("Y",1:nItemPerFactor))

  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
    }
  }
  catDat <- contDat
  for (i in 1:(3*nItemPerFactor)){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }
  catDat <- cbind(repno = repno, catDat)

  # MAR: Chung & Cai (2019)
  Z <- rowMeans(catDat[,c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor))])
  quartiles <- c(-Inf, quantile(Z)[2:4], Inf)

  propMiss20 <- as.numeric(as.character(cut(Z, quartiles, right=FALSE, labels= c(.50, .20, .075, .025) )))
  propMiss40 <- as.numeric(as.character(cut(Z, quartiles, right=FALSE, labels= c(1, .40, .15, .05) )))
  
  indMiss20 <- propMiss20 > runif(nrow(catDat), min = 0, max = 1)
  indMiss40 <- propMiss40 > runif(nrow(catDat), min = 0, max = 1)

  catDatMiss20 <- catDatMiss40 <- catDat
  
  catDatMiss20[indMiss20, paste0("Y",1:nItemPerFactor)] <- NA
  catDatMiss40[indMiss40, paste0("Y",1:nItemPerFactor)] <- NA

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}

### 2
genData_Method2 <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
  if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
  set.seed(seedused)
  # Conditions
  nCat <- nCat # Number of categories: C2, C5
  thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
  N <- N # Sample size of the observed data: N = 250, 500, 1000

  # Factor loadings and factor correlations
  factorLoading <- 0.8
  factorCorr <- 0.4

  # Generate observed data on the indicators 

  # Population loading matrix
  nItemPerFactor <- 6

  lambda <- matrix(0, nrow = nItemPerFactor*3, ncol=3)
  lambda[1:nItemPerFactor,1] <- factorLoading
  lambda[(nItemPerFactor+1):(nItemPerFactor+nItemPerFactor),2] <- factorLoading
  lambda[(nItemPerFactor*2+1):(nItemPerFactor*2+nItemPerFactor),3] <- factorLoading

  # Population factor covariance matrix
  psi <- matrix(c(1,factorCorr,factorCorr,
                  factorCorr,1,factorCorr,
                  factorCorr,factorCorr,1),
                nrow=3,ncol=3,byrow = TRUE)

  # Population covariance matrix of the latent responses of the ordinal items
  Sigma <- lambda %*% psi %*% t(lambda)
  diag(Sigma) <- 1
    
  # Generate data on the latent responses
  contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*nItemPerFactor), Sigma = Sigma)
  contDat <- as.data.frame(contDat)
  colnames(contDat) <- c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor), paste0("Y",1:nItemPerFactor))

  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
    }
  }
  catDat <- contDat
  for (i in 1:(3*nItemPerFactor)){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }
  catDat <- cbind(repno = repno, catDat)
  
  # MAR
  Z <- rowSums(catDat[,c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor))])
  
  propMiss20 <- sort(rep(c(.50, .20, .075, .025), length.out = N), decreasing = TRUE)[order(order(Z))]
  propMiss40 <- sort(rep(c(1, .40, .15, .05), length.out = N), decreasing = TRUE)[order(order(Z))]

  indMiss20 <- propMiss20 > runif(nrow(catDat), min = 0, max = 1)
  indMiss40 <- propMiss40 > runif(nrow(catDat), min = 0, max = 1)

  catDatMiss20 <- catDatMiss40 <- catDat
  
  catDatMiss20[indMiss20, paste0("Y",1:nItemPerFactor)] <- NA
  catDatMiss40[indMiss40, paste0("Y",1:nItemPerFactor)] <- NA

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}

### 3
genData_Method3 <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
  if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
  set.seed(seedused)
  # Conditions
  nCat <- nCat # Number of categories: C2, C5
  thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
  N <- N # Sample size of the observed data: N = 250, 500, 1000

  # Factor loadings and factor correlations
  factorLoading <- 0.8
  factorCorr <- 0.4

  # Generate observed data on the indicators 

  # Population loading matrix

  lambda <- matrix(0, ncol=3, nrow = 6*3)
  lambda[1:6,1] <- factorLoading
  lambda[(6+1):(6+6),2] <- factorLoading
  lambda[(6*2+1):(6*2+6),3] <- factorLoading

  # Population factor covariance matrix
  psi <- matrix(c(1,factorCorr,factorCorr,
                  factorCorr,1,factorCorr,
                  factorCorr,factorCorr,1),
                nrow=3,ncol=3,byrow = TRUE)

  # Population covariance matrix of the latent responses of the ordinal items
  Sigma <- lambda %*% psi %*% t(lambda)
  diag(Sigma) <- 1
    
  # Generate data on the latent responses
  contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*6), Sigma = Sigma)
  contDat <- as.data.frame(contDat)
  colnames(contDat) <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
      propMiss20 <- c(0, .4)
      propMiss40 <- c(0, .8)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
      propMiss20 <- c((20-15)/85, 1)
      propMiss40 <- c((40-15)/85, 1)
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
      propMiss20 <- c(0, 0, 0, (20-7)/24, 1)
      propMiss40 <- c(0, 0, (40-7-24)/38, 1, 1)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
      propMiss20 <- c(0, 0, 0, 1, 1)
      propMiss40 <- c(0, 0.4, 1, 1, 1)
    }
  }
  
  catDat <- contDat
  for (i in 1:(3*6)){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }

  Miss20 <- Miss40 <- matrix(0, nrow = N, ncol = ncol(catDat))
  colnames(Miss20) <- colnames(Miss40) <- colnames(catDat)
  # Items X1-X6 predict missingness of Y1-Y6
  for (i in 1:6){
    Miss20[,i+12] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss20)))
    Miss40[,i+12] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss40)))
  }
  IndMiss20 <- Miss20 > runif(max = 1, min = 0, n = ncol(Miss20)*N)
  IndMiss40 <- Miss40 > runif(max = 1, min = 0, n = ncol(Miss40)*N)

  catDatMiss20 <- catDatMiss40 <- catDat
  catDatMiss20[IndMiss20] <- NA
  catDatMiss40[IndMiss40] <- NA
  
  # Add repno
  catDat <- cbind(repno = repno, catDat)
  catDatMiss20 <- cbind(repno = repno, catDatMiss20)
  catDatMiss40 <- cbind(repno = repno, catDatMiss40)

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}

### 4
genData_Method4 <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
  if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
  set.seed(seedused)
  # Conditions
  nCat <- nCat # Number of categories: C2, C5
  thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
  N <- N # Sample size of the observed data: N = 250, 500, 1000

  # Factor loadings and factor correlations
  factorLoading <- 0.8
  factorCorr <- 0.4

  # Generate observed data on the indicators 

  # Population loading matrix

  lambda <- matrix(0, ncol=3, nrow = 6*3)
  lambda[1:6,1] <- factorLoading
  lambda[(6+1):(6+6),2] <- factorLoading
  lambda[(6*2+1):(6*2+6),3] <- factorLoading

  # Population factor covariance matrix
  psi <- matrix(c(1,factorCorr,factorCorr,
                  factorCorr,1,factorCorr,
                  factorCorr,factorCorr,1),
                nrow=3,ncol=3,byrow = TRUE)

  # Population covariance matrix of the latent responses of the ordinal items
  Sigma <- lambda %*% psi %*% t(lambda)
  diag(Sigma) <- 1
    
  # Generate data on the latent responses
  contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*6), Sigma = Sigma)
  contDat <- as.data.frame(contDat)
  colnames(contDat) <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
      propMiss20 <- c(0, .4)
      propMiss40 <- c(0, .8)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
      propMiss20 <- c((20-15)/85, 1)
      propMiss40 <- c((40-15)/85, 1)
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
      propMiss20 <- c(0, 0, 0, (20-7)/24, 1)
      propMiss40 <- c(0, 0, (40-7-24)/38, 1, 1)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
      propMiss20 <- c(0, 0, 0, 1, 1)
      propMiss40 <- c(0, 0.4, 1, 1, 1)
    }
  }
  
  catDat <- contDat
  for (i in 1:(3*6)){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }

  Miss20 <- Miss40 <- matrix(0, nrow = N, ncol = ncol(catDat))
  colnames(Miss20) <- colnames(Miss40) <- colnames(catDat)
  # Items X1-X6 predict missingness of M1-M6 and Y1-Y6
  for (i in 1:6){
    Miss20[,i+6] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss20)))
    Miss20[,i+12] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss20)))
    Miss40[,i+6] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss40)))
    Miss40[,i+12] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss40)))
  }
  IndMiss20 <- Miss20 > runif(max = 1, min = 0, n = ncol(Miss20)*N)
  IndMiss40 <- Miss40 > runif(max = 1, min = 0, n = ncol(Miss40)*N)

  catDatMiss20 <- catDatMiss40 <- catDat
  catDatMiss20[IndMiss20] <- NA
  catDatMiss40[IndMiss40] <- NA
  
  # Add repno
  catDat <- cbind(repno = repno, catDat)
  catDatMiss20 <- cbind(repno = repno, catDatMiss20)
  catDatMiss40 <- cbind(repno = repno, catDatMiss40)

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}

### 5
genData_Method5 <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
  if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
  set.seed(seedused)
  # Conditions
  nCat <- nCat # Number of categories: C2, C5
  thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
  N <- N # Sample size of the observed data: N = 250, 500, 1000

  # Factor loadings and factor correlations
  factorLoading <- 0.8
  factorCorr <- 0.4

  # Generate observed data on the indicators 

  # Population loading matrix

  lambda <- matrix(0, ncol=3, nrow = 6*3)
  lambda[1:6,1] <- factorLoading
  lambda[(6+1):(6+6),2] <- factorLoading
  lambda[(6*2+1):(6*2+6),3] <- factorLoading

  # Population factor covariance matrix
  psi <- matrix(c(1,factorCorr,factorCorr,
                  factorCorr,1,factorCorr,
                  factorCorr,factorCorr,1),
                nrow=3,ncol=3,byrow = TRUE)

  # Population covariance matrix of the latent responses of the ordinal items
  Sigma <- lambda %*% psi %*% t(lambda)
  diag(Sigma) <- 1
    
  # Generate data on the latent responses
  contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*6), Sigma = Sigma)
  contDat <- as.data.frame(contDat)
  colnames(contDat) <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
      propMiss20 <- c(0, .4)
      propMiss40 <- c(0, .8)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
      propMiss20 <- c((20-15)/85, 1)
      propMiss40 <- c((40-15)/85, 1)
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
      propMiss20 <- c(0, 0, 0, (20-7)/24, 1)
      propMiss40 <- c(0, 0, (40-7-24)/38, 1, 1)
    } else if(thDist == "asym") {
      Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
      propMiss20 <- c(0, 0, 0, 1, 1)
      propMiss40 <- c(0, 0.4, 1, 1, 1)
    }
  }
  
  catDat <- contDat
  for (i in 1:(3*6)){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }

  Miss20 <- Miss40 <- matrix(0, nrow = N, ncol = ncol(catDat))
  colnames(Miss20) <- colnames(Miss40) <- colnames(catDat)
  # Items 1-3 predict missingness of items 4-6 for each factor
  for (i in c(1:3, 7:9, 13:15)){
    Miss20[,i+3] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss20)))
    Miss40[,i+3] <- as.numeric(as.character(factor(catDat[,i], levels = 0:(nCat-1), labels = propMiss40)))
  }
  IndMiss20 <- Miss20 > runif(max = 1, min = 0, n = ncol(Miss20)*N)
  IndMiss40 <- Miss40 > runif(max = 1, min = 0, n = ncol(Miss40)*N)

  catDatMiss20 <- catDatMiss40 <- catDat
  catDatMiss20[IndMiss20] <- NA
  catDatMiss40[IndMiss40] <- NA
  
  # Add repno
  catDat <- cbind(repno = repno, catDat)
  catDatMiss20 <- cbind(repno = repno, catDatMiss20)
  catDatMiss40 <- cbind(repno = repno, catDatMiss40)

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}