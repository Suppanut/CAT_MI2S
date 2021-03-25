# CAT MI2S

###########################
# R functions and objects #
###########################

# General Purpose

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

# Data generation


genData_check_write <- function(maindir, nCat, thDist, N, repno, seed = NULL, writedat = TRUE, missflag = "9", sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  proper <- FALSE
  i <- 0    
  while(isFALSE(proper)) { # only generate data with proper number of categories in each variable
    if(!is.null(seed)) seed <- seed + i
    i <- i + 1
    dat_all <- genData(nCat, thDist, N, repno, seed)
    
    filename_comp <- genPath(maindir, nCat, thDist, N, repno, 0, createdir = TRUE)
    filename_miss20 <- genPath(maindir, nCat, thDist, N, repno, 20, createdir = TRUE)
    filename_miss40 <- genPath(maindir, nCat, thDist, N, repno, 40, createdir = TRUE)
    
    check_comp <- checkDat(dat_all$comp, nCat, comp = TRUE)
    check_miss20 <- checkDat(dat_all$miss20, nCat, comp = FALSE)
    check_miss40 <- checkDat(dat_all$miss40, nCat, comp = FALSE)
    
    proper <- all(c(check_comp, check_miss20, check_miss40))
    if (isFALSE(proper)) cat("Regenerate datasets:", basename(filename_comp), "\n")
    if (i > 20) stop("too many iterations")
  }
  
  if (isTRUE(writedat)) {
    write.table(dat_all$comp, filename_comp, sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(dat_all$miss20, filename_miss20, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
    write.table(dat_all$miss40, filename_miss40, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
    output <- c(filename_comp, filename_miss20, filename_miss40, dat_all$seed, proper, i)
  } 
  output
}

checkDat <- function(dat, nCat, comp) {
  check_nCat <- sapply(lapply(dat, unique), length)[-1] # exclude repno column
  if (isTRUE(comp)) {
    proper <- (check_nCat - nCat) == 0
  } else {
    proper <- (check_nCat - c(rep(nCat, 12), rep(nCat + 1, 6))) == 0 # missisng in last 6 out of 18 items
  }
  all(proper)
}

genData <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
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
  # Z is based on continuous latent responses, otherwise the proprotion of missing data will be off
  Z <- rowMeans(contDat[,c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor))])
  quartiles <- c(-Inf, quantile(Z)[2:4], Inf)

  propMiss20 <- as.numeric(as.character(cut(Z, quartiles, right=FALSE, labels= c(.50, .20, .075, .025) )))
  propMiss40 <- as.numeric(as.character(cut(Z, quartiles, right=FALSE, labels= c(1, .40, .15, .05) )))
  
  indMiss20 <- propMiss20 > runif(nrow(catDat), min = 0, max = 1)
  indMiss40 <- propMiss40 > runif(nrow(catDat), min = 0, max = 1)

  # # MAR: Shi et al. (2020)
  # # Z is based on continuous latent responses, otherwise the proprotion of missing data will be off
  # # first 6 items
  # Z <- rowMeans(contDat[,paste0("X",1:nItemPerFactor)]) # based on continuous latent responses
  
  # propMiss20 <- quantile(Z, .20)
  # propMiss40 <- quantile(Z, .40)

  # indMiss20 <- Z < propMiss20
  # indMiss40 <- Z < propMiss40

  catDatMiss20 <- catDatMiss40 <- catDat
  
  catDatMiss20[indMiss20, paste0("Y",1:nItemPerFactor)] <- NA
  catDatMiss40[indMiss40, paste0("Y",1:nItemPerFactor)] <- NA

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}


anaComp <- function(maindir, nCat, thDist, N, repno, propMiss = 0, anaModel, est, sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  cat_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  # Model CM1: 3-factor CFA for X, M, Y
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
                   estimator        = est))
  
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

anaFIML <- function(maindir, nCat, thDist, N, repno, propMiss = 0, anaModel, est, sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  cont_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  # Model CM1: 3-factor CFA for X, M, Y
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
  
  model <- get(anaModel) # CM1, ICM1, BM

  # Import data
  filename <- genPath(maindir, nCat, thDist, N, repno, propMiss)
  data <- data.frame(read.table(filename, sep = ",", header = FALSE, na.strings="9"))
  colnames(data) <- c("repno", cont_items)
  
  # Analysis of complete data
  fit <- Catch(sem(model            = model,
                   data             = data,
                   std.lv           = TRUE,
                   estimator        = est,
                   missing          = "fiml"))

  err_and_warn <- fit[2:3]
  
  fit <- fit[[1]]
  param <- try(lavaan::coef(fit), silent = TRUE)
  fitstat <- try(fitMeasures(fit), silent = TRUE)
  nobs <- try(inspect(fit, "nobs"))

  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, 
             anaModel = anaModel, est = est)

  RES <- list(conds     = conds,
              param     = param, # should save std.all
              fitstat   = fitstat,
              nobs      = nobs,
              err       = err_and_warn$err,
              warn      = err_and_warn$warn)
  RES
}


inp_template <- 
"DATA:
file=
'MISSING_DATA_FILE';

VARIABLE:
names = 
repno
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
usevariables =
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
categorical =
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
AUXILIARY = repno;
missing = ALL(MISSFLAG);

DATA IMPUTATION:
impute = (c) Y1 Y2 Y3 Y4 Y5 Y6;
ndatasets = NIMP;
save = SAVE_IMPUTED_DATA;
FORMAT = F5.0;
thin = BURNIN;

ANALYSIS:
TYPE = BASIC;
BITERATIONS = 100000 (BURNIN);
CHAINS = 2;
!BSEED = 0;
BCONVERGENCE = .025;

OUTPUT: TECH8;"

inp_template_convdiag <- 
"DATA:
file=
'MISSING_DATA_FILE';

VARIABLE:
names = 
repno
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
usevariables =
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
categorical =
X1 X2 X3 X4 X5 X6
M1 M2 M3 M4 M5 M6
Y1 Y2 Y3 Y4 Y5 Y6;
missing = ALL(MISSFLAG);

ANALYSIS:
estimator = BAYES;
BITERATIONS = 100000 (BURNIN);
CHAINS = 2;
!BSEED = 0;
BCONVERGENCE = .025;

MODEL:
X1-Y6 WITH X1-Y6;

OUTPUT: 
TECH8;

PLOT: 
TYPE = PLOT2;"

generateMplusSyntax <- function(inp_template, FILE, SAVE, nBurn, nImp, missflag = 9, BSEED = NULL) {
  syntax <- gsub("MISSING_DATA_FILE", FILE, inp_template)
  syntax <- gsub("SAVE_IMPUTED_DATA", SAVE, syntax)
  syntax <- gsub("BURNIN", nBurn, syntax)
  syntax <- gsub("NIMP", nImp, syntax)
  syntax <- gsub("MISSFLAG", missflag, syntax)
  if (!is.null(BSEED)) {
    if (BSEED == "RANDOM") BSEED <- sample(1:99999999, 1)
    syntax <- gsub("!BSEED = 0", paste0("BSEED = ", BSEED), syntax)
  }
  cat(writeLines(syntax))
  syntax
}

CheckPSR <- function(mplusoutfile, CheckLastPSR = TRUE) {
  impOut <- mplusoutfile$tech8$psr
  impOut <- cbind(impOut[,1]/2+1,impOut) # add FROM column
  colnames(impOut) <- c("FROM", "TO", "PSR", "LABEL")
  
  # check Max PSR
  below105 <- impOut$PSR < 1.05
  below110 <- impOut$PSR < 1.10
  
  res <- cbind(iteration = impOut$FROM-1, impOut, below105, below110)
  
  # find the last number of iteration that is beyond the criteria 
  nIter105 <- res$iteration[suppressWarnings(max(which(!below105)+1))]
  nIter110 <- res$iteration[suppressWarnings(max(which(!below110)+1))]
  
  nIter <- c(`PSR<1.05` = nIter105, `PSR<1.10` = nIter110)
  nIter[apply(cbind(below105,below110),2 ,all)] <- 50
  
  lastPSR <- as.numeric(impOut[nrow(impOut),c("TO","PSR")])
  # save the results
  if(CheckLastPSR == FALSE) {
    list("PSR" = res, "iterations" = nIter)
  } else {
    list("PSR" = res, "iterations" = nIter, "lastPSR" = lastPSR)
  }
}

runMplusMI <- function(nCat = NULL, thDist = NULL, N = NULL, repno = NULL, propMiss = NULL,
                       nBurn = NULL,
                       nImp = NULL,
                       maindir = NULL,
                       execute = TRUE,
                       BSEED = NULL,
                       sourcedir = NULL,
                       convdiag = FALSE,
                       savemaindir = NULL) {
  # Setup 
  suppressMessages(library(MplusAutomation))
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions
  
  if (isFALSE(convdiag)) {
    imptemp <- inp_template
  } else if (isTRUE(convdiag)) {
    imptemp <- inp_template_convdiag # Estimator = BAYES; MODEL: X1-Y6 WITH X1-Y6; PLOT: TYPE = PLOT2; 
  }

  # Replace text in Mplus imptemp
  nBurn <- nBurn # BURN and THIN iterations
  nImp <- nImp # m imputed datasets
  FILE <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss) # location of the .csv missing data file 
  SAVE <- paste0(sub("\\.csv$", "_imp*.dat", basename(FILE))) # save imputed data as ..._imp*.dat
  FILE <- paste0(dirname(FILE),"\n/",basename(FILE)) # seperate to multiple lines
  syntax <- generateMplusSyntax(imptemp, FILE, SAVE, as.character(nBurn), nImp = nImp, BSEED = BSEED) # generate Mplus syntax from inp_template  

  # Create Mplus .inp file
  if (!is.null(savemaindir)) maindir <- savemaindir
  inp_full_path <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss, MI = TRUE, createdir = TRUE) # absolute path 
  inp_full_path <- gsub(".dat",".inp",inp_full_path) # change .dat to .inp
  writeLines(syntax, inp_full_path) # write .inp file

  # Run Mplus via MplusAutomation (batch mode) 
  MplusAutomation::runModels(inp_full_path)

  # Check PSR
  out_full_path <- sub("\\.inp", "\\.out",inp_full_path) # absolute path to .out 
  tech8 <- MplusAutomation::readModels(out_full_path, what = "tech8") # extract tech8 from the .out 
  PSR <- CheckPSR(tech8)

  # R output
  c("inp" = inp_full_path, PSR)
}


Pool_ChungCai2019 <- function(H0.list, imp = NULL) {
  # Chung & Cai (2019, Eqs. 19, 21, 22)
  # Note: lavaan requires threshols
  
  if (!is.null(imp)) { 
    H0.list <- H0.list[imp]
  }

  m <- length(H0.list) # number of imputed datasets
  # Average polychoric correlation matrix
  pcorr.list <- lapply(H0.list, function(x) lavInspect(x, "cov.ov"))
  pcorr.avg <- Reduce('+', pcorr.list)/length(pcorr.list)
  pcorr.avg.vec <- lavaan::lav_matrix_vech(pcorr.avg, diagonal = FALSE)

  # Average threshold matrix
  threshold.list <- lapply(H0.list, function(x) lavInspect(x, "th"))
  threshold.avg <- Reduce('+', threshold.list)/length(threshold.list)
  threshold.avg <- as.numeric(threshold.avg) 
  attr(threshold.avg, "th.idx") <- lavInspect(H0.list[[1]], "th.idx") # add attribute to thresholds

  samstat.avg <- c(threshold.avg, pcorr.avg.vec) # Average thresholds and non-redundant correlations

  # Between-imputation variance of the thresholds and the polychoric correlation matrix
  samstat.list <- mapply(FUN = function(x,y) {
    c(x,lavaan::lav_matrix_vech(y, diagonal = FALSE))
  }, threshold.list, pcorr.list, SIMPLIFY=FALSE)
  acov.between.list <- lapply(samstat.list, function(x) (x-samstat.avg) %*% t(x-samstat.avg))
  acov.between <- Reduce('+', acov.between.list)/(m-1)

  # Within-imputation asymptotic covariance matrix of the thresholds and the polychoric correlation matrix
  acov.list <- lapply(H0.list, vcov)
  acov.within <- Reduce('+', acov.list)/length(acov.list)
  # CHECK: sort - thresholds come first
  thr.where <- grepl("\\|",colnames(acov.within)) 
  sort.acov.within <- c(which(thr.where), which(!thr.where))
  acov.within <- acov.within[sort.acov.within, sort.acov.within]

  # Corrected variance-covariance matrix of the thresholds and the polychoric correlation matrix 
  acov.total <- acov.within + acov.between + acov.between/m
  
  # lavaan uses nacov as input
  N <- H0.list[[1]]@SampleStats@nobs[[1]]
  nacov.total <- (N-1)*acov.total

  # OUTPUT
  OUTPUT <- list(nacov.total    = nacov.total,
                 pcorr.avg      = pcorr.avg,
                 threshold.avg  = threshold.avg,
                 m              = m,
                 N              = N)
  OUTPUT
}

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