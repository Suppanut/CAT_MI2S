# CAT MI2S

###########################
# R functions and objects #
###########################

# General Purpose
sortParam <- c("X=~X1","X=~X2","X=~X3","X=~X4","X=~X5","X=~X6",
"M=~M1","M=~M2","M=~M3","M=~M4","M=~M5","M=~M6",
"Y=~Y1","Y=~Y2","Y=~Y3","Y=~Y4","Y=~Y5","Y=~Y6",
"X~~M","X~~Y","M~~Y", #"Y~M","M~X",
"X1|t1","X1|t2","X1|t3","X1|t4","X2|t1","X2|t2","X2|t3","X2|t4",
"X3|t1","X3|t2","X3|t3","X3|t4","X4|t1","X4|t2","X4|t3","X4|t4",
"X5|t1","X5|t2","X5|t3","X5|t4","X6|t1","X6|t2","X6|t3","X6|t4",
"M1|t1","M1|t2","M1|t3","M1|t4","M2|t1","M2|t2","M2|t3","M2|t4",
"M3|t1","M3|t2","M3|t3","M3|t4","M4|t1","M4|t2","M4|t3","M4|t4",
"M5|t1","M5|t2","M5|t3","M5|t4","M6|t1","M6|t2","M6|t3","M6|t4",
"Y1|t1","Y1|t2","Y1|t3","Y1|t4","Y2|t1","Y2|t2","Y2|t3","Y2|t4",
"Y3|t1","Y3|t2","Y3|t3","Y3|t4","Y4|t1","Y4|t2","Y4|t3","Y4|t4",
"Y5|t1","Y5|t2","Y5|t3","Y5|t4","Y6|t1","Y6|t2","Y6|t3","Y6|t4")

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
  set.seed(NULL)
  proper <- FALSE
  i <- 0    
  while(isFALSE(proper)) { # only generate data with proper number of categories
    if(!is.null(seed)) seed <- seed + i
    i <- i + 1
    # Generate complete and incomplete datasets
    dat_all <- genData(nCat, thDist, N, repno, seed)
    
    # Generate absolute filepaths to datasets
    filename_comp <- genPath(maindir, nCat, thDist, N, repno, 0, createdir = TRUE)
    filename_miss20 <- genPath(maindir, nCat, thDist, N, repno, 20, createdir = TRUE)
    filename_miss40 <- genPath(maindir, nCat, thDist, N, repno, 40, createdir = TRUE)
    
    # Check whether the generated datasets are approiate
    check_comp <- checkDat(dat_all$comp, nCat, comp = TRUE)
    check_miss20 <- checkDat(dat_all$miss20, nCat, comp = FALSE)
    check_miss40 <- checkDat(dat_all$miss40, nCat, comp = FALSE)
    
    proper <- all(c(check_comp$proper, check_miss20$proper, check_miss40$proper)) # all must be proper
    if (isFALSE(proper)) cat("Regenerate datasets:", basename(filename_comp), "\n") # if not proper, regenerate the data
    if (i > 5) stop("too many iterations")
  }
  
  # Write the data on disk
  if (isTRUE(writedat)) {
    write.table(dat_all$comp, filename_comp, sep = ",", row.names = FALSE, col.names = FALSE)
    write.table(dat_all$miss20, filename_miss20, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
    write.table(dat_all$miss40, filename_miss40, sep = ",", row.names = FALSE, col.names = FALSE, na = missflag)
  } 

  # Output
  output <- list(comp = filename_comp, 
                 miss20 = filename_miss20, 
                 miss40 = filename_miss40, 
                 seed = dat_all$seed, 
                 proper = proper, 
                 missrate20 = check_miss20$missrate, 
                 missrate40 = check_miss40$missrate, 
                 iteration = i)
  return(output)
}

# MAR: Shi et al. (2020)
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

  # Specify threshold values and proportion of missing data 
  if(nCat == 2) {
    if(thDist == "sym") {
      p <- c(0, 0.5, 1) # diff(p) = .5, .5
      Th <- qnorm(p, mean = 0, sd = 1) 
    } else if(thDist == "asym") {
      p <- c(0, 0.85, 1) # diff(p) = .85, .15
      Th <- qnorm(p, mean = 0, sd = 1) 
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      p <- c(0, 0.07, 0.31, 0.69, 0.93, 1) # diff(p) = .07 .24 .38 .24 .07
      Th <- qnorm(p, mean = 0, sd = 1)
    } else if(thDist == "asym") {
      p <- c(0, 0.52, 0.67, 0.80, 0.91, 1) # diff(p) = .52 .15 .13 .11 .09
      Th <- qnorm(p, mean = 0, sd = 1) 
    }
  }
  
  # Generate categorical data
  catDat <- contDat
  for (i in 1:18){
    catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
  }
  
  # Generate missing data (Shi et al., 2020; EPM)
  Z <- rowSums(catDat[,c(paste0("X",1:6), paste0("M",1:6))])
  Z.rank <- rank(Z, ties.method = "first") # smaller Z values correspond to smaller Z.rank

  # Impose NA
  catDatMiss20 <- catDatMiss40 <- catDat
  catDatMiss20[Z.rank <= N*.2, paste0("Y",1:6)] <- NA
  catDatMiss40[Z.rank <= N*.4, paste0("Y",1:6)] <- NA
  
  # Add repno
  catDat <- cbind(repno = repno, catDat)
  catDatMiss20 <- cbind(repno = repno, catDatMiss20)
  catDatMiss40 <- cbind(repno = repno, catDatMiss40)

  list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
}

# # MAR: Modified Chung & Cai (2019)
# genData <- function (nCat, thDist, N, repno = 0, seed = NULL) { 
#   if (is.null(seed)) seedused <- sample(1000:.Machine$integer.max, 1) else seedused <- seed
#   set.seed(seedused)
#   # Conditions
#   nCat <- nCat # Number of categories: C2, C5
#   thDist <- thDist # Threshold distributions: symmetric (sym), asymmetric (asym)
#   N <- N # Sample size of the observed data: N = 250, 500, 1000

#   # Factor loadings and factor correlations
#   factorLoading <- 0.8
#   factorCorr <- 0.4

#   # Generate observed data on the indicators 

#   # Population loading matrix
#   nItemPerFactor <- 6

#   lambda <- matrix(0, nrow = nItemPerFactor*3, ncol=3)
#   lambda[1:nItemPerFactor,1] <- factorLoading
#   lambda[(nItemPerFactor+1):(nItemPerFactor+nItemPerFactor),2] <- factorLoading
#   lambda[(nItemPerFactor*2+1):(nItemPerFactor*2+nItemPerFactor),3] <- factorLoading

#   # Population factor covariance matrix
#   psi <- matrix(c(1,factorCorr,factorCorr,
#                   factorCorr,1,factorCorr,
#                   factorCorr,factorCorr,1),
#                 nrow=3,ncol=3,byrow = TRUE)

#   # Population covariance matrix of the latent responses of the ordinal items
#   Sigma <- lambda %*% psi %*% t(lambda)
#   diag(Sigma) <- 1
    
#   # Generate data on the latent responses
#   contDat <- MASS::mvrnorm(n=N, mu = rep(0, 3*nItemPerFactor), Sigma = Sigma)
#   contDat <- as.data.frame(contDat)
#   colnames(contDat) <- c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor), paste0("Y",1:nItemPerFactor))

#   if(nCat == 2) {
#     if(thDist == "sym") {
#       Th <- qnorm(c(0, 0.5, 1), mean = 0, sd = 1) # 50%,50%
#     } else if(thDist == "asym") {
#       Th <- qnorm(c(0, 0.85, 1), mean = 0, sd = 1) # 85%,15%
#     }
#   } else if(nCat == 5) {
#     if(thDist == "sym") {
#       Th <- qnorm(c(0, 0.07, 0.31, 0.69, 0.93, 1), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
#     } else if(thDist == "asym") {
#       Th <- qnorm(c(0, 0.52, 0.67, 0.80, 0.91, 1), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
#     }
#   }
#   catDat <- contDat
#   for (i in 1:(3*nItemPerFactor)){
#     catDat[,i] <- as.numeric(as.character(cut(contDat[,i], Th, right=FALSE, labels=c(0:(nCat-1)))))
#   }

#   # MAR: Similar to Chung & Cai (2019), but use proportion rather than percentile
#   # because percentile will not lead to 4 levels with equal size
  
#   # Proportion of missing data
#   propMiss20 <- c(.50, .20, .075, .025)
#   propMiss40 <- propMiss20*2

#   # Calculate the sum scores Z 
#   Z <- rowSums(catDat[,c(paste0("X",1:nItemPerFactor), paste0("M",1:nItemPerFactor))])
#   Z.rank <- rank(Z, ties.method="first") # smaller Z values correspond to smaller Z.rank

#   # Assign each propMiss to 25% of cases (or about 25% when N = 250, i.e., [62, 62, 63, 63])
#   # Lower sum scores correspond to more missing
#   indPropMiss20 <- sort(rep(rev(propMiss20), length.out = N), decreasing = TRUE)[Z.rank]
#   indPropMiss40 <- sort(rep(rev(propMiss40), length.out = N), decreasing = TRUE)[Z.rank]
  
#   # Cases with missing
#   indMiss20 <- indPropMiss20 > runif(N, min = 0, max = 1)
#   indMiss40 <- indPropMiss40 > runif(N, min = 0, max = 1)

#   # Impose NA
#   catDatMiss20 <- catDatMiss40 <- catDat
#   catDatMiss20[indMiss20, paste0("Y",1:nItemPerFactor)] <- NA
#   catDatMiss40[indMiss40, paste0("Y",1:nItemPerFactor)] <- NA

#   # Add repno
#   catDat <- cbind(repno = repno, catDat)
#   catDatMiss20 <- cbind(repno = repno, catDatMiss20)
#   catDatMiss40 <- cbind(repno = repno, catDatMiss40)

#   # # Some correlations
#   # missimpact <- round(cor(cbind(data.frame(Z, Z.rank, indPropMiss20, indMiss20, indMiss40), 
#   #                         catDat[,paste0("Y",1:nItemPerFactor)])), 2)[,1:5]
#   # print(missimpact)

#   list(comp = catDat, miss20 = catDatMiss20, miss40 = catDatMiss40, seed = seedused)
# }

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

anaComp <- function(maindir, nCat, thDist, N, repno, propMiss = 0, anaModel, est, sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  cat_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))
  
  model <- get(anaModel) # CM1, ICM1, BM

  # Import data
  filename <- genPath(maindir, nCat, thDist, N, repno, propMiss)
  data <- data.frame(read.table(filename, sep = ",", header = FALSE, na.strings="9"))
  colnames(data) <- c("repno", cat_items)
  
  # Analysis of complete data
  if (isTRUE(propMiss == 0)) {
      fit <- Catch(sem(model            = model,
                       data             = data,
                       ordered          = cat_items,
                       std.lv           = TRUE,
                       parameterization = "delta",
                       estimator        = est))
  } else {
      fit <- Catch(sem(model            = model,
                       data             = data,
                       ordered          = cat_items,
                       std.lv           = TRUE,
                       parameterization = "delta",
                       estimator        = est,
                       missing          = "pairwise"))
  }

  # fit is a list: 1) results from sem, 2) error msg, 3) warning msg
  err_and_warn <- fit[2:3]    
  
  fit <- fit[[1]]
  est_se <- lav_est_se(fit) # extract lavaan estimates and standard errors
  param <- est_se$est
  se <- est_se$se
  fitstat <- try(fitMeasures(fit), silent = TRUE)
  teststat <- try(calculateT(fit), silent = TRUE)
  nobs <- inspect(fit, "nobs")

  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss,
             anaModel = anaModel, est = est)

  output <- list(conds     = conds,
                 param     = param,
                 se        = se,
                 fitstat   = fitstat,
                 teststat  = teststat,
                 nobs      = nobs,
                 err       = err_and_warn$err,
                 warn      = err_and_warn$warn)
  return(output)
}

# Extract lavaan estimates and standard errors
# sum(cbind(lav_coef_se(fit, "all")$est,lav_coef_se(fit, "all")$se) - lavaan::parameterEstimates(fit)[,4:5]) # 0
# all.equal(names(lav_coef_se(fit, "all")$est), apply(lavaan::parameterEstimates(fit)[1:3], 1, paste, collapse="")) # TRUE
lav_est_se <- function(object, type = "free", add.labels = TRUE) {
  # Based on lavaan:::lav_object_inspect_coef
  if (type == "user" || type == "all") {
      type <- "user"
      idx <- 1:length(object@ParTable$lhs)
  }
  else if (type == "free") {
      idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
  }
  else {
      stop("lavaan ERROR: argument 'type' must be one of free or user")
  }
  EST <- object@Fit@est
  SE <- object@Fit@se
  cof <- EST[idx]
  se <- SE[idx]
  if (add.labels) {
      names(cof) <- names(se) <- lavaan:::lav_partable_labels(object@ParTable, type = type)
  }
  output <- list(est = cof, se = se)
  return(output)
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

## Cont (fix first loading to be 1)
# Model CM1: 3-factor CFA for X, M, Y with 6 items per factor
CM1_cont <- '
  X =~ X1 + X2 + X3 + X4 + X5 + X6
  M =~ M1 + M2 + M3 + M4 + M5 + M6
  Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6'

# Model ICM1: X --> M --> Y
ICM1_cont <- '
  X =~ X1 + X2 + X3 + X4 + X5 + X6
  M =~ M1 + M2 + M3 + M4 + M5 + M6
  Y =~ Y1 + Y2 + Y3 + Y4 + Y5 + Y6
  Y ~ M
  M ~ X
  Y ~ 0*X'

anaFIMLcont_lavaan <- function(maindir, nCat, thDist, N, repno, propMiss = 0, anaModel, sourcedir = NULL) {
  if (!is.null(sourcedir)) source(sourcedir)
  
  cont_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))
  
  model <- get(paste(anaModel,"cont",sep ="_")) # CM1, ICM1

  # Import data
  filename <- genPath(maindir, nCat, thDist, N, repno, propMiss)
  data <- data.frame(read.table(filename, sep = ",", header = FALSE, na.strings="9"))
  colnames(data) <- c("repno", cont_items)
  
  # Analysis of complete data
  fit <- Catch(sem(model            = model,
                   data             = data,
                   estimator        = "mlr",
                   missing          = "fiml"))

  # fit is a list: 1) results from sem, 2) error msg, 3) warning msg
  err_and_warn <- fit[2:3]    
  
  fit <- fit[[1]]
  paramEst <- parameterEstimates(fit, standardized = TRUE)
  param <- paramEst$std.all # std.all
  std.lv <- paramEst$std.lv # std.lv
  unstd <- paramEst$est # unstd
  se <- paramEst$se
  names(param) <- names(unstd) <- names(se) <- names(std.lv) <- apply(paramEst[, c("lhs", "op", "rhs")], 1, function(x) paste(x, collapse = ""))
  # To be consistent with cat
  exclude <- !grepl("X~~X|M~~M|Y~~Y|X~1|M~1|Y~1", names(param))
  param <- param[exclude]
  std.lv <- std.lv[exclude]
  unstd <- unstd[exclude]
  se <- se[exclude]

  fitstat <- try(fitMeasures(fit), silent = TRUE)
  nobs <- inspect(fit, "nobs")
    
  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss,
             anaModel = anaModel, est = "mlr")

  output <- list(conds     = conds,
                 unstd     = unstd,
                 se        = se,
                 param     = param,
                 std.lv    = std.lv,
                 fitstat   = fitstat,
                 nobs      = nobs,
                 err       = err_and_warn$err,
                 warn      = err_and_warn$warn)
  return(output)
}

# mplus
## fiml cont
## fiml cat (no fit stats)
## bayes cont (PPP + approx fit indices)
## bayes cat (PPP)
ana_mplus <- function(maindir, nCat, thDist, N, repno, propMiss, 
                      EST = "bayes", ITEM = "cat", std.lv = FALSE, fullmed = FALSE, suffix = NULL, missflag = 9, BSEED = NULL, sourcedir = sourcedir) {
  suppressMessages(library(MplusAutomation))
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions

  # Replace text in Mplus inp template
  FILE <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss) # location of the .csv data file 
  FILE <- paste0(dirname(FILE),"\n/",basename(FILE)) # seperate to multiple lines
  syntax <- gsub("MISSING_DATA_FILE", FILE, inp_template_ana)
  syntax <- gsub("ESTIMATOR_COMMAND", EST, syntax)
  syntax <- gsub("MISSFLAG", missflag, syntax)
  if (tolower(EST) == "bayes") {
    syntax <- gsub("!BITERATIONS", "BITERATIONS", syntax)
    syntax <- gsub("!CHAINS", "CHAINS", syntax)
    syntax <- gsub("!BCONVERGENCE", "BCONVERGENCE", syntax)
    syntax <- gsub("OUTPUT: ", "OUTPUT: TECH8 ", syntax)
  }
  if (isTRUE(std.lv)) {
    syntax <- gsub("X BY X1-X6;", "X BY X1* X2-X6; X@1;", syntax)
    syntax <- gsub("M BY M1-M6;", "M BY M1* M2-M6; M@1;", syntax)
    syntax <- gsub("Y BY Y1-Y6;", "Y BY Y1* Y2-Y6; Y@1;", syntax)
  }
  if (isTRUE(fullmed)) {
    syntax <- gsub("MODEL:", "MODEL:\nY ON M; M ON X; Y ON X@0;", syntax)
  }
  if (tolower(ITEM) == "cat") {
    syntax <- gsub("!categorical", "categorical", syntax)
  }
  if (!is.null(BSEED)) {
    if (BSEED == "RANDOM") BSEED <- sample(1:99999999, 1)
    syntax <- gsub("!BSEED = 0", paste0("BSEED = ", BSEED), syntax)
  }
  cat(writeLines(syntax))

  # Create Mplus .inp file
  inp_dir <- paste0(maindir,"/ana_mplus/",EST,ITEM,suffix,"/rep",repno)
  if (!dir.exists(inp_dir)) dir.create(inp_dir, recursive=TRUE)
  inp_basename <- sub("csv","inp",basename(FILE))
  inp_full_path <- paste(inp_dir, inp_basename, sep = "/")
  writeLines(syntax, inp_full_path) # write .inp file

  # Run Mplus via MplusAutomation (batch mode) 
  MplusAutomation::runModels(inp_full_path, logFile = NULL)

  out_full_path <- sub("inp","out",inp_full_path)
  out <- MplusAutomation::readModels(out_full_path)
  return(out)
}

inp_template_ana <- 
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
!categorical = X1-Y6;
missing = ALL(MISSFLAG);

! CM1
MODEL:
X BY X1-X6;
M BY M1-M6;
Y BY Y1-Y6;

ANALYSIS:
ESTIMATOR = ESTIMATOR_COMMAND;
!BITERATIONS = 100001 (100000);
!CHAINS = 2;
!BSEED = 0;
!BCONVERGENCE = .025;

OUTPUT: STDYX PATTERNS TECH1;"

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

inp_template_h0 <- 
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
ESTIMATOR = BAYES;
BITERATIONS = 100000 (BURNIN);
CHAINS = 2;
!BSEED = 0;
BCONVERGENCE = .025;

MODEL:
X BY X1* X2-X6; X@1;
M BY M1* M2-M6; M@1;
Y BY Y1* Y2-Y6; Y@1;
X WITH M Y;
M WITH Y;

OUTPUT: TECH8;"

inp_template_h0_cov <- 
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
ESTIMATOR = BAYES;
BITERATIONS = 100000 (BURNIN);
CHAINS = 2;
!BSEED = 0;
BCONVERGENCE = .025;

MODEL:
X1-Y6 WITH X1-Y6;

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
BITERATIONS = 100001 (BURNIN);
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
                       h0 = FALSE,
                       h0_cov = FALSE,
                       savemaindir = NULL) {
  # Setup 
  suppressMessages(library(MplusAutomation))
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions
  
  if (isFALSE(convdiag)) {
    imptemp <- inp_template
    if (isTRUE(h0)) imptemp <- inp_template_h0 # Estimator = BAYES; MODEL: "CM1";
    if (isTRUE(h0_cov)) imptemp <- inp_template_h0_cov # Estimator = BAYES; MODEL: "Saturated";
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
  inp_full_path <- gsub("_imp[0-9]+.dat",".inp",inp_full_path) # change .dat to .inp
  writeLines(syntax, inp_full_path) # write .inp file

  # Run Mplus via MplusAutomation (batch mode) 
  MplusAutomation::runModels(inp_full_path, logFile = NULL)

  # Check PSR
  out_full_path <- sub("\\.inp", "\\.out",inp_full_path) # absolute path to .out 
  tech8 <- MplusAutomation::readModels(out_full_path, what = "tech8") # extract tech8 from the .out 
  PSR <- CheckPSR(tech8)

  # R output
  c("inp" = inp_full_path, PSR)
}

# fitH0 <- function(maindir, nCat, thDist, N, repno, propMiss, sourcedir = NULL) {
  
#   if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions

#   out_full_path <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss, MI = TRUE, createdir = FALSE) # absolute path 
#   out_full_path <- sub("dat", "out", out_full_path)

#   out <- MplusAutomation::readModels(out_full_path)
#   impdat.list <- out$savedata
#   impdat.list <- impdat.list[stringr::str_order(names(impdat.list), numeric = TRUE)] # sort imp1, imp2, .. NOT imp1, imp10, ..
  
#   cat_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

#   lavCor.list <- lapply(impdat.list, function(x) {
#       Catch(lavCor(x[,cat_items], ordered=cat_items, output="full", se = "standard", baseline = FALSE))
#     })

#   err.list <- lapply(lavCor.list, "[[", "err")
#   warn.list <- lapply(lavCor.list, "[[", "warn")
#   lavCor.list <- lapply(lavCor.list, "[[", "res")
#   pcorr.list <- lapply(lavCor.list, function(x) lavInspect(x, "cov.ov"))
#   threshold.list <- lapply(lavCor.list, function(x) lavInspect(x, "th"))
#   th.idx <- lavInspect(lavCor.list[[1]], "th.idx") # add attribute to thresholds
#   N <- lavCor.list[[1]]@SampleStats@nobs[[1]]
#   acov.list <- lapply(lavCor.list, vcov)
  
#   conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss)

#   # OUTPUT
#   OUTPUT <- list(conds          = conds,
#                  N              = N,
#                  pcorr.list     = pcorr.list,
#                  threshold.list = threshold.list,
#                  th.idx         = th.idx,
#                  acov.list      = acov.list,
#                  err.list       = err.list,
#                  warn.list      = warn.list)
#   OUTPUT
# }

anaImp <- function(maindir, nCat, thDist, N, repno, propMiss, sourcedir = NULL) {
  
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions

  # Fit H0
  out_full_path <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss, MI = TRUE, createdir = FALSE) # absolute path 
  out_full_path <- gsub("_imp[0-9]+.dat",".out",out_full_path)

  out <- MplusAutomation::readModels(out_full_path)
  impdat.list <- out$savedata
  if (length(impdat.list) != 300) stop("m is not 300")
  impdat.list <- impdat.list[stringr::str_order(names(impdat.list), numeric = TRUE)] # sort imp1, imp2, .. NOT imp1, imp10, ..
  
  cat_items <- c(paste0("X",1:6), paste0("M",1:6), paste0("Y",1:6))

  lavCor.list <- lapply(impdat.list, function(x) {
    Catch(lavCor(x[,cat_items], ordered=cat_items, output="full", se = "standard", baseline = FALSE))
  })

  err.list <- lapply(lavCor.list, "[[", "err")
  warn.list <- lapply(lavCor.list, "[[", "warn")
  lavCor.list <- lapply(lavCor.list, "[[", "res")

  # Fit H1
  # Correct model
  fitH1_CM1_dwls_m20  <- fitH1(lavCor.list, impno =  1:20,  anaModel = "CM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_CM1_dwls_m50  <- fitH1(lavCor.list, impno = 21:70,  anaModel = "CM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_CM1_dwls_m100 <- fitH1(lavCor.list, impno = 71:170, anaModel = "CM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_CM1_dwls_m300 <- fitH1(lavCor.list, impno =  1:300, anaModel = "CM1", estimator = "wlsmv", cat_items = cat_items)

  fitH1_CM1_uls_m20   <- fitH1(lavCor.list, impno =  1:20,  anaModel = "CM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM1_uls_m50   <- fitH1(lavCor.list, impno = 21:70,  anaModel = "CM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM1_uls_m100  <- fitH1(lavCor.list, impno = 71:170, anaModel = "CM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM1_uls_m300  <- fitH1(lavCor.list, impno =  1:300, anaModel = "CM1", estimator = "ulsmv", cat_items = cat_items)

  # Incorrect model
  fitH1_ICM1_dwls_m20  <- fitH1(lavCor.list, impno =  1:20,  anaModel = "ICM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_ICM1_dwls_m50  <- fitH1(lavCor.list, impno = 21:70,  anaModel = "ICM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_ICM1_dwls_m100 <- fitH1(lavCor.list, impno = 71:170, anaModel = "ICM1", estimator = "wlsmv", cat_items = cat_items)
  fitH1_ICM1_dwls_m300 <- fitH1(lavCor.list, impno =  1:300, anaModel = "ICM1", estimator = "wlsmv", cat_items = cat_items)

  fitH1_ICM1_uls_m20   <- fitH1(lavCor.list, impno =  1:20,  anaModel = "ICM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_ICM1_uls_m50   <- fitH1(lavCor.list, impno = 21:70,  anaModel = "ICM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_ICM1_uls_m100  <- fitH1(lavCor.list, impno = 71:170, anaModel = "ICM1", estimator = "ulsmv", cat_items = cat_items)
  fitH1_ICM1_uls_m300  <- fitH1(lavCor.list, impno =  1:300, anaModel = "ICM1", estimator = "ulsmv", cat_items = cat_items)

  # Baseline model
  fitH1_BM_dwls_m20  <- fitH1(lavCor.list, impno =  1:20,  anaModel = "BM", estimator = "wlsmv", cat_items = cat_items)
  fitH1_BM_dwls_m50  <- fitH1(lavCor.list, impno = 21:70,  anaModel = "BM", estimator = "wlsmv", cat_items = cat_items)
  fitH1_BM_dwls_m100 <- fitH1(lavCor.list, impno = 71:170, anaModel = "BM", estimator = "wlsmv", cat_items = cat_items)
  fitH1_BM_dwls_m300 <- fitH1(lavCor.list, impno =  1:300, anaModel = "BM", estimator = "wlsmv", cat_items = cat_items)

  fitH1_BM_uls_m20   <- fitH1(lavCor.list, impno =  1:20,  anaModel = "BM", estimator = "ulsmv", cat_items = cat_items)
  fitH1_BM_uls_m50   <- fitH1(lavCor.list, impno = 21:70,  anaModel = "BM", estimator = "ulsmv", cat_items = cat_items)
  fitH1_BM_uls_m100  <- fitH1(lavCor.list, impno = 71:170, anaModel = "BM", estimator = "ulsmv", cat_items = cat_items)
  fitH1_BM_uls_m300  <- fitH1(lavCor.list, impno =  1:300, anaModel = "BM", estimator = "ulsmv", cat_items = cat_items)
  
  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss)
  
  # OUTPUT
  output <- list(conds          = conds,
                 CM1_dwls_m20   = fitH1_CM1_dwls_m20,
                 CM1_dwls_m50   = fitH1_CM1_dwls_m50,
                 CM1_dwls_m100  = fitH1_CM1_dwls_m100,
                 CM1_dwls_m300  = fitH1_CM1_dwls_m300,
                 CM1_uls_m20    = fitH1_CM1_uls_m20,
                 CM1_uls_m50    = fitH1_CM1_uls_m50,
                 CM1_uls_m100   = fitH1_CM1_uls_m100,
                 CM1_uls_m300   = fitH1_CM1_uls_m300,
                 ICM1_dwls_m20  = fitH1_ICM1_dwls_m20,
                 ICM1_dwls_m50  = fitH1_ICM1_dwls_m50,
                 ICM1_dwls_m100 = fitH1_ICM1_dwls_m100,
                 ICM1_dwls_m300 = fitH1_ICM1_dwls_m300,
                 ICM1_uls_m20   = fitH1_ICM1_uls_m20,
                 ICM1_uls_m50   = fitH1_ICM1_uls_m50,
                 ICM1_uls_m100  = fitH1_ICM1_uls_m100,
                 ICM1_uls_m300  = fitH1_ICM1_uls_m300,
                 BM_dwls_m20    = fitH1_BM_dwls_m20,
                 BM_dwls_m50    = fitH1_BM_dwls_m50,
                 BM_dwls_m100   = fitH1_BM_dwls_m100,
                 BM_dwls_m300   = fitH1_BM_dwls_m300,
                 BM_uls_m20     = fitH1_BM_uls_m20,
                 BM_uls_m50     = fitH1_BM_uls_m50,
                 BM_uls_m100    = fitH1_BM_uls_m100,
                 BM_uls_m300    = fitH1_BM_uls_m300,
                 err.list       = err.list,
                 warn.list      = warn.list)
  return(output)
}

fitH1 <- function(lavCor.list, impno, anaModel, estimator, cat_items) {
  sampstats_pooled <- pool_ChungCai2019(lavCor.list, impno)
  
  nacov.total <- sampstats_pooled$nacov.total
  pcorr.avg <- sampstats_pooled$pcorr.avg
  threshold.avg <- sampstats_pooled$threshold.avg
  N <- sampstats_pooled$N
  
  if (estimator %in% c("wlsm", "wlsmv")) {
    WLS.V <- solve(diag(diag(nacov.total))) # dwls
  } else if (estimator %in% c("ulsm", "ulsmv")) {
    WLS.V <- diag(nrow(nacov.total)) # uls
  } else stop("Only wlsm, wlsmv, ulsm, or ulsmv is supported")

  model <- get(anaModel)

  fit <- Catch(sem(model = model, 
                   std.lv = TRUE,
                   ordered = cat_items, 
                   sample.cov = pcorr.avg,
                   sample.th = threshold.avg,
                   sample.nobs = N,
                   WLS.V = WLS.V,
                   NACOV = nacov.total,
                   parameterization = "delta",
                   estimator = estimator))
  
  err_and_warn <- fit[2:3]
  
  fit <- fit[[1]]
  est_se <- lav_est_se(fit) # extract lavaan estimates and standard errors
  param <- est_se$est
  se <- est_se$se
  fitstat <- try(fitMeasures(fit), silent = TRUE)
  teststat <- try(calculateT(fit), silent = TRUE)

  output <- list(param     = param, 
                 se        = se,
                 fitstat   = fitstat,
                 teststat  = teststat,
                 err       = err_and_warn$err,
                 warn      = err_and_warn$warn)
  # OUTPUT
  return(output)
}

anaImp_CM2 <- function(maindir, nCat, thDist, N, repno, propMiss, sourcedir = NULL) {
  
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions

  # Fit H0
  out_full_path <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss, MI = TRUE, createdir = FALSE) # absolute path 
  out_full_path <- sub("dat", "out", out_full_path)

  out <- MplusAutomation::readModels(out_full_path)
  impdat.list <- out$savedata
  impdat.list <- impdat.list[stringr::str_order(names(impdat.list), numeric = TRUE)] # sort imp1, imp2, .. NOT imp1, imp10, ..
  
  cat_items <- c(paste0("X",1:5), paste0("M",1:5), paste0("Y",1:5)) # 5 items per factor

  lavCor.list <- lapply(impdat.list, function(x) {
    Catch(lavCor(x[,cat_items], ordered=cat_items, output="full", se = "standard", baseline = FALSE))
  })

  err.list <- lapply(lavCor.list, "[[", "err")
  warn.list <- lapply(lavCor.list, "[[", "warn")
  lavCor.list <- lapply(lavCor.list, "[[", "res")

  # Fit H1
  # Correct model
  fitH1_CM2_uls_m20   <- fitH1(lavCor.list, impno =  1:20,  anaModel = "CM2", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM2_uls_m50   <- fitH1(lavCor.list, impno = 21:70,  anaModel = "CM2", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM2_uls_m100  <- fitH1(lavCor.list, impno = 71:170, anaModel = "CM2", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM2_uls_m300  <- fitH1(lavCor.list, impno =  1:300, anaModel = "CM2", estimator = "ulsmv", cat_items = cat_items)

  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss)
  
  # OUTPUT
  RES <- list(conds          = conds,
              CM2_uls_m20    = fitH1_CM2_uls_m20,
              CM2_uls_m50    = fitH1_CM2_uls_m50,
              CM2_uls_m100   = fitH1_CM2_uls_m100,
              CM2_uls_m300   = fitH1_CM2_uls_m300,
              err.list       = err.list,
              warn.list      = warn.list)
  RES
}

anaImp_CM3 <- function(maindir, nCat, thDist, N, repno, propMiss, sourcedir = NULL) {
  
  if(!is.null(sourcedir)) source(sourcedir) # load R objects and functions

  # Fit H0
  out_full_path <- genPath(maindir = maindir, nCat, thDist, N, repno, propMiss, MI = TRUE, createdir = FALSE) # absolute path 
  out_full_path <- sub("dat", "out", out_full_path)

  out <- MplusAutomation::readModels(out_full_path)
  impdat.list <- out$savedata
  impdat.list <- impdat.list[stringr::str_order(names(impdat.list), numeric = TRUE)] # sort imp1, imp2, .. NOT imp1, imp10, ..
  
  cat_items <- c(paste0("X",1:3), paste0("M",1:3), paste0("Y",1:3)) # 3 items per factor

  lavCor.list <- lapply(impdat.list, function(x) {
    Catch(lavCor(x[,cat_items], ordered=cat_items, output="full", se = "standard", baseline = FALSE))
  })

  err.list <- lapply(lavCor.list, "[[", "err")
  warn.list <- lapply(lavCor.list, "[[", "warn")
  lavCor.list <- lapply(lavCor.list, "[[", "res")

  # Fit H1
  # Correct model
  fitH1_CM3_uls_m20   <- fitH1(lavCor.list, impno =  1:20,  anaModel = "CM3", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM3_uls_m50   <- fitH1(lavCor.list, impno = 21:70,  anaModel = "CM3", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM3_uls_m100  <- fitH1(lavCor.list, impno = 71:170, anaModel = "CM3", estimator = "ulsmv", cat_items = cat_items)
  fitH1_CM3_uls_m300  <- fitH1(lavCor.list, impno =  1:300, anaModel = "CM3", estimator = "ulsmv", cat_items = cat_items)

  conds <- c(nCat = nCat, thDist = thDist, N = N, repno = repno, propMiss = propMiss)
  
  # OUTPUT
  RES <- list(conds          = conds,
              CM3_uls_m20    = fitH1_CM3_uls_m20,
              CM3_uls_m50    = fitH1_CM3_uls_m50,
              CM3_uls_m100   = fitH1_CM3_uls_m100,
              CM3_uls_m300   = fitH1_CM3_uls_m300,
              err.list       = err.list,
              warn.list      = warn.list)
  RES
}


# pool_ChungCai2019 <- function(pcorr.list, threshold.list, th.idx, acov.list, N) {
#   # Chung & Cai (2019, Eqs. 19, 21, 22)
#   # Note: lavaan requires thresholds
  
#   m <- length(pcorr.list) # number of imputed datasets
#   # Average polychoric correlation matrix
#   pcorr.avg <- Reduce('+', pcorr.list)/length(pcorr.list)
#   pcorr.avg.vec <- lavaan::lav_matrix_vech(pcorr.avg, diagonal = FALSE)

#   # Average threshold matrix
#   threshold.avg <- Reduce('+', threshold.list)/length(threshold.list)
#   threshold.avg <- as.numeric(threshold.avg) 
#   attr(threshold.avg, "th.idx") <- th.idx # add attribute to thresholds

#   samstat.avg <- c(threshold.avg, pcorr.avg.vec) # Average thresholds and non-redundant correlations

#   # Between-imputation variance of the thresholds and the polychoric correlation matrix
#   samstat.list <- mapply(FUN = function(x,y) {
#     c(x,lavaan::lav_matrix_vech(y, diagonal = FALSE))
#   }, threshold.list, pcorr.list, SIMPLIFY=FALSE)
#   acov.between.list <- lapply(samstat.list, function(x) (x-samstat.avg) %*% t(x-samstat.avg))
#   acov.between <- Reduce('+', acov.between.list)/(m-1)

#   # Within-imputation asymptotic covariance matrix of the thresholds and the polychoric correlation matrix
#   acov.within <- Reduce('+', acov.list)/length(acov.list)
#   # CHECK: sort - thresholds come first
#   thr.where <- grepl("\\|",colnames(acov.within)) 
#   sort.acov.within <- c(which(thr.where), which(!thr.where))
#   acov.within <- acov.within[sort.acov.within, sort.acov.within]

#   # Corrected variance-covariance matrix of the thresholds and the polychoric correlation matrix 
#   acov.total <- acov.within + acov.between + acov.between/m
  
#   # lavaan uses nacov as input
#   nacov.total <- (N-1)*acov.total

#   # OUTPUT
#   OUTPUT <- list(nacov.total    = nacov.total,
#                  pcorr.avg      = pcorr.avg,
#                  threshold.avg  = threshold.avg,
#                  m              = m,
#                  N              = N)
#   OUTPUT
# }

pool_ChungCai2019 <- function(lavCor.list, impno = NULL) {
  # Chung & Cai (2019, Eqs. 19, 21, 22)
  # Note: lavaan requires thresholds
  
  if (!is.null(impno)) { 
    lavCor.list <- lavCor.list[impno]
  }

  m <- length(lavCor.list) # number of imputed datasets
  # Average polychoric correlation matrix
  pcorr.list <- lapply(lavCor.list, function(x) lavInspect(x, "cov.ov"))
  pcorr.avg <- Reduce('+', pcorr.list)/length(pcorr.list)
  pcorr.avg.vec <- lavaan::lav_matrix_vech(pcorr.avg, diagonal = FALSE)

  # Average threshold matrix
  threshold.list <- lapply(lavCor.list, function(x) lavInspect(x, "th"))
  threshold.avg <- Reduce('+', threshold.list)/length(threshold.list)
  threshold.avg <- as.numeric(threshold.avg) 
  attr(threshold.avg, "th.idx") <- lavInspect(lavCor.list[[1]], "th.idx") # add attribute to thresholds

  samstat.avg <- c(threshold.avg, pcorr.avg.vec) # Average thresholds and non-redundant correlations

  # Between-imputation variance of the thresholds and the polychoric correlation matrix
  samstat.list <- mapply(FUN = function(x,y) {
    c(x,lavaan::lav_matrix_vech(y, diagonal = FALSE))
  }, threshold.list, pcorr.list, SIMPLIFY=FALSE)
  acov.between.list <- lapply(samstat.list, function(x) (x-samstat.avg) %*% t(x-samstat.avg))
  acov.between <- Reduce('+', acov.between.list)/(m-1)

  # Within-imputation asymptotic covariance matrix of the thresholds and the polychoric correlation matrix
  acov.list <- lapply(lavCor.list, vcov)
  acov.within <- Reduce('+', acov.list)/length(acov.list)
  # CHECK: sort - thresholds come first
  thr.where <- grepl("\\|",colnames(acov.within)) 
  sort.acov.within <- c(which(thr.where), which(!thr.where))
  acov.within <- acov.within[sort.acov.within, sort.acov.within]

  # Corrected variance-covariance matrix of the thresholds and the polychoric correlation matrix 
  acov.total <- acov.within + acov.between + acov.between/m
  
  # lavaan uses nacov as input
  N <- lavCor.list[[1]]@SampleStats@nobs[[1]]
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
            "TYB" = TYB, "pvalue.TYB" = 1-pchisq(TYB,df),
            "a" = a, "b" = b, "c" = c)
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


## getBias()
getBias <- function(param.onerow, param.SE.onerow = NULL) {
  
  cond <- param.onerow %>% select(-any_of(sortParam)) # nonparam columns 
  nCat <- param.onerow$nCat
  thDist <- param.onerow$thDist
  anaModel <- param.onerow$anaModel
  nItemPerFactor <- ifelse(anaModel == "CM1", 6, NA)
  nFactors <- 3

  ## need to be in a one row format
  ## remove irrelevant columns if any
  param <- param.onerow[1,!is.na(param.onerow)] %>% select(any_of(sortParam))

  factorLoading <- 0.8 # factor loadings
  lambda <- rep(factorLoading, nItemPerFactor*nFactors)

  factorCorr <- 0.4 # factor correlations
  psi <- rep(factorCorr, nFactors*(nFactors-1)/2)
  
  # Thresholds
  if(nCat == 2) {
    if(thDist == "sym") {
      Th <- qnorm(c(0.5), mean = 0, sd = 1) # 50%,50%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0.85), mean = 0, sd = 1) # 85%,15%
    }
  } else if(nCat == 5) {
    if(thDist == "sym") {
      Th <- qnorm(c(0.07, 0.31, 0.69, 0.93), mean = 0, sd = 1) # 7%,24%,38%,24%,7%
    } else if(thDist == "asym") {
      Th <- qnorm(c(0.52, 0.67, 0.80, 0.91), mean = 0, sd = 1) # 52%,15%,13%,11%,9%
    }
  }
  thresholds <- rep(Th, nItemPerFactor*nFactors)

  if (isTRUE(length(c(lambda, psi, thresholds)) != length(param))) stop("unequal param length")

  # Population param
  pop.param <- param
  pop.param[grep("=~", colnames(pop.param))] <- lambda
  pop.param[grep("~~", colnames(pop.param))] <- psi
  pop.param[grep("\\|", colnames(pop.param))] <- thresholds

  # Bias
  bias <- param-pop.param
  bias.rel <- bias/pop.param
  bias.rel[pop.param == 0] <- NA

  bias.std <- NA
  if(!is.null(param.SE.onerow)) {
    param.SE <- param.SE.onerow[1,!is.na(param.SE.onerow)] %>% select(any_of(sortParam))
    bias.std <- bias/param.SE
    bias.std[param.SE == 0] <- 0
  }

  list("cond"       = cond,
       "param"      = param, 
       "pop.param"  = pop.param, 
       "bias"       = bias, 
       "bias.rel"   = bias.rel,
       "bias.std"   = bias.std) 
}

extractResults <- function(RES, where) {
  where <- where
  cond_temp1 <- as.data.frame(t(sapply(RES, function(x) x$conds)), stringsAsFactors = FALSE)
  cond_temp2 <- as.data.frame(str_split(where, "_", simplify = TRUE), stringsAsFactors = FALSE) 
  colnames(cond_temp2) <- c("anaModel","estimator","m")
  cond_temp2$m <- gsub("m", "",cond_temp2$m)
  cond <- bind_cols(cond_temp1, cond_temp2)
  cond <- mutate(cond, across(c(nCat,N,repno,propMiss,m), as.numeric))
  err <- sapply(RES, function(x) ifelse(is.null(x[[where]]$err), 0, 1))
  warn <- sapply(RES, function(x) ifelse(is.null(x[[where]]$warn), 0, 1))

  teststat <- t(sapply(RES, function(x) x[[where]]$teststat))
  Full_teststat <- cbind(cond, err, warn, teststat)

  param.list <- lapply(RES, function(x) data.frame(rbind(x[[where]]$param), check.names = FALSE))
  Full_param <- data.table::rbindlist(param.list, fill = TRUE, use.names=TRUE) %>% relocate(any_of(sortParam))
  Full_param <- cbind(cond, err, warn, Full_param)

  list(Full_param = Full_param, Full_teststat = Full_teststat)
}

## abs_max() is used to find max(|x|) and its corresponding name
abs_max <- function(x) { 
    if( all(is.na(x)) ) {
       abs_max <- NA
       which_max <- NA
    } else {
       abs_max <- max(abs(x), na.rm = TRUE)
       which_max <- names(which.max(abs(x))) 
    }
    data.frame(max = abs_max, which_max = which_max)
}

## A list containing functions (written in tidy-style), mean var min max nReps
getMeanVarMinMaxN <- list(mean = ~mean(.x, na.rm = TRUE),
                          var  = ~var(.x, na.rm = TRUE),
                          min  = ~min(.x, na.rm = TRUE),
                          max  = ~max(.x, na.rm = TRUE),
                          nReps = ~sum(!is.na(.x)))