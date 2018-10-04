#' 1. Calculating estimated ages on supported models with given parameters
#' 2. Returning standard diviation, mean, meadian and median absolute error(MAE) of difference in age between estimated age and actual age
#' @param   bs       a BSseq object (must have assays named `M` and `Cov`)
#' @param   actual_ages    a data.frame which must be given with two columnas, sid and age
#' @param   padding how many bases +/- to pad the target CpG by (default is 15)
#' @param   minCovg minimum regional read coverage desired to estimate 5mC (5)
#' @param   minSamp minimum number of non-NA samples to perform imputation (5)
#' @param   genome  genome to use as reference, if no genome(x) is set (NULL)

getOptAgeModel <- function(bs, actual_ages, imputed = NULL, padding = 15, minCovg = 5, minSamp = 5, genome=NULL) {
    
    if(!is(actual_ages, "data.frame")) {
        stop("actual_ages is not given in data.frame")
    }
    
    if(ncol(actual_ages) != 2) {
        stop("actual_ages must provide with two columns, sid and age")
    }
    
    if(!all(c("sid", "age") == colnames(actual_ages))) {
        stop("actual_ages must provide with two columns, sid and age")
    }
    if(is.null(imputed)) {
        imputed = c(TRUE, FALSE)
    }
    
    if(typeof(imputed) != "logical") {
        stop("imputed must be give in logical values(TRUE and FALSE) or NULL")
    }
    
    param <- list()
    param$padding <- padding
    param$minCovg <- minCovg
    param$minSamp <- minSamp
    
    samples <- data.frame(sid=bs$sampleNames)
    param$actual_ages <- merge(actual_ages, samples, by="sid")
    if(nrow(param$actual_ages) == 0) {
        message("The number of matched sample ids between actual_ages and bsseq object is 0")
        stop("Need to check bs$sampleNames or actual_ages")
    }
    
    models <- c("horvath","horvathshrunk","hannum","skinandblood")
    
    # sort out assemblies
    param$genome <- unique(genome(bs))
    if (is.null(param$genome)) param$genome <- genome
    if (!param$genome %in% c("hg19","GRCh37","hg38","GRCh38")) {
        message("genome(bs) is set to ",unique(genome(bs)),", which is unsupported.")
        stop("Provide a `genome` argument, or set genome(bs) manually, to proceed.")
    }
    
    
    
    res <- NA
    
    for(i in 1:length(models)) {
        param$model <- models[i]
        clocks <- .getClockAllConds(param)
        message("Calculating WGBSage on data set with ", param$model)
        for(j in 1:length(clocks)) {
            param$clock_name <- names(clocks)[j]
            param$clock <- clocks[[param$clock_name]]
            
            for(k in 1:length(imputed)) {
                param$imputed <- imputed[k]
                
                param$age <- .calcWGBSage(bs=bs, param)
                
                res <- .mergeRes(res=res, param)
            }
        }
        message("Done.")

    }
    
    ## sort result by columns, abs(diff_mean) and abs(diff_median)
    res <- res[order(abs(res[,9]), abs(res[,10]) ),]
    rownames(res) <- NULL
    
    return(res)
}

.mergeRes <- function(res, param) {
    model <- param$model
    clock_name <- param$clock_name
    age <- param$age
    padding <- param$padding
    minCovg <- param$minCovg
    minSamp <- param$minSamp
    actual_ages <- param$actual_ages
    imputed <- param$imputed
    
    if(clock_name == "N") {
        useENSR <- FALSE; useHMMI <- FALSE;
    } else if(clock_name == "EN") {
        useENSR <- TRUE; useHMMI <- FALSE;
    } else if(clock_name == "HM") {
        useENSR <- FALSE; useHMMI <- TRUE;
    } else if(clock_name == "ALL") {
        useENSR <- TRUE; useHMMI <- TRUE;
    }
    
    temp_res <- data.frame(model=model, padding=padding, minCovg=minCovg, minSamp=minSamp, useENSR=useENSR, useHMMI=useHMMI, IMPUTED=imputed)
    
    m <- merge(age, actual_ages, by="sid")
    diff <- m$est_age - m$age
    diff <- diff[!is.na(diff)]
    
    temp_res <- data.frame(temp_res, diff_sd=sd(diff), diff_mean=mean(diff), diff_median=median(diff), MAE=median(abs(diff)))
    
    if(is.na(res)) {
        return(temp_res)
    } else {
        return(rbind(res, temp_res))
    }
}

.getClockAllConds <- function(param) {
    model <- param$model; padding <- param$padding; g <- param$genome;
    
    return(list(N = getClock(model=model, padding=padding, genome=g, useENSR=F, useHMMI=F), EN = getClock(model=model, padding=padding, genome=g, useENSR=T, useHMMI=F), HM = getClock(model=model, padding=padding, genome=g, useENSR=F, useHMMI=T), ALL = getClock(model=model, padding=padding, genome=g, useENSR=T, useHMMI=T)))
}

.calcWGBSage <- function(bs, param) {
    clock <- param$clock
    minCovg <- param$minCovg
    minSamp <- param$minSamp
    imputed <- param$imputed
    
    covgWGBSage <- getCoverage(bs, regions=clock$gr, what="perRegionTotal")
    rownames(covgWGBSage) <- names(clock$gr)
    
    methWGBSage <- getMeth(bs, regions=clock$gr, type="raw", what="perRegion")
    rownames(methWGBSage) <- as.character(clock$gr)
    methWGBSage[covgWGBSage < minCovg] <- NA
    methWGBSage <- as(methWGBSage, "matrix") # 353 x ncol(x) is not too huge
    
    if(imputed) { methWGBSage <- impute.knn(methWGBSage, k=minSamp)$data }
    
    keep <- (rowSums2(is.na(methWGBSage)) < 1)
    
    if (!all(keep)) methWGBSage <- methWGBSage[which(keep), ]
    
    names(clock$gr) <- as.character(granges(clock$gr))
    
    coefs <- clock$gr[rownames(methWGBSage)]$score
    
    names(coefs) <- rownames(methWGBSage)
    
    agePredRaw <- (clock$intercept + (t(methWGBSage) %*% coefs))
    agePredRaw <- clock$cleanup(agePredRaw)
    agePredRaw <- data.frame(sid=rownames(agePredRaw), est_age=agePredRaw[,1], row.names = NULL)
    
    return(agePredRaw)
}
