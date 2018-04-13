

#' Measures for Karyotype Heterogeneity
#' 
#' Computes measures for karyotype heterogeneity. See the Details section for how these measures are defined.
#'
#' We define \eqn{x} as the vector of copy number states for each position. The number of HMMs is \eqn{S}. The measures are computed for each bin as follows:
#' \describe{
#' \item{Aneuploidy:}{\eqn{D = mean( abs(x-P) )}, where P is the physiological number of chromosomes at that position.}
#' \item{Heterogeneity:}{\eqn{H = sum( table(x) * 0:(length(table(x))-1) ) / S}}
#' }
#'
#' @param hmms A list with \code{\link{aneuHMM}} objects or a list of files that contain such objects.
#' @param normalChromosomeNumbers A named integer vector or matrix with physiological copy numbers, where each element (vector) or column (matrix) corresponds to a chromosome. This is useful to specify male or female samples, e.g. \code{c('X'=2)} for female samples or \code{c('X'=1,'Y'=1)} for male samples. Specify a vector if all your \code{hmms} have the same physiological copy numbers. Specify a matrix if your \code{hmms} have different physiological copy numbers (e.g. a mix of male and female samples). If not specified otherwise, '2' will be assumed for all chromosomes.
#' @param regions A \code{\link{GRanges-class}} object containing ranges for which the karyotype measures will be computed.
#' @param exclude.regions A \code{\link{GRanges-class}} with regions that will be excluded from the computation of the karyotype measures. This can be useful to exclude regions with artifacts.
#' @return A \code{list} with two \code{data.frame}s, containing the karyotype measures $genomewide and $per.chromosome. If \code{region} was specified, a third list entry $regions will contain the regions with karyotype measures.
#' @author Aaron Taudt
#' @importFrom stats weighted.mean
#' @export
#'@examples
#'### Example 1 ###
#'## Get results from a small-cell-lung-cancer
#'lung.folder <- system.file("extdata", "primary-lung", "hmms", package="AneuFinderData")
#'lung.files <- list.files(lung.folder, full.names=TRUE)
#'## Get results from the liver metastasis of the same patient
#'liver.folder <- system.file("extdata", "metastasis-liver", "hmms", package="AneuFinderData")
#'liver.files <- list.files(liver.folder, full.names=TRUE)
#'## Compare karyotype measures between the two cancers
#'normal.chrom.numbers <- rep(2, 23)
#'names(normal.chrom.numbers) <- c(1:22,'X')
#'lung <- karyotypeMeasures(lung.files, normalChromosomeNumbers=normal.chrom.numbers)
#'liver <- karyotypeMeasures(liver.files, normalChromosomeNumbers=normal.chrom.numbers)
#'print(lung$genomewide)
#'print(liver$genomewide)
#'
#'### Example 2 ###
#'## Construct a matrix with physiological copy numbers for a mix of 5 male and 5 female samples
#'normal.chrom.numbers <- matrix(2, nrow=10, ncol=24,
#'                               dimnames=list(sample=c(paste('male', 1:5), paste('female', 6:10)),
#'                                             chromosome=c(1:22,'X','Y')))
#'normal.chrom.numbers[1:5,c('X','Y')] <- 1
#'normal.chrom.numbers[6:10,c('Y')] <- 0
#'print(normal.chrom.numbers)
#'
#'### Example 3 ###
#'## Exclude artifact regions with high variance
#'consensus <- consensusSegments(c(lung.files, liver.files))
#'variance <- apply(consensus$copy.number, 1, var)
#'exclude.regions <- consensus[variance > quantile(variance, 0.999)]
#'## Compare karyotype measures between the two cancers
#'normal.chrom.numbers <- rep(2, 23)
#'names(normal.chrom.numbers) <- c(1:22,'X')
#'lung <- karyotypeMeasures(lung.files, normalChromosomeNumbers=normal.chrom.numbers,
#'                          exclude.regions = exclude.regions)
#'liver <- karyotypeMeasures(liver.files, normalChromosomeNumbers=normal.chrom.numbers,
#'                           exclude.regions = exclude.regions)
#'print(lung$genomewide)
#'print(liver$genomewide)
karyotypeMeasures <- function(hmms, normalChromosomeNumbers=NULL, regions=NULL, exclude.regions=NULL) {

    ## Check user input
    if (is.matrix(normalChromosomeNumbers)) {
        if (nrow(normalChromosomeNumbers) != length(hmms)) {
            stop("nrow(normalChromosomeNumbers) must be equal to length(hmms)")
        }
    }
    ## Load the files
    hmms <- loadFromFiles(hmms, check.class=c("aneuHMM", "aneuBiHMM"))
  
    ## If all binsizes are the same the consensus template can be chosen equal to the bins
    ptm <- startTimedMessage("Making consensus template ...")
    binsizes <- unlist(lapply(hmms, function(x) { width(x$bins)[1] }))
    # Filter out HMMs where segments of bins$state are NULL
    mask <- !sapply(hmms, function(hmm) { is.null(hmm$segments) | is.null(hmm$bins$state) })
    hmms <- hmms[mask]
    if (all(binsizes==binsizes[1])) {
        consensus <- hmms[[1]]$bins
        mcols(consensus) <- NULL
        constates <- array(NA, dim=c(length(consensus), length(hmms)), dimnames=list(chromosome=as.vector(seqnames(consensus)), sample=sapply(hmms, '[[', 'ID')))
        for (i1 in 1:length(hmms)) {
            hmm <- hmms[[i1]]
            constates[,i1] <- hmm$bins$copy.number
        }
    } else { # binsizes differ
        consensus <- consensusSegments(hmms)
        constates <- consensus$copy.number
    }
    stopTimedMessage(ptm)
  
    ### Exclude regions ###
    if (!is.null(exclude.regions)) {
        ind <- findOverlaps(consensus, exclude.regions)@from
    		constates <- constates[-ind,]
    		consensus <- consensus[-ind]
    }

    ### Karyotype measures ###
    ptm <- startTimedMessage("Karyotype measures ...")
    result <- list()
    S <- ncol(constates)
    ## Genomewide
    physioState <- array(2, dim=dim(constates), dimnames=list(chromosome=rownames(constates), sample=NULL))
    if (!is.null(normalChromosomeNumbers)) {
        if (is.vector(normalChromosomeNumbers)) {
            mask <- rownames(physioState) %in% names(normalChromosomeNumbers)
            physioState[mask,] <- normalChromosomeNumbers[rownames(physioState)[mask]]
        } else if (is.matrix(normalChromosomeNumbers)) {
            mask <- rownames(physioState) %in% colnames(normalChromosomeNumbers)
            physioState[mask, ] <- t(normalChromosomeNumbers[, rownames(physioState)[mask]])
        }
    }
    consensus$Aneuploidy <- rowMeans(abs(constates - physioState))
    tabs <- apply(constates - physioState, 1, function(x) { sort(table(x), decreasing=TRUE) }) # Heterogeneity score in reference to the physiological state
#     tabs <- apply(constates, 1, function(x) { sort(table(x), decreasing=TRUE) })
    if (is.list(tabs) | is.vector(tabs)) {
        consensus$Heterogeneity <- unlist(lapply(tabs, function(x) { sum(x * 0:(length(x)-1)) })) / S
    } else if (is.matrix(tabs)) {
        consensus$Heterogeneity <- colSums( tabs * 0:(nrow(tabs)-1) ) / S
    }
    weights <- as.numeric(width(consensus))
    result[['genomewide']] <- data.frame(Aneuploidy = stats::weighted.mean(consensus$Aneuploidy, weights),
                                        Heterogeneity = stats::weighted.mean(consensus$Heterogeneity, weights))
    ## Chromosomes
    consensus.split <- split(consensus, seqnames(consensus))
    weights.split <- split(weights, seqnames(consensus))
    result[['per.chromosome']] <- data.frame(Aneuploidy = unlist(mapply(function(x,y) { stats::weighted.mean(x$Aneuploidy, y) }, consensus.split, weights.split)),
                                              Heterogeneity = unlist(mapply(function(x,y) { stats::weighted.mean(x$Heterogeneity, y) }, consensus.split, weights.split)))
    
    ## Region
    if (!is.null(regions)) {
        regions <- subsetByOverlaps(regions, consensus)
        # Find split vector
        ind <- findOverlaps(consensus, regions, select='first')
        consensus.split <- split(consensus, ind)
        weights.split <- split(weights, ind)
        regions$Aneuploidy <- unlist(mapply(function(x,y) { stats::weighted.mean(x$Aneuploidy, y) }, consensus.split, weights.split))
        regions$Heterogeneity <- unlist(mapply(function(x,y) { stats::weighted.mean(x$Heterogeneity, y) }, consensus.split, weights.split))
        result[['regions']] <- regions
    }
    stopTimedMessage(ptm)
    
    return(result)

}


