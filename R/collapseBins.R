#' Collapse consecutive bins
#'
#' The function will collapse consecutive bins which have, for example, the same combinatorial state.
#'
#' The following tables illustrate the principle of the collapsing:
#'
#' Input data:
#' \tabular{rrrrrr}{
#' seqnames \tab start \tab end \tab column2collapseBy \tab moreColumns \tab columns2sumUp \cr
#' chr1     \tab     0 \tab 199 \tab                 2 \tab        1 10 \tab           1 3 \cr
#' chr1     \tab   200 \tab 399 \tab                 2 \tab        2 11 \tab           0 3 \cr
#' chr1     \tab   400 \tab 599 \tab                 2 \tab        3 12 \tab           1 3 \cr
#' chr1     \tab   600 \tab 799 \tab                 1 \tab        4 13 \tab           0 3 \cr
#' chr1     \tab   800 \tab 999 \tab                 1 \tab        5 14 \tab           1 3 \cr
#' }
#' Output data:
#' \tabular{rrrrrr}{
#' seqnames \tab start \tab end \tab column2collapseBy \tab moreColumns \tab columns2sumUp \cr
#' chr1     \tab     0 \tab 599 \tab                 2 \tab        1 10 \tab           2 9 \cr
#' chr1     \tab   600 \tab 999 \tab                 1 \tab        4 13 \tab           1 6 \cr
#' }
#' 
#' @param data A data.frame containing the genomic coordinates in the first three columns.
#' @param column2collapseBy The number of the column which will be used to collapse all other inputs. If a set of consecutive bins has the same value in this column, they will be aggregated into one bin with adjusted genomic coordinates. If \code{NULL} directly adjacent bins will be collapsed.
#' @param columns2sumUp Column numbers that will be summed during the aggregation process.
#' @param columns2average Column numbers that will be averaged during the aggregation process.
#' @param columns2getMax Column numbers where the maximum will be chosen during the aggregation process.
#' @param columns2drop Column numbers that will be dropped after the aggregation process.
#' @return A data.frame.
#' @author Aaron Taudt
#' @export
#' @examples
#'## Get an example BED file with single-cell-sequencing reads
#'bedfile <- system.file("extdata", "KK150311_VI_07.bam.bed.gz", package="AneuFinderData")
#'## Bin the BAM file into bin size 1Mp
#'binned <- binReads(bedfile, assembly='mm10', binsize=1e6,
#'                   chromosomes=c(1:19,'X','Y'))
#'## Collapse the bins by chromosome and get average, summed and maximum read count
#'df <- as.data.frame(binned[[1]])
#'# Remove one bin for illustration purposes
#'df <- df[-3,]
#'head(df)
#'collapseBins(df, column2collapseBy='seqnames', columns2sumUp=c('width','counts'),
#'                        columns2average='counts', columns2getMax='counts',
#'                        columns2drop=c('mcounts','pcounts'))
#'collapseBins(df, column2collapseBy=NULL, columns2sumUp=c('width','counts'),
#'                        columns2average='counts', columns2getMax='counts',
#'                        columns2drop=c('mcounts','pcounts'))
#'
collapseBins = function(data, column2collapseBy=NULL, columns2sumUp=NULL, columns2average=NULL, columns2getMax=NULL, columns2drop=NULL) {

	## Name to index
	if (is.character(column2collapseBy)) {
		column2collapseBy <- which(column2collapseBy == names(data))
	}
	if (is.character(columns2sumUp)) {
		columns2sumUp <- unlist(lapply(columns2sumUp, function(x) { which(x == names(data)) }))
	}
	if (is.character(columns2average)) {
		columns2average <- unlist(lapply(columns2average, function(x) { which(x == names(data)) }))
	}
	if (is.character(columns2getMax)) {
		columns2getMax <- unlist(lapply(columns2getMax, function(x) { which(x == names(data)) }))
	}
	if (is.character(columns2drop)) {
		columns2drop <- unlist(lapply(columns2drop, function(x) { which(x == names(data)) }))
	}
	## Indices
	ind_coords <- 1:3
	ind_morecols <- setdiff(1:ncol(data), c(ind_coords, columns2sumUp, columns2average, columns2getMax, columns2drop))
	ind_sumcols <- columns2sumUp
	ind_meancols <- columns2average
	ind_maxcols <- columns2getMax

	## Make the comparison vector
	ptm <- startTimedMessage('Making comparison vector ...')
	if (is.null(column2collapseBy)) {
		c <- data$start
		cShift1 <- rep(NA,length(c))
		cShift1[2:length(cShift1)] <- data$end[-length(c)] + 1
	} else {
		if (is(data[,column2collapseBy], "factor")) {
			c <- as.integer(data[,column2collapseBy])
		} else {
			c <- data[,column2collapseBy]
		}
		cShift1 <- rep(NA,length(c))
		cShift1[-1] <- c[-length(c)]
	}
	compare_custom <- c != cShift1
	## Make the comparison vector to separate chromosomes
	c <- as.integer(data[,1])
	cShift1 <- rep(NA,length(c))
	cShift1[-1] <- c[-length(c)]
	compare_chrom <- c != cShift1
	## Combine the vectors
	compare <- compare_custom | compare_chrom
	compare[1] <- TRUE
	numcollapsedbins <- length(which(compare==TRUE))
	numbins <- nrow(data)
	stopTimedMessage(ptm)
	if (any(is.na(compare))) {
		stop("NAs in vector 'compare'")
	}

	## Select the collapsed rows
	ptm <- startTimedMessage('Selecting rows ...')
	collapsed.bins <- list()
	collapsed.bins[[names(data)[1]]] <- data[which(compare),1] #which to remove NAs which shouldn't be there in the first place
	collapsed.bins[[names(data)[2]]] <- data[which(compare),2]
	collapsed.bins[[names(data)[3]]] <- data[c((which(compare)-1)[-1],numbins), 3]
	if (length(ind_morecols)==1) {
		collapsed.bins[[names(data)[ind_morecols]]] <- data[which(compare), ind_morecols]
	} else if (length(ind_morecols)>1) {
		lcb <- length(collapsed.bins)
		lmc <- length(ind_morecols)
		collapsed.bins[(lcb+1):(lcb+lmc)] <- data[which(compare), ind_morecols]
		names(collapsed.bins)[(lcb+1):(lcb+lmc)] <- names(data)[ind_morecols]
	}
	stopTimedMessage(ptm)

	## Sum up columns
	xfuns <- list(sum, mean, max)
	xstrings <- list('sum', 'mean', 'max')
	columns2xs <- list(columns2sumUp, columns2average, columns2getMax)
	inds_xcols <- list(ind_sumcols, ind_meancols, ind_maxcols)
	for (ix in 1:length(xfuns)) {
		xfun <- xfuns[[ix]]
		xstring <- xstrings[[ix]]
		columns2x <- columns2xs[[ix]]
		ind_xcols <- inds_xcols[[ix]]
		if (!is.null(columns2x)) {
			ptm <- startTimedMessage('Calculating ',xstring,' ...')
			xcols <- as.matrix(data[,columns2x])
			collapsed.xcols <- matrix(NA, nrow=numcollapsedbins, ncol=length(columns2x))
			icount <- 1
			i1_lasttrue <- 1
			for (i1 in 1:length(compare)) {
				if (compare[i1]==TRUE) {
					if (length(columns2x)==1) {
						collapsed.xcols[icount-1] <- xfun(xcols[i1_lasttrue:(i1-1),])
					} else if (length(columns2x) > 1) {
						if (i1_lasttrue==i1-1 | i1==1) {
							collapsed.xcols[icount-1,] <- as.numeric(xcols[i1_lasttrue,])
						} else {
							collapsed.xcols[icount-1,] <- apply(xcols[i1_lasttrue:(i1-1),],2,xfun)
						}
					}
					icount <- icount+1
					i1_lasttrue <- i1
				}
			}
			i1 = i1+1
			if (length(columns2x)==1) {
				collapsed.xcols[icount-1] <- xfun(xcols[i1_lasttrue:(i1-1),])
			} else if (length(columns2x) > 1) {
				if (i1_lasttrue==i1-1 | i1==1) {
					collapsed.xcols[icount-1,] <- as.numeric(xcols[i1_lasttrue,])
				} else {
					collapsed.xcols[icount-1,] <- apply(xcols[i1_lasttrue:(i1-1),],2,xfun)
				}
			}
			if (length(ind_xcols) > 0) {
				lcb <- length(collapsed.bins)
				lsc <- length(ind_xcols)
				collapsed.bins[(lcb+1):(lcb+lsc)] <- as.data.frame(collapsed.xcols)
				names(collapsed.bins)[(lcb+1):(lcb+lsc)] <- paste(xstring, names(data)[ind_xcols], sep='.')
			}
			stopTimedMessage(ptm)
		}
	}

	return(as.data.frame(collapsed.bins))

}

