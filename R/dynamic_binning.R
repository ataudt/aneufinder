dynamic.binning <- function(file, index, chromosomes=NULL, min.mapq=10, win_length=20000) {

	library(GenomicRanges)
	library(Rsamtools)
	file.header <- Rsamtools::scanBamHeader(file)[[1]]
	chrom.lengths <- file.header$targets
	chroms.in.data <- names(chrom.lengths)
	if (is.null(chromosomes)) {
		chromosomes <- chroms.in.data
	}

	chroms2use <- intersect(chromosomes, chroms.in.data)
	gr <- GenomicRanges::GRanges(seqnames=S4Vectors::Rle(chroms2use), ranges=IRanges::IRanges(start=rep(1, length(chroms2use)), end=chrom.lengths[chroms2use]))

	message("Reading data")
	data <- GenomicAlignments::readGAlignments(file, index=index, param=ScanBamParam(which=range(gr), what='mapq', flag=scanBamFlag(isDuplicate=F)))

	if (!is.null(min.mapq)) {
		data <- data[mcols(data)$mapq >= min.mapq]
	}

	ranges.map <- GenomicRanges::GRangesList()
	for (chromosome in chroms2use) {
		message(chromosome)
		chr_data <- data[seqnames(data) == chromosome]
		ranges <- ranges(chr_data)
		reduced.ranges <- reduce(ranges)
		mappable.pos <- coverage(reduced.ranges)
		map.pos.cs <- cumsum(mappable.pos)

		bins.starts <- seq(from = 1, to = max(map.pos.cs)-win_length, by = win_length)
		bins.ends <- seq(from = win_length, to = max(map.pos.cs), by = win_length)
		bins <- c(rbind(bins.starts, bins.ends))

		bins.pos <- findInterval(bins, as.vector(map.pos.cs), rightmost.closed = T)

		## from vector to data frame
		df.bins <- reformat(bins.pos) 

		numbin <- nrow(df.bins)
		chrom.binned.data <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(chromosome, numbin),
							ranges = IRanges::IRanges(start = df.bins$start, end = df.bins$end),
							strand = S4Vectors::Rle(strand("*"), numbin)
							)

		ranges.map[[chromosome]] <- chrom.binned.data
	}
	seqlengths(ranges.map) <- seqlengths(data)
	return(ranges.map)
}

reformat <- function(x) {
    out_list <- list() 
    for ( i in seq(1, length(x), 2) ) {
        out_list[[i]] <- c(x[i], x[i+1])
    }
    mt <- do.call("rbind",out_list)
    df <- data.frame(mt)
    colnames(df) <- c("start", "end")
    df
}
