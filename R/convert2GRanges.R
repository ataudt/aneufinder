hmm2GRanges <- function(hmm, reduce=TRUE) {

# 	library(GenomicRanges)
	### Check user input ###
	if (check.univariate.model(hmm)!=0) stop("argument 'hmm' expects a univariate hmm object (type ?uni.hmm for help)")
	if (check.logical(reduce)!=0) stop("argument 'reduce' expects TRUE or FALSE")

	### Create GRanges ###
	# Transfer coordinates
	gr <- GenomicRanges::GRanges(
			seqnames = Rle(hmm$coordinates$chrom),
			ranges = IRanges(start=hmm$coordinates$start, end=hmm$coordinates$end),
			strand = Rle(strand("*"), nrow(hmm$coordinates))
			)
	seqlengths(gr) <- hmm$seqlengths[names(seqlengths(gr))]
	# Reorder seqlevels
	gr <- GenomicRanges::keepSeqlevels(gr, names(hmm$seqlengths))

	if (reduce) {
		# Reduce state by state
		red.gr.list <- GenomicRanges::GRangesList()
		for (state in unique(hmm$states)) {
			red.gr <- GenomicRanges::reduce(gr[hmm$states==state])
			mcols(red.gr)$states <- rep(as.factor(state),length(red.gr))
			red.gr.list[[length(red.gr.list)+1]] <- red.gr
		}
		# Merge and sort
		red.gr <- GenomicRanges::sort(GenomicRanges::unlist(red.gr.list))
		remove(red.gr.list)
		return(red.gr)
	} else {
		mcols(gr)$reads <- hmm$reads
		mcols(gr)$posteriors <- hmm$posteriors
		mcols(gr)$states <- hmm$states
		return(gr)
	}

}

