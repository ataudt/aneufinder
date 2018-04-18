

#' Simulate reads from genome
#' 
#' Simulate single or paired end reads from any \pkg{\link[BSgenome]{BSgenome-class}} object. These simulated reads can be mapped to the reference genome using any aligner to produce BAM files that can be used for mappability correction.
#' 
#' Reads are simulated by splitting the genome into reads with the specified \code{readLength}.
#' 
#' @param bsgenome A \pkg{\link[BSgenome]{BSgenome-class}} object containing the sequence of the reference genome.
#' @param readLength The length in base pairs of the simulated reads that are written to file.
#' @param bamfile A BAM file. This file is used to estimate the distribution of Phred quality scores.
#' @param file The filename that is written to disk. The ending .fastq.gz will be appended.
#' @param pairedEndFragmentLength If this option is specified, paired end reads with length \code{readLength} will be simulated coming from both ends of fragments of this size. NOT IMPLEMENTED YET.
#' @param every.X.bp Stepsize for simulating reads. A read fragment will be simulated every X bp.
#' @return A fastq.gz file is written to disk.
#' @author Aaron Taudt
#' @importFrom Biostrings DNAStringSet BString BStringSet
#' @export
#'@examples
#'## Get an example BAM file with single-cell-sequencing reads
#'bamfile <- system.file("extdata", "BB150803_IV_074.bam", package="AneuFinderData")
#'## Simulate 51bp reads for at a distance of every 5000bp
#'if (require(BSgenome.Mmusculus.UCSC.mm10)) {
#'simulateReads(BSgenome.Mmusculus.UCSC.mm10, bamfile=bamfile, readLength=51,
#'              file=tempfile(), every.X.bp=5000)
#'}
#'
simulateReads <- function(bsgenome, readLength, bamfile, file, pairedEndFragmentLength=NULL, every.X.bp=500) {
	
	if (!is.null(pairedEndFragmentLength)) {
		pairedEndReads <- TRUE
	} else {
		pairedEndReads <- FALSE
	}
	
	## Estimate quality score distribution
	ptm <- startTimedMessage("Estimating quality score distribution ...")
	gr <- suppressMessages( bam2GRanges(bamfile, pairedEndReads=pairedEndReads, remove.duplicate.reads=FALSE, min.mapq=NA, what=c('mapq','seq','qual')) )
	seqbp <- unlist(strsplit(as.character(gr$seq),''))
	mapbp <- unlist(strsplit(as.character(gr$qual),''))
	qualdistr <- table(mapbp, seqbp)
	qualdistr <- sweep(qualdistr, MARGIN=2, STATS=colSums(qualdistr), FUN='/')
	remove(seqbp, mapbp)
	stopTimedMessage(ptm)
	
	bases <- colnames(qualdistr)
	missing.bases <- setdiff(bases, c('A','T','C','G','N'))
	if (length(missing.bases)>0) {
		stop("'bamfile' does not contain bases ", paste(missing.bases, collapse=', '))
	}
	
	## Simulate the reads
	file.gz <- gzfile(paste0(file,'.fastq.gz'), 'w')
	cat('', file=file.gz)
	identifier.string <- attributes(bsgenome)$pkgname
	itotal <- 1
	for (chrom in seqlevels(bsgenome)) {
		message("Simulating reads for chromosome ",chrom)
		## Simulate the reads
		ptm <- startTimedMessage("  Simulating reads ...")
		starts <- seq(1, seqlengths(bsgenome)[chrom]-readLength+1L, by=every.X.bp)
		reads <- Biostrings::DNAStringSet(bsgenome[[chrom]], start=starts, width=readLength)
		mask.N <- reads == paste0(rep('N', readLength), collapse='')
		reads <- reads[!mask.N]
		stopTimedMessage(ptm)
		## Constructing quality fields
		ptm <- startTimedMessage("  Constructing quality fields ...")
		reads.full <- unlist(reads)
		reads.split <- Biostrings::DNAStringSet(reads.full, start=1:length(reads.full), width=1)
		quality.char <- character(length=length(reads.split))
		for (base in bases) {
			idx <- which(reads.split == base)
			quality.char[idx] <- sample(rownames(qualdistr), size=length(idx), replace=TRUE, prob=qualdistr[,base])
		}
		quality <- Biostrings::BString(paste0(quality.char, collapse=''))
		starts <- seq(1, length(reads.split)-readLength+1L, by=every.X.bp)
		quality <- Biostrings::BStringSet(quality, start=starts, width=readLength)
		remove(reads.full, reads.split)
		stopTimedMessage(ptm)
		## Construct ID fields
		ptm <- startTimedMessage("  Constructing ID fields ...")
		id.strings <- paste0(identifier.string, '_', itotal:(itotal+length(reads)-1), ' simulated length=',readLength)
		itotal <- itotal + length(reads)
		stopTimedMessage(ptm)
		## Rearrange into @ID\nSeq\n+ID\nQual
		ptm <- startTimedMessage("  Rearranging for writing ...")
		data <- paste0('@',id.strings,'\n', as.character(reads),'\n', '+',id.strings,'\n', as.character(quality),'\n')
		stopTimedMessage(ptm)
		ptm <- startTimedMessage("  Writing to file ...")
		cat(data, file=file.gz, sep='', append=TRUE)
		stopTimedMessage(ptm)
	}
	close(file.gz)

}
