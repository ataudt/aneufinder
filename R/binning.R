bedGraph2bin <- function(bedGraphfile, chrom.length.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE) {
	return(align2bin(bedGraphfile, format="bedGraph", chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bam2bin <- function(bamfile, bamindex=bamfile, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE) {
	return(align2bin(bamfile, format="bam", index=bamindex, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

bed2bin <- function(bedfile, chrom.length.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE) {
	return(align2bin(bedfile, format="bed", chrom.length.file=chrom.length.file, outputfolder=outputfolder, binsizes=binsizes, chromosomes=chromosomes, separate.chroms=separate.chroms, save.as.RData=save.as.RData))
}

align2bin <- function(file, format, index=file, chrom.length.file, outputfolder="binned_data", binsizes=200, chromosomes=NULL, separate.chroms=TRUE, save.as.RData=TRUE) {

	## Load libraries
# 	library(GenomicRanges)

	## Create outputfolder if not exists
	if (!file.exists(outputfolder) & save.as.RData==TRUE) {
		dir.create(outputfolder)
	}

	## Read in the data
	if (format == "bed") {
		cat("Reading file",basename(file),"...")
		# File with chromosome lengths
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		# File with reads, determine classes first for faster import
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:3] <- classes.in.bed[1:3]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2], end=data[,3]), strand=Rle(strand("*"), nrow(data)))
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		chroms.in.data <- seqlevels(data)
	} else if (format == "bam") {
		library(Rsamtools)
		cat("Reading header of",basename(file),"...")
		file.header <- scanBamHeader(file)[[1]]
		chrom.lengths <- file.header$targets
		chroms.in.data <- names(chrom.lengths)
	} else if (format == "bedGraph") {
		cat("Reading file",basename(file),"...")
		# File with chromosome lengths
		chrom.lengths.df <- read.table(chrom.length.file)
		chrom.lengths <- chrom.lengths.df[,2]
		names(chrom.lengths) <- chrom.lengths.df[,1]
		# File with reads, determine classes first for faster import
		tab5rows <- read.table(file, nrows=5)
		classes.in.bed <- sapply(tab5rows, class)
		classes <- rep("NULL",length(classes.in.bed))
		classes[1:4] <- classes.in.bed[1:4]
		data <- read.table(file, colClasses=classes)
		# Convert to GRanges object
		data <- GRanges(seqnames=Rle(data[,1]), ranges=IRanges(start=data[,2], end=data[,3]), strand=Rle(strand("*"), nrow(data)), signal=data[,4])
		seqlengths(data) <- as.integer(chrom.lengths[names(seqlengths(data))])
		chroms.in.data <- seqlevels(data)
	}
	cat(" done\n")

 
	### Do the loop for all binsizes
	for (binsize in binsizes) {
		cat("Binning into binsize",binsize,"\n")

		### Iterate over all chromosomes
		binned.data.allchroms <- NULL
		if (is.null(chromosomes)) {
			chromosomes <- chroms.in.data
		}
		for (chromosome in chromosomes) {
			## Check if chromosome exists in data
			if ( !(chromosome %in% chroms.in.data) ) {
				warning("Skipped chromosome ",chromosome,", not in the data!")
				next
			} else if ( !(chromosome %in% names(chrom.lengths)) ) {
				warning("Skipped chromosome ",chromosome,", no length found!")
				next
			}
			cat(chromosome,"                              \n")
			## Check last incomplete bin
			incomplete.bin <- chrom.lengths[chromosome] %% binsize > 0
			if (incomplete.bin) {
				numbins <- ceiling(chrom.lengths[chromosome]/binsize)
			} else {
				numbins <- chrom.lengths[chromosome]/binsize
			}
			## Initialize vectors
			cat("initialize vectors...\r")
			chroms <- rep(chromosome,numbins)
			reads <- rep(0,numbins)
			start <- seq(from=0, by=binsize, length.out=numbins)
			end <- seq(from=binsize-1, by=binsize, length.out=numbins)
			end[length(end)] <- chrom.lengths[chromosome]

			## Create binned chromosome as GRanges object
			cat("creating GRanges container...            \r")
			ichrom <- GRanges(seqnames = Rle(chromosome, numbins),
							ranges = IRanges(start=start, end=end),
							strand = Rle(strand("*"), numbins)
							)

			if (format=="bam") {
				cat("reading reads from file...               \r")
				data <- readGAlignmentsFromBam(file, index=index, param=ScanBamParam(what=c("pos"),which=range(ichrom)))
			}

			## Count overlaps
			cat("counting overlaps...                     \r")
			if (format=="bam" | format=="bed") {
				reads <- countOverlaps(ichrom, data[seqnames(data)==chromosome])
			} else if (format=="bedGraph") {
				# Take the max value from all regions that fall into / overlap a given bin as read count
				midx <- as.matrix(findOverlaps(ichrom, data[seqnames(data)==chromosome]))
				reads <- rep(0,length(ichrom))
				signal <- mcols(data)$signal
				rle <- rle(midx[,1])
				read.idx <- rle$values
				max.idx <- cumsum(rle$lengths)
				maxvalues <- rep(NA, length(read.idx))
				maxvalues[1] <- max(signal[midx[1:(max.idx[1]),2]])
				for (i1 in 2:length(read.idx)) {
					maxvalues[i1] <- max(signal[midx[(max.idx[i1-1]+1):(max.idx[i1]),2]])
				}
				reads[read.idx] <- maxvalues
			}
			
			## Concatenate
			cat("concatenate...                           \r")
			binned.data <- data.frame(chroms,start,end,reads)
			names(binned.data) <- c("chrom","start","end","reads")

			if (separate.chroms==TRUE) {
				if (save.as.RData==TRUE) {
					## Print to file
					filename <- paste(basename(file),"_binsize_",binsize,"_",chromosome,".RData", sep="")
					cat("save...                                  \r")
					save(binned.data, file=file.path(outputfolder,filename) )
				} else {
					cat("                                         \r")
					return(binned.data)
				}
			} else {
				binned.data.allchroms[[length(binned.data.allchroms)+1]] <- binned.data
			}
			cat("                                         \r")

		}
		if (separate.chroms!=TRUE) {
			cat("Concatenating chromosomes ...")
			binned.data.allchroms <- do.call("rbind",binned.data.allchroms)
			cat(" done\n")
			if (save.as.RData==TRUE) {
				# Print to file
				filename <- paste(basename(file),"_binsize_",binsize,".RData", sep="")
				cat("Saving to file ...")
				save(binned.data.allchroms, file=file.path(outputfolder,filename) )
				cat(" done\n")
			} else {
				return(binned.data.allchroms)
			}
		}

	}

}
