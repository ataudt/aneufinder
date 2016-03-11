

#' Wrapper function for the \code{\link{AneuFinder}} package
#'
#' This function is an easy-to-use wrapper to \link[AneuFinder:binning]{bin the data}, \link[AneuFinder:findCNVs]{find copy-number-variations}, \link[AneuFinder:findSCEs]{find sister-chromatid-exchange} events, plot \link[AneuFinder:heatmapGenomewide]{genomewide heatmaps}, \link[AneuFinder:plot.aneuHMM]{distributions, profiles and karyograms}.
#'
#' @param inputfolder Folder with either BAM or BED files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param format Either 'bam' or 'bed', depending if your \code{inputfolder} contains files in BAM or BED format.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing files in \code{outputfolder} should be reused.
#' @param stepsize Fraction of the binsize that the sliding window is offset at each step. Example: If \code{stepsize=0.1} and \code{binsizes=c(200000,500000)}, the actual stepsize in basepairs is 20000 and 50000, respectively.
#' @inheritParams bam2GRanges
#' @inheritParams bed2GRanges
#' @inheritParams binReads
#' @param correction.method Correction methods to be used for the binned read counts. Currently any combination of \code{c('GC','mappability')}.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
#' @param mappability.reference A file that serves as reference for mappability correction. Has to be the same format as specified by \code{format}.
#' @param method Any combination of \code{c('univariate','bivariate')}. Option \code{'univariate'} treats both strands as one, while option \code{'bivariate'} treats both strands separately. NOTE: SCEs can only be called when \code{method='bivariate'}.
#' @inheritParams univariate.findCNVs
#' @param most.frequent.state.univariate One of the states that were given in \code{states}. The specified state is assumed to be the most frequent one when running the univariate HMM. This can help the fitting procedure to converge into the correct fit. Default is '2-somy'.
#' @param most.frequent.state.bivariate One of the states that were given in \code{states}. The specified state is assumed to be the most frequent one when running the bivariate HMM. This can help the fitting procedure to converge into the correct fit. Default is '1-somy'.
#' @inheritParams getSCEcoordinates
#' @param bw Bandwidth for SCE hotspot detection (see \code{\link{hotspotter}} for further details).
#' @param pval P-value for SCE hotspot detection (see \code{\link{hotspotter}} for further details).
#' @param cluster.plots A logical indicating whether plots should be clustered by similarity.
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom utils read.table
#' @importFrom cowplot plot_grid
#' @export
#'
#'@examples
#'\dontrun{
#'## The following call produces plots and genome browser files for all BAM files in "my-data-folder"
#'Aneufinder(inputfolder="my-data-folder", outputfolder="my-output-folder", format='bam')}
#'
Aneufinder <- function(inputfolder, outputfolder, format, configfile=NULL, numCPU=1, reuse.existing.files=TRUE, binsizes=1e6, variable.width.reference=NULL, reads.per.bin=NULL, pairedEndReads=FALSE, stepsize=NULL, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, blacklist=NULL, correction.method=NULL, GC.BSgenome=NULL, mappability.reference=NULL, method='univariate', eps=0.1, max.time=60, max.iter=5000, num.trials=15, states=c('zero-inflation',paste0(0:10,'-somy')), most.frequent.state.univariate='2-somy', most.frequent.state.bivariate='1-somy', resolution=c(3,6), min.segwidth=2, min.reads=50, bw=4*binsizes[1], pval=1e-8, cluster.plots=TRUE) {

#=======================
### Helper functions ###
#=======================
as.object <- function(x) {
	return(eval(parse(text=x)))
}

#========================
### General variables ###
#========================
conf <- NULL
if (is.character(configfile)) {
	## Read config file ##
	errstring <- tryCatch({
		conf <- readConfig(configfile)
		errstring <- ''
	}, error = function(err) {
		errstring <- paste0("Could not read configuration file ",configfile)
	})
	if (errstring!='') {
		stop(errstring)
	}
}
total.time <- proc.time()

## Convert GC.BSgenome to string if necessary
if (class(GC.BSgenome)=='BSgenome') {
	GC.BSgenome <- attributes(GC.BSgenome)$pkgname
}

## Put options into list and merge with conf
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, binsizes=binsizes, variable.width.reference=variable.width.reference, reads.per.bin=reads.per.bin, pairedEndReads=pairedEndReads, stepsize=stepsize, format=format, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, blacklist=blacklist, correction.method=correction.method, GC.BSgenome=GC.BSgenome, mappability.reference=mappability.reference, method=method, eps=eps, max.time=max.time, max.iter=max.iter, num.trials=num.trials, states=states, most.frequent.state.univariate=most.frequent.state.univariate, most.frequent.state.bivariate=most.frequent.state.bivariate, resolution=resolution, min.segwidth=min.segwidth, min.reads=min.reads, bw=bw, pval=pval, cluster.plots=cluster.plots)
conf <- c(conf, params[setdiff(names(params),names(conf))])

## Input checks
if (! conf[['format']] %in% c('bam','bed')) {
	stop("Unknown file format ",conf[['format']])
}
if (conf[['format']] == 'bed' & is.null(conf[['assembly']])) {
	stop("Please specify 'assembly' if format=\"bed\"")
}

## Helpers
binsizes <- conf[['binsizes']]
reads.per.bins <- conf[['reads.per.bin']]
patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=TRUE, trim=TRUE),'_'))
patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
pattern <- NULL #ease R CMD check
numcpu <- conf[['numCPU']]

## Set up the directory structure ##
readspath <- file.path(outputfolder,'data')
readsbrowserpath <- file.path(outputfolder,'browserfiles_data')
binpath.uncorrected <- file.path(outputfolder,'binned')
CNVpath <- file.path(outputfolder,'hmms')
CNVplotpath <- file.path(outputfolder,'plots')
CNVbrowserpath <- file.path(outputfolder,'browserfiles')
SCEpath <- file.path(outputfolder,'hmms_bivariate')
SCEplotpath <- file.path(outputfolder,'plots_bivariate')
SCEbrowserpath <- file.path(outputfolder,'browserfiles_bivariate')
## Delete old directory if desired ##
if (conf[['reuse.existing.files']]==FALSE) {
	if (file.exists(outputfolder)) {
		message("Deleting old directory ",outputfolder)
		unlink(outputfolder, recursive=TRUE)
	}
}
if (!file.exists(outputfolder)) {
	dir.create(outputfolder)
}
## Make a copy of the conf file
writeConfig(conf, configfile=file.path(outputfolder, 'AneuFinder.config'))

## Parallelization ##
if (numcpu > 1) {
	ptm <- startTimedMessage("Setting up parallel execution with ", numcpu, " CPUs ...")
	cl <- parallel::makeCluster(numcpu)
	doParallel::registerDoParallel(cl)
	on.exit(
		if (conf[['numCPU']] > 1) {
			parallel::stopCluster(cl)
		}
	)
	stopTimedMessage(ptm)
}


#==============
### Binning ###
#==============
if (!file.exists(binpath.uncorrected)) { dir.create(binpath.uncorrected) }
files <- list.files(inputfolder, full.names=TRUE, recursive=TRUE)
chrom.lengths.df <- NULL
if (conf[['format']]=='bam') {
	files <- grep('\\.bam$', files, value=TRUE)
} else if (conf[['format']]=='bed') {
	files <- grep('\\.bed$|\\.bed\\.gz$', files, value=TRUE)
	## Get chromosome lengths
	ptm <- startTimedMessage("Fetching chromosome lengths from UCSC ...")
	chrom.lengths.df <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(conf[['assembly']])
	stopTimedMessage(ptm)
	## Determining chromosome format
	df <- utils::read.table(files[1], header=FALSE, nrows=5)
	if (grepl('^chr',as.character(df[nrow(df),1]))) {
		chromosome.format <- 'UCSC'
	} else {
		chromosome.format <- 'NCBI'
	}
	remove(df)
} else {
	stop("Unknown format ", conf[['format']])
}

### Make bins ###
message("==> Making bins:")
if (!is.null(conf[['variable.width.reference']])) {
	if (conf[['format']] == 'bam') {
		reads <- bam2GRanges(conf[['variable.width.reference']], chromosomes=conf[['chromosomes']], pairedEndReads=conf[['pairedEndReads']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']])
	} else if (conf[['format']]=='bed') {
		reads <- bed2GRanges(conf[['variable.width.reference']], assembly=chrom.lengths.df, chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']])
	}
	bins <- variableWidthBins(reads, binsizes=conf[['binsizes']], chromosomes=conf[['chromosomes']])
} else {
	if (conf[['format']] == 'bam') {
		bins <- fixedWidthBins(bamfile=files[1], binsizes=conf[['binsizes']], chromosomes=conf[['chromosomes']])
	} else if (conf[['format']] == 'bed') {
		bins <- fixedWidthBins(assembly=chrom.lengths.df, chromosome.format=chromosome.format, binsizes=conf[['binsizes']], chromosomes=conf[['chromosomes']])
	}
}
message("==| Finished making bins.")

### Binning ###
parallel.helper <- function(file) {
	existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value=TRUE)
	existing.binsizes <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 2)))
	existing.rpbin <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 3)))
	binsizes.todo <- setdiff(binsizes, existing.binsizes)
	rpbin.todo <- setdiff(reads.per.bins, existing.rpbin)
	if (length(c(binsizes.todo,rpbin.todo)) > 0) {
		tC <- tryCatch({
			binReads(file=file, format=conf[['format']], assembly=chrom.lengths.df, pairedEndReads=conf[['pairedEndReads']], binsizes=NULL, variable.width.reference=NULL, reads.per.bin=rpbin.todo, bins=bins[as.character(binsizes.todo)], stepsize=conf[['stepsize']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']], outputfolder.binned=binpath.uncorrected, save.as.RData=TRUE, reads.store=TRUE, outputfolder.reads=readspath)
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
}
if (numcpu > 1) {
	ptm <- startTimedMessage("Binning the data ...")
	temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
		parallel.helper(file)
	}
	stopTimedMessage(ptm)
} else {
	temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
		parallel.helper(file)
	}
}
	
### Read fragments that are not produced yet ###
parallel.helper <- function(file) {
	savename <- file.path(readspath,paste0(basename(file),'.RData'))
	if (!file.exists(savename)) {
		tC <- tryCatch({
			binReads(file=file, format=conf[['format']], assembly=chrom.lengths.df, pairedEndReads=conf[['pairedEndReads']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']], calc.complexity=FALSE, reads.store=TRUE, outputfolder.reads=readspath, reads.only=TRUE)
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
}

if (numcpu > 1) {
	ptm <- startTimedMessage("Saving reads as .RData ...")
	temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
		parallel.helper(file)
	}
	stopTimedMessage(ptm)
} else {
	temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
		parallel.helper(file)
	}
}

### Export read fragments as browser file ###
if (!file.exists(readsbrowserpath)) { dir.create(readsbrowserpath) }
readfiles <- list.files(readspath,pattern='.RData$',full.names=TRUE)

parallel.helper <- function(file) {
	savename <- file.path(readsbrowserpath,sub('.RData','',basename(file)))
	if (!file.exists(paste0(savename,'.bed.gz'))) {
		tC <- tryCatch({
			gr <- loadGRangesFromFiles(file)[[1]]
			exportGRanges(gr, filename=savename, trackname=basename(savename), score=gr$mapq)
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
}

if (numcpu > 1) {
	ptm <- startTimedMessage("Exporting data as browser files ...")
	temp <- foreach (file = readfiles, .packages=c("AneuFinder")) %dopar% {
		parallel.helper(file)
	}
	stopTimedMessage(ptm)
} else {
	temp <- foreach (file = readfiles, .packages=c("AneuFinder")) %do% {
		parallel.helper(file)
	}
}

#=================
### Correction ###
#=================
if (!is.null(conf[['correction.method']])) {

	binpath.corrected <- binpath.uncorrected
	for (correction.method in conf[['correction.method']]) {
		binpath.corrected <- paste0(binpath.corrected, '-', correction.method)
		if (!file.exists(binpath.corrected)) { dir.create(binpath.corrected) }

		if (correction.method=='GC') {
			## Load BSgenome
			if (class(conf[['GC.BSgenome']])!='BSgenome') {
				if (is.character(conf[['GC.BSgenome']])) {
					suppressPackageStartupMessages(library(conf[['GC.BSgenome']], character.only=TRUE))
					conf[['GC.BSgenome']] <- as.object(conf[['GC.BSgenome']]) # replacing string by object
				}
			}

			## Go through patterns
			parallel.helper <- function(pattern) {
				binfiles <- list.files(binpath.uncorrected, pattern='RData$', full.names=TRUE)
				binfiles <- grep(gsub('\\+','\\\\+',pattern), binfiles, value=TRUE)
				binfiles.corrected <- list.files(binpath.corrected, pattern='RData$', full.names=TRUE)
				binfiles.corrected <- grep(gsub('\\+','\\\\+',pattern), binfiles.corrected, value=TRUE)
				binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
				if (length(binfiles.todo)>0) {
					binfiles.todo <- paste0(binpath.uncorrected,.Platform$file.sep,binfiles.todo)
					if (grepl('binsize',gsub('\\+','\\\\+',pattern))) {
						binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf[['GC.BSgenome']], same.binsize=TRUE))
					} else {
						binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf[['GC.BSgenome']], same.binsize=FALSE))
					}
					for (i1 in 1:length(binned.data.list)) {
						binned.data <- binned.data.list[[i1]]
						savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
						save(binned.data, file=savename)
					}
				}
			}
			if (numcpu > 1) {
				ptm <- startTimedMessage(paste0(correction.method," correction ..."))
				temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
					parallel.helper(pattern)
				}
				stopTimedMessage(ptm)
			} else {
				temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
					parallel.helper(pattern)
				}
			}
		}

		if (correction.method=='mappability') {

			## Go through patterns
			parallel.helper <- function(pattern) {
				binfiles <- list.files(binpath.uncorrected, pattern='RData$', full.names=TRUE)
				binfiles <- grep(gsub('\\+','\\\\+',pattern), binfiles, value=TRUE)
				binfiles.corrected <- list.files(binpath.corrected, pattern='RData$', full.names=TRUE)
				binfiles.corrected <- grep(gsub('\\+','\\\\+',pattern), binfiles.corrected, value=TRUE)
				binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
				if (length(binfiles.todo)>0) {
					binfiles.todo <- paste0(binpath.uncorrected,.Platform$file.sep,binfiles.todo)
					if (grepl('binsize',gsub('\\+','\\\\+',pattern))) {
						binned.data.list <- suppressMessages(correctMappability(binfiles.todo, reference=conf[['mappability.reference']], format=conf[['format']], assembly=chrom.lengths.df, pairedEndReads = conf[['pairedEndReads']], min.mapq = conf[['min.mapq']], remove.duplicate.reads = conf[['remove.duplicate.reads']], same.binsize=TRUE))
					} else {
						binned.data.list <- suppressMessages(correctMappability(binfiles.todo, reference=conf[['mappability.reference']], format=conf[['format']], assembly=chrom.lengths.df, pairedEndReads = conf[['pairedEndReads']], min.mapq = conf[['min.mapq']], remove.duplicate.reads = conf[['remove.duplicate.reads']], same.binsize=FALSE))
					}
					for (i1 in 1:length(binned.data.list)) {
						binned.data <- binned.data.list[[i1]]
						savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
						save(binned.data, file=savename)
					}
				}
			}
			if (numcpu > 1) {
				ptm <- startTimedMessage(paste0(correction.method," correction ..."))
				temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
					parallel.helper(pattern)
				}
				stopTimedMessage(ptm)
			} else {
				temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
					parallel.helper(pattern)
				}
			}
		}

	}
	binpath <- binpath.corrected

} else {
	binpath <- binpath.uncorrected
}

#===============
### findCNVs ###
#===============
if ('univariate' %in% conf[['method']]) {

	if (!file.exists(CNVpath)) { dir.create(CNVpath) }

	files <- list.files(binpath, full.names=TRUE, recursive=TRUE, pattern='.RData$')

	parallel.helper <- function(file) {
		tC <- tryCatch({
			savename <- file.path(CNVpath,basename(file))
			if (!file.exists(savename)) {
				model <- findCNVs(file, eps=conf[['eps']], max.time=conf[['max.time']], max.iter=conf[['max.iter']], num.trials=conf[['num.trials']], states=conf[['states']], most.frequent.state=conf[['most.frequent.state.univariate']]) 
				save(model, file=savename)
			}
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Running univariate HMMs ...")
		temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(file)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
			parallel.helper(file)
		}
	}

	#===================
	### Plotting CNV ###
	#===================
	if (!file.exists(CNVplotpath)) { dir.create(CNVplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=TRUE, trim=TRUE),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(CNVpath, full.names=TRUE, recursive=TRUE, pattern='.RData$')

	#------------------
	## Plot heatmaps ##
	#------------------
	parallel.helper <- function(pattern) {
		ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
		if (length(ifiles)>0) {
			savename=file.path(CNVplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=FALSE, cluster=conf[['cluster.plots']]))
			}
		} else {
			warning("Plotting genomewide heatmaps: No files for pattern ",pattern," found.")
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Plotting genomewide heatmaps ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}
	parallel.helper <- function(pattern) {
		ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
		if (length(ifiles)>0) {
			savename=file.path(CNVplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf[['cluster.plots']]))
				grDevices::pdf(savename, width=30, height=0.3*length(ifiles))
				print(ggplt)
				d <- grDevices::dev.off()
			}
		} else {
			warning("Plotting chromosome heatmaps: No files for pattern ",pattern," found.")
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Plotting chromosome heatmaps ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	#------------------------------------
	## Plot profiles and distributions ##
	#------------------------------------
	parallel.helper <- function(pattern) {
		savename <- file.path(CNVplotpath,paste0('profiles_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			grDevices::pdf(file=savename, width=20, height=10)
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					p1 <- graphics::plot(model, type='profile')
					p2 <- graphics::plot(model, type='histogram')
					cowplt <- cowplot::plot_grid(p1, p2, nrow=2, rel_heights=c(1.2,1))
					print(cowplt)
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- grDevices::dev.off()
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Making profile and distribution plots ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	#-------------------------
	## Export browser files ##
	#-------------------------
	if (!file.exists(CNVbrowserpath)) { dir.create(CNVbrowserpath) }
	parallel.helper <- function(pattern) {
		savename <- file.path(CNVbrowserpath,sub('_$','',pattern))
		if (!file.exists(paste0(savename,'_CNV.bed.gz'))) {
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
			exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.SCE=FALSE)
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Exporting browser files ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

}

#===============
### findSCEs ###
#===============
if ('bivariate' %in% conf[['method']]) {

	if (!file.exists(SCEpath)) { dir.create(SCEpath) }

	files <- list.files(binpath, full.names=TRUE, recursive=TRUE, pattern='.RData$')
	parallel.helper <- function(file) {
		tC <- tryCatch({
			savename <- file.path(SCEpath,basename(file))
			if (!file.exists(savename)) {
				model <- findSCEs(file, eps=conf[['eps']], max.time=conf[['max.time']], max.iter=conf[['max.iter']], num.trials=conf[['num.trials']], states=conf[['states']], most.frequent.state=conf[['most.frequent.state.bivariate']]) 
				## Add SCE coordinates to model
				reads.file <- file.path(readspath, paste0(model$ID,'.RData'))
				model$sce <- suppressMessages( getSCEcoordinates(model, resolution=conf[['resolution']], min.segwidth=conf[['min.segwidth']], fragments=reads.file, min.reads=conf[['min.reads']]) )
				ptm <- startTimedMessage("Saving to file ",savename," ...")
				save(model, file=savename)
				stopTimedMessage(ptm)
# 			} else {
# 				model <- get(load(savename))
# 				model$sce <- suppressMessages( getSCEcoordinates(model, resolution=conf[['resolution']], min.segwidth=conf[['min.segwidth']]) )
# 				ptm <- startTimedMessage("Saving to file ",savename," ...")
# 				save(model, file=savename)
# 				stopTimedMessage(ptm)
			}
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Running bivariate HMMs ...")
		temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(file)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
			parallel.helper(file)
		}
	}

	### Finding hotspots ###
	parallel.helper <- function(pattern) {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
		sces <- list()
		for (file in ifiles) {
			hmm <- loadHmmsFromFiles(file)[[1]]
			sces[[file]] <- hmm$sce
		}
		hotspot <- hotspotter(sces, bw=conf[['bw']], pval=conf[['pval']])
		return(hotspot)
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Finding SCE hotspots ...")
		hotspots <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		hotspots <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}
	names(hotspots) <- patterns

	#===================
	### Plotting SCE ###
	#===================
	if (!file.exists(SCEplotpath)) { dir.create(SCEplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=TRUE, trim=TRUE),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(SCEpath, full.names=TRUE, recursive=TRUE, pattern='.RData$')

	#------------------
	## Plot heatmaps ##
	#------------------
	parallel.helper <- function(pattern) {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
		if (length(ifiles)>0) {
			savename=file.path(SCEplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=TRUE, hotspots=hotspots[[pattern]], cluster=conf[['cluster.plots']]))
			}
		} else {
			warning("Plotting genomewide heatmaps: No files for pattern ",pattern," found.")
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Plotting genomewide heatmaps ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	parallel.helper <- function(pattern) {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
		if (length(ifiles)>0) {
			savename=file.path(SCEplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				grDevices::pdf(savename, width=30, height=0.3*length(ifiles))
				ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf[['cluster.plots']]))
				print(ggplt)
				d <- grDevices::dev.off()
			}
		} else {
			warning("Plotting chromosome heatmaps: No files for pattern ",pattern," found.")
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Plotting chromosome heatmaps ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	#------------------
	## Plot profiles ##
	#------------------
	parallel.helper <- function(pattern) {
		savename <- file.path(SCEplotpath,paste0('profiles_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			grDevices::pdf(file=savename, width=20, height=10)
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					p1 <- graphics::plot(model, type='profile')
					p2 <- graphics::plot(model, type='histogram')
					cowplt <- cowplot::plot_grid(p1, p2, nrow=2, rel_heights=c(1.2,1))
					print(cowplt)
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- grDevices::dev.off()
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Making profile and distribution plots ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	#--------------------
	## Plot karyograms ##
	#--------------------
	parallel.helper <- function(pattern) {
		savename <- file.path(SCEplotpath,paste0('karyograms_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			grDevices::pdf(file=savename, width=12*1.4, height=2*4.6)
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(graphics::plot(model, type='karyogram'))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- grDevices::dev.off()
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Plotting karyograms ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

	#-------------------------
	## Export browser files ##
	#-------------------------
	if (!file.exists(SCEbrowserpath)) { dir.create(SCEbrowserpath) }
	parallel.helper <- function(pattern) {
		savename <- file.path(SCEbrowserpath,sub('_$','',pattern))
		if (!file.exists(paste0(savename,'_CNV.bed.gz'))) {
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
			exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.SCE=TRUE)
		}
		savename <- file.path(SCEbrowserpath,paste0(pattern,'SCE-hotspots'))
		if (!file.exists(paste0(savename,'.bed.gz'))) {
			exportGRanges(hotspots[[pattern]], filename=savename, trackname=basename(savename), score=hotspots[[pattern]]$num.events)
		}
	}
	if (numcpu > 1) {
		ptm <- startTimedMessage("Exporting browser files ...")
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
			parallel.helper(pattern)
		}
		stopTimedMessage(ptm)
	} else {
		temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
			parallel.helper(pattern)
		}
	}

}

total.time <- proc.time() - total.time
message("==> Total time spent: ", round(total.time[3]), "s <==")

}
