

#' Wrapper function for the aneufinder package
#'
#' This function is an easy-to-use wrapper to \link[aneufinder:binning]{bin the data}, \link[aneufinder:findCNVs]{find copy-number-variations}, \link[aneufinder:findSCEs]{find sister-chromatid-exchange} events and plot the results.
#'
#' @param inputfolder Folder with either BAM or BED files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing files in \code{outputfolder} should be reused.

#' @param stepsize Fraction of the binsize that the sliding window is offset at each step. Example: If \code{stepsize=0.1} and \code{binsizes=c(200000,500000)}, the actual stepsize in basepairs is 20000 and 50000, respectively.
#' @inheritParams align2binned

#' @param correction.method Correction methods to be used for the binned read counts. Currently any combination of \code{c('GC')}.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
#' @param method Any combination of \code{c('univariate','bivariate')}. Option \code{'univariate'} treats both strands as one, while option \code{'bivariate'} treats both strands separately. NOTE: SCEs can only be called when \code{method='bivariate'}.

#' @inheritParams univariate.findCNVs
#' @param most.frequent.state.univariate One of the states that were given in \code{states}. The specified state is assumed to be the most frequent one when running the univariate HMM. This can help the fitting procedure to converge into the correct fit. Default is 'disomy'.
#' @param most.frequent.state.bivariate One of the states that were given in \code{states}. The specified state is assumed to be the most frequent one when running the bivariate HMM. This can help the fitting procedure to converge into the correct fit. Default is 'monosomy'.

#' @inheritParams getSCEcoordinates

#' @param bw Bandwidth for SCE hotspot detection (see \code{\link{hotspotter}} for further details).
#' @param pval P-value for SCE hotspot detection (see \code{\link{hotspotter}} for further details).
#' @param cluster.plots A logical indicating whether plots should be clustered by similarity.
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @export
Aneufinder <- function(inputfolder, outputfolder, configfile=NULL, numCPU=1, reuse.existing.files=TRUE, binsizes=500000, reads.per.bin=NULL, pairedEndReads=FALSE, stepsize=NULL, format='bam', chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, correction.method=NULL, GC.BSgenome=NULL, method='univariate', eps=0.1, max.time=60, max.iter=5000, num.trials=15, states=c('zero-inflation','nullsomy','monosomy','disomy','trisomy','tetrasomy','multisomy'), most.frequent.state.univariate='disomy', most.frequent.state.bivariate='monosomy', resolution=c(3,6), min.segwidth=2, min.reads=50, bw=4*binsizes[1], pval=0.05, cluster.plots=TRUE) {

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

## Convert GC.BSgenome to string if necessary
if (class(GC.BSgenome)=='BSgenome') {
	GC.BSgenome <- attributes(GC.BSgenome)$pkgname
}

## Put options into list and merge with conf
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, binsizes=binsizes, reads.per.bin=reads.per.bin, pairedEndReads=pairedEndReads, stepsize=stepsize, format=format, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, correction.method=correction.method, GC.BSgenome=GC.BSgenome, method=method, eps=eps, max.time=max.time, max.iter=max.iter, num.trials=num.trials, states=states, most.frequent.state.univariate=most.frequent.state.univariate, most.frequent.state.bivariate=most.frequent.state.bivariate, resolution=resolution, min.segwidth=min.segwidth, min.reads=min.reads, bw=bw, pval=pval, cluster.plots=cluster.plots)
conf <- c(conf, params[setdiff(names(params),names(conf))])

## Input checks
if (! conf[['format']] %in% c('bam')) {
	stop("Unknown file format ",conf[['format']])
}

## Helpers
binsizes <- conf[['binsizes']]
reads.per.bins <- conf[['reads.per.bin']]
patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=T, trim=T),'_'))
patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
pattern <- NULL #ease R CMD check

## Set up the directory structure ##
readspath <- file.path(outputfolder,'data')
binpath.uncorrected <- file.path(outputfolder,'binned')
CNVpath <- file.path(outputfolder,'results_univariate')
CNVplotpath <- file.path(outputfolder,'plots_univariate')
CNVbrowserpath <- file.path(outputfolder,'browserfiles_univariate')
SCEpath <- file.path(outputfolder,'results_bivariate')
SCEplotpath <- file.path(outputfolder,'plots_bivariate')
SCEbrowserpath <- file.path(outputfolder,'browserfiles_bivariate')
## Delete old directory if desired ##
if (conf[['reuse.existing.files']]==FALSE) {
	if (file.exists(outputfolder)) {
		message("Deleting old directory ",outputfolder)
		unlink(outputfolder, recursive=T)
	}
}
if (!file.exists(outputfolder)) {
	dir.create(outputfolder)
}
## Make a copy of the conf file
writeConfig(conf, configfile=file.path(outputfolder, 'aneufinder.config'))

## Parallelization ##
message("Using ",conf[['numCPU']]," CPUs")
cl <- makeCluster(conf[['numCPU']])
doParallel::registerDoParallel(cl)


#==============
### Binning ###
#==============
message("Binning the data ...", appendLF=F); ptm <- proc.time()
if (!file.exists(binpath.uncorrected)) { dir.create(binpath.uncorrected) }

files <- list.files(inputfolder, full.names=TRUE, recursive=T, pattern=paste0('.',conf[['format']],'$'))
### Binning ###
temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
	existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value=T)
	existing.binsizes <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 2)))
	existing.rpbin <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 3)))
	binsizes.todo <- setdiff(binsizes, existing.binsizes)
	rpbin.todo <- setdiff(reads.per.bins, existing.rpbin)
	if (length(c(binsizes.todo,rpbin.todo)) > 0) {
		tC <- tryCatch({
			if (conf[['format']]=='bam') {
				bam2binned(bamfile=file, pairedEndReads=conf[['pairedEndReads']], binsizes=binsizes.todo, reads.per.bin=rpbin.todo, stepsize=conf[['stepsize']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], outputfolder.binned=binpath.uncorrected, save.as.RData=TRUE, reads.store=TRUE, outputfolder.reads=readspath)
			}
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
}
### Read fragments that are not produced yet ###
temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
	savename <- file.path(readspath,paste0(basename(file),'.RData'))
	if (!file.exists(savename)) {
		tC <- tryCatch({
			if (conf[['format']]=='bam') {
				bam2binned(bamfile=file, pairedEndReads=conf[['pairedEndReads']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], calc.complexity=FALSE, reads.store=TRUE, outputfolder.reads=readspath, reads.only=TRUE)
			}
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
}
time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

#=================
### Correction ###
#=================
if (!is.null(conf[['correction.method']])) {

	message(paste0(conf[['correction.method']]," correction ..."), appendLF=F); ptm <- proc.time()
	if (conf[['correction.method']]=='GC') {
		binpath.corrected <- file.path(outputfolder,'binned_GC-corrected')
		if (!file.exists(binpath.corrected)) { dir.create(binpath.corrected) }
		## Load BSgenome
		if (class(conf[['GC.BSgenome']])!='BSgenome') {
			if (is.character(conf[['GC.BSgenome']])) {
				suppressPackageStartupMessages(library(conf[['GC.BSgenome']], character.only=T))
				conf[['GC.BSgenome']] <- as.object(conf[['GC.BSgenome']]) # replacing string by object
			}
		}

		## Go through patterns
		temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
			binfiles <- list.files(binpath.uncorrected, pattern='RData$', full.names=TRUE)
			binfiles <- grep(gsub('\\+','\\\\+',pattern), binfiles, value=T)
			binfiles.corrected <- list.files(binpath.corrected, pattern='RData$', full.names=TRUE)
			binfiles.corrected <- grep(gsub('\\+','\\\\+',pattern), binfiles.corrected, value=T)
			binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
			if (length(binfiles.todo)>0) {
				binfiles.todo <- paste0(binpath.uncorrected,.Platform$file.sep,binfiles.todo)
				if (grepl('binsize',gsub('\\+','\\\\+',pattern))) {
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf[['GC.BSgenome']], same.GC.content=TRUE))
				} else {
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf[['GC.BSgenome']]))
				}
				for (i1 in 1:length(binned.data.list)) {
					binned.data <- binned.data.list[[i1]]
					savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
					save(binned.data, file=savename)
				}
			}
		}
		binpath <- binpath.corrected
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

} else {
	binpath <- binpath.uncorrected
}

#===============
### findCNVs ###
#===============
if ('univariate' %in% conf[['method']]) {

	message("Running univariate HMMs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(CNVpath)) { dir.create(CNVpath) }

	files <- list.files(binpath, full.names=TRUE, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
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
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	#===================
	### Plotting CNV ###
	#===================
	if (!file.exists(CNVplotpath)) { dir.create(CNVplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=T, trim=T),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(CNVpath, full.names=TRUE, recursive=T, pattern='.RData$')

	#-----------------------
	## Plot distributions ##
	#-----------------------
	message("Plotting distributions ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVplotpath,paste0('distributions_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=10, height=7)
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='histogram') + theme_bw(base_size=18) + theme(legend.position=c(1,1), legend.justification=c(1,1)))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#------------------
	## Plot heatmaps ##
	#------------------
	message("Plotting genomewide heatmaps ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
		if (length(ifiles)>0) {
			savename=file.path(CNVplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=F, cluster=conf[['cluster.plots']]))
			}
		} else {
			warning("Plotting genomewide heatmaps: No files for pattern ",pattern," found.")
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	message("Plotting chromosome heatmaps ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
		if (length(ifiles)>0) {
			savename=file.path(CNVplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf[['cluster.plots']]))
				pdf(savename, width=30, height=0.3*length(ifiles))
				print(ggplt)
				d <- dev.off()
			}
		} else {
			warning("Plotting chromosome heatmaps: No files for pattern ",pattern," found.")
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#------------------
	## Plot arrayCGH ##
	#------------------
	message("Making arrayCGH plots ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVplotpath,paste0('arrayCGH_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=20, height=5)
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='arrayCGH'))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#--------------------
	## Plot karyograms ##
	#--------------------
	message("Plotting karyograms ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVplotpath,paste0('karyograms_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=12*1.4, height=2*4.6)
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='karyogram'))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#-------------------------
	## Export browser files ##
	#-------------------------
	message("Exporting browser files ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(CNVbrowserpath)) { dir.create(CNVbrowserpath) }
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVbrowserpath,sub('_$','',pattern))
		if (!file.exists(paste0(savename,'.bed.gz'))) {
			ifiles <- list.files(CNVpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.SCE=FALSE)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

#===============
### findSCEs ###
#===============
if ('bivariate' %in% conf[['method']]) {

	message("Running bivariate HMMs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(SCEpath)) { dir.create(SCEpath) }

	files <- list.files(binpath, full.names=TRUE, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
		tC <- tryCatch({
			savename <- file.path(SCEpath,basename(file))
			if (!file.exists(savename)) {
				model <- findSCEs(file, eps=conf[['eps']], max.time=conf[['max.time']], max.iter=conf[['max.iter']], num.trials=conf[['num.trials']], states=conf[['states']], most.frequent.state=conf[['most.frequent.state.bivariate']]) 
				## Add SCE coordinates to model
				reads.file <- file.path(readspath, paste0(model$ID,'.RData'))
				model$sce <- suppressMessages( getSCEcoordinates(model, resolution=conf[['resolution']], min.segwidth=conf[['min.segwidth']], fragments=reads.file, min.reads=conf[['min.reads']]) )
				message("Saving to file ",savename," ...", appendLF=F); ptm <- proc.time()
				save(model, file=savename)
				time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
# 			} else {
# 				model <- get(load(savename))
# 				model$sce <- suppressMessages( getSCEcoordinates(model, resolution=conf[['resolution']], min.segwidth=conf[['min.segwidth']]) )
# 				message("Saving to file ",savename," ...", appendLF=F); ptm <- proc.time()
# 				save(model, file=savename)
# 				time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
			}
		}, error = function(err) {
			stop(file,'\n',err)
		})
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	### Finding hotspots ###
	message("Finding SCE hotspots ...", appendLF=F); ptm <- proc.time()
	hotspots <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
		sces <- list()
		for (file in ifiles) {
			hmm <- loadHmmsFromFiles(file)[[1]]
			sces[[file]] <- hmm$sce
		}
		hotspot <- hotspotter(sces, bw=conf[['bw']], pval=conf[['pval']])
		return(hotspot)
	}
	names(hotspots) <- patterns
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	#===================
	### Plotting SCE ###
	#===================
	if (!file.exists(SCEplotpath)) { dir.create(SCEplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=T, trim=T),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(SCEpath, full.names=TRUE, recursive=T, pattern='.RData$')

	#-----------------------
	## Plot distributions ##
	#-----------------------
	message("Plotting distributions ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEplotpath,paste0('distributions_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=10, height=7)
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='histogram') + theme_bw(base_size=18) + theme(legend.position=c(1,1), legend.justification=c(1,1)))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#------------------
	## Plot heatmaps ##
	#------------------
	message("Plotting heatmaps ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
		if (length(ifiles)>0) {
			savename=file.path(SCEplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=T, hotspots=hotspots[[pattern]], cluster=conf[['cluster.plots']]))
			}
		} else {
			warning("Plotting genomewide heatmaps: No files for pattern ",pattern," found.")
		}
	}
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
		ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
		if (length(ifiles)>0) {
			savename=file.path(SCEplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
			if (!file.exists(savename)) {
				pdf(savename, width=30, height=0.3*length(ifiles))
				ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf[['cluster.plots']]))
				print(ggplt)
				d <- dev.off()
			}
		} else {
			warning("Plotting chromosome heatmaps: No files for pattern ",pattern," found.")
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#------------------
	## Plot arrayCGH ##
	#------------------
	message("Making arrayCGH plots ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEplotpath,paste0('arrayCGH_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=20, height=5)
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='arrayCGH'))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#--------------------
	## Plot karyograms ##
	#--------------------
	message("Plotting karyograms ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEplotpath,paste0('karyograms_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=12*1.4, height=2*4.6)
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			for (ifile in ifiles) {
				tC <- tryCatch({
					model <- get(load(ifile))
					print(plot(model, type='karyogram'))
				}, error = function(err) {
					stop(ifile,'\n',err)
				})
			}
			d <- dev.off()
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	#-------------------------
	## Export browser files ##
	#-------------------------
	message("Exporting browser files ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(SCEbrowserpath)) { dir.create(SCEbrowserpath) }
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEbrowserpath,sub('_$','',pattern))
		if (!file.exists(paste0(savename,'.bed.gz'))) {
			ifiles <- list.files(SCEpath, pattern='RData$', full.names=TRUE)
			ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=T)
			exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.SCE=TRUE)
		}
		savename <- file.path(SCEbrowserpath,paste0(pattern,'SCE-hotspots'))
		if (!file.exists(paste0(savename,'.bed.gz'))) {
			exportGRanges(hotspots[[pattern]], filename=savename, trackname=basename(savename), score=hotspots[[pattern]]$num.events)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

stopCluster(cl)

}
