# aneufinder - An R-package for CNV detection in whole-genome single cell sequencing data
# Copyright (C) 2015  Aaron Taudt
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


#' Wrapper function for the aneufinder package
#'
#' This function is an easy-to-use wrapper to \link[aneufinder:binning]{bin the data}, \link[aneufinder:findCNVs]{find copy-number-variations}, \link[aneufinder:findSCEs]{find sister-chromatid-exchange} events and plot the results.
#'
#' @param inputfolder Folder with either BAM or BED files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param config A file with parameters that are used by the various subfunctions. Alternatively a list with parameters created by \code{\link{makeConfig}}.
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @export
Aneufinder <- function(inputfolder, outputfolder, config=makeConfig()) {

#=======================
### Helper functions ###
#=======================
as.object <- function(x) {
	return(eval(parse(text=x)))
}

#========================
### General variables ###
#========================
if (is.character(config)) {
	## Read config file ##
	errstring <- tryCatch({
		conf <- readConfig(config)
		errstring <- ''
	}, error = function(err) {
		errstring <- paste0("Could not read configuration file ",config)
	})
	if (errstring!='') {
		stop(errstring)
	}
	## Make a copy of the conf file ##
	temp <- file.copy(from=config, to=file.path(outputfolder, basename(config)), overwrite=TRUE)
} else {
	conf <- config
}

## Helpers
binsizes <- conf$Binning$binsizes
reads.per.bins <- conf$Binning$reads.per.bin
patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=F, trim=T),'_'))
patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))

## Set up the directory structure ##
binpath.uncorrected <- file.path(outputfolder,'binned')
CNVpath <- file.path(outputfolder,'CNV_results')
CNVplotpath <- file.path(outputfolder,'CNV_plots')
CNVbrowserpath <- file.path(outputfolder,'CNV_browser_files')
SCEpath <- file.path(outputfolder,'SCE_results')
SCEplotpath <- file.path(outputfolder,'SCE_plots')
SCEbrowserpath <- file.path(outputfolder,'SCE_browser_files')
## Delete old directory if desired ##
if (conf$General$reuse.existing.files==FALSE) {
	if (file.exists(outputfolder)) {
		message("Deleting old directory ",outputfolder)
		unlink(outputfolder, recursive=T)
	}
}
if (!file.exists(outputfolder)) {
	dir.create(outputfolder)
}

## Parallelization ##
message("Using ",conf$General$numCPU," CPUs")
cl <- makeCluster(conf$General$numCPU)
doParallel::registerDoParallel(cl)

#==============
### Binning ###
#==============
message("Binning the data ...", appendLF=F); ptm <- proc.time()
if (!file.exists(binpath.uncorrected)) { dir.create(binpath.uncorrected) }
## Prepare parameter list ##
paramlist <- conf$Binning
paramlist$outputfolder <- binpath.uncorrected
paramlist$save.as.RData <- TRUE
paramlist$format <- NULL
paramlist$GC.correction.bsgenome.string <- NULL

files <- list.files(inputfolder, full=T, recursive=T, pattern=paste0('.',conf$Binning$format,'$'))
temp <- foreach (file = files, .packages=c('aneufinder',conf$Binning$GC.correction.bsgenome.string)) %dopar% {
	message(file)
	## Check for existing files ##
	existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value=T)
	existing.binsizes <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 2)))
	existing.rpbin <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 3)))
	binsizes.todo <- setdiff(binsizes, existing.binsizes)
	rpbin.todo <- setdiff(reads.per.bins, existing.rpbin)
	if (length(c(binsizes.todo,rpbin.todo)) > 0) {
		# Make paramter list
		paramlist$binsizes <- binsizes.todo
		paramlist$reads.per.bin <- rpbin.todo
		if (conf$Binning$format=='bam') {
			paramlist$bamfile <- file
			do.call('bam2binned', paramlist)
		} else {
			warning("Unknown file format ",conf$Binning$format)
		}
	}
}
time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

#=================
### Correction ###
#=================
if (!is.null(conf$Correction$method)) {

	message(paste0(conf$Correction$method," correction ..."), appendLF=F); ptm <- proc.time()
	if (conf$Correction$method=='GC') {
		binpath.corrected <- file.path(outputfolder,'binned_GC-corrected')
		if (!file.exists(binpath.corrected)) { dir.create(binpath.corrected) }
		## Load BSgenome
		suppressPackageStartupMessages(library(conf$Correction$GC.bsgenome, character.only=T))
		conf$Correction$GC.bsgenome.string <- conf$Correction$GC.bsgenome
		conf$Correction$GC.bsgenome <- as.object(conf$Correction$GC.bsgenome.string) # replacing string by object

		## Go through patterns
		temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
			binfiles <- list.files(binpath.uncorrected, pattern='RData$', full=T)
			binfiles <- grep(pattern, binfiles, value=T)
			binfiles.corrected <- list.files(binpath.corrected, pattern='RData$', full=T)
			binfiles.corrected <- grep(pattern, binfiles.corrected, value=T)
			binfiles.todo <- setdiff(basename(binfiles), basename(binfiles.corrected))
			if (length(binfiles.todo)>0) {
				binfiles.todo <- paste0(binpath.uncorrected,.Platform$file.sep,binfiles.todo)
				if (grepl('binsize',pattern)) {
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf$Correction$GC.bsgenome, same.GC.content=TRUE))
				} else {
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,conf$Correction$GC.bsgenome))
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
if (conf$CNV$findCNVs==TRUE) {

	message("Finding CNVs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(CNVpath)) { dir.create(CNVpath) }
	## Prepare parameter list ##
	paramlist <- conf$CNV
	paramlist$findCNVs <- NULL

	files <- list.files(binpath, full=T, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVpath,basename(file))
		if (!file.exists(savename)) {
			# Make parameter list
			paramlist$binned.data <- file
			model <- do.call('findCNVs', paramlist)
			save(model, file=savename)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	#===================
	### Plotting CNV ###
	#===================
	if (!file.exists(CNVplotpath)) { dir.create(CNVplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=F, trim=T),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(CNVpath, full=T, recursive=T, pattern='.RData$')

	#-----------------------
	## Plot distributions ##
	#-----------------------
	message("Plotting distributions ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVplotpath,paste0('distributions_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=10, height=7)
			ifiles <- list.files(CNVpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='histogram') + theme_bw(base_size=18) + theme(legend.position=c(1,1), legend.justification=c(1,1)))
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
		ifiles <- list.files(CNVpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(CNVplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=F, cluster=conf$Plotting$cluster))
		}
	}
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(CNVpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(CNVplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(savename, width=30, height=0.3*length(ifiles))
			ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf$Plotting$cluster))
			print(ggplt)
			d <- dev.off()
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
			ifiles <- list.files(CNVpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='arrayCGH'))
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
			ifiles <- list.files(CNVpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='karyogram'))
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
		savename <- file.path(CNVbrowserpath,paste0(pattern))
		if (!file.exists(paste0(savename,'.bed.gz'))) {
			ifiles <- list.files(CNVpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			exportCNVs(ifiles, filename=savename, cluster=conf$Plotting$cluster, export.CNV=TRUE, export.SCE=FALSE)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

#===============
### findSCEs ###
#===============
if (conf$SCE$findSCEs==TRUE) {

	message("Finding SCEs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(SCEpath)) { dir.create(SCEpath) }
	## Prepare parameter list ##
	paramlist <- conf$SCE
	paramlist$findSCEs <- NULL

	files <- list.files(binpath, full=T, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEpath,basename(file))
		if (!file.exists(savename)) {
			# Make parameter list
			paramlist$binned.data <- file
			model <- do.call('findSCEs', paramlist)
			save(model, file=savename)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

	#===================
	### Plotting SCE ###
	#===================
	if (!file.exists(SCEplotpath)) { dir.create(SCEplotpath) }
	patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=F, trim=T),'_'))
	patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
	files <- list.files(SCEpath, full=T, recursive=T, pattern='.RData$')

	#-----------------------
	## Plot distributions ##
	#-----------------------
	message("Plotting distributions ...", appendLF=F); ptm <- proc.time()
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEplotpath,paste0('distributions_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(file=savename, width=10, height=7)
			ifiles <- list.files(SCEpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='histogram') + theme_bw(base_size=18) + theme(legend.position=c(1,1), legend.justification=c(1,1)))
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
		ifiles <- list.files(SCEpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(SCEplotpath,paste0('genomeHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=T, cluster=conf$Plotting$cluster))
		}
	}
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(SCEpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(SCEplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(savename, width=30, height=0.3*length(ifiles))
			ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=conf$Plotting$cluster))
			print(ggplt)
			d <- dev.off()
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
			ifiles <- list.files(SCEpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='arrayCGH'))
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
			ifiles <- list.files(SCEpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			for (ifile in ifiles) {
				model <- get(load(ifile))
				print(plot(model, type='karyogram'))
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
			ifiles <- list.files(SCEpath, pattern='RData$', full=T)
			ifiles <- grep(pattern, ifiles, value=T)
			exportCNVs(ifiles, filename=savename, cluster=conf$Plotting$cluster, export.CNV=TRUE, export.SCE=TRUE)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

stopCluster(cl)

}
