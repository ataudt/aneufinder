#' Wrapper function for the aneufinder package
#'
#' This function is an easy-to-use wrapper to \link[aneufinder:binning]{bin the data}, \link[aneufinder:findCNVs]{find copy-number-variations}, \link[aneufinder:findSCEs]{find sister-chromatid-exchange} events and plot the results.
#'
#' @param inputfolder Folder with BAM files.
#' @param configfile A file with parameters that are used by the various functions.
#' @param outputfolder Folder to put the results. If it does not exist it will be created.
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @export
easyAneufinder <- function(inputfolder, configfile, outputfolder) {

#=======================
### Helper functions ###
#=======================
as.object <- function(x) {
	return(eval(parse(text=x)))
}

#========================
### General variables ###
#========================
## Read config file ##
tC <- tryCatch({
	config <- readConfig(configfile)
}, error = function(err) {
	message("Error: Could not read configuration file ",configfile)
})

## Helpers
binsizes <- config$Binning$binsizes
reads.per.bins <- config$Binning$reads.per.bin
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
if (config$General$reuse.existing.files==FALSE) {
	if (file.exists(outputfolder)) {
		message("Deleting old directory ",outputfolder)
		unlink(outputfolder, recursive=T)
	}
}
if (!file.exists(outputfolder)) {
	dir.create(outputfolder)
}
## Make a copy of the config file ##
temp <- file.copy(configfile, file.path(outputfolder, basename(configfile)))

## Parallelization ##
message("Using ",config$General$numCPU," CPUs")
cl <- makeCluster(config$General$numCPU)
doParallel::registerDoParallel(cl)

#==============
### Binning ###
#==============
message("Binning the data ...", appendLF=F); ptm <- proc.time()
if (!file.exists(binpath.uncorrected)) { dir.create(binpath.uncorrected) }
## Prepare parameter list ##
paramlist <- config$Binning
paramlist$outputfolder <- binpath.uncorrected
paramlist$save.as.RData <- TRUE
paramlist$format <- NULL
paramlist$GC.correction.bsgenome.string <- NULL

files <- list.files(inputfolder, full=T, recursive=T, pattern='.bam$')
temp <- foreach (file = files, .packages=c('aneufinder',config$Binning$GC.correction.bsgenome.string)) %dopar% {
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
		if (config$Binning$format=='bam') {
			paramlist$bamfile <- file
			do.call('bam2binned', paramlist)
		} else {
			warning("Unknown file format ",config$Binning$format)
		}
	}
}
time <- proc.time() - ptm; message(" ",round(time[3],2),"s")

#=================
### Correction ###
#=================
if (!is.null(config$Correction$method)) {

	message(paste0(config$Correction$method," correction ..."), appendLF=F); ptm <- proc.time()
	if (config$Correction$method=='GC') {
		binpath.corrected <- file.path(outputfolder,'binned_GC-corrected')
		if (!file.exists(binpath.corrected)) { dir.create(binpath.corrected) }
		## Load BSgenome
		suppressPackageStartupMessages(library(config$Correction$GC.bsgenome, character.only=T))
		config$Correction$GC.bsgenome.string <- config$Correction$GC.bsgenome
		config$Correction$GC.bsgenome <- as.object(config$Correction$GC.bsgenome.string) # replacing string by object

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
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,config$Correction$GC.bsgenome, same.GC.content=TRUE))
				} else {
					binned.data.list <- suppressMessages(correctGC(binfiles.todo,config$Correction$GC.bsgenome))
				}
				for (i1 in 1:length(binned.data.list)) {
					binned.data <- binned.data.list[[i1]]
					savename <- file.path(binpath.corrected, basename(names(binned.data.list)[i1]))
					save(binned.data, file=savename)
				}
			}
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
	binpath <- binpath.corrected

} else {
	binpath <- binpath.uncorrected
}

#===============
### findCNVs ###
#===============
if (config$CNV$findCNVs==TRUE) {

	message("Finding CNVs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(CNVpath)) { dir.create(CNVpath) }
	## Prepare parameter list ##
	paramlist <- config$CNV
	paramlist$findCNVs <- NULL

	files <- list.files(binpath, full=T, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(CNVpath,basename(file))
		if (!file.exists(savename)) {
			# Make parameter list
			paramlist$binned.data <- file
			paramlist$ID <- basename(file)
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
			suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=F, cluster=config$Plotting$cluster))
		}
	}
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(CNVpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(CNVplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(savename, width=30, height=0.3*length(ifiles))
			ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=config$Plotting$cluster))
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
			exportCNVs(ifiles, filename=savename, cluster=config$Plotting$cluster, export.CNV=TRUE, export.SCE=FALSE)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

#===============
### findSCEs ###
#===============
if (config$SCE$findSCEs==TRUE) {

	message("Finding SCEs ...", appendLF=F); ptm <- proc.time()
	if (!file.exists(SCEpath)) { dir.create(SCEpath) }
	## Prepare parameter list ##
	paramlist <- config$SCE
	paramlist$findSCEs <- NULL

	files <- list.files(binpath, full=T, recursive=T, pattern='.RData$')
	temp <- foreach (file = files, .packages=c('aneufinder')) %dopar% {
		savename <- file.path(SCEpath,basename(file))
		if (!file.exists(savename)) {
			# Make parameter list
			paramlist$binned.data <- file
			paramlist$ID <- basename(file)
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
			suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.SCE=T, cluster=config$Plotting$cluster))
		}
	}
	temp <- foreach (pattern = patterns, .packages=c('aneufinder')) %dopar% {
		ifiles <- list.files(SCEpath, pattern='RData$', full=T)
		ifiles <- grep(pattern, ifiles, value=T)
		savename=file.path(SCEplotpath,paste0('aneuploidyHeatmap_',sub('_$','',pattern),'.pdf'))
		if (!file.exists(savename)) {
			pdf(savename, width=30, height=0.3*length(ifiles))
			ggplt <- suppressMessages(heatmapAneuploidies(ifiles, cluster=config$Plotting$cluster))
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
			exportCNVs(ifiles, filename=savename, cluster=config$Plotting$cluster, export.CNV=TRUE, export.SCE=TRUE)
		}
	}
	time <- proc.time() - ptm; message(" ",round(time[3],2),"s")
}

stopCluster(cl)

}
