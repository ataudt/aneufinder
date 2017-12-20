

#' Wrapper function for the \code{\link{AneuFinder}} package
#'
#' This function is an easy-to-use wrapper to \link[AneuFinder:binning]{bin the data}, \link[AneuFinder:findCNVs]{find copy-number-variations}, \link[AneuFinder:getBreakpoints]{locate breakpoints}, plot \link[AneuFinder:heatmapGenomewide]{genomewide heatmaps}, \link[AneuFinder:plot.aneuHMM]{distributions, profiles and karyograms}.
#'
#' @param inputfolder Folder with either BAM or BED files.
#' @param outputfolder Folder to output the results. If it does not exist it will be created.
#' @param configfile A file specifying the parameters of this function (without \code{inputfolder}, \code{outputfolder} and \code{configfile}). Having the parameters in a file can be handy if many samples with the same parameter settings are to be run. If a \code{configfile} is specified, it will take priority over the command line parameters.
#' @param numCPU The numbers of CPUs that are used. Should not be more than available on your machine.
#' @param reuse.existing.files A logical indicating whether or not existing files in \code{outputfolder} should be reused.
#' @inheritParams bam2GRanges
#' @inheritParams bed2GRanges
#' @inheritParams binReads
#' @param reads.store Set \code{reads.store=TRUE} to store read fragments as RData in folder 'data' and as BED files in 'BROWSERFILES/data'. This option will force \code{use.bamsignals=FALSE}.
#' @param correction.method Correction methods to be used for the binned read counts. Currently only \code{'GC'}.
#' @param GC.BSgenome A \code{BSgenome} object which contains the DNA sequence that is used for the GC correction.
# #' @param mappability.reference A file that serves as reference for mappability correction.
#' @param strandseq A logical indicating whether the data comes from Strand-seq experiments. If \code{TRUE}, both strands carry information and are treated separately.
#' @inheritParams edivisive.findCNVs
#' @inheritParams HMM.findCNVs
#' @inheritParams findCNVs
#' @param confint Desired confidence interval for breakpoints. Set \code{confint=NULL} to disable confidence interval estimation. Confidence interval estimation will force \code{reads.store=TRUE}.
#' @param refine.breakpoints A logical indicating whether breakpoints from the HMM should be refined with read-level information. \code{refine.breakpoints=TRUE} will force \code{reads.store=TRUE}.
#' @param hotspot.bandwidth A vector the same length as \code{binsizes} with bandwidths for breakpoint hotspot detection (see \code{\link{hotspotter}} for further details). If \code{NULL}, the bandwidth will be chosen automatically as the average distance between reads.
#' @param hotspot.pval P-value for breakpoint hotspot detection (see \code{\link{hotspotter}} for further details). Set \code{hotspot.pval = NULL} to skip hotspot detection.
#' @param cluster.plots A logical indicating whether plots should be clustered by similarity.
#' @return \code{NULL}
#' @author Aaron Taudt
#' @import foreach
#' @import doParallel
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot
#' @importFrom utils read.table write.table
#' @importFrom cowplot plot_grid
#' @export
#'
#'@examples
#'\dontrun{
#'## The following call produces plots and genome browser files for all BAM files in "my-data-folder"
#'Aneufinder(inputfolder="my-data-folder", outputfolder="my-output-folder")}
#'
Aneufinder <- function(inputfolder, outputfolder, configfile=NULL, numCPU=1, reuse.existing.files=TRUE, binsizes=1e6, stepsizes=binsizes, variable.width.reference=NULL, reads.per.bin=NULL, pairedEndReads=FALSE, assembly=NULL, chromosomes=NULL, remove.duplicate.reads=TRUE, min.mapq=10, blacklist=NULL, use.bamsignals=FALSE, reads.store=FALSE, correction.method=NULL, GC.BSgenome=NULL, method=c('edivisive'), strandseq=FALSE, R=10, sig.lvl=0.1, eps=0.01, max.time=60, max.iter=5000, num.trials=15, states=c('zero-inflation',paste0(0:10,'-somy')), confint=NULL, refine.breakpoints=FALSE, hotspot.bandwidth=NULL, hotspot.pval=5e-2, cluster.plots=TRUE) {

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

## reads.store
if (!is.null(conf[['confint']])) {
  confint <- conf[['confint']]
}
if (!is.null(conf[['refine.breakpoints']])) {
  refine.breakpoints <- conf[['refine.breakpoints']]
}
if (refine.breakpoints & is.null(confint)) {
  confint <- 0.99
  conf[['confint']] <- confint
  warning("Changed 'confint=NULL' to 'confint=0.99' because 'refine.breakpoints=TRUE'.")
}
if (!is.null(confint) | refine.breakpoints) {
  if (reads.store == FALSE) {
    warning("Changed 'reads.store=TRUE' because we need it for confidence intervals and breakpoint refinement.")
  }
  reads.store <- TRUE
  conf[['reads.store']] <- TRUE
}
if (!is.null(conf[['reads.store']])) {
  reads.store <- conf[['reads.store']]
}
if (reads.store) {
  use.bamsignals <- FALSE
  conf[['use.bamsignals']] <- FALSE
}

if (length(hotspot.bandwidth) != length(binsizes) & !is.null(hotspot.bandwidth)) {
  hotspot.bandwidth <- rep(hotspot.bandwidth, length(binsizes))[1:length(binsizes)]
}

## Convert numCPU to numeric
numCPU <- as.numeric(numCPU)

## Put options into list and merge with conf
params <- list(numCPU=numCPU, reuse.existing.files=reuse.existing.files, binsizes=binsizes, stepsizes=stepsizes, variable.width.reference=variable.width.reference, reads.per.bin=reads.per.bin, pairedEndReads=pairedEndReads, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, blacklist=blacklist, reads.store=reads.store, use.bamsignals=use.bamsignals, correction.method=correction.method, GC.BSgenome=GC.BSgenome, method=method, strandseq=strandseq, eps=eps, max.time=max.time, max.iter=max.iter, num.trials=num.trials, states=states, R=R, sig.lvl=sig.lvl, confint=confint, refine.breakpoints=refine.breakpoints, hotspot.bandwidth=hotspot.bandwidth, hotspot.pval=hotspot.pval, cluster.plots=cluster.plots)
conf <- c(conf, params[setdiff(names(params),names(conf))])

## Check user input
if ('GC' %in% conf[['correction.method']] & is.null(conf[['GC.BSgenome']])) {
    stop("Option 'GC.bsgenome' has to be given if correction.method='GC'.")
}

## Determine format
files <- list.files(inputfolder, full.names=TRUE)
files.clean <- sub('\\.gz$','', files)
formats <- sapply(strsplit(files.clean, '\\.'), function(x) { rev(x)[1] })
datafiles <- files[formats %in% c('bam','bed')]
files.clean <- sub('\\.gz$','', datafiles)
formats <- sapply(strsplit(files.clean, '\\.'), function(x) { rev(x)[1] })
if (any(formats == 'bed') & is.null(conf[['assembly']])) {
  stop("Please specify 'assembly' if you have BED files in your inputfolder.")
}

## Helpers
binsizes <- conf[['binsizes']]
stepsizes <- conf[['stepsizes']]
reads.per.bins <- conf[['reads.per.bin']]
patterns <- c(paste0('reads.per.bin_',reads.per.bins,'_'), paste0('binsize_',format(binsizes, scientific=TRUE, trim=TRUE),'_stepsize_',format(stepsizes, scientific=TRUE, trim=TRUE),'_'))
patterns <- setdiff(patterns, c('reads.per.bin__','binsize__'))
pattern <- NULL #ease R CMD check
numcpu <- conf[['numCPU']]
if (!is.null(conf[['hotspot.bandwidth']])) {
  names(conf[['hotspot.bandwidth']]) <- patterns
}

## Set up the directory structure ##
readspath <- file.path(outputfolder,'data')
binpath.uncorrected <- file.path(outputfolder,'binned')
modelpath <- file.path(outputfolder, 'MODELS')
if (strandseq) {
    modelpath <- paste0(modelpath, '-StrandSeq')
}
refinedmodelpath <- paste0(modelpath, '_refined')
plotpath <- file.path(outputfolder, 'PLOTS')
browserpath <- file.path(outputfolder, 'BROWSERFILES')
readsbrowserpath <- file.path(browserpath,'data')
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


#===================
### Write README ###
#===================
savename <- file.path(outputfolder, 'README.txt')
cat("", file=savename)
cat("This folder contains the following files:\n", file=savename, append=TRUE)
cat("-----------------------------------------\n", file=savename, append=TRUE)
cat("- chrominfo.tsv: A tab-separated file with chromosome lengths.\n", file=savename, append=TRUE)
cat("- AneuFinder.config: A text file with all the parameters that were used to run Aneufinder().\n", append=TRUE, file=savename)

cat("\n", file=savename, append=TRUE)
cat("This folder contains the following folders. Some folders are [optional] depending on the parameters:\n", file=savename, append=TRUE)
cat("-------------------------------------------\n", file=savename, append=TRUE)
cat("- binned: RData files with the results of the binnig step. Contains GRanges objects with binned genomic coordinates and read counts.\n", file=savename, append=TRUE)
cat("- [binned-GC]: RData files with the results of the GC-correction step. Depends on option 'correction.method=\"GC\"'. Contains GRanges objects with binned genomic coordinates and read counts.\n", file=savename, append=TRUE)
cat("- BROWSERFILES: Bed files for upload to the UCSC genome browser.\n", file=savename, append=TRUE)
cat("- [data]: RData files with GRanges containing read-fragments. Depends on option 'reads.store=TRUE'.\n", file=savename, append=TRUE)
cat("- MODELS [or MODELS-StrandSeq]: RData files with aneuHMM objects. Result of the copy-number and breakpoint estimation step. Depends on option 'strandseq=TRUE'.\n", file=savename, append=TRUE)
cat("- MODELS_refined [or MODELS-StrandSeq_refined]: RData files with aneuHMM objects after the breakpoint refinement step. Depends on option 'strandseq=TRUE'.\n", file=savename, append=TRUE)
cat("- PLOTS: Several plots that are produced by default.\n", file=savename, append=TRUE)


#======================
### Parallelization ###
#======================
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
### Get chromosome lengths ###
## Get first bam file
bamfile <- grep('bam$', datafiles, value=TRUE)[1]
if (!is.na(bamfile)) {
    ptm <- startTimedMessage("Obtaining chromosome length information from file ", bamfile, " ...")
    chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
    stopTimedMessage(ptm)
} else {
    ## Read chromosome length information
    if (is.character(conf[['assembly']])) {
        if (file.exists(conf[['assembly']])) {
            ptm <- startTimedMessage("Obtaining chromosome length information from file ", conf[['assembly']], " ...")
            df <- utils::read.table(conf[['assembly']], sep='\t', header=TRUE)
            stopTimedMessage(ptm)
        } else {
            ptm <- startTimedMessage("Obtaining chromosome length information from UCSC ...")
            df.chroms <- GenomeInfoDb::fetchExtendedChromInfoFromUCSC(conf[['assembly']])
            ## Get first bed file
            bedfile <- grep('bed$|bed.gz$', datafiles, value=TRUE)[1]
            if (!is.na(bedfile)) {
                firstline <- read.table(bedfile, nrows=1)
                if (grepl('^chr',firstline[1,1])) {
                    df <- df.chroms[,c('UCSC_seqlevel','UCSC_seqlength')]
                } else {
                    df <- df.chroms[,c('NCBI_seqlevel','UCSC_seqlength')]
                }
            }
            stopTimedMessage(ptm)
        }
    } else if (is.data.frame(conf[['assembly']])) {
        df <- conf[['assembly']]
    } else {
        stop("'assembly' must be either a data.frame with columns 'chromosome' and 'length' or a character specifying the assembly.")
    }
    chrom.lengths <- df[,2]
    names(chrom.lengths) <- df[,1]
    chrom.lengths <- chrom.lengths[!is.na(chrom.lengths) & !is.na(names(chrom.lengths))]
}
chrom.lengths.df <- data.frame(chromosome=names(chrom.lengths), length=chrom.lengths)
## Write chromosome length information to file
utils::write.table(chrom.lengths.df, file=file.path(outputfolder, 'chrominfo.tsv'), sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    
### Make bins ###
message("==> Making bins:")
if (!is.null(conf[['variable.width.reference']])) {
  ## Determine format
  file <- conf[['variable.width.reference']]
  file.clean <- sub('\\.gz$','', file)
  format <- rev(strsplit(file.clean, '\\.')[[1]])[1]
  if (format == 'bam') {
    reads <- bam2GRanges(conf[['variable.width.reference']], chromosomes=conf[['chromosomes']], pairedEndReads=conf[['pairedEndReads']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']])
  } else if (format == 'bed') {
    reads <- bed2GRanges(conf[['variable.width.reference']], assembly=chrom.lengths.df, chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']])
  }
  bins <- variableWidthBins(reads, binsizes=conf[['binsizes']], stepsizes=conf[['stepsizes']], chromosomes=conf[['chromosomes']])
} else {
  bins <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=conf[['chromosomes']], binsizes=conf[['binsizes']], stepsizes=conf[['stepsizes']])
}
message("==| Finished making bins.")

### Binning ###
parallel.helper <- function(file) {
  existing.binfiles <- grep(basename(file), list.files(binpath.uncorrected), value=TRUE)
  existing.binsizes.stepsizes <- sub('_reads.per.bin.*', '', sub('.*_binsize', 'binsize', existing.binfiles))
  existing.rpbin <- as.numeric(unlist(lapply(strsplit(existing.binfiles, split='binsize_|_reads.per.bin_|_\\.RData'), '[[', 3)))
  binsizes.stepsizes.todo <- setdiff(sub('_$','',patterns), existing.binsizes.stepsizes)
  rpbin.todo <- setdiff(reads.per.bins, existing.rpbin)
  if (length(c(binsizes.stepsizes.todo,rpbin.todo)) > 0) {
    tC <- tryCatch({
      binReads(file=file, assembly=chrom.lengths.df, pairedEndReads=conf[['pairedEndReads']], binsizes=NULL, variable.width.reference=NULL, reads.per.bin=rpbin.todo, bins=bins[as.character(binsizes.stepsizes.todo)], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']], outputfolder.binned=binpath.uncorrected, save.as.RData=TRUE, reads.store=conf[['reads.store']], outputfolder.reads=readspath, use.bamsignals=conf[['use.bamsignals']])
    }, error = function(err) {
      stop(file,'\n',err)
    })
  }
}

## Bin the files
if (!file.exists(binpath.uncorrected)) { dir.create(binpath.uncorrected) }
files <- list.files(inputfolder, full.names=TRUE, pattern='\\.bam$|\\.bed$|\\.bed\\.gz$')
if (length(files) == 1 & conf[['cluster.plots']] == TRUE) {
  conf[['cluster.plots']] <- FALSE
  warning("Need more than one file for cluster.plots=TRUE. Continuing with cluster.plots=FALSE.")
}
if (numcpu > 1) {
  ptm <- startTimedMessage("Binning the data ...")
  temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
    parallel.helper(file)
  }
  stopTimedMessage(ptm)
} else {
  # temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
  for (file in files) {
    parallel.helper(file)
  }
}
  
### Read fragments that are not produced yet ###
if (!conf[['use.bamsignals']] & conf[['reads.store']]) {
  parallel.helper <- function(file) {
    savename <- file.path(readspath,paste0(basename(file),'.RData'))
    if (!file.exists(savename)) {
      tC <- tryCatch({
        binReads(file=file, assembly=chrom.lengths.df, pairedEndReads=conf[['pairedEndReads']], chromosomes=conf[['chromosomes']], remove.duplicate.reads=conf[['remove.duplicate.reads']], min.mapq=conf[['min.mapq']], blacklist=conf[['blacklist']], calc.complexity=FALSE, reads.store=TRUE, outputfolder.reads=readspath, reads.only=TRUE)
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
    # temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
    for (file in files) {
      parallel.helper(file)
    }
  }
  
  ### Export read fragments as browser file ###
  if (!file.exists(readsbrowserpath)) { dir.create(readsbrowserpath, recursive=TRUE) }
  readfiles <- list.files(readspath,pattern='.RData$',full.names=TRUE)
  
  parallel.helper <- function(file) {
    savename <- file.path(readsbrowserpath,sub('.RData','',basename(file)))
    if (!file.exists(paste0(savename,'.bed.gz'))) {
      tC <- tryCatch({
        gr <- loadFromFiles(file, check.class='GRanges')[[1]]
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
    # temp <- foreach (file = readfiles, .packages=c("AneuFinder")) %do% {
    for (file in readfiles) {
      parallel.helper(file)
    }
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
        ptm <- startTimedMessage(paste0(correction.method," correction ..."))
        # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
        for (pattern in patterns) {
          parallel.helper(pattern)
        }
        stopTimedMessage(ptm)
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
for (method in conf[['method']]) {

    modeldir <- file.path(modelpath, paste0('method-', method))
    if (!file.exists(modeldir)) { dir.create(modeldir, recursive=TRUE) }
    if (conf[['refine.breakpoints']]) {
        refinedmodeldir <- file.path(refinedmodelpath, paste0('method-', method))
        if (!file.exists(refinedmodeldir)) { dir.create(refinedmodeldir, recursive=TRUE) }
    } else {
        refinedmodeldir <- modeldir
    }
    plotdir <- file.path(plotpath, paste0('method-', method))
    if (!file.exists(plotdir)) { dir.create(plotdir, recursive=TRUE) }
    browserdir <- file.path(browserpath, paste0('method-', method))
    if (!file.exists(browserdir)) { dir.create(browserdir, recursive=TRUE) }
  
    files <- list.files(binpath, full.names=TRUE, pattern='.RData$')
    files <- grep(paste(gsub('\\+','\\\\+',patterns), collapse = '|'), files, value=TRUE)
    
    parallel.helper <- function(file) {
        tC <- tryCatch({
            savename <- file.path(modeldir,basename(file))
            if (!file.exists(savename)) {
                if (conf[['strandseq']]) {
                    findCNV <- findCNVs.strandseq
                } else {
                    findCNV <- findCNVs
                }
                if (method == 'dnacopy') {
                    model <- findCNV(file, method='dnacopy') 
                } else if (method == 'HMM') {
                    model <- findCNV(file, method='HMM', eps=conf[['eps']], max.time=conf[['max.time']], max.iter=conf[['max.iter']], num.trials=conf[['num.trials']], states=conf[['states']]) 
                } else if (method == 'edivisive') {
                    model <- findCNV(file, method='edivisive', R=conf[['R']], sig.lvl=conf[['sig.lvl']]) 
                }
                # Breakpoints and confidence intervals
                if (is.null(conf[['confint']])) {
                    reads.file <- NULL
                } else {
                    reads.file <- file.path(readspath, paste0(model$ID,'.RData'))
                }
                model$breakpoints <- getBreakpoints(model, fragments=reads.file, confint = conf[['confint']])
                ptm <- startTimedMessage("Saving to file ",savename," ...")
                save(model, file=savename)
                stopTimedMessage(ptm)
            }
        }, error = function(err) {
          stop(file,'\n',err)
        })
    }
    if (numcpu > 1) {
        if (method == 'dnacopy') {
            ptm <- startTimedMessage("Running DNAcopy ...")
        } else if (method == 'HMM') {
            ptm <- startTimedMessage("Running HMMs ...")
        } else if (method == 'edivisive') {
            ptm <- startTimedMessage("Running edivisive ...")
        }
        temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
            parallel.helper(file)
        }
        stopTimedMessage(ptm)
      } else {
          # temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
          for (file in files) {
              parallel.helper(file)
          }
    }
  
    #========================
    ### refineBreakpoints ###
    #========================
    if (conf[['refine.breakpoints']]) {
        files <- list.files(modeldir, full.names=TRUE, pattern='.RData$')
        parallel.helper <- function(file) {
            tC <- tryCatch({
                savename <- file.path(refinedmodeldir,basename(file))
                if (!file.exists(savename)) {
                    model <- loadFromFiles(file)[[1]]
                    reads.file <- file.path(readspath, paste0(model$ID,'.RData'))
                    ## Refining breakpoints
                    if (!is.null(model$breakpoints)) {
                        message("Breakpoint refinement for ", basename(file))
                        model <- refineBreakpoints(model, fragments=reads.file, breakpoints=model$breakpoints, confint = conf[['confint']])
                    }
                    ptm <- startTimedMessage("Saving breakpoints to file ",savename," ...")
                    save(model, file=savename)
                    stopTimedMessage(ptm)
                }
            }, error = function(err) {
                stop(file,'\n',err)
            })
        }
        if (numcpu > 1) {
            ptm <- startTimedMessage("Refining breakpoints ...")
            temp <- foreach (file = files, .packages=c("AneuFinder")) %dopar% {
                parallel.helper(file)
            }
            stopTimedMessage(ptm)
        } else {
            # temp <- foreach (file = files, .packages=c("AneuFinder")) %do% {
            for (file in files) {
                parallel.helper(file)
            }
        }
    }
  
    #=======================
    ### Finding hotspots ###
    #=======================
    hslist <- list()
    if (is.null(hotspot.pval)) {
        message("Skipping breakpoint hotspot detection because 'hotspot.pval=NULL'.")
    } else {
        parallel.helper <- function(pattern) {
            ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
            ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
            breakpoints <- list()
            total.read.count <- numeric()
            for (file in ifiles) {
                hmm <- suppressMessages( loadFromFiles(file)[[1]] )
                breakpoints[[basename(file)]] <- hmm$breakpoints
                total.read.count[basename(file)] <- hmm$qualityInfo$total.read.count
            }
            if (is.null(conf[['hotspot.bandwidth']])) {
                bw <- sum(as.numeric(seqlengths(hmm$bins))) / mean(total.read.count)
            } else {
                bw <- conf[['hotspot.bandwidth']][pattern]
            }
            hslist <- suppressMessages( hotspotter(breakpoints, bw=bw, pval=conf[['hotspot.pval']], spacing.bp = bw) )
            return(hslist)
        }
        if (numcpu > 1) {
            ptm <- startTimedMessage("Finding breakpoint hotspots ...")
            hslist <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
                parallel.helper(pattern)
            }
            stopTimedMessage(ptm)
        } else {
            ptm <- startTimedMessage("Finding breakpoint hotspots ...")
            hslist <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
                parallel.helper(pattern)
            }
            stopTimedMessage(ptm)
        }
        names(hslist) <- patterns
    }
  
    #===============
    ### Plotting ###
    #===============
    if (!file.exists(plotdir)) { dir.create(plotdir) }
    files <- list.files(refinedmodeldir, full.names=TRUE, pattern='.RData$')
    if (strandseq) {
        strandseq.string <- '_StrandSeq'
    } else {
        strandseq.string <- ''
    }
    if (conf[['refine.breakpoints']]) {
        strandseq.string <- paste0(strandseq.string, '_refined')
    }
  
    #------------------
    ## Plot heatmaps ##
    #------------------
    parallel.helper <- function(pattern) {
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        if (length(ifiles)>0) {
            savename=file.path(plotdir,paste0('genomeHeatmap_',sub('_$','',pattern), strandseq.string, '.pdf'))
            if (!file.exists(savename)) {
                suppressMessages(heatmapGenomewide(ifiles, file=savename, plot.breakpoints=FALSE, hotspots=NULL, cluster=conf[['cluster.plots']]))
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
        # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
        for (pattern in patterns) {
            parallel.helper(pattern)
        }
    }
  
    parallel.helper <- function(pattern) {
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        if (length(ifiles)>0) {
            savename=file.path(plotdir,paste0('aneuploidyHeatmap_',sub('_$','',pattern), strandseq.string,'.pdf'))
            if (!file.exists(savename)) {
                grDevices::pdf(savename, width=30, height=max(0.3*length(ifiles), 2/2.54))
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
        ptm <- startTimedMessage("Plotting chromosome heatmaps ...")
        # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
        for (pattern in patterns) {
            parallel.helper(pattern)
        }
        stopTimedMessage(ptm)
    }
  
    #-------------------------
    ## Export browser files ##
    #-------------------------
    if (!file.exists(browserdir)) { dir.create(browserdir) }
    parallel.helper <- function(pattern) {
        ## Export CNV and breakpoints
        savename <- paste0(file.path(browserdir,sub('_$','',pattern)), strandseq.string)
        if (!file.exists(paste0(savename,'_CNV.bed.gz'))) {
            ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
            ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
            exportCNVs(ifiles, filename=savename, cluster=conf[['cluster.plots']], export.CNV=TRUE, export.breakpoints=TRUE)
        }
        ## Breakpoint hotspots
        savename <- file.path(browserdir,paste0(sub('_$','',pattern), strandseq.string, '_breakpoint-hotspots'))
        if (!file.exists(paste0(savename,'.bed.gz'))) {
            hotspots <- hslist[[pattern]]$hotspots
            if (!is.null(hotspots)) {
                exportGRanges(hotspots, filename=savename, trackname=basename(savename), score=hotspots$num.events, thickStart = hotspots$start.max, thickEnd = hotspots$end.max, priority=41)
            }
        }
        ## Hotspot densities
        savename <- file.path(browserdir,paste0(sub('_$','',pattern), strandseq.string, '_breakpoint-hotspot-densities'))
        if (!file.exists(paste0(savename,'.wig.gz'))) {
            densities <- hslist[[pattern]]$densities
            if (!is.null(hotspots)) {
                exportGRanges(densities, filename=savename, trackname=basename(savename), as.wiggle = TRUE, wiggle.val = densities$kde, priority=40)
            }
        }
    }
    if (numcpu > 1) {
      ptm <- startTimedMessage("Exporting browser files ...")
      temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %dopar% {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    } else {
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
    }
  
    #------------------
    ## Plot profiles ##
    #------------------
    parallel.helper <- function(pattern) {
      savename <- file.path(plotdir,paste0('profiles_',sub('_$','',pattern), strandseq.string,'.pdf'))
      if (!file.exists(savename)) {
        grDevices::pdf(file=savename, width=20, height=10)
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        for (ifile in ifiles) {
          tC <- tryCatch({
            model <- get(load(ifile))
            p1 <- graphics::plot(model, type='profile', plot.breakpoints=FALSE)
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
      ptm <- startTimedMessage("Making profile and distribution plots ...")
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    }
  
    #--------------------
    ## Plot karyograms ##
    #--------------------
    parallel.helper <- function(pattern) {
      savename <- file.path(plotdir,paste0('karyograms_',sub('_$','',pattern), strandseq.string,'.pdf'))
      if (!file.exists(savename)) {
        grDevices::pdf(file=savename, width=12*1.4, height=2*4.6)
        ifiles <- list.files(refinedmodeldir, pattern='RData$', full.names=TRUE)
        ifiles <- grep(gsub('\\+','\\\\+',pattern), ifiles, value=TRUE)
        for (ifile in ifiles) {
          tC <- tryCatch({
            model <- get(load(ifile))
            print(graphics::plot(model, type='karyogram', plot.breakpoints=TRUE))
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
      ptm <- startTimedMessage("Plotting karyograms ...")
      # temp <- foreach (pattern = patterns, .packages=c("AneuFinder")) %do% {
      for (pattern in patterns) {
        parallel.helper(pattern)
      }
      stopTimedMessage(ptm)
    }

}

total.time <- proc.time() - total.time
message("==> Total time spent: ", round(total.time[3]), "s <==")

}
