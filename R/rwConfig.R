

#' Read AneuFinder configuration file
#'
#' Read an AneuFinder configuration file into a list structure. The configuration file has to be specified in INI format. R expressions can be used and will be evaluated.
#'
#' @param configfile Path to the configuration file
#' @return A \code{list} with one entry for each element in \code{configfile}.
#' @author Aaron Taudt
#' @importFrom utils read.table
readConfig <- function(configfile) {

	connection <- file(configfile) 
  Lines  <- readLines(connection) 
  close(connection) 

  Lines <- chartr("[]", "==", Lines) # change section headers 
	Lines <- gsub(" ", "", Lines) # no spaces

  connection <- textConnection(Lines) 
  data <- utils::read.table(connection, as.is = TRUE, sep = "=", fill = TRUE, quote="") 
  close(connection) 
	names(data) <- c('argument','value','section')

  L <- data$argument == "" # location of section breaks 
  data$section <- data$value[which(L)[cumsum(L)]]
  data <- data[data$argument!="",]

  configlist <- list() 
	ToParse <- paste0("configlist$", data$argument, " <- ", data$value)
#   ToParse  <- paste0("configlist$", data$section, "$",  data$argument, " <- ", data$value) # with sections

  eval(parse(text=ToParse)) 

  return(configlist) 
} 

#' Write AneuFinder configuration file
#'
#' Write an AneuFinder configuration file from a list structure.
#'
#' @param conf A list structure with parameter values. Each entry will be written in one line.
#' @param configfile Filename of the outputfile.
#' @return \code{NULL}
#' @author Aaron Taudt
writeConfig <- function(conf, configfile) {

	## Printing function
	formatstring <- function(string) {
		if (is.character(string) & length(string)>1) {
			string <- paste0("c('",paste0(string,collapse="','"),"')")
		} else if (is.character(string) & length(string)==1) {
			string <- paste0("'",string,"'")
		} else if (is.numeric(string) & length(string)>1) {
			string <- paste0("c(",paste0(string,collapse=','),")")
		} else if (is.numeric(string) & length(string)==1) {
			string <- string
		} else if (is.null(string)) {
			string <- "NULL"
		} else if (is.data.frame(string)) {
        if (all(names(string) %in% c('chromosome','length'))) {
            string <- paste0("'", file.path(dirname(configfile), 'chrominfo.tsv'), "'")
        }
    }
		return(string)
	}
		
	f <- file(configfile, open='w')
	cat("#============== Aneufinder configuration file ===============#\n", file=f)
	cat("\n[General]\n", file=f)
	for (i1 in c('numCPU','reuse.existing.files')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[Binning]\n", file=f)
	for (i1 in c('binsizes', 'stepsizes', 'variable.width.reference', 'reads.per.bin', 'pairedEndReads', 'assembly', 'chromosomes', 'remove.duplicate.reads', 'min.mapq', 'blacklist', 'reads.store', 'use.bamsignals')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[Correction]\n", file=f)
	for (i1 in c('correction.method', 'GC.BSgenome')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[CopyNumberCalling]\n", file=f)
	for (i1 in c('method', 'strandseq')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[CopyNumberCalling_HMM]\n", file=f)
	for (i1 in c('eps', 'max.time', 'max.iter', 'num.trials', 'states', 'most.frequent.state', 'most.frequent.state.strandseq')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[CopyNumberCalling_edivisive]\n", file=f)
	for (i1 in c('R', 'sig.lvl')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[Breakpoint_Detection]\n", file=f)
	for (i1 in c('confint','refine.breakpoints','hotspot.bandwidth','hotspot.pval')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	cat("\n[Plotting]\n", file=f)
	for (i1 in c('cluster.plots')) {
		cat(i1," = ",formatstring(conf[[i1]]),"\n", file=f)
	}
	close(f, type='w')
}
