#' Read aneufinder configuration file
#'
#' Read an aneufinder configuration file into a list structure. The configuration file has to be specified in INI format. R expressions can be used and will be evaluated.
#'
#' @param configfile Path to the configuration file
#' @author Aaron Taudt
#' @export
readConfig <- function(configfile) {

	connection <- file(configfile) 
  Lines  <- readLines(connection) 
  close(connection) 

  Lines <- chartr("[]", "==", Lines) # change section headers 
	Lines <- gsub(" ", "", Lines) # no spaces

  connection <- textConnection(Lines) 
  data <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE, quote="") 
  close(connection) 
	names(data) <- c('argument','value','section')

  L <- data$argument == "" # location of section breaks 
  data <- subset(transform(data, section = value[which(L)[cumsum(L)]])[1:3], argument != "") 

  configlist <- list() 
  ToParse  <- paste0("configlist$", data$section, "$",  data$argument, " <- ", data$value) 

  eval(parse(text=ToParse)) 

  return(configlist) 
} 
