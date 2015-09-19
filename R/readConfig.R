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
  data$section <- data$value[which(L)[cumsum(L)]]
  data <- data[data$argument!="",]

  configlist <- list() 
  ToParse  <- paste0("configlist$", data$section, "$",  data$argument, " <- ", data$value) 

  eval(parse(text=ToParse)) 

  return(configlist) 
} 
