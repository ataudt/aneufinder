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


#' Aneufinder configuration
#'
#' Make a list of configuration parameters for the \code{\link{Aneufinder}} function.
#'
#' @param General General options.
#' @param Binning Parameters for binning. See \code{\link{binning}} for which options are available.
#' @param Correction Parameters for correction methods.
#' @param CNV Parameters for CNV detection. See \code{\link{findCNVs}} for which options are available.
#' @param SCE Parameters for SCE detection. See \code{\link{findSCEs}} for which options are available.
#' @param Plotting Plotting parameters.
#' @author Aaron Taudt
#' @export
makeConfig <- function(General=list(reuse.existing.files=TRUE, numCPU=1), Binning=list(format='bam', binsizes=500000), Correction=list(), CNV=list(findCNVs=TRUE), SCE=list(findSCEs=FALSE), Plotting=list()) {
	return(list(General=General, Binning=Binning, Correction=Correction, CNV=CNV, SCE=SCE, Plotting=Plotting))
}
