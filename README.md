!!!THIS PROJECT IS CURRENTLY BEING DEVELOPED!!!

The easiest way to install it is currently the following:
1) Install a recent version of R from https://www.r-project.org/
2) Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3) Open R and install all dependencies
	source("http://bioconductor.org/biocLite.R")
	biocLite("GenomicRanges")
	biocLite("GenomicAlignments")
	biocLite("Rsamtools")
	install.packages("devtools")
	install_github("chakalakka/aneufinder")

