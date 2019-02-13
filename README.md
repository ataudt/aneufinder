[![bioc](http://www.bioconductor.org/shields/years-in-bioc/AneuFinder.svg)](http://bioconductor.org/packages/devel/bioc/html/AneuFinder.html)
[![bioc](http://www.bioconductor.org/shields/downloads/AneuFinder.svg)](https://bioconductor.org/packages/stats/bioc/AneuFinder/) 


AneuFinder
==========

Copy-number detection in single-cell whole genome sequencing (scWGS) and Strand-seq data using a Hidden Markov Model or binary bisection method. The package implements copy-number detection, estimation of breakpoints with confidence intervals, hotspot detection, commonly used plotting functions, export to BED format for upload to genome browsers, and measures for assessment of karyotype heterogeneity and quality metrics.

Installation
------------

### Stable release version from Bioconductor
To install the *current stable* version from Bioconductor, please visit http://bioconductor.org/packages/AneuFinder/ and follow the provided instructions.

### Development version from Github
To install the *development* version from Github, follow the steps given below. The installation has only been tested on Ubuntu so far, if you need to install on Windows or Mac additional steps might be necessary (e.g. installation of Rtools from https://cran.r-project.org/bin/windows/Rtools/)

1. Install a recent version of R (>3.3) from https://www.r-project.org/
2. Optional: For ease of use, install Rstudio from https://www.rstudio.com/
3. Open R and install all dependencies. Please ensure that you have writing permissions to install packages. Execute the following lines one by one:

   install.packages("devtools")  
	 if (!requireNamespace("BiocManager", quietly=TRUE))
    	 install.packages("BiocManager")
	 BiocManager::install("GenomicRanges")  
	 BiocManager::install("GenomicAlignments")  
	 library(devtools)  
	 install_github("ataudt/aneufinderData")  
	 install_github("ataudt/aneufinder")  
	 #### Or alternatively if the above line doesn't work:  
	 install_git("git://github.com/ataudt/aneufinderData.git", branch = "master")  
	 install_git("git://github.com/ataudt/aneufinder.git", branch = "master")

How to use AneuFinder
--------------------

Please refer to the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/AneuFinder/inst/doc/AneuFinder.pdf) on Bioconductor for tutorials on how to use AneuFinder.

Report Errors
-------------

If you encounter errors of any kind, please file an [issue here](https://github.com/ataudt/aneufinder/issues/new). I will try to react within a week.
