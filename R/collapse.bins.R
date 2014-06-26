collapse.bins = function(data, column2collapseBy=NULL, columns2sumUp=NULL) {

	# Indexing stuff
	ind_coords = 1:3
	ind_morecols = setdiff(1:ncol(data), c(ind_coords,columns2sumUp))
	ind_sumcols = columns2sumUp

	# Make the comparison vector
	if (is.null(column2collapseBy)) {
		c = data$start
		cShift1 = rep(NA,length(c))
		cShift1[2:length(cShift1)] = data$end[-length(c)] + 1
	} else {
		if (is(data[,column2collapseBy], "factor")) {
			c <- as.integer(data[,column2collapseBy])
		} else {
			c = data[,column2collapseBy]
		}
		cShift1 = rep(NA,length(c))
		cShift1[-1] = c[-length(c)]
	}
	compare_custom = c != cShift1
	# Make the comparison vector to separate chromosomes
	c <- as.integer(data[,1])
	cShift1 = rep(NA,length(c))
	cShift1[-1] = c[-length(c)]
	compare_chrom = c != cShift1
	# Combine the vectors
	compare <- compare_custom | compare_chrom
	compare[1] = TRUE
	numcollapsedbins = length(which(compare==TRUE))
	numbins = nrow(data)

	# Select the collapsed rows
	collapsed.bins = NULL
	collapsed.bins$chrom = data[compare,1]
	collapsed.bins$start = data[compare,2]
	collapsed.bins$end = data[c((which(compare==TRUE)-1)[-1],numbins), 3]
	if (length(ind_morecols)==1) {
		collapsed.bins[[4]] = data[compare, ind_morecols]
	} else if (length(ind_morecols)>1) {
		lcb = length(collapsed.bins)
		lmc = length(ind_morecols)
		collapsed.bins[(lcb+1):(lcb+lmc)] = data[compare, ind_morecols]
	}

	# Sum up columns
	if (!is.null(columns2sumUp)) {
		sumcols = as.matrix(data[,columns2sumUp])
		collapsed_sumcols = matrix(rep(NA,numcollapsedbins*length(columns2sumUp)), ncol=length(columns2sumUp))
		pb = txtProgressBar(min=1, max=length(compare), style=3)
		icount = 1
		i1_lasttrue = 1
		for (i1 in 1:length(compare)) {
			if (compare[i1]==TRUE) {
				if (length(columns2sumUp)==1) {
					collapsed_sumcols[icount-1] = sum(sumcols[i1_lasttrue:(i1-1),])
				} else {
					if (i1_lasttrue==i1-1 | i1==1) {
						collapsed_sumcols[icount-1,] = as.numeric(sumcols[i1_lasttrue,])
					} else {
						collapsed_sumcols[icount-1,] = apply(sumcols[i1_lasttrue:(i1-1),],2,sum)
					}
				}
					
				icount = icount+1
				i1_lasttrue = i1
				setTxtProgressBar(pb, i1)
			}
		}
		i1 = i1+1
		if (length(columns2sumUp)==1) {
			collapsed_sumcols[icount-1] = sum(sumcols[i1_lasttrue:(i1-1),])
		} else {
			if (i1_lasttrue==i1-1 | i1==1) {
				collapsed_sumcols[icount-1,] = as.numeric(sumcols[i1_lasttrue,])
			} else {
				collapsed_sumcols[icount-1,] = apply(sumcols[i1_lasttrue:(i1-1),],2,sum)
			}
		}
		setTxtProgressBar(pb, length(compare))
		close(pb)

		lcb = length(collapsed.bins)
		lsc = length(ind_sumcols)
		collapsed.bins[(lcb+1):(lcb+lsc)] = as.data.frame(collapsed_sumcols)
	}

	# Create the return data.frame with same order as input data.frame
	collapsed.bins.reordered = NULL
	collapsed.bins.reordered[c(ind_coords,ind_morecols,ind_sumcols)] = collapsed.bins
	names(collapsed.bins.reordered) = names(data)

	return(as.data.frame(collapsed.bins.reordered))

}

