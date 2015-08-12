bam2binned(binfile, binsize=500000, min.mapq=35) -> binned.data
findCNVs(binned.data, ID=1) -> model
attr(model$bins, 'fragments') -> gr
dgr <- disjoin(gr)
dgr$overlap <- countOverlaps(dgr, gr)
dgr$state <- model$segments$state[findOverlaps(dgr, model$segments, select='first')]
dgr <- dgr[seqnames(dgr) %in% c(1:22,'X')]
dgr <- dgr[!is.na(dgr$state)]
dgr$multiplicity <- aneufinder:::initializeStates(states=levels(dgr$state))$multiplicity[as.character(dgr$state)]
dgr[dgr$overlap>dgr$multiplicity & dgr$state!='multisomy' & seqnames(dgr)=='X']
