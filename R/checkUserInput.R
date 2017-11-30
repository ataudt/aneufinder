

check.positive.integer = function(testvar) {
	if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
	if (length(testvar)>1) return(2)
	if (testvar <= 0) return(3)
	if (testvar != as.integer(testvar)) return(4)
	return(0)
}
check.positive.integer.vector = function(testvec) {
	if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
	for (elem in testvec) {
		if (elem <= 0) return(2)
		if (elem != as.integer(elem)) return(3)
	}
	return(0)
}
check.nonnegative.integer.vector = function(testvec) {
	if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
	for (elem in testvec) {
		if (elem < 0) return(2)
		if (elem != as.integer(elem)) return(3)
	}
	return(0)
}
check.positive = function(testvar) {
	if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
	if (length(testvar)>1) return(2)
	if (testvar <= 0) return(3)
	return(0)
}
check.positive.vector = function(testvec) {
	if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
	for (elem in testvec) {
		if (elem <= 0) return(2)
	}
	return(0)
}
check.nonnegative.vector = function(testvec) {
	if (!is(testvec,"numeric") & !is(testvec,"integer")) return(1)
	for (elem in testvec) {
		if (elem < 0) return(2)
	}
	return(0)
}
check.integer = function(testvar) {
	if (!is(testvar,"numeric") & !is(testvar,"integer")) return(1)
	if (length(testvar)>1) return(2)
	if (testvar != as.integer(testvar)) return(3)
	return(0)
}

check.univariate.modellist = function(modellist) {
	if (!is(modellist,"list")) return(1)
	for (model in modellist) {
		if (!is(model,"aneuHMM")) return(2)
	}
	return(0)
}
check.univariate.model = function(model) {
	if (!is(model,"aneuHMM")) return(1)
	return(0)
}
check.multivariate.model = function(model) {
	if (!is(model,"aneuMultiHMM")) return(1)
	return(0)
}

check.logical = function(testbool) {
	if (!is(testbool,"logical")) return(1)
	if (length(testbool)>1) return(2)
	return(0)
}

check.strand = function(teststrand) {
	if (teststrand!='+' & teststrand!='-' & teststrand!='*') return(1)
	return(0)
}
