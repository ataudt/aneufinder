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
		if (!is(model,class.univariate.hmm)) return(2)
	}
	return(0)
}
check.univariate.model = function(model) {
	if (!is(model,class.univariate.hmm)) return(1)
	return(0)
}
check.multivariate.model = function(model) {
	if (!is(model,class.multivariate.hmm)) return(1)
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
