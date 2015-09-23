

#' Chromosome lengths for different assemblies
#'
#' Get the chromosome lengths of different common assemblies. Only the standard chromosomes (no chrN_random) will be returned.
#'
#' @param assembly One of the following: \code{c('hg19','hg18','mm10','mm9','rn4')}.
#' @return A named vector with chromosome lengths.
#' @export
getChromLengths <- function(assembly='hg19') {

	if (assembly=='hg19') {
		chrom.lengths <- c(
			chr1	=	249250621,
			chr2	=	243199373,
			chr3	=	198022430,
			chr4	=	191154276,
			chr5	=	180915260,
			chr6	=	171115067,
			chr7	=	159138663,
			chr8	=	146364022,
			chr9	=	141213431,
			chr10	=	135534747,
			chr11	=	135006516,
			chr12	=	133851895,
			chr13	=	115169878,
			chr14	=	107349540,
			chr15	=	102531392,
			chr16	=	90354753,
			chr17	=	81195210,
			chr18	=	78077248,
			chr19	=	59128983,
			chr20	=	63025520,
			chr21	=	48129895,
			chr22	=	51304566,
			chrX	=	155270560,
			chrY	=	59373566,
			chrM	=	16571
		)
	} else if (assembly=='hg18') {
		chrom.lengths <- c(
			chr1	=	247249719,
			chr2	=	242951149,
			chr3	=	199501827,
			chr4	=	191273063,
			chr5	=	180857866,
			chr6	=	170899992,
			chr7	=	158821424,
			chr8	=	146274826,
			chr9	=	140273252,
			chr10	=	135374737,
			chr11	=	134452384,
			chr12	=	132349534,
			chr13	=	114142980,
			chr14	=	106368585,
			chr15	=	100338915,
			chr16	=	88827254 ,
			chr17	=	78774742 ,
			chr18	=	76117153 ,
			chr19	=	63811651 ,
			chr20	=	62435964,
			chr21	=	46944323,
			chr22	=	49691432,
			chrX	=	154913754,
			chrY	=	57772954,
			chrM	=	16571
		)
	} else if (assembly=='mm10') {
		chrom.lengths <- c(
			chr1    = 195471971,
			chr2    = 182113224,
			chr3    = 160039680,
			chr4    = 156508116,
			chr5    = 151834684,
			chr6    = 149736546,
			chr7    = 145441459,
			chr8    = 129401213,
			chr9    = 124595110,
			chr10   = 130694993,
			chr11   = 122082543,
			chr12   = 120129022,
			chr13   = 120421639,
			chr14   = 124902244,
			chr15   = 104043685,
			chr16   = 98207768,
			chr17   = 94987271,
			chr18   = 90702639,
			chr19   = 61431566,
			chrX    = 171031299,
			chrY    = 91744698,
			chrM   = 16299
		)
	} else if (assembly=='mm9') {
		chrom.lengths <- c(
			chr1 = 197195432,
			chr2 = 181748087,
			chr3 = 159599783,
			chr4 = 155630120,
			chr5 = 152537259,
			chr6 = 149517037,
			chr7 = 152524553,
			chr8 = 131738871,
			chr9 = 124076172,
			chr10 = 129993255,
			chr11 = 121843856,
			chr12 = 121257530,
			chr13 = 120284312,
			chr14 = 125194864,
			chr15 = 103494974,
			chr16 = 98319150,
			chr17 = 95272651,
			chr18 = 90772031,
			chr19 = 61342430,
			chrX = 166650296,
			chrY = 15902555,
			chrM = 16299
		)
	} else if (assembly=='rn4') {
		chrom.lengths <- c(
			chr1 = 267910886,
			chr2 = 258207540,
			chr3 = 171063335,
			chr4 = 187126005,
			chr5 = 173096209,
			chr6 = 147636619,
			chr7 = 143002779,
			chr8 = 129041809,
			chr9 = 113440463,
			chr10 = 110718848,
			chr11 = 87759784,
			chr12 = 46782294,
			chr13 = 111154910,
			chr14 = 112194335,
			chr15 = 109758846,
			chr16 = 90238779,
			chr17 = 97296363,
			chr18 = 87265094,
			chr19 = 59218465,
			chr20 = 55268282,
			chrX = 160699376,
			chrM = 16300
		)
	}


	return(chrom.lengths)
}
