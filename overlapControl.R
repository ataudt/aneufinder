## Blacklist
dac <- bed2GRanges('inst/extdata//wgEncodeDacMapabilityConsensusExcludable.bed.gz')
duke <- bed2GRanges('inst/extdata/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
blacklist <- union(dac,duke)

binfiles <- list.files('/media/aaron/seagate/DATA/work_ERIBA/aneuploidy/DATA_ANALYSIS/HILDA/cellscan_HB150623_1/', full=T, pattern='bam$')
dgrs <- list()
for (binfile in binfiles) {
  bam2binned(binfile, binsize=500000, min.mapq=42, return.fragments=T) -> fragments
  gr <- setdiff(fragments, subsetByOverlaps(fragments, blacklist))
  dgr <- disjoin(gr)
  dgr$overlap <- countOverlaps(dgr, gr)
  dgrs[[basename(binfile)]] <- dgr
}

## Ensembl
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', host='grch37.ensembl.org', dataset='hsapiens_gene_ensembl')
filters <- listFilters(ensembl)
attributes <- listAttributes(ensembl)
