bait.loadSCE <- function(filename) {
  
  data <- read.table(filename, header=FALSE)
  sce <- GRanges(seqnames=data[,1], ranges=IRanges(start=data[,2], end=data[,3]), ID=data[,4])
	sce.split <- split(sce, sce$ID)
  return(sce.split)
  
}

bait.loadSWITCH <- function(filename) {
  
  data <- read.table(filename, header=TRUE)
  sce <- GRanges(seqnames=data[,1], ranges=IRanges(start=data[,2], end=data[,3]), ID=data[,4])
  sce.split <- split(sce, sce$ID)
  return(sce.split)
  
}


# ------------------------------------------------------------
# Plot state categorization for all chromosomes
# ------------------------------------------------------------
bait.plot <- function(model, both.strands=FALSE, percentages=TRUE, file=NULL, bait.scecoord=NULL) {
  
  ## Convert to GRanges
  gr <- model$bins
  grl <- split(gr, seqnames(gr))
  
  ## Get some variables
  num.chroms <- length(levels(seqnames(gr)))
  maxseqlength <- max(seqlengths(gr))
  # 	if (!is.null(model$distributions)) {
  # 		custom.xlim <- model$distributions['disomy','mu'] * 4
  # 	} else {
  tab <- table(gr$reads)
  tab <- tab[names(tab)!='0']
  custom.xlim <- as.numeric(names(tab)[which.max(tab)]) * 4
  # 	}
  if (both.strands) {
    custom.xlim <- custom.xlim / 2
  }
  
  ## Setup page
  fs_title <- 20
  fs_x <- 13
  nrows <- 2	# rows for plotting chromosomes
  nrows.text <- 2	# two additional rows for displaying ID and qualityInfo
  nrows.total <- nrows + nrows.text
  ncols <- ceiling(num.chroms/nrows)
  if (!is.null(file)) {
    pdf(file=file, width=ncols*1.4, height=nrows*4.6)
  }
  grid.newpage()
  layout <- matrix(1:((nrows.total)*ncols), ncol=ncols, nrow=nrows.total, byrow=T)
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights=c(1,1,21,21))))
  # Main title
  grid.text(model$ID, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_title))
  # Quality info
  quality.string <- paste0('complexity = ',round(model$qualityInfo$complexity),',  spikyness = ',round(model$qualityInfo$spikyness,2),',  entropy = ',round(model$qualityInfo$shannon.entropy,2))
  grid.text(quality.string, vp = viewport(layout.pos.row = 2, layout.pos.col = 1:ncols), gp=gpar(fontsize=fs_x))
  
  ## Get SCE coordinates
  if (both.strands) {
    scecoords <- getSCEcoordinates(model)
  }
  
  ## Go through chromosomes and plot
  for (i1 in 1:num.chroms) {
    # Get the i,j matrix positions of the regions that contain this subplot
    matchidx <- as.data.frame(which(layout == i1+nrows.text*ncols, arr.ind = TRUE))
    if (!is.null(grl[[i1]]$state) & percentages) {
      # Percentage of chromosome in state
      tstate <- table(mcols(grl[[i1]])$state)
      pstate.all <- tstate / sum(tstate)
      pstate <- round(pstate.all*100)[-1]	# without 'nullsomy / unmapped' state
      pstring <- apply(pstate, 1, function(x) { paste0(": ", x, "%") })
      pstring <- paste0(names(pstring), pstring)
      pstring <- paste(pstring[which.max(pstate)], collapse="\n")
      pstring2 <- round(pstate.all*100)[1]	# only 'nullsomy / unmapped'
      pstring2 <- paste0(names(pstring2), ": ", pstring2, "%")
    } else {
      pstring <- ''
      pstring2 <- ''
    }
    
    # Plot the read counts
    dfplot <- as.data.frame(grl[[i1]])
    # Set values too big for plotting to limit
    dfplot$reads[dfplot$reads>=custom.xlim] <- custom.xlim
    dfplot.points <- dfplot[dfplot$reads>=custom.xlim,]
    dfplot.points$reads <- rep(custom.xlim, nrow(dfplot.points))
    
    if (both.strands) {
      dfplot$mreads <- - dfplot$mreads	# negative minus reads
      dfplot$preads[dfplot$preads>=custom.xlim] <- custom.xlim
      dfplot$mreads[dfplot$mreads<=-custom.xlim] <- -custom.xlim
      dfplot.points.plus <- dfplot[dfplot$preads>=custom.xlim,]
      dfplot.points.plus$reads <- rep(custom.xlim, nrow(dfplot.points.plus))
      dfplot.points.minus <- dfplot[dfplot$mreads<=-custom.xlim,]
      dfplot.points.minus$reads <- rep(-custom.xlim, nrow(dfplot.points.minus))
    }
    
    empty_theme <- theme(axis.line=element_blank(),
                         axis.text.x=element_blank(),
                         axis.text.y=element_blank(),
                         axis.ticks=element_blank(),
                         axis.title.x=element_text(size=fs_x),
                         axis.title.y=element_blank(),
                         legend.position="none",
                         panel.background=element_blank(),
                         panel.border=element_blank(),
                         panel.grid.major=element_blank(),
                         panel.grid.minor=element_blank(),
                         plot.background=element_blank())
    ggplt <- ggplot(dfplot, aes_string(x='start', y='reads'))	# data
    if (!is.null(grl[[i1]]$state)) {
      if (both.strands) {
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='preads', col='pstate'), size=0.2)	# read count
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mreads', col='mstate'), size=0.2)	# read count
        ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='reads', col='pstate'), size=5, shape=21)	# outliers
        ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='reads', col='mstate'), size=5, shape=21)	# outliers
      } else {
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='reads', col='state'), size=0.2)	# read count
        ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='reads', col='state'), size=2, shape=21)	# outliers
      }
      ggplt <- ggplt + scale_color_manual(values=state.colors, drop=F)	# do not drop levels if not present
    } else {
      if (both.strands) {
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='preads'), size=0.2, col='gray20')	# read count
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='mreads'), size=0.2, col='gray20')	# read count
        ggplt <- ggplt + geom_point(data=dfplot.points.plus, mapping=aes_string(x='start', y='reads'), size=5, shape=21, col='gray20')	# outliers
        ggplt <- ggplt + geom_point(data=dfplot.points.minus, mapping=aes_string(x='start', y='reads'), size=5, shape=21, col='gray20')	# outliers
      } else {
        ggplt <- ggplt + geom_linerange(aes_string(ymin=0, ymax='reads'), size=0.2, col='gray20')	# read count
        ggplt <- ggplt + geom_point(data=dfplot.points, mapping=aes_string(x='start', y='reads'), size=2, shape=21, col='gray20')	# outliers
      }
    }
    if (both.strands) {
      ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim, ymax=0.05*custom.xlim, xmin=0, mapping=aes(xmax=max(start)), col='white', fill='gray20')	# chromosome backbone as simple rectangle
    } else {
      ggplt <- ggplt + geom_rect(ymin=-0.05*custom.xlim-0.1*custom.xlim, ymax=-0.05*custom.xlim, xmin=0, mapping=aes(xmax=max(start)), col='white', fill='gray20')	# chromosome backbone as simple rectangle
    }
    if (both.strands) {
      dfsce <- as.data.frame(scecoords[seqnames(scecoords)==names(grl)[i1]])
      if (nrow(dfsce)>0) {
        ggplt <- ggplt + geom_segment(data=dfsce, aes(x=start, xend=start), y=-custom.xlim, yend=-0.5*custom.xlim, arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
      }
    }
		if (!is.null(bait.scecoord)) {
			chrom.string <- if (!grepl('chr', names(grl)[i1])) { paste0('chr', names(grl)[i1]) }
      dfsce <- as.data.frame(bait.scecoord[seqnames(bait.scecoord)==chrom.string])
      if (nrow(dfsce)>0) {
        ggplt <- ggplt + geom_segment(data=dfsce, aes(x=start, xend=start), y=custom.xlim, yend=0.5*custom.xlim, color='gray30', arrow=arrow(length=unit(0.5, 'cm'), type='closed'))
      }
		}
    ggplt <- ggplt + empty_theme	# no axes whatsoever
    ggplt <- ggplt + ylab(paste0(seqnames(grl[[i1]])[1], "\n", pstring, "\n", pstring2))	# chromosome names
    if (both.strands) {
      ggplt <- ggplt + xlim(0,maxseqlength) + ylim(-custom.xlim,custom.xlim)	# set x- and y-limits
    } else {
      ggplt <- ggplt + xlim(0,maxseqlength) + ylim(-0.6*custom.xlim,custom.xlim)	# set x- and y-limits
    }
    ggplt <- suppressMessages( ggplt + coord_flip() + scale_x_reverse() )
    suppressWarnings(
      print(ggplt, vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
    )
    
  }
  if (!is.null(file)) {
    d <- dev.off()
  }
}
