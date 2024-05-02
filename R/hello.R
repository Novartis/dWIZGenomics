
#' @import AnnotationHub rtracklayer biovizBase cowplot ggbio ggExtra ggplot2 ggpubr ggrepel magrittr `org.Hs.eg.db` scales GenomicRanges IRanges
#' @export
plotRNAVolcano <- function(plotData){
  plotData <- plotData |>
    dplyr::filter( baseMean > 20 ) |>
    dplyr::mutate(
      padj=p.adjust(pvalue, method="BH"),
      logpvalue=-log10(pvalue),
      cl=ifelse( padj < 0.1 & abs(log2FoldChange) > 1, "yes", "no"),
      lab=ifelse(SYMBOL %in% c("HBG1", "HBG2", "HBB", "HBD"), SYMBOL, NA))
  yline <-
    min(plotData$logpvalue[which(plotData$padj < 0.1)])
  plotData |>
    ggplot( aes( log2FoldChange, logpvalue, label=lab, col=cl ) ) +
    geom_point( size=1.1 ) +
    geom_hline(yintercept=yline, color="black", linetype="dashed", alpha=0.3) +
    geom_vline(xintercept=c(-1, 1), color="black", size=0.4, linetype="dashed", alpha=0.6) +
    theme_cowplot(font_family="Helvetica") +
    theme(legend.position="none") +
    coord_cartesian(x=c(-3.3, 3.3), ylim=c(0, 30)) +
    scale_color_manual(values=c(no="#D3D3D390", yes="#52525250")) +
    geom_text_repel(data=plotData[!is.na(plotData$lab),],
                    nudge_x = c(5, 1.5, -1.5, -1.5), nudge_y=10,
                    max.overlaps = Inf, col="#ef3b2c", size=4) +
    geom_point( data=plotData[!is.na(plotData$lab),], size=0.8, color="#ef3b2c") +
    labs(x=expression(log[2]~"("~frac("dWIZ-2", DMSO)~")"),
         y=expression(-log[10]~"(p-value)"))
}

#' @importMethodsFrom BiocGenerics mean
#' @export
getPlottingData <- function( bigwigFile, plottingRegion, resolution=100, smoothIter=1 ){
  smoothWindowWidth <- resolution * 3
  plottingRegion <- resize(plottingRegion, width(plottingRegion)*1.5,fix="center")
  plottingRegion <- resize(plottingRegion, width(plottingRegion) + (smoothWindowWidth - (width(plottingRegion) %% smoothWindowWidth)), fix="center")
  plottingRegionTiles <- tile(plottingRegion, width=resolution)[[1]]
  smoothingTiles <- slidingWindows(plottingRegion, width=smoothWindowWidth, step=resolution)[[1]]
  covData <- import( bigwigFile, which=plottingRegion )
  seqlevels(covData) <- seqlevelsInUse(covData)
  binnedCov <- binnedAverage(plottingRegionTiles, coverage( covData, weight=covData$score ), varname="cov")
  while( smoothIter > 0 ){
    ovl <- findOverlaps(binnedCov, smoothingTiles)
    smoothingTiles$cov <- mean(NumericList(split(binnedCov$cov[queryHits(ovl)], subjectHits(ovl))))
    smoothedTiles <- resize( smoothingTiles, width=resolution, fix="center" )
    ovl2 <- findOverlaps( binnedCov, smoothedTiles, type="equal" )
    binnedCov$cov[queryHits(ovl2)] <- smoothedTiles$cov[subjectHits(ovl2)]
    binnedCov$cov[1] <- binnedCov$cov[2]
    binnedCov$cov[length(binnedCov)-1] <- binnedCov$cov[length(binnedCov)]
    smoothIter <- smoothIter - 1
  }
  binnedCov
}

#' @export
importAndAverage <- function( inputDF, plottingRegion, resolution=100, smoothIter=2 ){
  covData <- lapply( inputDF$bigwigfile, function(xx){
    covData <- getPlottingData( xx, plottingRegion, resolution=resolution, smoothIter=smoothIter )
    covData <- as.data.frame(covData)
    covData
  })
  names(covData) <- inputDF$sample_alias
  for( i in names(covData) ){
    covData[[i]]$sample_alias <- i
  }
  covData <- dplyr::bind_rows( covData )
  covData$replicate_group <- gsub("_R\\d+$", "", covData$sample_alias)
  covData <- covData %>%
    dplyr::group_by( seqnames, start, end, replicate_group ) %>%
    dplyr::summarize( cov=mean(cov) )
  covData
}

#' @export
plotCovFromDF <- function( covData, wh="WIZ D2", plottingRegion, ylim=NULL, cl="black", yexpand=TRUE ){
  covPlot <- covData %>%
    dplyr::filter( replicate_group %in% wh ) %>%
    ggplot( ) +
    geom_rect(aes(xmin = start-1, xmax = end+1, ymax = cov, ymin=0 ), fill=cl) +
    scale_x_reverse( limits=c(end(plottingRegion), start(plottingRegion)), expand = c(0, 0)) +
    coord_cartesian( ylim=ylim, expand=yexpand ) +
    theme_cowplot() +
    labs(y=wh)
  covPlot
}

