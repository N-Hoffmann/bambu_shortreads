prepareShortReadsFromBam <- function(shortreadBam, fasta, minWidth = 70, maxWidth = 100000){
    shortReads <- prepareDataFromBam(shortreadBam, verbose = verbose)
    uniqueShortReadsIntrons <- unique(unlistIntrons(shortReads, use.ids = FALSE))
    rm(shortReads)

    #Filtering short read introns below minimum reference width
    uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(uniqueShortReadsIntrons@ranges@width < minWidth)),]

    #Filtering short read introns above a maximum intron width
    uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(uniqueShortReadsIntrons@ranges@width >= maxWidth)),]

    #Filtering short read introns with low short read support
    rawbam <- GenomicAlignments::readGAlignments(shortreadBam, param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE)))
    raw_juncs <- data.frame(GenomicAlignments::summarizeJunctions(rawbam))
    rm(rawbam)
    sr_exclusive_df <- data.frame(anti_join(raw_juncs, data.frame(uniqueAnnotatedIntrons)))
    x <- left_join(data.frame(uniqueShortReadsIntrons), sr_exclusive_df[,-5])
    uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-which(x$score < 1 ),] #score set to 1, can be changed

    #Filtering introns with non canonical splice sites
    uniqueShortReadsIntrons <- get_splice(uniqueShortReadsIntrons, fasta = fasta)
    uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(mcols(uniqueShortReadsIntrons)$splice != c("GTAG", "CTAC"))),]
    mcols(uniqueShortReadsIntrons) <- NULL

    return(uniqueShortReadsIntrons)
}

#'update shortReadsGranges object with splice sites
#' to metadata
get_splice <- function(shortReadsGranges, fasta){
    fa <- FaFile(fasta)
    junction_start <- getSeq(fa,
                         resize(shortReadsGranges, width = 2, fix='start'))
    junction_end <- getSeq(fa,
                         resize(shortReadsGranges, width = 2, fix='end'))
    mcols(shortReadsGranges)$splice <- paste(junction_start, junction_end, sep = "")
    return(shortReadsGranges)
}