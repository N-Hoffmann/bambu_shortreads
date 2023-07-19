#' Isoform reconstruction using genomic alignments
#' @param unlisted_junctions unlisted_junctions
#' @param annotations annotations
#' @param genomeSequence genomeSequence
#' @param stranded stranded
#' @param verbose verbose
#' @inheritParams bambu
#' @importFrom GenomicRanges match
#' @importFrom dplyr tibble %>% mutate select
#' @noRd
isore.constructJunctionTables <- function(unlisted_junctions, annotations, shortReads,
    genomeSequence, stranded = FALSE, verbose = FALSE, combined = FALSE,
    juncDist = 10,
    intron_limit = NULL) {
    start.ptm <- proc.time()
    if(length(unlisted_junctions)==0) return(NULL)
    #summarise junction counts and strand for all reads
    uniqueJunctions <- createJunctionTable(unlisted_junctions,
        genomeSequence = genomeSequence)
    end.ptm <- proc.time()
    if (verbose) message("Finished creating junction list with splice motif ",
        "in ", round((end.ptm - start.ptm)[3] / 60, 1), " mins.")

    uniqueAnnotatedIntrons <- unique(unlistIntrons(annotations, 
        use.ids = FALSE))

    #Processing short reads
    if (length(shortReads) > 0){
        message("Using shortReads for correction")
        uniqueShortReadsIntrons <- unique(unlistIntrons(shortReads, use.ids = FALSE))
        #Removing strand information and duplicates in uniqueShortReadsIntrons
        z <- data.frame(uniqueShortReadsIntrons)
        z$strand <- NULL
        uniqueShortReadsIntrons <- unique(makeGRangesFromDataFrame(z))

        #Filtering short read introns below minimum reference width
        uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(uniqueShortReadsIntrons@ranges@width < 70)),]
        #Filtering short read introns with low short read support
        #rawbam <- GenomicAlignments::readGAlignments("chrIS_shortreads.bam", param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE)))
        rawbam <- GenomicAlignments::readGAlignments("SGNex_Hct116_Illumina_replicate3_run1.bam", param = ScanBamParam(flag = scanBamFlag(isSecondaryAlignment = FALSE)))
        raw_juncs <- data.frame(GenomicAlignments::summarizeJunctions(rawbam))
        rm(rawbam)
        sr_exclusive_df <- data.frame(anti_join(raw_juncs, data.frame(uniqueAnnotatedIntrons)))
        x <- left_join(data.frame(uniqueShortReadsIntrons), sr_exclusive_df[,-5])
        uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-which(x$score < 5 ),] #score set to 5, can be changed
        #Filtering short read introns above a maximum intron width
        if (!is.null(intron_limit)){
            uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(uniqueShortReadsIntrons@ranges@width >= intron_limit)),]
        }
        #Filtering introns with non canonical splice sites
        # uniqueShortReadsIntrons <- get_splice(uniqueShortReadsIntrons)
        # uniqueShortReadsIntrons <- uniqueShortReadsIntrons[-(which(mcols(uniqueShortReadsIntrons)$splice != "GTAG")),]
        # mcols(uniqueShortReadsIntrons) <- NULL

        if(combined == TRUE){
            message("Using combination of shortReads and Annotations")
            #uniqueAnnotatedIntrons <- SparseSummarizedExperiment::combine(uniqueShortReadsIntrons, uniqueAnnotatedIntrons)
            #uniqueAnnotatedIntrons <- unique(c(uniqueShortReadsIntrons, uniqueAnnotatedIntrons))

            #Assigning origin of junction when using combination of shortreads and annotation
            combined_df <- data.frame(unique(c(uniqueShortReadsIntrons, uniqueAnnotatedIntrons)))
            anno_exclusive <- anti_join(combined_df, data.frame(uniqueShortReadsIntrons), by=c("seqnames","start","end","width")) #SR are unstranded, remove strand info
            if(nrow(anno_exclusive) > 0){
                anno_exclusive$source = "Annotation"
            }
            sr_exclusive <- anti_join(combined_df, data.frame(uniqueAnnotatedIntrons), by=c("seqnames","start","end","width"))
            if(nrow(sr_exclusive) > 0){
                sr_exclusive$source = "ShortReads"
            }
            #both <- semi_join(data.frame(uniqueShortReadsIntrons), data.frame(uniqueAnnotatedIntrons))
            both <- semi_join(data.frame(uniqueAnnotatedIntrons), data.frame(uniqueShortReadsIntrons), by=c("seqnames","start","end","width"))
            if(nrow(both) > 0){
                both$source = "Both"
            }
            sources <- c(anno_exclusive$source, sr_exclusive$source, both$source)
            df_list <- list(anno_exclusive[,1:5], sr_exclusive[,1:5], both[,1:5])
            rebuild_df <- (df_list %>% purrr::reduce(full_join))
            uniqueAnnotatedIntrons <- makeGRangesFromDataFrame(rebuild_df)
            mcols(uniqueAnnotatedIntrons)$source = sources
        } else{
        message("Using shortReads only as annotations")
        uniqueAnnotatedIntrons <- uniqueShortReadsIntrons
        }
    }
    # end of short read processing

    # ### Splitting Annotation into half, make bad annotations
    # indexes <- sample(c(TRUE,FALSE), length(uniqueAnnotatedIntrons), replace = TRUE, prob = c(0.5,0.5))
    # correct_annot <- uniqueAnnotatedIntrons[indexes]
    # wrong_annot <- uniqueAnnotatedIntrons[!indexes]
    # wrong_annot <- shift(wrong_annot, 10)
    # mcols(correct_annot)$wrong <- "FALSE"
    # mcols(wrong_annot)$wrong <- "TRUE"
    # uniqueAnnotatedIntrons <- unique(c(correct_annot, wrong_annot))

    # correct strand of junctions based on (inferred) strand of reads
    strand(uniqueJunctions) <- junctionStrandCorrection(uniqueJunctions,
        unlisted_junctions, uniqueAnnotatedIntrons,
        stranded = stranded, verbose = verbose)
    # add annotation labels to junctions
    mcols(uniqueJunctions) <- tibble(as.data.frame(uniqueJunctions)) %>% 
        mutate(annotatedJunction = (!is.na(match(uniqueJunctions,
            uniqueAnnotatedIntrons)))) %>% group_by(seqnames) %>% 
        mutate(annotatedStart = start %in% start[annotatedJunction],
            annotatedEnd = end %in% end[annotatedJunction]) %>% ungroup() %>%
        dplyr::select(score, spliceMotif, spliceStrand, junctionStartName, 
            junctionEndName, startScore, endScore, id, annotatedJunction,
            annotatedStart, annotatedEnd)
    # correct junction coordinates using logistic regression classifier
    uniqueJunctions <- junctionErrorCorrection(uniqueJunctions, verbose, juncDist)
    return(uniqueJunctions)
}



#' Get unlisted intron ranges from exon ranges list
#' @importFrom GenomicRanges GRanges
#' @noRd
unlistIntrons <- function(x, use.ids = TRUE, use.names = FALSE) {
    # License note: This function is adopted from the GenomicAlignments 
    # package (Author: Hervé Pagès, Valerie Obenchain, Martin Morgan)
    # License Artistic-2.0
    # https://doi.org/doi:10.18129/B9.bioc.GenomicAlignments
    
    flat <- unlist(x, use.names = FALSE)
    gaps <- gaps(ranges(x))
    
    firstseg <- start(PartitioningByWidth(x))
    seqnms <- rep(seqnames(flat)[firstseg], elementNROWS(gaps))
    strand <- rep(strand(flat)[firstseg], elementNROWS(gaps))
    
    gr <- GRanges(seqnms, unlist(gaps, use.names = use.names), strand)
    if (use.ids & !is.null(mcols(x, use.names = FALSE)$id)) 
        mcols(gr)$id <- rep(mcols(x)$id, elementNROWS(gaps))
    return(gr)
}



#' Create Junction tables from unlisted junction granges
#' @importFrom BiocGenerics unstrand
#' @importFrom IRanges shift
#' @importFrom dplyr tibble group_by %>% mutate ungroup select
#' @noRd
createJunctionTable <- function(unlisted_junctions,
    genomeSequence = NULL) {
    # License note: This function is adopted from the GenomicAlignments package 
    uniqueJunctions <- sort(unique(unstrand(unlisted_junctions)))
    names(uniqueJunctions) <- paste("junc", seq_along(uniqueJunctions),
        sep = ".")
    plus_score <- countMatches(uniqueJunctions,
        unlisted_junctions[strand(unlisted_junctions) == '+'], 
        ignore.strand = TRUE)
    minus_score <- countMatches(uniqueJunctions,
        unlisted_junctions[strand(unlisted_junctions) == '-'], 
        ignore.strand = TRUE)

    junctionSeqStart <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2), 2))#shift from IRanges
    junctionSeqEnd <- BSgenome::getSeq(genomeSequence,
        IRanges::shift(flank(uniqueJunctions,width = 2, start = FALSE), -2))
    
    mcols(uniqueJunctions) <- DataFrame(tibble(
        chr = as.factor(seqnames(uniqueJunctions)), 
        start = start(uniqueJunctions),
        end = end(uniqueJunctions),
        score = plus_score + minus_score,
        plus_score = plus_score,
        minus_score = minus_score,
        spliceMotif = paste(junctionSeqStart, junctionSeqEnd, sep = "-"),
        spliceStrand = spliceStrand(spliceMotif),
        junctionStartName = paste(chr, start, sep = ":"),
        junctionEndName = paste(chr, end, sep = ":"),
        id = seq_along(uniqueJunctions)) %>%
        group_by(chr, start) %>% 
        mutate(startScore = sum(score)) %>% 
        group_by(chr, end) %>%  
        mutate(endScore = sum(score)) %>%
        ungroup() %>%
        dplyr::select(score, plus_score, minus_score, spliceMotif, spliceStrand,
            junctionStartName, junctionEndName, startScore, endScore, id))
    strand(uniqueJunctions) <- uniqueJunctions$spliceStrand
    return(uniqueJunctions)
}

#' @param motif motif
#' @noRd
spliceStrand <- function(motif) {
    NATURAL_INTRON_MOTIFS_RC <- as.character(Biostrings::reverseComplement(
        Biostrings::DNAStringSet(GenomicAlignments::NATURAL_INTRON_MOTIFS)))
    
    motifStrand <- ifelse(motif %in% GenomicAlignments::NATURAL_INTRON_MOTIFS,
        "+", "*")
    motifStrand[motif %in% NATURAL_INTRON_MOTIFS_RC] <- "-"
    return(motifStrand)
}


#' JUNCTIONSTRANDCORRECTION
#' @noRd
junctionStrandCorrection <- function(uniqueJunctions, unlisted_junctions,
    uniqueAnnotatedIntrons, stranded, verbose = FALSE) {
    # note: strand sometimes incorrectly infered based on motifs, might 
    # introduce systematic errors due to alignment (biased to splice motifs)
    
    uniqueJunctionsUpdate <- uniqueJunctions
    # make a copy to revert to if strand correction does not improve results
    annotatedIntronNumber <- evalAnnotationOverlap(uniqueJunctions,
        uniqueAnnotatedIntrons,ignore.strand = FALSE)["TRUE"]
    if (verbose) {
        message("before strand correction, annotated introns: ", 
        annotatedIntronNumber, " (", annotatedIntronNumber / length(uniqueJunctions),"%)")
    }
    # infer strand for each read based on strand of junctions
    strandStep <- TRUE
    while (strandStep) { # iterate twice to improve strand prediction w.t.
        # mean junction counts, annotate junction strand with read strand
        if (!stranded) { # update junction strand score
            strandScoreByRead <- updateStrandScoreByRead(unlisted_junctions,
                uniqueJunctionsUpdate)
        } else {# just use strand from reads for stranded data
            strandScoreByRead <- uniqueJunctionsUpdate$minus_score -
                uniqueJunctionsUpdate$plus_score
        }
        # overwrite info from motif which increases overlap with known junc
        strand(uniqueJunctionsUpdate[strandScoreByRead < 0]) <- "+"
        strand(uniqueJunctionsUpdate[strandScoreByRead > 0]) <- "-"
        updatedList <- updateJunctionwimprove(annotatedIntronNumber,
            uniqueJunctions, uniqueJunctionsUpdate, uniqueAnnotatedIntrons,
            strandStep, verbose)
        annotatedIntronNumber <- updatedList$annotatedIntronNumber
        uniqueJunctions <- updatedList$uniqueJunctions
        strandStep <- updatedList$strandStep
    }
    return(strand(uniqueJunctions))
}


#' Evaluate annoation overlap
#' @importFrom GenomicRanges match
#' @noRd
evalAnnotationOverlap <- function(intronRanges, uniqueAnnotatedIntrons,
    ignore.strand = FALSE) {
    return(table(!is.na(GenomicRanges::match(intronRanges,
        uniqueAnnotatedIntrons, ignore.strand = ignore.strand ))))
}


#' This function assigns a strand to each read based on the majority of 
#' junctions. The strand of the junctions is infered by the sequence in 
#' createJunctionTables
#' @noRd
updateStrandScoreByRead <- function(unlisted_junctions, uniqueJunctions){
    allJunctionToUniqueJunctionMatch <- match(unlisted_junctions,
        uniqueJunctions, ignore.strand = TRUE)

    unlisted_junction_granges_strandList <-
        splitAsList(strand(uniqueJunctions)[
        allJunctionToUniqueJunctionMatch],
        mcols(unlisted_junctions)$id)

    strandJunctionSum <-
        as.integer(sum(unlisted_junction_granges_strandList == "-") -
        sum(unlisted_junction_granges_strandList == "+"))

    readStrand <- factor(rep("*", length(unlisted_junction_granges_strandList)),
        levels = c('+','-','*'))
    readStrand[strandJunctionSum < 0] <- "+"
    readStrand[strandJunctionSum > 0] <- "-"

    strand_unlisted_junctions <-
        readStrand[match(mcols(unlisted_junctions)$id,
        as.integer(names(unlisted_junction_granges_strandList)))]
    plus_score <- countMatches(uniqueJunctions,
        unlisted_junctions[strand_unlisted_junctions == '+'], 
        ignore.strand = TRUE)
    minus_score <- countMatches(uniqueJunctions,
        unlisted_junctions[strand_unlisted_junctions == '-'], 
        ignore.strand = TRUE)
    return(minus_score - plus_score)
}



#' update junctions object if strand prediction improves overlap 
#' with annotations
#' @param annotatedIntronNumber annotatedIntronNumber
#' @param uniqueJunctions uniqueJunctions
#' @param uniqueJunctionsUpdate uniqueJunctionsUpdate
#' @param uniqueAnnotatedIntrons uniqueAnnotatedIntrons
#' @param strandStep strandStep
#' @param verbose verbose
#' @noRd
updateJunctionwimprove <- function(annotatedIntronNumber, uniqueJunctions,
    uniqueJunctionsUpdate, uniqueAnnotatedIntrons, strandStep, verbose) {
    annotatedIntronNumberNew <- evalAnnotationOverlap(uniqueJunctionsUpdate,
        uniqueAnnotatedIntrons, ignore.strand = FALSE)["TRUE"]
    if (annotatedIntronNumberNew > annotatedIntronNumber & !is.na(
        annotatedIntronNumber)) {
        # update junctions object if strand prediction improves overlap
        # with annotations
        if (verbose) {
            message("after strand correction, annotated introns: ", annotatedIntronNumberNew,
            " (", annotatedIntronNumberNew / length(uniqueJunctionsUpdate), ")")
        }
        annotatedIntronNumber <- annotatedIntronNumberNew
        uniqueJunctions <- uniqueJunctionsUpdate
    } else {
        strandStep <- FALSE
    }
    outputList <- list(
        "strandStep" = strandStep,
        "annotatedIntronNumber" = annotatedIntronNumber,
        "uniqueJunctions" = uniqueJunctions
    )
    return(outputList)    
}

#'update shortReadsGranges object with splice sites
#' to metadata
get_splice <- function(shortReadsGranges){
    #TODO add variable path to .fa
    #fa <- FaFile("hg38_sequins_SIRV_ERCCs_longSIRVs.fa")
    fa <- FaFile("chrIS_fasta.fa")
    junction_start <- getSeq(fa,
                         resize(shortReadsGranges, width = 2, fix='start'))
    junction_end <- getSeq(fa,
                         resize(shortReadsGranges, width = 2, fix='end'))
    mcols(shortReadsGranges)$splice <- paste(junction_start, junction_end, sep = "")
    return(shortReadsGranges)
}