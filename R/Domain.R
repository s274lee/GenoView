# Establish PFAM mapping objects for initial PFAM file
makePFAMObjs <- function() {
    x <- PFAM.db::PFAMDE
    pfam.desc <- AnnotationDbi::contents(x)
    y <- PFAM.db::PFAMID
    pfam.ids <- unlist(AnnotationDbi::contents(y))
    return(list(desc = pfam.desc, ids = pfam.ids))
}

# Processing the raw ucscGenePfam table to obtain domain information
updatePFAM <- function(raw.file = "ucscGenePfam.txt", 
                       search.dir = getwd(), 
                       file.name = "", ID2DE = FALSE) {
    
    # Save original directory
    first.dir <- getwd()
    setwd(search.dir)
    
    if (file.exists(raw.file)){
        
        # Extract the relevant columns
        raw.pfam <- read.delim(raw.file, sep = '\t', header = FALSE)
        pfam.df <- raw.pfam[,c(5, 2, 3, 4, 7)]
        colnames(pfam.df) <- c("PFAMID", "chr", "start", "end", "strand")
        
        if (ID2DE) {
            
            pfam.objs <- makePFAMObjs()
            pfam.desc <- pfam.objs$desc
            pfam.ids <- pfam.objs$ids
            
            # Make PFAM GRanges object with domain IDs as values 
            pfam.gr <- GRanges(seqnames = as.character(pfam.df[,"chr"]), 
                               IRanges(start = pfam.df[,"start"], 
                                       end = pfam.df[,"end"]), 
                               strand = pfam.df[,"strand"])
            values(pfam.gr) <- list(domain = 
                                        as.character(as.matrix(pfam.df[,1])))
            
            # Map IDs to desciptions
            pfam.df.upd <- biovizBase::mold(PFAMIDE(pfam.gr, pfam.desc, 
                                                    pfam.ids))
            
            if (file.name == "") {
                file.name = paste0("Updated", raw.file)
            }
            
            write.table(pfam.df.upd, file = file.name, sep = '\t', 
                        quote = FALSE, row.names = FALSE)
            
        } else {
            
            if (file.name == "") {
                file.name = paste0("Initial", raw.file)
            }
            
            write.table(pfam.df, file = file.name, 
                        sep = '\t', quote = FALSE, row.names = FALSE)
            
        }
        
    } else {
        warning("raw.file could not be found in search.dir")
        setwd(first.dir)
        return(FALSE)
    }
    
    setwd(first.dir)
    return(file.name)
}

# Change PFAM domain ID to PFAM domain description in GRanges transcripts
PFAMIDE <- function(transcripts, desc, ids) {
    domain <- values(transcripts)[["domain"]]
    pfam.ids <- apply(matrix(domain), 1, function(row) {
        names(ids)[match(row, ids)]
    })
    na.ind <- which(is.na(pfam.ids))
    if(length(na.ind) > 0) {
        # Don't translate values which mapped to NA
        good.ind <- which(!(1:length(domain) %in% na.ind))
        dom.desc <- domain
        dom.desc[good.ind] <- unlist(desc[pfam.ids[good.ind]])
    } else {
        dom.desc <- unlist(desc[pfam.ids])
    }
    values(transcripts)[["domain"]] <- dom.desc
    return(transcripts)
}
