# Main function which selects mutation data and interface for human
mepHuman <- function(m.data, gui = FALSE) {
    
    genesymbol <- NULL
    hg19Ideogram <- NULL
    
    # Determine default published genomic datasets
    tx.db <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    seqs.hg19 <- seqlengths(tx.db)
    data(genesymbol, package = "biovizBase", envir = environment())
    
    # Create PFAM objects
    pfam.objs <- makePFAMObjs()
    pfam.desc <- pfam.objs$desc
    pfam.ids <- pfam.objs$ids
    
    # Load the pfam table
    data(pfam.df, package = "covert", envir = environment())
    pfam.gr <- GRanges(seqnames = as.character(pfam.df[,"chr"]), 
                       IRanges(start = pfam.df[,"start"], 
                               end = pfam.df[,"end"]), 
                       strand = pfam.df[,"strand"])
    values(pfam.gr) <- list(domain = as.character(as.matrix(pfam.df[,1])))
    
    if (gui) {
        mepGUI(m.data = m.data, tx.db = tx.db, gene.loc = genesymbol, 
               seq.lens = seqs.hg19, pfam.desc = pfam.desc, 
               pfam.ids = pfam.ids, pfam.gr = pfam.gr)
    } else {
        mepTxtInt(m.data = m.data, tx.db = tx.db, gene.loc = genesymbol, 
                  seq.lens = seqs.hg19, pfam.desc = pfam.desc, 
                  pfam.ids = pfam.ids, pfam.gr = pfam.gr)
    }
}

# Truncate an interval of an existing GRanges object
intSel <- function(gr, new.start = NULL, new.end = NULL) {
    if (is.null(new.start) && is.null(new.end)) {
        # Determine characteristics of that gene name
        len <- biovizBase::mold(gr)[2:3]
        cat(sprintf("Please enter the start location between %s: \n", 
                    paste(as.character(len), collapse = "-")))
        new.start <- scan(n=1)
        cat(sprintf("Please enter the end location between %s: \n", 
                    paste(as.character(c(new.start, len[2])), collapse = "-")))
        new.end <- scan(n=1)
    }
    
    # Replace old range
    ranges(gr) <- IRanges(new.start, new.end)
    return(gr)
}

# Make a filtered GRanges object from input dataFrame
makeGR <- function(dataFrame, chr.col, s.p.col, e.p.col, str.col, 
                   id.col = NULL, plot.int, des.seqs, show.legend = FALSE, 
                   l.height = 0) {
    
    # Discard all entries which are incomplete across the column criteria
    dataFrame <- dataFrame[complete.cases(dataFrame[, c(chr.col, s.p.col, 
                                                        e.p.col, str.col)]), 
                           c(chr.col, s.p.col, e.p.col, str.col, id.col)]
    
    # Truncate the input dataFrame based on selected interval
    if (length(which(dataFrame[,s.p.col] > end(plot.int) | 
                         dataFrame[,e.p.col] < start(plot.int))) > 0) {
        dataFrame <- dataFrame[-which(dataFrame[,s.p.col] > 
                                          end(plot.int) | 
                                          dataFrame[,e.p.col] <
                                          start(plot.int)), ]
    }
    
    # Create GRanges object from filtered dataFrame
    mut.gr <- GRanges(seqnames = as.character(dataFrame[, chr.col]), 
                      IRanges(start = dataFrame[, s.p.col], 
                              end = dataFrame[, e.p.col]), 
                      strand = as.character(dataFrame[, str.col]))
    
    if (show.legend | l.height > 0){
        values(mut.gr)$value <- dataFrame[, id.col]
    }
    
    # Determine seqlengths
    chr.sub <- as.character(unique(biovizBase::mold(mut.gr)[, "seqnames"]))
    seqlevels(mut.gr) <- chr.sub
    mut.gr <- keepSeqlevels(mut.gr, chr.sub)
    seqs.sub <- des.seqs[chr.sub]
    
    # remove intervals which extend past the end of the chr
    bidx <- end(mut.gr) <= seqs.sub[match(as.character(seqnames(mut.gr)),
                                          names(seqs.sub))]
    mut.gr <- mut.gr[which(bidx)]
    
    # Assign seqlengths
    seqlengths(mut.gr) <- seqs.sub
    
    return(mut.gr)
}

# Plot exon components and color according to a subset
plotExonRect <- function(plt, component = c("cds", "utr"), comp.df,
                         col.name = "type", filt.by.comp = FALSE, 
                         size = c(0.4, 0.15), subset.df, 
                         comp.col = "dark grey", 
                         sub.col = RColorBrewer::brewer.pal(3, "Set1")[2], 
                         aesCst) {
    
    # Add component layer with appropriate sizing
    component <- match.arg(component)
    if (filt.by.comp) {
        comp.df <- comp.df[which(comp.df[,col.name] == component),]
        subset.df <- subset.df[which(subset.df[,col.name] == component),]
    }
    if (nrow(comp.df) > 0 ) {
        plt <- plt + 
            ggplot2::geom_rect(data = comp.df, 
                               aesCst(xmin = comp.df[,"start"], 
                                      xmax = comp.df[,"end"], 
                                      ymin = comp.df[,"stepping"]-size,
                                      ymax = comp.df[,"stepping"]+size), 
                               fill = comp.col, colour = comp.col)
    }
    # Superimpose coloured layer in mutated subset
    if (nrow(subset.df) > 0) {
        plt <- plt + 
            ggplot2::geom_rect(data = subset.df, 
                               aesCst(xmin = subset.df[,"start"], 
                                      xmax = subset.df[,"end"], 
                                      ymin = subset.df[,"stepping"]-size, 
                                      ymax = subset.df[,"stepping"]+size), 
                               fill = sub.col, colour = sub.col)
    }
    return(plt)
}

# Declares warnings when there is a plotting problem
plotFailed <- function(txt = "", gui = FALSE, stat.bar = NULL) {
    
    if (gui) {
        # Create a popup warning
        gWidgets::gmessage(message = txt, title = "Plotting component failed",
                           icon = "error")
        
        gWidgets::svalue(stat.bar) <- paste(gWidgets::svalue(stat.bar), txt)
        
    } else {
        warning(txt, call. = FALSE, immediate. = TRUE)
    }
}

# Changes GRanges m.data to data.frame for column name selection
checkCols <- function(m.data) {
    output <- NULL
    if (class(m.data) == "GRanges") {
        output = biovizBase::mold(m.data)
    } else if (class(m.data) != "data.frame") {
        stop(sprintf("Input m.data of invalid class: %s", class(m.data)))
    } else {
        output = m.data
    }
    return(output)
}

# Checks input m.data and produces a GRanges object if applicable
checkMData <- function(m.data, ...) {
    mut.gr <- NULL
    if (class(m.data) == "data.frame") {
        mut.gr <- makeGR(dataFrame = m.data, ...)
    } else if (class(m.data) == "GRanges") {
        mut.gr <- m.data
        inp.lst <- list(...)
        if (inp.lst$l.height > 0){
            values(mut.gr)$value = values(m.data)[,inp.lst$id.col]
        }
    } else {
        stop(sprintf("Input m.data of invalid class: %s", class(m.data)))
    }
    return(mut.gr)
}

# Helper function which determines how each interface handles 
# failed domain and legend plotting
intPlotFailed <- function(ml, d.height, l.height, gui, ...) {
    
    # Declare errors
    domain.error = "No domains found within plotting interval.
    Domain plotting disabled. "
    legend.error = "No mutations overlapped with exons. 
    Legend plotting disabled. "
    
    ml.sd <- ml$show.domain
    ml.sl <- ml$show.legend
    # Domain is allocated but none found within the plotting interval
    if (d.height > 0 & !ml.sd) {
        plotFailed(txt = domain.error, gui = gui, ...)
    }
    
    # Legend is allocated but no overlap found within the plotting 
    # interval
    if (l.height > 0 & !ml.sl) {
        plotFailed(txt = legend.error, gui = gui, ...)
    }
}

# Check that m.data is in the right format before plotting
checkAndPlot <- function(m.data, Chromosome, Start.position, End.position, 
                         Strand, Fill, plot.int, l.height, seq.lens, 
                         exon.int, disp.track, p.height, 
                         d.height, plt.title, tx.db, gr, pfam.gr, 
                         pfam.desc, pfam.ids, ...){
    
    # Filter the m.data for relevant mutations and plotting data 
    mut.gr <- checkMData(m.data = m.data, chr.col = Chromosome, 
                         s.p.col = Start.position, e.p.col = End.position, 
                         str.col = Strand, id.col = Fill, 
                         plot.int = plot.int, l.height = l.height, 
                         des.seqs = seq.lens)
    
    # Create desired mutations exon plot
    ml <- mutExonPlot(mut.gr = mut.gr, exon.int = exon.int, 
                      plot.int = plot.int, disp.track = disp.track,
                      p.height = p.height, d.height = d.height, 
                      l.height = l.height, plt.title = plt.title, id.col = Fill, 
                      tx.db = tx.db, seq.lens = seq.lens, gr = gr, pfam.gr = pfam.gr, 
                      pfam.desc = pfam.desc, pfam.ids = pfam.ids, ...)
    return(ml)
}

# Create a plot using all inputs from the selection interfaces
mutExonPlot <- function(mut.gr, exon.int, plot.int = exon.int, disp.track = 2, 
                        p.height = 1/4, d.height = 0, l.height = 1/2, 
                        plt.title = "mutExonPlot", id.col = "Fill", 
                        tx.db = NULL, pfam.gr = NULL, pfam.desc = NULL, 
                        pfam.ids = NULL, gr = NULL, ...) {
    
    midpoint <- NULL
    Value <- NULL
    name <- NULL
    r <- NULL
    stepping <- NULL
    domain <- NULL
    
    show.domain = d.height > 0
    show.legend = l.height > 0
    
    # Using an aes function which evaluates arguments locally without parsing
    aesCst <- function(...) {
        structure(list(...),  class = "uneval")
    }
    
    # Determine all tracks of RNA
    if (is.null(gr)) {
        suppressWarnings(gr <- biovizBase::crunch(tx.db, which=exon.int, type = "all"))
    }
    
    # Determine relevant, reduced exons in plot interval
    found.ex <- reduce(gr[gr$type == "exon"])
    red.exon <- found.ex[found.ex %over% exon.int]
    if (show.legend) {
        genestrand <- as.character(strand(exon.int))
        # Reverse exon list if strand is negative, counting from
        # the correct direction
        if (genestrand == "-") {
            red.exon <- rev(red.exon)
        }
    }
    
    # Determine location of the gene
    exon.xlim <- c(start(ranges(plot.int)), end(ranges(plot.int)))
    
    if (disp.track == 1) {
        
        combo <- gr[gr$type %in% c("cds", "utr")]
        
        plot.ids <- sort(unique(values(combo)$tx_id[combo %over% plot.int]))
        combo <- combo[values(combo)$tx_id %in% plot.ids]
        
        # Stepping is based on tx_id
        values(combo)$stepping <-  match(values(combo)$tx_id, plot.ids)
        
        # Determine gaps between exon regions (utr and cds) for each track
        gr.cds <- combo[combo$type == "cds"]
        gr.utr <- combo[combo$type == "utr"]
        df.gaps <- biovizBase::getGaps(c(gr.cds, gr.utr), group.name = "tx_id")
        
    } else {
        
        # Stepping is the same accross reduced track
        gr.cds <- reduce(gr[gr$type == "cds"])
        values(gr.cds) <- list(type = "cds", stepping = 1)
        gr.utr <- setdiff(red.exon, gr.cds)
        values(gr.utr) <- list(type = "utr", stepping = 1)
        gr.blocks <- c(gr.cds, gr.utr)
        
        # Determine gaps between exon regions (utr and cds) for one track
        gr.rr <- reduce(ranges(gr.blocks))
        df.gaps <- gaps(gr.rr, start = min(start(gr.rr)), 
                        end = max(end(gr.rr)))
        chrs <- unique(as.character(seqnames(gr.blocks)))
        df.gaps <- GRanges(chrs, df.gaps)
        
    }
    
    if (show.domain) {
        
        # Obtain the domain transcripts within exon.int
        suppressWarnings(pfam.ovl <- subsetByOverlaps(pfam.gr, plot.int, 
                                                      type = "within"))
        
        if (length(pfam.ovl) > 0) {
            values(pfam.ovl) <- list(values(pfam.ovl), 
                                     stepping = 1:length(pfam.ovl))
            
            if (disp.track == 2) {
                # Display only longest domain in plot
                pfam.ovl <- pfam.ovl[which.max(width(pfam.ovl))]
            }
            
            # Rename pfam.ids to domain names
            pfam.ovl <- PFAMIDE(pfam.ovl, pfam.desc, pfam.ids)
            
        } else {
            show.domain = FALSE
        }
    }
    
    # Perform check to determine relevance of plot
    if (length(which(mut.gr %over% plot.int)) == 0) {
        checked <- FALSE
    } else {
        
        if (show.legend){
            # Match each mutation with a corresponding exon if applicable
            ex.col <- rep(NA, length = length(mut.gr))
            for (index in 1:length(mut.gr)) {
                a <- which(red.exon %over% mut.gr[index])
                if(length(a)) {
                    ex.col[index] <- a
                }
            }
            # Select exons which have mutations based on matched indices
            match.ind <- sort(unique(ex.col))
        } else {
            match.ind <- red.exon %over% mut.gr
        }
        match.exons <- red.exon[match.ind]
        
        # Create the exon annotations in the plot for the legend
        if (show.legend && length(which(!is.na(match.ind))) > 0){
            
            exon.df <- biovizBase::mold(match.exons)
            exon.df <- cbind(exon.df, Value = rep(toupper(letters), 
                                                  length = nrow(exon.df)), 
                             Exon = paste("Exon", match.ind))
            # Assign exon labels to mutations if applicable
            lab.col <- apply(matrix(ex.col), 1, function(row, ref.vec) {
                # Strip any leading spaces
                myval <- as.character(exon.df[match(paste("Exon", 
                                                          gsub("\\s", "", row),
                                                          sep = " "), 
                                                    ref.vec),"Value"])
                if (!(is.na(myval))) {
                    return(myval)
                } else {
                    return("NA")
                }
            }, as.character(exon.df[,"Exon"]))
            
            # Make legend data with exon annotation columns
            leg.data <- cbind(id.col = as.character(values(mut.gr)$value), 
                              Exon = ex.col, Exonlabel = lab.col)
            colnames(leg.data)[1] <- id.col
        } else {
            show.legend = FALSE
        }
        
        if (length(mut.gr) == 0) {
            checked <- FALSE
        } else {
            checked <- TRUE
        }
    } 
    
    if (checked){
        
        cds.df <- biovizBase::mold(gr.cds)
        utr.df <- biovizBase::mold(gr.utr)
        
        # Determine which exons should be colored in plot
        cds.color <- biovizBase::mold(gr.cds[gr.cds %over% match.exons])
        utr.color <- biovizBase::mold(gr.utr[gr.utr %over% match.exons])
        
        if (disp.track == 1) {
            cds.size = 0.4
            utr.size = 0.16
        } else {
            cds.size = 1.5
            utr.size = 0.6
        }
        
        if (show.domain | disp.track == 1) {
            # Determine rate of arrow indicating the strand
            gap.arrow.rate <- 0.03 * (end(range(ranges(plot.int))) - 
                                          start(range(ranges(plot.int))))/
                (end(range(ranges(df.gaps))) - start(range(ranges(df.gaps))))
        }
        
        p1 <- ggplot()
        if (disp.track == 1) {
            # Plot the arrows on the axis of the track(s) for a unique strand
            p1 <- p1 + ggbio::geom_arrow(data = df.gaps, stat = "identity", 
                                         aes(y=stepping), arrow.rate = gap.arrow.rate,
                                         angle = 50, color = "dark grey")+
                coord_cartesian(xlim = exon.xlim)
        } else {
            # Plot introns and gaps
            p1 <- p1 + geom_abline(intercept = 1, slope = 0, 
                                   aes(xmin = -1, xmax = 1, ymin = 0.5, 
                                       ymax = 1.5), 
                                   color = "dark grey")+
                coord_cartesian(xlim = exon.xlim, ylim = c(-0.6, 2.6))
        }
        # Plot the cds regions with appropriate size
        p1 <- plotExonRect(plt = p1, component = "cds", comp.df = cds.df, 
                           size = cds.size, subset.df = cds.color, 
                           aesCst = aesCst)
        # Plot the utr regions with appropriate size
        p1 <- plotExonRect(plt = p1, component = "utr", comp.df = utr.df, 
                           size = utr.size, subset.df = utr.color, 
                           aesCst = aesCst)
        # Impose x-axis scale and labels, limits
        if (exon.xlim[2] > 1e6) {
            xunits = "Mb"
        } else if (exon.xlim[2] > 1e3) {
            xunits = "kb"
        } else {
            xunits = "bp"
        }
        p1 <- p1 + theme(axis.text.y = element_blank(), 
                         axis.title.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.text.x = element_text(size = 8)) + 
            ggbio::scale_x_sequnit(xunits)
        
        if (show.domain) {
            p1 <- p1 + 
                theme(plot.margin=unit(c(0.5,1,0.5,1), "cm"))
        } else {
            p1 <- p1 + 
                theme(plot.margin=unit(c(0.5,1,0.25,1), "cm"))
        }
        
        p2 <- ggplot()
        # Plots karyogram of mutations and exon annotations
        if (disp.track == 1) {
            p2 <- p2 + ggbio::layout_karyogram(mut.gr, color = "red", fill = "red", 
                                               geom = "rect", ylab = NULL, 
                                               xlab = NULL) +
                coord_cartesian(xlim = exon.xlim)
        } else {
            # Vary height of mutation bars to represent a tally
            height.col = countOverlaps(mut.gr)
            max.height = max(height.col)
            values(mut.gr) <- list(values(mut.gr), height = height.col)
            mut.df <- biovizBase::mold(mut.gr)
            p2 <- p2 + 
                ggplot2::geom_rect(data = mut.df, 
                                   aesCst(xmin = mut.df[,"start"], 
                                          xmax = mut.df[,"end"], 
                                          ymin = -1.5/max.height*
                                              mut.df[,"height"],
                                          ymax = 1.5/max.height*
                                              mut.df[,"height"]), 
                                   fill = "red", colour = "red")+
                coord_cartesian(xlim = exon.xlim, ylim=c(-1.6, 1.6))
        }
        
        # Add labels to Exons to be referenced in plot legend
        if (show.legend){
            p2 <- p2 + geom_text(data = exon.df, 
                                 aes(x = midpoint, y = 0, label = Value), 
                                 vjust = 0, hjust = 0.4, colour = "black", 
                                 fontface = "bold", size = 6)
        }
        
        # Impose x-axis scale and labels, limits, theme options
        p2 <- p2 + theme(panel.background = element_rect(fill = NA), 
                         axis.text.y = element_blank(), 
                         axis.title.y = element_blank(), 
                         axis.ticks.y = element_blank(), 
                         axis.text.x = element_text(size = 8), 
                         plot.margin=unit(c(0.5,1,0.5,1), "cm"), 
                         panel.grid = element_blank())
        p2 <- p2 + ggbio::scale_x_sequnit(xunits)
        
        if(show.domain & disp.track == 2) {
            max.df <- biovizBase::mold(pfam.ovl)
            domainy = -1.7
            domain.arr = -1.8
            p2ylim <- c(-1.9, 1.6)
            p1ylim <- c(-0.9, 2.6)
            if(show.legend) {
                domain.arr = domain.arr - 0.5
                domainy = domainy - 0.2
                p2ylim[1] = p2ylim[1] - 0.5
                p1ylim[1] = p1ylim[1] - 0.5
            }
            # Add the longest domain track to the reduced track plot
            p2 <- p2 + 
                geom_text(data = max.df, 
                          aesCst(x = max.df[,"midpoint"], 
                                 y = domainy, label = max.df[,"domain"]), 
                          colour = "black", fontface = "bold", size = 4)
            p2 <- p2 + ggbio::geom_arrow(data = pfam.ovl, stat = "identity", 
                                         aesCst(y = domain.arr), angle = 50, 
                                         color = "black")
            p2 <- p2 + coord_cartesian(xlim = exon.xlim, ylim = p2ylim)
            p1 <- p1 + coord_cartesian(xlim = exon.xlim, ylim = p1ylim)
        }
        
        # Draw the two plots on top of each other
        # Supplied by http://rpubs.com/kohske/dual_axis_in_ggplot2
        g1 <- ggplotGrob(p1+theme(legend.position = "none"))
        g2 <- ggplotGrob(p2+theme(legend.position = "none"))
        pp <- c(subset(g1$layout, name == "panel", se = t:r))
        me <- gtable::gtable_add_grob(g1, 
                                      g2$grobs[[which(g2$layout$name == 
                                                          "panel")]], 
                                      pp$t, pp$l, pp$b, pp$l)
        
        # Plot all domains in a seperate plot
        if (show.domain && disp.track == 1){
            
            domains <- ggplot()
            # Plot exons in background
            red.exon.df <- biovizBase::mold(red.exon)
            domains <- domains + 
                ggplot2::geom_rect(data = red.exon.df, 
                                   aesCst(xmin = red.exon.df[,"start"], 
                                          xmax = red.exon.df[,"end"], 
                                          ymin = 1, ymax = length(pfam.ovl)), 
                                   colour = "grey", fill = "grey") +
                coord_cartesian(xlim = exon.xlim)
            domains <- domains + 
                geom_text(data = biovizBase::mold(pfam.ovl), 
                          aes(x = midpoint, y = stepping+0.5, label = domain), 
                          colour = "black", fontface = "bold", size = 4)
            # Plot domain transcripts and names
            domains <- domains + 
                ggbio::geom_arrow(data = pfam.ovl, stat = "stepping", 
                                  angle = 50, 
                                  color = "black") 
            domains <- domains + 
                theme(axis.text.y = element_blank(), 
                      axis.title.y = element_blank(), 
                      axis.ticks.y = element_blank(), 
                      axis.text.x = element_blank(), 
                      axis.ticks.x = element_blank(),
                      plot.margin=unit(c(1,1,1,1), "cm"))
            # Make domains grob object
            doms <- ggplotGrob(domains + 
                                   theme(legend.position = "none"))
            
            # specify width in output plot
            maxWidth = unit.pmax(me$widths[2:5], doms$widths[2:5])
            me$widths[2:5] <- as.list(maxWidth)
            doms$widths[2:5] <- as.list(maxWidth)
        }
        
        if (show.legend){
            # Creates a fixed-scale legend grob from the mutation data
            grob1 = gridExtra::tableGrob(leg.data[,c(id.col, "Exon", 
                                                     "Exonlabel")], 
                                         gpar.coretext = gpar(fontsize = 10), 
                                         gpar.coltext  = gpar(fontsize = 12),            
                                         show.rownames = FALSE, 
                                         equal.height = FALSE, 
                                         padding.h=unit(2, "mm"),
                                         padding.v=unit(2, "mm"), 
                                         h.even.alpha=0.15, 
                                         h.odd.alpha=0.35,  
                                         v.even.alpha=0.5, 
                                         v.odd.alpha=0.5,
                                         show.hlines = FALSE,
                                         show.vlines = FALSE,
                                         show.box = TRUE,
                                         separator = "black",
                                         gpar.corefill = gpar(fill = "grey80", 
                                                              col = "black"),
                                         gpar.colfill = gpar(fill = "white", 
                                                             col = "black"))
        }
        
        # Combine all plot components
        if (show.domain & disp.track == 1) {
            if (show.legend) {
                mep <- gridExtra::arrangeGrob(me, doms, grob1, ncol=1, 
                                              heights=c(p.height, d.height, 
                                                        l.height),
                                              # added main title
                                              main = textGrob(plt.title, 
                                                              vjust = 0.9, 
                                                              gp = gpar(fontface = "bold", 
                                                                        cex = 1.5)))
                returns <- list(mep = mep, plot = me, domain = doms, 
                                legend = grob1, show.domain = show.domain, 
                                show.legend = show.legend, gr = gr, 
                                exon.int = exon.int)
            } else {
                mep <- gridExtra::arrangeGrob(me, doms, ncol=1, 
                                              heights=c(p.height, d.height),
                                              # added main title
                                              main = textGrob(plt.title, 
                                                              vjust = 0.9, 
                                                              gp = gpar(fontface = "bold", 
                                                                        cex = 1.5)))
                returns <- list(mep = mep, plot = me, domain = doms, 
                                show.domain = show.domain, 
                                show.legend = show.legend, gr = gr, 
                                exon.int = exon.int)
            }
        } else if (show.legend) {
            mep <- gridExtra::arrangeGrob(me, grob1, ncol=1, heights=c(p.height, l.height),
                                          # added main title
                                          main = textGrob(plt.title, 
                                                          vjust = 0.9, 
                                                          gp = gpar(fontface = "bold", 
                                                                    cex = 1.5)))
            returns <- list(mep = mep, plot = me, legend = grob1, 
                            show.domain = show.domain, 
                            show.legend = show.legend, gr = gr, 
                            exon.int = exon.int)
        } else {
            mep <- gridExtra::arrangeGrob(me, main = textGrob(plt.title, 
                                                              vjust = 0.9, 
                                                              gp = gpar(fontface = "bold",
                                                                        cex = 1.5)))
            returns <- list(mep = mep, plot = me, show.domain = show.domain, 
                            show.legend = show.legend, gr = gr, 
                            exon.int = exon.int)
        }
        
    } else {
        returns <- checked
    }
    return(returns)
}
