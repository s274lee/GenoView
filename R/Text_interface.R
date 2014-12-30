# Mutations Exon plot wrapper for one gene using text interface
mepTxtInt <- function(m.data, tx.db, gene.loc, seq.lens, pfam.gr = NULL, pfam.desc = NULL, 
                      pfam.ids = NULL, ...){
    
    # Test if using hg19 data
    meta <- metadata(tx.db)
    genome <- meta[meta[,"name"] == "Genome", "value"]
    prot.dom <- toupper(genome) == "HG19"
    
    # Create options for text interface
    int.opt <- c("Whole Gene", "Custom Interval")
    disp.opt <- c("Display all RNA tracks", "Display condensed track")
    
    # Error text if plotting components fail
    err.txt <- "No overlap between exon and mutation locations. 
    No plot created. "
    
    # Declare boolean which signals that plotting is complete
    done <- FALSE
    # Declare initial null gr
    gr <- NULL
    
    # Create temporary data.frame from m.data, mainly for column selection
    temp.mut <- checkCols(m.data)
    
    while (!done) {
        
        # Assign gene.name of the plot
        cat("Please enter the gene name: \n")
        if (exists("gene.name")) {
            new.name <- toupper(readline())
            if (!(new.name == gene.name)) {
                gene.name <- toupper(readline())
                exon.int <- gene.loc[gene.name]
                gr <- NULL
            }
        } else {
            gene.name <- toupper(readline())
            exon.int <- gene.loc[gene.name]
        }
        
        # Assign the interval option of the plot
        exon.method <- printOption(int.opt, desc = "the interval", 
                                   match.ind = TRUE)
        
        # Assign custom start and end locations if necessary
        if (exon.method == 2) {
            plot.int <- intSel(gr = exon.int)
        } else {
            plot.int = exon.int
        }
        
        # Assigning necessary plotting criteria
        Chromosome <- colSel("Chromosome", temp.mut)
        Start.position <- colSel("Start Position", temp.mut)
        End.position <- colSel("End Position", temp.mut)
        Strand <- colSel("Strand", temp.mut)
        
        # Selection of exon plot type
        disp.track <- printOption(disp.opt, match.ind = TRUE)
        
        # Assign the title and height allocations for the plot
        cat("Please enter a title for the plot: \n")
        plt.title <- readline()
        cat("Please enter a height allocation for each component", 
            "(0 to hide component). \n")
        cat("Please enter a height for the plot (recommended: 1/4): \n")
        p.height <- eval(parse(text = readline()))
        
        if (prot.dom) {
            cat("Please enter a height for the domain (recommended: 1/4): \n")
            d.height <- eval(parse(text = readline()))
        } else {
            d.height = 0
        }
        
        cat("Please enter a height for the legend (recommended: 1/2): \n")
        l.height <- eval(parse(text = readline()))
        if (l.height > 0) {
            Fill <- colSel("Fill", temp.mut)
        } else {
            Fill = NULL
        }
        
        # Create desired mutations exon plot
        ml <- checkAndPlot(m.data=m.data, Chromosome = Chromosome, 
                           Start.position = Start.position, 
                           End.position = End.position, Strand = Strand, 
                           Fill = Fill, plot.int = plot.int, 
                           l.height = l.height, seq.lens = seq.lens, 
                           exon.int = exon.int, disp.track = disp.track, 
                           p.height = p.height, d.height = d.height, 
                           plt.title = plt.title, tx.db = tx.db, gr = gr, 
                           pfam.gr, pfam.desc, pfam.ids, ...)
        
        if (class(ml) == "list") {
            
            # Check domain and legend plotting
            intPlotFailed(ml = ml, d.height = d.height, 
                          l.height = l.height, gui = FALSE)
            
            mep <- ml$mep
            gr <- ml$gr
            oldex.gr <- ml$exon.int
            
            # Display plot
            grid.newpage()
            grid.draw(mep)
            
            # Saving plot options
            cat("Would you like to save your plot? Y/N \n")
            save.plt <- toupper(readline())
            if (save.plt == "Y") {
                cat("Please enter a filename with suffix for the plot. \n")
                cat("Please note: \".jpg\" is not supported but \".jpeg\" is. \n")
                save.loc <- readline()
                save.loc <- gsub("/", "\\", save.loc, fixed=TRUE)
                graphics.off()
                
                # Alter saving function to allow grobs
                body(ggsave) <- body(ggsave)[-2]
                ggsave(filename = save.loc, plot = mep, scale = 1)
                grid.draw(mep)
            }
        } else {
            plotFailed(txt = err.txt, gui = FALSE)
        }
        cat("Are you finished plotting? Y/N \n")
        done.bool <- toupper(readline())
        if (done.bool == "Y"){
            done <- TRUE
        }
    }
    graphics.off()
}

# Prints a numbered list of options
printOption <- function(opt.vec, inst = TRUE, desc = "an option", 
                        match.opt = FALSE, match.ind = FALSE) {
    
    # If instructions should be printed
    if (inst) {
        cat(sprintf("Please pick %s from the following: \n", desc));
    }
    
    holder <- apply(matrix(opt.vec), 1, function(x) {
        cat(match(x, opt.vec), ". ", x, '\n', sep = "")
    })
    
    # Determines returned option
    while (match.opt) {
        input <- scan(n=1)
        tryCatch({
            # Returns a valid selected option
            sel <- opt.vec[input]
            match.opt = FALSE
            return(sel)
        }, error = function(e) print("Input is not valid"))
    }
    
    # Determines returned input
    while (match.ind) {
        input <- scan(n=1)
        if (input <= length(opt.vec)) {
            match.ind = FALSE
            return(input)
        } else{
            print("Input is not valid")
        }
    }
}

# Chooses column name for criteria based on user selection
colSel <- function(col.crit, dataFrame) {
    name <- printOption(colnames(dataFrame), 
                        desc = sprintf("the %s criteria", col.crit), 
                        match.opt = TRUE)
    return(name)
}
