# Mutations Exon plot wrapper for one gene using GUI
mepGUI <- function(m.data, tx.db, gene.loc, seq.lens, pfam.gr = NULL, 
                   pfam.desc = NULL, pfam.ids = NULL, ...){
    
    # Declare GUI function's environment
    mep.env <- environment()
    
    pass <- FALSE
    gui.lst <- NULL
    exon.xlim <- NULL
    gene.name <- NULL
    exon.int <- NULL
    disp.track <- NULL
    exon.method <- NULL
    
    # Test if using hg19 data
    meta <- metadata(tx.db)
    genome <- meta[meta[,"name"] == "Genome", "value"]
    prot.dom <- toupper(genome) == "HG19"
    
    # Text for plot description notebook page
    plt.txt <- "The taller blocks represent coding sequences (cds). 
    The shorter blocks represent untranslated regions (utr). 
    The lines in between blocks represent introns. \n
    The exons which are blue and labelled overlap one or more input mutations.
    The exons which are grey do not overlap any input mutations. \n
    The red lines denote the intervals between start and end columns. \n
    The legend is not scalable. 
    Please resize the table to make the necessary accommodations. \n
    Right-click the plot to save. 
    Please include a suffix in your filename (e.g. \".png\", \".jpeg\").
    \".jpg\" is not supported, \".jpeg\" is. \n
    The plot needs to be uncovered to export properly, 
    verify that the export is successful before closing the tab."
    
    # Error text if plotting components fail
    err.txt = "No overlap between exon and mutation locations. 
    No plot created. "
    
    # Create options for gui domain selection
    int.opt <- c("Whole Gene", "Custom Interval")
    disp.opt <- c("Display all RNA tracks", "Display condensed track")
    dom.opt <- c("Plots all domains in a seperate section", 
                 "Longest domain included in plot")
    inv.name <- "Gene name invalid. Interval selection failed. "
    inv.int <- "Interval selection invalid. "
    
    # Create temporary data.frame from m.data, mainly for column selection
    temp.mut <- checkCols(m.data)
    
    # Incremental suffix for different plots
    name.suff <- 1
    # Declare initial null gr and oldex.gr, which are assigned after plotting
    gr <- NULL
    oldex.gr <- NULL
    
    # List of widgets used in gformlayout
    gfl.input <- list(type = "ggroup", 
                      horizontal = FALSE, 
                      children = list(
                          list(type = "fieldset",
                               columns = 1, 
                               label = "Column selection", 
                               label.font = c(weight = "bold"), 
                               
                               # Creates comboboxes for each of the required 
                               # four input criteria
                               children = list(
                                   list(name = "chrom.box", 
                                        label = "Chromosome column:", 
                                        type = "gcombobox", 
                                        items = colnames(temp.mut)
                                   ), 
                                   list(name = "start.box", 
                                        label = "Start column:", 
                                        type = "gcombobox", 
                                        items = colnames(temp.mut)
                                   ), 
                                   list(name = "end.box", 
                                        label = "End column:", 
                                        type = "gcombobox", 
                                        items = colnames(temp.mut)
                                   ), 
                                   list(name = "strand.box", 
                                        label = "Strand column:", 
                                        type = "gcombobox", 
                                        items = colnames(temp.mut)
                                   )
                               )
                          )
                      )
    )
    
    # Creating GUI interface
    # Creates main window to house all components
    w <- gWidgets::gwindow("Mutations Over Exons Plot")
    
    # Creates group which holds plotting settings
    g <- gWidgets::ggroup(horizontal = FALSE, container = w)
    g.lay <- gWidgets::glayout(horizontal = FALSE, container = g)
    
    # Creates notebook to display temp.mut data.frame and plotted results
    my.nb <- gWidgets::gnotebook(container = w, expand = TRUE, closebuttons = TRUE, 
                                 dontCloseThese = c(1, 2))
    gWidgets::size(my.nb) <- c(250, 250)
    g.data <- gWidgets::gtable(temp.mut, expand = TRUE, container = my.nb,
                               label="Datatable")
    g.plt.txt <- gWidgets::glabel(plt.txt, expand = TRUE, container = my.nb, 
                                  label="Plot Description")
    gWidgets::svalue(my.nb) <- 1
    
    # Creates the status bar which displays progress
    stat.bar <-gWidgets::gstatusbar(container = w)
    
    # Creates frame which contains the current gene name
    geneF <- gWidgets::gframe("Gene")
    gene.lab <- gWidgets::glabel("Gene: ", container = geneF)
    gWidgets::font(gene.lab) <- list(weight = "bold")
    gene.box <- gWidgets::gedit("", container = geneF)
    g.lay[1,1] <- geneF
    
    # Updates the interval selection data based on gene intervals
    b1 <- gWidgets::gbutton("Update Gene", container = geneF, 
                            handler = function(h, ...){
                                
                                # Check that provided gene name is valid
                                tryCatch({assign("gene.name", gWidgets::svalue(gene.box), 
                                                 envir = mep.env)
                                          assign("exon.int", gene.loc[gene.name], 
                                                 envir = mep.env)
                                          gWidgets::svalue(stat.bar) <- ""
                                          assign("exon.xlim", ranges(exon.int), 
                                                 envir = mep.env)
                                          assign("pass", TRUE, mep.env)
                                }, error = function(e) {
                                    gWidgets::svalue(stat.bar) <- inv.name
                                    assign("pass", FALSE, mep.env)})
                                
                                if (pass){
                                    
                                    # Assign start and end values of gene 
                                    # to interval selection
                                    # Must be updated because gene name is 
                                    # subject to change
                                    gWidgets::svalue(start.lab) <- paste("Gene start:", 
                                                                         start(exon.xlim))
                                    gWidgets::svalue(end.lab) <- paste("Gene end:", 
                                                                       end(exon.xlim))
                                    gWidgets::svalue(cust.start) <- start(exon.xlim)
                                    gWidgets::svalue(cust.end) <- end(exon.xlim)
                                    
                                    handInt(int.sel, int.opt, cust.start, cust.end, 
                                            stat.bar, inv.int)
                                }
                            })
    
    # Create the droplist which controls interval selection availability
    intOptF <- gWidgets::gframe("Interval Selection")
    int.sel.lab <- gWidgets::glabel("Interval selection: ", container = intOptF)
    gWidgets::font(int.sel.lab) <- list(weight = "bold")
    int.sel <- gWidgets::gcombobox(c("", int.opt), selected = 1, 
                                   container = intOptF)
    gWidgets::addHandlerClicked(int.sel, function(h, ...){
        try(handInt(int.sel, int.opt, cust.start, cust.end, stat.bar, inv.int))})
    g.lay[2,1] <- intOptF
    
    # Creates all interval selection menu components
    intF <- gWidgets::gframe("Interval", horizontal = FALSE)
    intT <- gWidgets::ggroup(container = intF)
    intB <- gWidgets::ggroup(container = intF)
    start.lab <- gWidgets::glabel("Gene start:", container = intT)
    gWidgets::font(start.lab) <- list(weight = "bold")
    end.lab <- gWidgets::glabel("Gene end:", container = intT)
    gWidgets::font(end.lab) <- list(weight = "bold")
    cust.start.lab <- gWidgets::glabel("Custom \nstart:", container = intB)
    gWidgets::font(cust.start.lab) <- list(weight = "bold")
    cust.start <- gWidgets::gedit(container = intB)
    cust.end.lab <- gWidgets::glabel("Custom \nend:", container = intB)
    gWidgets::font(cust.end.lab) <- list(weight = "bold")
    cust.end <- gWidgets::gedit(container = intB)
    g.lay[3,1] <- intF
    
    # Display the gformlayout widgets
    fl2 <- gWidgets::gformlayout(gfl.input, container = g, expand = FALSE)
    
    # Track selection
    traF <- gWidgets::gframe("Tracks", container = g)
    tra.lab <- gWidgets::glabel("Tracks: ", container = traF)
    gWidgets::font(tra.lab) <- list(weight = "bold")
    track.sel <- gWidgets::gradio(items = disp.opt, container = traF) 
    
    # Plot title and height selection
    ploF <- gWidgets::gframe("Plot", container = g, horizontal = FALSE)
    ploT <- gWidgets::ggroup(container = ploF)
    ploB <- gWidgets::ggroup(container = ploF)
    plot.lab <- gWidgets::glabel("Plot Title: ", container = ploT)
    gWidgets::font(plot.lab) <- list(weight = "bold")
    plot.title <- gWidgets::gedit("", container = ploT)
    p.height.lab <- gWidgets::glabel("Plot Height: ", container = ploB)
    gWidgets::font(p.height.lab) <- list(weight = "bold")
    plot.height <- gWidgets::gedit("1/4", container = ploB)
    
    if (prot.dom) {
        
        # Selection of domain inputs
        domF <- gWidgets::gframe("Domain", container = g, horizontal = FALSE)
        domT <- gWidgets::ggroup(container = domF)
        domM <- gWidgets::ggroup(container = domF)
        domB <- gWidgets::ggroup(container = domF)
        dom.lab <- gWidgets::glabel("Show Domain", container = domT)
        gWidgets::font(dom.lab) <- list(weight = "bold")
        show.domain <- gWidgets::gcheckbox(checked = FALSE, container = domT, 
                                           handler = function(h, ...) {
                                               handDom(show.domain, 
                                                       domain.height, track.sel, 
                                                       disp.opt, dom.track.lab, 
                                                       dom.opt)
                                           })
        dom.track.lab <- gWidgets::glabel("", container = domM)
        d.height.lab <- gWidgets::glabel("Domain Height: ", container = domB)
        gWidgets::font(d.height.lab) <- list(weight = "bold")
        domain.height <- gWidgets::gedit("1/4", container = domB)
        gWidgets::enabled(domain.height) <- FALSE
        gWidgets::addHandlerClicked(track.sel, function(h, ...) {
            handDom(show.domain, domain.height, track.sel, disp.opt, 
                    dom.track.lab, dom.opt)
        })
    }
    
    # Selection of legend inputs
    legF <- gWidgets::gframe("Legend", container = g, horizontal = FALSE)
    legT <- gWidgets::ggroup(container = legF)
    legM <- gWidgets::ggroup(container = legF)
    legB <- gWidgets::ggroup(container = legF)
    leg.lab <- gWidgets::glabel("Show Legend", container = legT)
    gWidgets::font(leg.lab) <- list(weight = "bold")
    show.legend <- gWidgets::gcheckbox(checked = FALSE, container = legT, 
                                       handler = function(h, ...) {
                                           if (gWidgets::svalue(show.legend)) {
                                               gWidgets::enabled(legend.height) <- TRUE
                                               gWidgets::enabled(fill.box) <- TRUE
                                           } else {
                                               gWidgets::enabled(legend.height) <- FALSE
                                               gWidgets::enabled(fill.box) <- FALSE
                                           }
                                       })
    
    l.heightlab <- gWidgets::glabel("Legend Height: ", container = legM)
    gWidgets::font(l.heightlab) <- list(weight = "bold")
    legend.height <- gWidgets::gedit("1/2", container = legM)
    gWidgets::enabled(legend.height) <- FALSE
    fill.lab <- gWidgets::glabel("Fill column:", container = legB)
    gWidgets::font(fill.lab) <- list(weight = "bold")
    fill.box <- gWidgets::gcombobox(items = colnames(temp.mut), container = legB)
    gWidgets::enabled(fill.box) <- FALSE
    
    b2 <- gWidgets::gbutton(text = "Plot", container = g, expand = FALSE,
                            handler = function(h, ...) {
                                
                                environment(mep.env)
                                
                                # Creates objects containing all plot settings after
                                # verifying their presence
                                
                                if (exists("pass")){
                                    if (pass) {
                                        if (gWidgets::svalue(int.sel) != "") {
                                            if (length(grep("Plotting", gWidgets::svalue(stat.bar), 
                                                            fixed = TRUE)) == 0) {
                                                
                                                gui.lst = c(exon.int = exon.int, 
                                                            gene.name = gene.name, 
                                                            int.sel = gWidgets::svalue(int.sel), 
                                                            gWidgets::svalue(fl2), 
                                                            start.loc = gWidgets::svalue(cust.start), 
                                                            end.loc = gWidgets::svalue(cust.end), 
                                                            track.sel = gWidgets::svalue(track.sel), 
                                                            plot.title = gWidgets::svalue(plot.title), 
                                                            plot.height = gWidgets::svalue(plot.height), 
                                                            show.legend = gWidgets::svalue(show.legend), 
                                                            legend.height = gWidgets::svalue(legend.height), 
                                                            fill.box = gWidgets::svalue(fill.box))
                                                
                                                # Add domain data if applicable
                                                if (prot.dom) {
                                                    gui.lst = c(gui.lst, 
                                                                show.domain = gWidgets::svalue(show.domain), 
                                                                domain.height = gWidgets::svalue(domain.height))
                                                }
                                                
                                                assign("gui.lst", gui.lst, envir = mep.env)
                                                gWidgets::svalue(stat.bar) <- paste0("Plotting (", name.suff ,")...")
                                                
                                                # Assign the interval option of the plot
                                                exon.method <- match(gui.lst$int.sel, int.opt)
                                                # Determine if we need to calculate gr through gene.name
                                                if (!(identical(oldex.gr, gui.lst$exon.int))) {
                                                    gene.name <- gui.lst$gene.name
                                                    exon.int <- gui.lst$exon.int
                                                    gr <- NULL
                                                }
                                                
                                                # Assign custom start and end locations if necessary
                                                if (exon.method == 2) {
                                                    start.loc <- eval(parse(text = gui.lst$start.loc))
                                                    end.loc <- eval(parse(text = gui.lst$end.loc))
                                                    plot.int <- intSel(exon.int, start.loc, end.loc)
                                                } else {
                                                    plot.int = exon.int
                                                }
                                                
                                                # Assigning necessary plotting criteria
                                                Chromosome <- gui.lst$chrom.box
                                                Start.position <- gui.lst$start.box
                                                End.position <- gui.lst$end.box
                                                Strand <- gui.lst$strand.box
                                                
                                                # Selection of exon plot type
                                                disp.track <- match(gui.lst$track.sel, disp.opt)
                                                
                                                # Assign the title and height allocation for the plot
                                                plt.title <- gui.lst$plot.title
                                                p.height <- eval(parse(text = gui.lst$plot.height))
                                                
                                                # Assign domain height allocation if applicable
                                                if (prot.dom) {
                                                    if (gui.lst$show.domain) {
                                                        d.height <- eval(parse(text = gui.lst$domain.height))
                                                        if(disp.track == 2) {
                                                            d.height = 1/4
                                                        }
                                                    } else {
                                                        d.height <- 0
                                                    }
                                                } else {
                                                    d.height <- 0
                                                }
                                                
                                                # Assign legend height allocation if applicable
                                                if (gui.lst$show.legend) {
                                                    Fill <- gui.lst$fill.box
                                                    l.height <- eval(parse(text = gui.lst$legend.height))
                                                } else {
                                                    Fill = NULL
                                                    l.height <- 0
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
                                                
                                                gWidgets::svalue(stat.bar) <- paste0("Finished plotting (", name.suff, "). ")
                                                
                                                if (class(ml) == "list") {
                                                    
                                                    # Check domain and legend plotting
                                                    intPlotFailed(ml = ml, d.height = d.height, 
                                                                  l.height = l.height, gui = TRUE, stat.bar = stat.bar)
                                                    
                                                    assign("gr", ml$gr, envir = mep.env)
                                                    assign("oldex.gr", ml$exon.int, envir = mep.env)
                                                    
                                                    # Display the plot and provide saving options
                                                    holder = notebookVis(nb = my.nb, 
                                                                         plt = ml$mep, 
                                                                         type = "Mutations over Exons", 
                                                                         suff = name.suff)
                                                    if (ml$show.legend) {
                                                        holder2 = notebookVis(nb = my.nb, 
                                                                              plt = ml$legend, 
                                                                              type = "Legend", 
                                                                              suff = name.suff)
                                                    }
                                                    
                                                    # Increase name.suff for the next round of plotting
                                                    assign("name.suff", name.suff + 1, envir = mep.env)
                                                } else {
                                                    plotFailed(txt = err.txt, gui = TRUE, stat.bar = stat.bar)
                                                }
                                                
                                            } else {
                                                gWidgets::gmessage(message = "Plotting in progress", 
                                                                   title = "Please wait...",
                                                                   icon = "error")
                                            }
                                        } else {
                                            gWidgets::svalue(stat.bar) <- inv.int
                                        }
                                    } else {
                                        gWidgets::svalue(stat.bar) <- inv.name
                                    } 
                                } else {
                                    gWidgets::svalue(stat.bar) <- inv.name
                                }
                            })
    gWidgets::addSpring(g)
    b3 <- gWidgets::gbutton(text = "Done plotting", container = g, expand = FALSE, 
                            handler = function(h, ...) {
                                gWidgets::dispose(w)
                                graphics.off()
                            })
    # Additional handler to stop plotting once window is closed
    gWidgets::addHandlerUnrealize(w, handler = function(h, ...) {
        conf <- gWidgets::gconfirm("Exit/stop plotting?", title = "Confirm action") 
        if (conf) {
            graphics.off()
            return(FALSE)
        } else {
            return(TRUE)
        }
    })
}

# GUI handler to adjust interval selection components' availability/blur
handInt <- function(int.sel, int.opt, cust.start, cust.end, stat.bar, 
                    inv.int){
    
    if (match(gWidgets::svalue(int.sel), c("", int.opt))-1 == 1) {
        gWidgets::enabled(cust.start) <- FALSE
        gWidgets::enabled(cust.end) <- FALSE
        gWidgets::svalue(stat.bar) <- ""
    } else if (match(gWidgets::svalue(int.sel), c("", int.opt))-1 == 2) {
        gWidgets::enabled(cust.start) <- TRUE
        gWidgets::enabled(cust.end) <- TRUE
        gWidgets::svalue(stat.bar) <- ""
    } else {
        gWidgets::svalue(stat.bar) <- inv.int
        gWidgets::enabled(cust.start) <- FALSE
        gWidgets::enabled(cust.end) <- FALSE
    }
}

# GUI handler to adjust domain selection components' availability/blur
handDom <- function(show.domain, domain.height, track.sel, disp.opt, 
                    dom.track.lab, dom.opt){
    
    if (gWidgets::svalue(show.domain)) {
        gWidgets::enabled(domain.height) <- TRUE
        
        # Description changes with the selected track option
        if (gWidgets::svalue(track.sel) == disp.opt[1]) {
            gWidgets::svalue(dom.track.lab) <- dom.opt[1]
        } else {
            gWidgets::svalue(dom.track.lab) <- dom.opt[2]
            gWidgets::enabled(domain.height) <- FALSE
        }
    } else {
        gWidgets::enabled(domain.height) <- FALSE
        gWidgets::svalue(dom.track.lab) <- ""
    }
}

# Visualize a grob plt in a gnotebook object
notebookVis <- function(nb, plt, type, suff) {
    
    plt.lab <- paste(type, suff)
    # Create a new ggraphics gwidget
    myplot <- gWidgets::ggraphics(container = nb, label = plt.lab, 
                                  expand = TRUE, visible = TRUE)
    grid.draw(plt)
    return(list(plt, plt.lab))
}
