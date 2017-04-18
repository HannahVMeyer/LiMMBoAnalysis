# adapted from Turner, S.D. qqman: an R package for visualizing GWAS results 
# using Q-Q and manhattan plots. biorXiv DOI: 10.1101/005165 (2014).

manhattan <- function(dataframe, title=NULL, max.y="max", min.y="min", 
                     suggestiveline=0, genomewideline=-log10(5e-8), 
                     size.x.labels=6, size.y.labels=10, annotate=FALSE, 
                     mtvsst=FALSE, trial=FALSE, xscale=FALSE, SNPlist=NULL, 
                     cols=c("#67a9cf", "#016c59"), a=0.5, 
                     colGenomewideline='red', colSuggestiveline='blue') {
    
    if (annotate & is.null(SNPlist)) {
        stop("You requested annotation but provided no SNPlist!")
    }
    
    d <- dataframe
    names(d) <- toupper(names(d))
    
    d <- d[d$CHR %in% 1:22, ]
    
    if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
        
        d <- na.omit(d)
        d <- d[d$P>0 & d$P<=1, ]
        d$logp <- -log10(d$P)
        
        d$pos <- NA
        ticks <- NULL
        lastbase <- 0
        
        d <- d[order(d$CHR),]
        
        numchroms <- length(unique(d$CHR))
        if (numchroms == 1) {
            d$pos <- d$BP
        } else {
            for (i in unique(d$CHR)) {
                if (i == 1) {
                    d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP
                }	else {
                    lastbase <- lastbase + max(subset(d,CHR==i-1)$BP)
                    d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP + lastbase
                }
                ticks <- c(ticks, 
                        d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
            }
            ticklim <- c(min(d$pos),max(d$pos))
            
        }
        
        mycols <- rep(c("gray10","gray60"), max(d$CHR))
        
        if (max.y=="max") {
            maxy <- ceiling(max(d$logp)) 
        } else {
            maxy <- max.y 
        } 
        if (min.y=="min") {
            miny <- floor(min(d$logp))
        } else {
            miny <-min.y 
        } 
        if (maxy < 8) {
            maxy <- 8
        }
        if (annotate) {
            d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
        }
        
        if (numchroms == 1) {
            p <- ggplot(data=d, aes(x=pos, y=logp))
            p <- p + geom_point()
            p <- p + ylab(expression(-log[10](italic(p)))) + 
                xlab(paste("Chromosome",unique(d$CHR),"position"))
        }	else {
            p <- ggplot(data=d, aes(x=pos, y=logp))
            p <- p + ylab(expression(-log[10](italic(p))))
            if (!xscale) p  <- p + 
                scale_x_continuous(name="Chromosome", breaks=ticks, 
                                   limits=ticklim, labels=(unique(d$CHR)))
            if (xscale)  p  <- p + 
                scale_x_continuous(name="Chromosome", breaks=ticks, 
                                   expand=c(0,0), limits=ticklim, 
                                   labels=(unique(d$CHR)))
            p <- p + scale_y_continuous(limits=c(miny,maxy))
        }
        
        if (trial) {
            p <- p + 
                geom_point(aes(color=factor(SETUP), shape=factor(ANALYSIS)), 
                           alpha=a)
            p <- p + facet_grid(THRESHOLD~TYPE)
        } else if (mtvsst) {
            p <- p + geom_point(aes(color=TYPE, shape=MARKER), alpha = a)
            p <- p + scale_colour_manual(values=cols)            
            p <- p + scale_shape_manual(values=c(20,15), guide=FALSE)
        } else {
            p <- p + geom_point(aes(color=as.factor(CHR)))
            p <- p + scale_colour_manual(values=mycols)
            p <- p + theme(legend.position = "none") 
        }
        if (annotate) 	p <- p + geom_point(data=d.annotate, 
                                              colour=I("green3")) 
        
        
        p <- p + theme(title=title)
        p <- p + theme_bw()
        p <- p + theme(
            axis.text.x=element_text(size=size.x.labels, colour="grey50"), 
            axis.text.y=element_text(size=size.y.labels, colour="grey50"), 
            axis.ticks=element_blank()
        )
        
        if (suggestiveline) {
            p <- p + 
                geom_hline(yintercept=suggestiveline, 
                           colour=colSuggestiveline, alpha=I(1/3))
        }
        if (genomewideline) {
            p <- p + geom_hline(yintercept=genomewideline,
                                      colour=colGenomewideline)
        }
        p
        
    } else {
        stop("Make sure your data frame contains columns CHR, BP, and P")
    }
}