################################################################################
## R script for ResistoXplorer
## Description: R code for miscellaneous functions and R packages 
###################################################

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
    path <- system("which bash", intern=TRUE);
    if((length(path) == 0) && (typeof(path) == "character")){
        print("Could not find bash in the PATH!");
        return("NA");
    }
    return(path);
}
###############################################################################
GetUniqueEntries <- function(db.path, statement){
    mir.db <- dbConnect(SQLite(), db.path);
    query <- dbSendQuery(mir.db, statement);
    res <- fetch(query, n=-1);
    dbDisconnect(mir.db);
    res <- sort(unique(as.character(res[,1])));
    return (res);
}
###############################################################################
# new range [a, b]
rescale2NewRange <- function(qvec, a, b){
    q.min <- min(qvec);
    q.max <- max(qvec);
    if(length(qvec) < 50){
        a <- a*2;
    }
    if(q.max == q.min){
        new.vec <- rep(8, length(qvec));
    }else{
        coef.a <- (b-a)/(q.max-q.min);
        const.b <- b - coef.a*q.max;
        new.vec <- coef.a*qvec + const.b;
    }
    return(new.vec);
}
###############################################################################
`%fin%` <- function(x, table) {
  fmatch(x, table, nomatch = 0L) > 0L
}
###############################################################################
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
###############################################################################
# shorthand
ShowMemoryUse <- function(..., n=30) {
    print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
}
# note, try to use the fread, however, it has issues with 
# some windows 10 files "Line ending is \r\r\n. .... appears to add the extra \r in text mode on Windows"
# in such as, use the slower read.table method
.readDataTable <- function(fileName){
    dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE,na.strings=c("","NA")));
    if(class(dat) == "try-error"){
        # try to use "tr" to remove double return characters
        trFileName <- paste("tr -d \'\\r\' <", fileName);
        dat <- try(data.table::fread(trFileName, header=TRUE, check.names=FALSE, data.table=FALSE,na.strings=c("","NA")));
        if(class(dat) == "try-error"){
            print("Using slower file reader ...");
            formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
            if(formatStr == "txt"){
                dat <-try(read.table(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T,na.strings=c("","NA")));
            }else{ # note, read.csv is more than read.table with sep=","
                dat <-try(read.csv(fileName,header=TRUE,comment.char = "", check.names=F, as.is=T,na.strings=c("","NA")));
            }  
        }
    }
    return(dat);
}
###############################################################################
GetExtendRange <-function(vec, unit){
    var.max <- max(vec, na.rm=T);
    var.min <- min(vec, na.rm=T);
    exts <- (var.max - var.min)/unit;
    c(var.min-exts, var.max+exts);
}
###############################################################################
GetSampleNamesaftNorm <-function(){
    return(sample_names(dataSet$norm.phyobj));  
}
###############################################################################
GetAllSampleGrpInfo <- function(){
    #sample variable having more than one group will be selected as default
    sam_var<-which(sapply(sample_data(dataSet$norm.phyobj)[,sapply(sample_data(dataSet$norm.phyobj), is.factor)], nlevels)>1);
    if(length(sam_var)>0){
        return(names(sam_var[1]));
    } else {
        return(NULL);
    }
}
###############################################################################
GetMetaInfo <- function(){
    dataSet$sample_data<-dataSet$sample_data[sapply(dataSet$sample_data, function(x) length(unique(na.omit(x)))) > 1];
    return(colnames(dataSet$sample_data));
}
###############################################################################
GetGeneLvlInfo <- function(){
    return(c(colnames(dataSet$taxa_table)));  
}
###############################################################################
GetTaxaLvlInfo <- function(){
    return(c(colnames(dataSet$m.taxa_table)));  
}
###############################################################################
GetSampleGrpInfo <- function(class){
    sa <-sample_data(dataSet$norm.phyobj);
    return(levels(get_variable(dataSet$norm.phyobj, class)));
}
###############################################################################
GetSampleGrpNo <- function(class){
    return(length(levels(get_variable(dataSet$norm.phyobj, class))));
}
###############################################################################
GetXYaxisCompCCA <-function (ncomp){
    return(c(1:ncomp));
}
###############################################################################
GetTaxaFeatName <- function(taxlvl){
    if(taxlvl=="OTU"){
        return(taxa_names(dataSet$norm.phyobj));
    }else {
        taxa_table <-tax_table(dataSet$proc.phyobj);
        data <-merge_phyloseq(dataSet$norm.phyobj,taxa_table);
        nm <-unique(as.character(tax_table(data)[,taxlvl]));
        indx <-which(is.na(nm)==TRUE);
        nm[indx] <-"Unassigned";
        return(nm);
    }
}
###############################################################################
GetLowerTaxaLvlNm <- function(taxrank){
    indx<- which(colnames(tax_table(dataSet$proc.phyobj))==taxrank);
    rem<-ncol(tax_table(dataSet$proc.phyobj))-indx;
    return(colnames(tax_table(dataSet$proc.phyobj))[indx+1:rem]);
}
###############################################################################
GetHighTaxaLvlNm <- function(taxrank){
    if(taxrank=="OTU"){
        return(colnames(tax_table(dataSet$proc.phyobj))[1:length(colnames(tax_table(dataSet$proc.phyobj)))]);
    }else{ 
    indx <- which(colnames(tax_table(dataSet$proc.phyobj))==taxrank);
    rem <-(indx):1;
    return(colnames(tax_table(dataSet$proc.phyobj))[rev(rem)]);
    }
}
###############################################################################
GetSampleGrpUser <- function(class){ 
    return(levels(get_variable(dataSet$proc.phyobj, class)));
}
###############################################################################
GetNameMapCol <-function(colInx){
    return(analSet$resTable[,colInx]);
}
###############################################################################
GetResRowNames <- function(){
    return(rownames(analSet$resTable));
}
###############################################################################
GetMirResColDE <-function(colInx){
    return(analSet$resTable[,colInx]);
}
###############################################################################
GetResColNames <- function(){
    return(colnames(analSet$resTable));
}
###############################################################################
GetResMat <- function(){
    return(as.matrix(analSet$resTable));
}
###############################################################################
fast_tax_glom_first <- function(physeq, taxrank){

    # setup data. we are going to glob the OTU table based on the Class Taxonomy
    CN  <- which( rank_names(physeq) %in% taxrank); 
    tax <- as(access(physeq, "tax_table"), "matrix")[, 1:CN, drop=FALSE];
    tax <- apply(tax, 1, function(i){paste(i, sep=";_;", collapse=";_;")});    
    # using Map-Reduce/vectorized
    otab2 <- data.frame(otu_table(physeq))
    taxdf <- data.frame(tax);
    otab2 <- merge(otab2, taxdf, by = "row.names");
    row.names(otab2) <- otab2$Row.names
    otab2 <- otab2[ , 2:ncol(otab2)]
    otab2 <- condenseOTUs(otab2,"tax");
    otab2 <-otu_table(otab2,taxa_are_rows = T);
    colnames(otab2)<-sample_names(physeq);
    phy_data <-merge_phyloseq(otab2,tax_table(physeq),sample_data(physeq));    
    return(phy_data)
}
###############################################################################
condenseOTUs <- function(otutable, splitcol) {
     #'  Helper Function that will take the OTU Table (as a data.frame()) 
     #'  to which column representing taxa has been added.
     #'  There is some munging to handle rownames and
     #'  the tax column but the basic idea is there.  
      #convert tax col to numeric
      # needed to allow colSums to work
      otutable[[splitcol]] <- as.numeric(factor(otutable[[splitcol]], labels=c(1:length(unique(otutable[[splitcol]])))))

      # split apply, combine
      splits   <- split(otutable, otutable[[splitcol]])  
      summed   <- Map(colSums, splits)
      summeddf <- Reduce(rbind, summed, init = NULL)

      #get rownames.
      # this needs to be changed to get most abundant OTU if thats what is used elsewhere
      # rightnow it returns the first OTU in the group
      newrownames <- Map(function(x){rownames(x)[[1]]}, splits)

      #add back rowname and remove tax column
      rownames(summeddf) <- newrownames
      summeddf[, !colnames(summeddf) %in% c(splitcol)]
}

#######################################
###########Color-palette#############
#######################################
custom_col42 <-  c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4",
                "#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278",
                "#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5",
                "#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B",
                "#EA6C81","#F098A7");

custom_final28 <-c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777", 
                "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455", 
                "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E", 
                "#DDAA77");

custom_col21 <-  c("#771155","#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", 
                "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455", "#DD7788");
custom_col74 <-c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#A6CEE3","#6A3D9A",
              "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB",
              "#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",
              "#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#FFFF99","#FBB4AE","#B3CDE3",
              "#984EA3","#FF7F00","#FFFF33","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2",
              "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#666666","#1B9E77","#D95F02",
              "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#CCCCCC","#E41A1C","#377EB8","#4DAF4A","#D9D9D9",
              "#BC80BD","#CCEBC5","#FFED6F","#B2DF8A","#FF7F00","#B15928","#CAB2D6","#1F78B4","#A65628");
###############################################################################
GetColorSchema <- function(grayscale=F){
    # test if total group number is over 9
    claslbl <- as.factor(sample_data(dataSet$norm.phyobj)[[variable]]);
    grp.num <- length(levels(claslbl));
     if(grayscale){
        dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
        lvs <- levels(claslbl);
        colors <- vector(mode="character", length=length(claslbl));
        for(i in 1:length(lvs)){
            colors[analSet$cls == lvs[i]] <- dist.cols[i];
        }
    }else if(grp.num > 9){
        pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                    "#FFFF99", "#B15928");
        dist.cols <- colorRampPalette(pal12)(grp.num);
        lvs <- levels(claslbl);
        colors <- vector(mode="character", length=length(analSet$cls));
        for(i in 1:length(lvs)){
            colors[claslbl == lvs[i]] <- dist.cols[i];
        }
    }else{
        if(exists("colVec") && !any(colVec =="#NA") ){
            cols <- vector(mode="character", length=length(claslbl));
            clsVec <- as.character(claslbl);
            grpnms <- names(colVec);
            for(i in 1:length(grpnms)){
                cols[clsVec == grpnms[i]] <- colVec[i];
            }
            colors <- cols;
        }else{
            colors <- as.numeric(claslbl)+1;
        }
    }
    return (colors);
}
###############################################################################
ComputeColorGradient <- function(nd.vec, centered=TRUE){
    require("RColorBrewer");
    if(sum(nd.vec<0, na.rm=TRUE) > 0){ 
        centered <- T;
    }else{
        centered <- F;
    }
    color <- colorRampPalette(c("green", "yellow", "red"))(100);
    breaks <- generate_breaks(nd.vec, length(color), center = centered);
    return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}
###############################################################################
GetBetaDiversityStats <-function(){
    return(analSet$stat.info);
}
###############################################################################
GetProtestStats <-function(){
    return(dataSet$pro.stat.info);
}
###############################################################################
GetCIAtestStats <-function(){
    return(dataSet$cia.stat.info);
}
###############################################################################
generate_breaks = function(x, n, center = F){
    if(center){
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
        res = seq(-m, m, length.out = n + 1)
    }
    else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
    }
    return(res)
}
###############################################################################
scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    breaks <- sort(unique(breaks));
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}
###############################################################################
# based on phyloseq post: https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr_transform <- function(x, base=2){
    x <- log((x / gm_mean(x)), base)
    x[!is.finite(x) | is.na(x)] <- 0.0
    return(x)
}
###############################################################################
gm_mean <- function(x, na.rm=TRUE){
    # The geometric mean, with some error-protection bits.
    exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}
###############################################################################
#CLR transformation from robComposition R package
cenLR_rcom <- function(x, base = exp(1)){
    #if(dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
    if(dim(x)[2] == 1){
        res <- list(x.clr=x, gm=rep(1,dim(x)[1]))	    	
    } else{
        geometricmean <- function (x) {
            if (any(na.omit(x == 0)))
                0
            else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
        }
        gm <- apply(x, 1, geometricmean)
        x.clr <- log(x/gm, base);
        x.clr[!is.finite(x.clr) | is.na(x.clr)] <- 0;
        res <- as.matrix(t(x.clr));
    }
    return(res)  
}
###############################################################################
#ALR transformation from robComposition R package (Note: one feature reduced)
addLR_rcom <- function (x, ivar=ncol(x), base = exp(1)){
    if(dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
    x.alr <- log(x/x[, ivar], base);
    x.alr[!is.finite(x.alr) | is.na(x.alr)] <- 0;
    res <- as.matrix(t((x.alr[,-ivar])));
    return(res)
}
###############################################################################
cpm.default <- function(x,log,lib.size=NULL, prior.count=0.25, ...){
    x <- as.matrix(x)
    if(is.null(lib.size)) lib.size <- colSums(x)
    if(log) {
        prior.count.scaled <- lib.size/mean(lib.size)*prior.count
        lib.size <- lib.size+2*prior.count.scaled
    }
    lib.size <- 1e-6*lib.size
    if(log)
        log2(t( (t(x)+prior.count.scaled) / lib.size ))
    else
        t(t(x)/lib.size)
}
###############################################################################
GetRandomNumbers <- function(){
    rm(.Random.seed);
    runif(1);
    return(.Random.seed[3:626]);
}
###############################################################################
#################################################
##########Cross-correlation (InterOmics: associate R function)
#################################################
associate <- function(x, y=NULL, method, p.adj.threshold, p.adj.method,mode="table",filter.self.correlations=TRUE,cth=NULL, order=FALSE, n.signif=0) {
    if (is.null(y)) {
    message("Cross-correlating the data with itself")
    y <- x
        if (filter.self.correlations) {
          # Ignore self-correlations in filtering
          n.signif <- n.signif + 1
        }
    }
    x <- as.data.frame(x)  # numeric or discrete
    y <- y  # numeric
    if (is.null(colnames(y))) {
    colnames(y) <- paste("column-", seq_len(ncol(y)), sep="")
    }
    xnames <- colnames(x)
    ynames <- colnames(y)
    qv <- NULL
    numeric.methods <- c("spearman", "pearson")
    categorical.methods <- c("categorical")

    # Rows paired.
    if (method %in% numeric.methods) {
    inds <- vapply(x, is.numeric, TRUE)
    if (any(!inds)) {
      warning("Considering only numeric annotations for \n       
              pearson/spearman")
    }
    inds <- names(which(inds))
    } else if (method %in% categorical.methods) {
      inds <- vapply(x, is.factor, TRUE)
      if (any(!inds)) {
        warning("Considering only categorical annotations for factors")
      }
      inds <- names(which(inds))
    }
    xnames <- inds
    if (!is.vector(x)) {
    x <- suppressWarnings(as.matrix(x[, inds], ncol=length(inds)))
    } else {
    x <- as.matrix(x[inds], ncol=length(inds))
    }
    colnames(x) <- xnames
    Pc <- matrix(NA, ncol(x), ncol(y))
    Cc <- matrix(NA, ncol(x), ncol(y))
    rownames(Cc) <- colnames(x)
    colnames(Cc) <- colnames(y)
    rownames(Pc) <- colnames(x)
    colnames(Pc) <- colnames(y)
    if (method %in% c("pearson", "spearman")) {
    minobs <- 8
    for (j in seq_len(ncol(y))) {
      jc <- apply(x, 2, function(xi) {
        if (sum(!is.na(xi)) >= minobs) {
          res <- suppressWarnings(
            cor.test(xi, unlist(y[, j], use.names=FALSE), 
                     method=method, use="pairwise.complete.obs"))
          res <- c(res$estimate, res$p.value)
        } else {
          warning(paste("Not enough observations (",
                        minobs, "required); \n   
                        (", 
                        sum(!is.na(xi)), ") \n \n 
                        - skipping correlation estimation"))
          res <- c(NA, NA)
        }
        res
      })
      Cc[, j] <- jc[1, ]
      Pc[, j] <- jc[2, ]
    }
    
    #} else if (method == "bicor") {
    #    
    #    t1 <- suppressWarnings(
    #    bicorAndPvalue(x, y, use="pairwise.complete.obs"))
    #    Pc <- t1$p
    #    Cc <- t1$bicor
    #    
    #    
    } else if (method == "categorical") {
      Cc <- matrix(NA, nrow=ncol(x), ncol=ncol(y))
      rownames(Cc) <- colnames(x)
      colnames(Cc) <- colnames(y)
      for (varname in colnames(x)) {
        for (lev in colnames(y)) {
          xvec <- x[, varname]
          yvec <- y[, lev]
          keep <- rowSums(is.na(cbind(xvec, yvec))) == 0
          xvec <- xvec[keep]
          yvec <- yvec[keep]
          
          # Number of data-annotation samples for
          # calculating the correlations
          n <- sum(keep)
          Cc[varname, lev] <- gktau(xvec, yvec) 
        }
      }
    }
  if (!all(is.na(Pc))) {
    rownames(Pc) <- xnames
    colnames(Pc) <- ynames
    rownames(Cc) <- xnames
    colnames(Cc) <- ynames
    
    # Corrected p-values
    qv <- array(NA, dim=dim(Pc))
    qv <- matrix(p.adjust(Pc, method=p.adj.method), nrow=nrow(Pc))
    dimnames(qv) <- dimnames(Pc)
  }
  # Filter
  if (!is.null(p.adj.threshold) || !is.null(cth)) {
    
    # Replace NAs with extreme values for filtering purposes
    qv[is.na(qv)] <- 1
    Pc[is.na(qv)] <- 1
    Cc[is.na(Cc)] <- 0
    
    # Filter by adjusted pvalues and correlations
    inds1.q <- inds2.q <- inds1.c <- inds2.c <- NULL
    if (!is.null(p.adj.threshold)) {
      inds1.q <- apply(qv, 1, function(x) {
        sum(x < p.adj.threshold) >= n.signif
      })
      inds2.q <- apply(qv, 2, function(x) {
        sum(x < p.adj.threshold) >= n.signif
      })
    }
    if (!is.null(cth)) {
      inds1.c <- apply(abs(Cc), 1, function(x) {
        sum(x > cth) >= n.signif
      })
      inds2.c <- apply(abs(Cc), 2, function(x) {
        sum(x > cth) >= n.signif
      })
    }
    if (!is.null(p.adj.threshold) && !is.null(cth)) {
      
      inds1 <- inds1.q & inds1.c
      inds2 <- inds2.q & inds2.c
    } else if (is.null(p.adj.threshold) && !is.null(cth)) {
      inds1 <- inds1.c
      inds2 <- inds2.c
    } else if (!is.null(p.adj.threshold) && is.null(cth)) {
      inds1 <- inds1.q
      inds2 <- inds2.q
    }
    Cmat <- as.matrix(0)
    
    # TODO: add also correlation filter, not only significance
    # Require each has at least n.signif. correlations
    if (sum(inds1) >= n.signif && sum(inds2) >= n.signif) {
      rnams <- rownames(Cc)[inds1]
      cnams <- colnames(Cc)[inds2]
      Cc <- matrix(Cc[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      Pc <- matrix(Pc[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      qv <- matrix(qv[inds1, inds2, drop=FALSE], nrow=sum(inds1))
      rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
      colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams
      if (order && sum(inds1) >= 2 && sum(inds2) >= 2) {
        
        # Order in visually appealing order
        tmp <- Cc
        rownames(tmp) <- NULL
        colnames(tmp) <- NULL
        rind <- hclust(as.dist(1 - cor(t(tmp),
                                       use="pairwise.complete.obs")))$order
        cind <- hclust(as.dist(1 - cor(tmp,
                                       use="pairwise.complete.obs")))$order
        rnams <- rownames(Cc)[rind]
        cnams <- colnames(Cc)[cind]
        Cc <- Cc[rind, cind]
        Pc <- Pc[rind, cind]
        qv <- qv[rind, cind]
        rownames(qv) <- rownames(Pc) <- rownames(Cc) <- rnams
        colnames(qv) <- colnames(Pc) <- colnames(Cc) <- cnams
      }
    } else {
      message("No significant correlations with the given criteria\n")
      Cc <- Pc <- qv <- NULL
    }
  }
  res <- list(cor=Cc, pval=Pc, p.adj=qv)
  
  # message('Ignore self-correlations in filtering')
  if (nrow(x) == nrow(y) && ncol(x) == ncol(y) && filter.self.correlations) {
    diag(res$cor) <- diag(res$pval) <- diag(res$p.adj) <- NA
  }
  if (mode == "table") {
    res <- cmat2table(res)
  } 
  res
}
#################################################
#' @title Convert Correlation Matrix into a Table
#' @description Arrange correlation matrices from associate into a table format.
#' @param res Output from associate
#' @return Correlation table
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
################################################
cmat2table <- function(res, verbose=FALSE) {
    ctab <- ID <- NULL
    if (!is.null(res$cor)) {
        ctab <- as.data.frame(res$cor)
        ctab$ID <- rownames(res$cor)
        suppressMessages(require(reshape));
        ctab <- melt(ctab, "ID")
        colnames(ctab) <- c("X1", "X2", "Correlation")
        ctab$Correlation <- as.numeric(as.character(ctab$Correlation))
    }
    correlation <- NULL  # circumwent warning on globabl vars
    if (!is.null(res$p.adj)) {
        if (verbose) {
            message("Arranging the table")
        }
        ctab2 <- as.data.frame(res$p.adj)
        ctab2$ID <- rownames(res$p.adj)
        ctab2 <- melt(ctab2, "ID")
        colnames(ctab2) <- c("X1", "X2", "p.adj")
        ctab2$p.adj <- as.numeric(as.character(ctab2$p.adj))
        ctab <- cbind(ctab, ctab2$p.adj)
        colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
        ctab <- ctab[order(ctab$p.adj), ]
        colnames(ctab) <- c("X1", "X2", "Correlation", "p.adj")
    } else {
        message("No significant adjusted p-values")
        if (!is.null(ctab)) {
            ctab2 <- as.data.frame(res$pval)
            ctab2$ID <- rownames(res$pval)
            ctab2 <- melt(ctab2, "ID")
            colnames(ctab2) <- c("X1", "X2", "value")
            ctab2$value <- as.numeric(as.character(ctab2$value))
            ctab <- cbind(ctab, ctab2$value)
            ctab <- ctab[order(-abs(ctab$Correlation)), ]
            colnames(ctab) <- c("X1", "X2", "Correlation", "pvalue")
        }
    }
    ctab$X1 <- as.character(ctab$X1)
    ctab$X2 <- as.character(ctab$X2)
    
    # Keep the original order of factor levels
    ctab$X1 <- factor(as.character(ctab$X1), levels=rownames(res$cor))
    ctab$X2 <- factor(as.character(ctab$X2), levels=colnames(res$cor))
    
    # Remove NAs
    ctab <- ctab[!is.na(ctab$Correlation), ]
    
    # Order the table by p-value
    if ("p.adj" %in% colnames(ctab)) {
        ctab <- ctab[order(ctab$p.adj), ]
    } else if ("pvalue" %in% colnames(ctab)) {
        ctab <- ctab[order(ctab$pvalue), ]
    }
    ctab
}
###############################################################################
ANCOM <- function (OTUdat, sig = sig, multcorr = multcorr, tau = 0.02, theta = 0.1, 
                   repeated = FALSE, ncore = 2) {
    num_col <- ncol(OTUdat)
    if (repeated == FALSE) {
        colnames(OTUdat)[num_col] <- "Group"
        num_OTU <- ncol(OTUdat) - 1
        sub_drop <- data.frame(nm_drop = "N/A")
        sub_keep <- data.frame(nm_keep = "All subjects")
        colnames(sub_drop) <- "Subjects removed"
        colnames(sub_keep) <- "Subjects retained"
        n_summary <- paste0("No subjects entirely removed (not a repeated-measures design)")
    }else {
        colnames(OTUdat)[num_col - 1] <- "Group"
        colnames(OTUdat)[num_col] <- "ID"
        OTUdat$ID <- factor(OTUdat$ID)
        num_OTU <- ncol(OTUdat) - 2
        crossTab <- table(OTUdat$Group, OTUdat$ID) == 0
        id_drop <- apply(crossTab, 2, FUN = function(x) any(x))
        nm_drop <- names(which(id_drop))
        idx_drop <- OTUdat$ID %in% nm_drop
        OTUdat <- OTUdat[idx_drop == FALSE, ]
        if (nrow(OTUdat) == 0) {
          stop("Too many missing values in data, all subjects dropped")
        }
        OTUdat$ID <- droplevels(OTUdat$ID)
        num_dropped <- sum(id_drop)
        num_retain <- length(id_drop) - num_dropped
        sub_drop <- as.data.frame(nm_drop = paste(nm_drop, collapse = ", "))
        sub_keep <- as.data.frame(nm_keep = paste(levels(OTUdat$ID), 
                                                  collapse = ", "))
        colnames(sub_drop) <- "Subjects removed"
        colnames(sub_keep) <- "Subjects retained"
        n_summary <- paste0("Analysis used ", num_retain, " subjects (", 
                            num_dropped, " were removed due to incomplete data)")
    }
    OTUdat$Group <- factor(OTUdat$Group)
    OTUdat <- as.data.frame(OTUdat[which(is.na(OTUdat$Group) == 
                                         FALSE), ], row.names = NULL)
    W.detected <- ancom.detect(OTUdat, num_OTU, sig, multcorr, 
                             ncore = ncore);
    W_stat <- W.detected
    if (num_OTU < 10) {
    detected <- colnames(OTUdat)[which(W.detected > num_OTU - 
                                         1)]
    }else {
    if (max(W.detected)/num_OTU >= theta) {
            c.start <- max(W.detected)/num_OTU
            cutoff <- c.start - c(0.05, 0.1, 0.15, 0.2, 0.25)
            prop_cut <- rep(0, length(cutoff))
            for (cut in 1:length(cutoff)) {
            prop_cut[cut] <- length(which(W.detected >= num_OTU * 
                                            cutoff[cut]))/length(W.detected)
            }
            del <- rep(0, length(cutoff) - 1)
            for (ii in 1:(length(cutoff) - 1)) {
            del[ii] <- abs(prop_cut[ii] - prop_cut[ii + 1])
            }
            if (del[1] < tau & del[2] < tau & del[3] < tau) {
            nu = cutoff[1]
            }
            else if (del[1] >= tau & del[2] < tau & del[3] < 
                   tau) {
            nu = cutoff[2]
            }
            else if (del[2] >= tau & del[3] < tau & del[4] < 
                   tau) {
            nu = cutoff[3]
            }
            else {
            nu = cutoff[4]
            }
            up_point <- min(W.detected[which(W.detected >= nu * 
                                             num_OTU)])
            W.detected[W.detected >= up_point] <- 99999
            W.detected[W.detected < up_point] <- 0
            W.detected[W.detected == 99999] <- 1
            detected <- colnames(OTUdat)[which(W.detected == 
                                               1)]
    }
    else {
        W.detected <- 0
        detected <- "No significant features detected"
    }
    }
    results <- list(W = W_stat, detected = detected, dframe = OTUdat, 
                  repeated = repeated, n_summary = n_summary, sub_drop = sub_drop, 
                  sub_keep = sub_keep)
    class(results) <- "ancom"
    return(results)
    }
    ancom.detect <- function (otu_data, n_otu, alpha, multcorr, ncore) {
    if (ncol(otu_data) == n_otu + 1) {
    Group <- otu_data[, ncol(otu_data)]
    ID <- rep(1, nrow(otu_data))
    repeated <- FALSE
    fformula <- formula("lr ~ Group")
    }
    else if (ncol(otu_data) == n_otu + 2) {
    Group <- otu_data[, ncol(otu_data) - 1]
    ID <- otu_data[, ncol(otu_data)]
    repeated <- TRUE
    fformula <- formula("lr ~ Group | ID")
    }
    else {
    stop("Problem with data. Dataset should contain OTU abundances, groups, \n 
        and optionally an ID for repeated measures.")
    }
    if (repeated == FALSE) {
    if (length(unique(Group)) == 2) {
      tfun <- stats::wilcox.test
    }
    else {
      tfun <- stats::kruskal.test
    }
    }
    else {
    tfun <- stats::friedman.test
    }
    if (FALSE) {
    registerDoParallel(cores = ncore)
    aa <- bb <- NULL
    logratio.mat <- foreach(bb = 1:n_otu, .combine = "rbind", 
                    .packages = "foreach") %:% foreach(aa = 1:n_otu, 
                       .combine = "c", .packages = "foreach") %dopar% {
                         if (aa == bb) {
                           p_out <- NA
                         }
                         else {
                           data.pair <- otu_data[, c(aa, bb)]
                           lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                                                                         as.numeric(data.pair[, 2])))
                           lr_dat <- data.frame(lr = lr, Group = Group, 
                                                ID = ID)
                           p_out <- tfun(formula = fformula, data = lr_dat)$p.value
                         }
                         p_out
                       }
    rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
    }else {
    logratio.mat <- matrix(NA, nrow = n_otu, ncol = n_otu)
    for (ii in 1:(n_otu - 1)) {
        for (jj in (ii + 1):n_otu) {
        data.pair <- otu_data[, c(ii, jj)]
        lr <- log((1 + as.numeric(data.pair[, 1]))/(1 + 
                                                      as.numeric(data.pair[, 2])))
        lr_dat <- data.frame(lr = lr, Group = Group, 
                             ID = ID)
        inx1 <- which(lr_dat$Group =="1");
        inx2 <- which(lr_dat$Group =="2");
        logratio.mat[ii, jj] <- tfun(lr_dat[inx1,"lr"],lr_dat[inx2,"lr"])$p.value
      }
    }
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
    }
    logratio.mat[which(is.finite(logratio.mat) == FALSE)] <- 1
    mc.pval <- t(apply(logratio.mat, 1, function(x) {
    s <- p.adjust(x, method = "BH")
    return(s)
    }))
    a <- logratio.mat[upper.tri(logratio.mat, diag = FALSE) == 
                      TRUE]
    b <- matrix(0, ncol = n_otu, nrow = n_otu)
    b[upper.tri(b) == T] <- p.adjust(a, method = "BH")
    diag(b) <- NA
    ind.1 <- lower.tri(b)
    b[ind.1] <- t(b)[ind.1]
    if (multcorr == 1) {
    W <- apply(b, 1, function(x) {
      subp <- length(which(x < alpha))
    })
    }
    else if (multcorr == 2) {
    W <- apply(mc.pval, 1, function(x) {
      subp <- length(which(x < alpha))
    })
    }
    else if (multcorr == 3) {
    W <- apply(logratio.mat, 1, function(x) {
      subp <- length(which(x < alpha))
    })
    }
    return(W)
}
###############################################################################
#Alpha diversity eveness based indexes
# Source: microbiome R package (Lahti L, Shetty S)
##############################################
evenness <- function(x, index="all", zeroes=TRUE, detection = 0) { 
    # Only include accepted indices
    index <- tolower(index)    
    accepted <- c("camargo", "pielou", "simpson", "evar", "bulla")    
    accepted <- tolower(accepted)
    # Return all indices
    if (length(index) == 1 && index == "all") {
    index <- accepted
    }
    if (!is.null(index)) {
    index <- intersect(index, accepted)
    }
    if (!is.null(index) && length(index) == 0) {
    return(NULL)
    }
    if (detection > 0) {
    x[x <= detection] <- 0
    }
    tab <- evenness_help(x, index, zeroes)
    if (is.vector(tab)) {
    tab <- as.matrix(tab, ncol=1)
    colnames(tab) <- index        
    }
    as.data.frame(tab)
}
###############################################################################
evenness_help <- function(x, index="all", zeroes=TRUE) { 
    # Only include accepted indices
    accepted <- c("camargo", "pielou", "simpson", "evar", "bulla")

    # Return all indices
    if ("all" %in% index) {
    index <- accepted
    }
    index <- intersect(index, accepted)
    if (length(index) == 0) {
    return(NULL)
    }
    if (length(index) > 1) {
    ev <- NULL
    for (idx in index) {
      ev <- cbind(ev,
                  evenness_help(x, index=idx, zeroes))
    }
    colnames(ev) <- index
    return(as.data.frame(ev))
    }

    # Pick data
    otu <- abundances(x)
    if (index == "camargo") {
    ev <- apply(otu, 2, function(x) {
      camargo(x, zeroes=zeroes)
    })
    } else if (index == "simpson") {
    ev <- apply(otu, 2, function(x) {
      simpson_evenness(x)
    })
    } else if (index == "pielou") {
    ev <- apply(otu, 2, function(x) {
      pielou(x)
    })
    } else if (index == "evar") {
    ev <- apply(otu, 2, function(x) {
      evar(x, zeroes=zeroes)
    })
    } else if (index == "bulla") {
    ev <- apply(otu, 2, function(x) {
      bulla(x, zeroes=zeroes)
    })
    }
    names(ev) <- colnames(otu)
    ev
}
###############################################################################
bulla <- function(x, zeroes=TRUE) {
    if (!zeroes) {
    x[x > 0]
    }
    # Species richness (number of species)
    S <- sum(x>0, na.rm = TRUE)

    # Relative abundances
    p <- x/sum(x)
    O <- sum(pmin(p, 1/S))  
    # Bulla's Evenness
    (O - 1/S)/(1 - 1/S)
}
###############################################################################
# Camargo's eveness x: species counts zeroes: include zeros Inspired
# by code from Pepijn de Vries and Zhou Xiang at
# researchgate.net/post/How_can_we_calculate_the_Camargo_evenness_index_in_R
# but rewritten here
camargo <- function(x, zeroes=TRUE) {
    if (!zeroes) {
    x[x > 0]
    }
    N <- sum(x > 0, na.rm = TRUE)
    xx <- 0
    for (i in seq_len(N - 1)) {
    xx <- xx + sum(abs(x[(i + 1):N] - x[i]))
    }
    # Return
    1 - xx/(sum(x) * N)
}
###############################################################################
# x: Species count vector
pielou <- function(x) {
    # Remove zeroes
    x <- x[x > 0]
    # Species richness (number of detected species)
    S <- sum(x > 0, na.rm = TRUE)
    # Relative abundances
    p <- x/sum(x)
    # Shannon index
    H <- (-sum(p * log(p)))
    # Simpson evenness
    H/log(S)
}
###############################################################################
# Smith and Wilson’s Evar index
evar <- function(x, zeroes=TRUE) {
    if (!zeroes) {
    x[x > 0]
    }
    n <- sum(x, na.rm = TRUE)
    d <- rep(NA, n)

    # Log abundance
    a <- ifelse(x != 0, log(x), 0)

    # Richness
    S <- sum(x > 0)
    b <- a/S
    c <- sum(b)
    d <- ifelse(x != 0, (a - c)^2/S, 0)
    f <- sum(d)
    (1 - 2/pi * atan(f))
}
###################################################
################################################################################
# Source: cclasso.R (Author : Fang Huaying (Peking University))
#
###############################################################################
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = k_max, n_boot = n_boot) {
    n <- nrow(x);
    p <- ncol(x);
    if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
    }
    x <- log(x);
    vx2 <- stats::var(x);

    # Diagonal weight for loss function
    rmean_vx2 <- rowMeans(vx2);
    wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
    wd2 <- sqrt(wd);

    # Some global parameters for optimization with single lambda
    rho <- 1;
    u_f <- eigen(diag(p) - 1/p)$vectors;
    wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
    wd_u_eig <- eigen(wd_u);
    d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
    (2 * rho) + 1 );
    u0_wd <- wd_u_eig$vectors;

    # Golden section method for the selection of lambda (log10 scale)
    sigma <- vx2;
    lam_int2 <- log10(range(lam_int));
    a1 <- lam_int2[1]; 
    b1 <- lam_int2[2];

    # Store lambda and corresponding cross validation's loss
    lams <- NULL; 
    fvals <- NULL;

    # Two trial points in first 
    a2 <- a1 + 0.382 * (b1 - a1); 
    b2 <- a1 + 0.618 * (b1 - a1);
    fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
    sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
    lams <- c(lams, b2); 
    fvals <- c(fvals, fb2$cv_loss);
    fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
    sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
    lams <- c(lams, a2); 
    fvals <- c(fvals, fa2$cv_loss);

    # Error tolerance for convergence
    err_lam2 <- 1e-1 * max(1, lam_int2);
    err_fval <- 1e-4;

    err <- b1 - a1;
    k <- 0;
    while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);

    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
        sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);
      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);
      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);  
    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
    }
    info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
    lam_int = 10^c(a1, b1)); 
    if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
      10^b1, "])\n", sep = "");
    }
    lambda <- 10^((a2 + b2)/2);

    # Bootstrap for cclasso
    lambda2 <- lambda / rho;
    info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
    n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
    p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}
###############################################################################
# Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);

  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);

  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
    ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
      lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);
  
  #----------------------------------------  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  #----------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
    lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
    (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  #----------------------------------------

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
###############################################################################
# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);
  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);
    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
      wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }
  return(list(cv_loss = cv_loss, sigma = sigma));
}
###############################################################################
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;
  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
      d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);
    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
      abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }
  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
      "&& Relative error:", err, "!\n");
  }
  return(sigma);
}
###############################################################################
#############
#Plot variable source R code (mixOmics: R package; version: 3.8) 
#############
#########################################################
plotVar2 <- function(object,
    comp = NULL,
    comp.select = comp,
    plot = TRUE,
    var.names = NULL,
    blocks = NULL,
    X.label = NULL,
    Y.label = NULL,
    Z.label = NULL,
    abline = TRUE,
    col,
    cex,
    pch,
    font,
    cutoff = 0,
    rad.in = 0.5,
    title = "Correlation Circle Plots",
    legend = FALSE,
    legend.title = "Block",
    style="ggplot2", # can choose between graphics,3d, lattice or ggplot2,
    overlap = TRUE,
    axes.box = "all",
    label.axes.box = "both"){
        class.object = class(object)
        object.pls=c("mixo_pls","mixo_spls","mixo_mlspls","mixo_mlsplsda","rcc")
        object.pca=c("ipca","sipca","pca","spca")
        object.blocks=c("sgcca","rgcca")

        #-- check that the user did not enter extra arguments
        arg.call = match.call()
        user.arg = names(arg.call)[-1]
        err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
        error = function(e) e)
        if ("simpleError" %in% class(err))
        stop(err[[1]], ".", call. = FALSE)

        #-- style
        if (!style %in% c("ggplot2", "lattice", "graphics","3d"))
        stop("'style' must be one of 'ggplot2', '3d' , lattice' or 'graphics'.", call. = FALSE)

        #-- plot
        if (length(plot) > 1)
        stop("'plot' must be single logical value.", call. = FALSE)
        else if (!is.logical(plot))
        stop("'plot' must be logical.", call. = FALSE)
        if(!plot)
        {
            style="N"}

        ### Start: Validation of arguments
        ncomp = object$ncomp
        if (any(class.object %in% object.blocks))
        {
            #-- legend: logical only
            if (length(legend) != 1 || !is.logical(legend))
            stop("'legend' must be a logical value.", call. = FALSE)
            if (is.null(blocks))
            {
                blocks = names(object$X)#names$blocks
                if (any(class.object == "DA"))
                blocks = names(object$X)#blocks[-object$indY]
            } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <= length(object$names$blocks)) {
                blocks = object$names$blocks[blocks]
            } else if (is.character(blocks)) {
                if (!any(blocks %in% object$names$blocks))
                stop("One element of 'blocks' does not match with the names of the blocks")
            } else {
                stop("Incorrect value for 'blocks'", call. = FALSE)
            }
            object$variates = object$variates[names(object$variates) %in% blocks]
            object$names$colnames = object$names$colnames[names(object$names$colnames) %in% blocks]
            object$blocks = object$X[names(object$X) %in% blocks]
            if (any(object$ncomp[blocks] == 1))
            {
                stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "),"' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
            }
            ncomp = object$ncomp[blocks]
        } else if (any(class.object %in% c("rcc", "mixo_pls", "mixo_spls", "mixo_mlspls")) & all(class.object !="DA")) {
            #-- legend: logical or a name for X and Y
            if ( ! (length(legend) == 1 & is.logical(legend) || (length(legend)==2)))
            stop("'legend' must be a logical value or a vector of 2 names for X and Y.", call. = FALSE)
            if(length(legend)==2){
                blocks=legend
                legend=TRUE
            } else{
                blocks = c("X","Y")
            }
        } else {
            #-- legend: logical or a name for X
            if (length(legend) != 1 )
            stop("'legend' must be a logical value or a vector of 1 name for X.", call. = FALSE)
            if(is.logical(legend)){
                blocks = "X"
            }else{
                blocks = legend
                legend=TRUE
            }
        }
        #-- legend.title
        if (length(legend.title)>1)
        stop("'legend.title' needs to be a single value (length 1)")

        #-- ellipse.level
        if (!is.numeric(rad.in) | (rad.in > 1) | (rad.in < 0))
        stop("The value taken by 'rad.in' must be between 0 and 1", call. = FALSE)

        #-- cutoff correlation
        if (!is.numeric(cutoff) | (cutoff > 1) | (cutoff < 0))
        stop("The value taken by 'cutoff' must be between 0 and 1", call. = FALSE)

        #-- comp
        if(is.null(comp))
        {
            if (style=="3d")
            {
                comp = seq_len(3)
            } else {
                comp = seq_len(2)
            }
        }
        if (length(comp) != 2 && !(style=="3d"))
        {
            stop("'comp' must be a numeric vector of length 2.", call. = FALSE)
        } else if(length(comp) != 3 && (style=="3d")) {
            stop("'comp' must be a numeric vector of length 3.", call. = FALSE)
        }

        if (!is.numeric(comp))
        stop("Invalid vector for 'comp'.")

        if (any(ncomp < max(comp)) || min(comp) <= 0)
        stop("Each element of 'comp' must be positive smaller or equal than ", min(object$ncomp), ".", call. = FALSE)

        comp1 = round(comp[1])
        comp2 = round(comp[2])
        if (style=="3d")
        comp3 = round(comp[3])

        #-- comp.select
        if (!is.null(comp.select))
        {
            if (!is.numeric(comp.select))
            stop("Invalid vector for 'comp'.", call. = FALSE)

            if (any(ncomp < max(comp.select)) || min(comp.select) <= 0)
            stop("Each element of 'comp.select' must be positive and smaller or equal than ", max(object$ncomp), ".", call. = FALSE)
        } else {
            comp.select = comp
        }

        #-- abline
        if (length(abline) > 1)
        {
            stop("'abline' must be single logical value.", call. = FALSE)
        }else if (!is.logical(abline)) {
            stop("'abline' must be logical.", call. = FALSE)
        }

        #-- Start: Retrieve variates from object
        cord.X = sample.X = ind.var.sel = list()
            if (any(class.object %in%  c(object.pls, object.blocks)))
            {
                if (any(class.object == "rcc"))
                {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                } else if (any(class.object %in% "mixo_plsda")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                } else if (any(class.object %in%  "mixo_pls")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, if(object$mode ==  "canonical"){object$variates$Y[, c(comp1, comp2)]} else {object$variates$X[, c(comp1, comp2)]}, use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                } else if (any(class.object %in%  c("mixo_splsda", "mixo_mlsplsda"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$name}))), drop = FALSE],
                    object$variates$X[, unique(c(comp1, comp2))], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                    
                    #if (!is.null(comp.select)) {
                    #   cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                } else if (any(class.object %in%  c("mixo_spls", "mixo_mlspls"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$X$name}))), drop = FALSE],
                    object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y[, colnames(object$Y) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$Y$name}))), drop = FALSE],
                    if(object$mode ==  "canonical")
                    {
                        object$variates$Y[, c(comp1, comp2)]
                    } else {
                        object$variates$X[, c(comp1, comp2)]
                    }, use = "pairwise")
                    
                    #ind.var.sel[[1]] =
                    sample.X[[1]] = seq_len(length(colnames(object$X)))
                    
                    #ind.var.sel[[2]] =
                    sample.X[[2]] = seq_len(length(colnames(object$Y)))
                    
                    #if (!is.null(comp.select)) {
                    #   cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$X$name}))), ,drop = FALSE]
                    #   cord.X[[2]] = cord.X[[2]][row.names(cord.X[[2]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$Y$name}))), , drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                    ind.var.sel[[2]] = which(colnames(object$Y) %in% rownames(cord.X[[2]]))
                } else { #block object
                    cord.X = lapply(blocks, function(x){cor(object$blocks[[x]], object$variates[[x]][, c(comp1, comp2)], use = "pairwise")})
                    ind.var.sel = sample.X = lapply(object$blocks, function(x){seq_len(ncol(x))})
                    if (!is.null(comp.select))
                    {
                        cord.X = lapply(seq_len(length(cord.X)), function(z){cord.X[[z]][row.names(cord.X[[z]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, block = blocks[z], comp = x)[[1]]$name}))), ,drop = FALSE]})
                    }
                    for (i in seq_len(length(cord.X)))
                    {
                        ind.var.sel[[i]] = which(colnames(object$blocks[[i]]) %in% rownames(cord.X[[i]]))
                    }
                }
            } else if (any(class.object %in%  object.pca)) {
                if (any(class.object %in%  c("sipca", "spca"))){
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$name}))), drop = FALSE],
                    object$x[, c(comp1, comp2)], use = "pairwise")
                    
                    #ind.var.sel[[1]] =
                    sample.X[[1]] = seq_len(length(colnames(object$X)))
                    
                    #if (!is.null(comp.select)) {
                    #    cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                } else {
                    cord.X[[1]] = cor(object$X, object$x[, c(comp1, comp2)], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                }
            }

        # output a message if some variates are anti correlated among blocks
        if (any(class.object %in%  object.blocks))
        {
            VarX = lapply(seq_len(2), function(j){do.call(cbind, lapply(object$variates, function(i) i[, comp[j]]))})
            corX = lapply(VarX, cor)
            if(any(sapply(corX, function(j){any(j < 0)})))
            warning("We detected negative correlation between the variates of some blocks, which means that some clusters of variables observed on the correlation circle plot are not necessarily positively correlated.")
        }
        if (any(sapply(cord.X, nrow) == 0))
        stop("No variable selected on at least one block")

        #-- End: Retrieve variates from object
        #-- Names of labels X and Y
        if (is.null(X.label)) X.label = paste("Component ",comp1)
        if (is.null(Y.label)) Y.label = paste("Component ",comp2)
        if (is.null(Z.label) && style=="3d") Z.label = paste("Component ",comp3)
        if (!is.character(X.label))
        stop("'X.label' must be a character.", call. = FALSE)
        if (!is.character(Y.label))
        stop("'Y.label' must be a character.", call. = FALSE)

        #-- pch argument
        missing.pch = FALSE
        if (missing(pch))
        {
            missing.pch = TRUE
            if(style=="3d")
            {
                pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], sum(sapply(cord.X[x], nrow)))}))
            } else {
                pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(seq_len(20)[x], sum(sapply(cord.X[x], nrow)))}))
            }
        } else if (((is.vector(pch, mode = "double") || is.vector(pch, mode = "integer")) && !(style=="3d"))
        || (is.vector(pch, mode = "character") && style=="3d")) {
            if (length(pch) != length(sample.X))
            stop.message('pch', sample.X)
            pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(pch[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(pch)) {
            if (length(pch) != length(sample.X) || length(unlist(pch)) != sum(sapply(sample.X, length)))
            stop.message('pch', sample.X)
            if (length(ind.var.sel) != 0)
            pch = lapply(seq_len(length(pch)), function(x){pch[[x]][ind.var.sel[[x]]]})
            pch = unlist(pch)
        } else if (style=="3d") {
            if (!all(pch %in% c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")) && style=="3d")
            stop("pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
            call. = FALSE)
        }
        else {
            stop.message('pch', sample.X)
        }

        #-- col argument
        if (missing(col)) {
            if (length(cord.X) < 10) {
                col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(color.mixo(x), sum(sapply(cord.X[x], nrow)))}))
            } else {
                col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(color.jet(length(cord.X))[x], sum(sapply(cord.X[x], nrow)))}))
            }
        } else if (is.vector(col, mode = "double") | is.vector(col, mode = "character")) {
            if (length(col) != length(sample.X))
            stop.message('col', sample.X)
            col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(col[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(col)) {
            if (length(col) != length(sample.X) || length(unlist(col)) != sum(sapply(sample.X, length)))
            stop.message('col', sample.X)
            if (length(ind.var.sel) != 0)
            col = lapply(seq_len(length(col)), function(x){col[[x]][ind.var.sel[[x]]]})
            col = unlist(col)
        } else {
            stop.message('col', sample.X)
        }

        #-- cex argument
        if (missing(cex)){
            if (style == "ggplot2"){
                cex = rep(5, sum(sapply(cord.X, nrow)))
            } else {
                cex = rep(1, sum(sapply(cord.X, nrow)))
            }
        } else if (is.vector(cex, mode = "double")) {
            if (length(cex) != length(cord.X))
            stop.message('cex', sample.X)
            cex = unlist(lapply(seq_len(length(cord.X)), function(x){rep(cex[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(cex)) {
            if (length(cex) != length(sample.X) || length(unlist(cex)) != sum(sapply(sample.X, length)))
            stop.message('cex', sample.X)
            if (length(ind.var.sel) != 0)
            cex = lapply(seq_len(length(cex)), function(x){cex[[x]][ind.var.sel[[x]]]})
            cex = unlist(cex)
        } else {
            stop.message('cex', sample.X)
        }

        #-- font argument
        if (missing(font)) {
            font = rep(1, sum(sapply(cord.X, nrow)))
        } else if (is.vector(font, mode = "numeric")) {
            if (length(font) != length(cord.X))
            stop.message('font', sample.X)
            font = unlist(lapply(seq_len(length(cord.X)), function(x){rep(font[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(font)) {
            if (length(font) != length(sample.X) || length(unlist(font)) != sum(sapply(sample.X, length)))
            stop.message('font', sample.X)
            if (length(ind.var.sel) != 0)
            font = lapply(seq_len(length(font)), function(x){font[[x]][ind.var.sel[[x]]]})
            font = unlist(font)
        } else {
            stop.message('font', sample.X)
        }

        #-- var.names
        ind.group = cumsum(c(0, sapply(cord.X, nrow)))
        if (is.null(var.names)){
            var.names.list = unlist(sapply(cord.X, rownames))
            if (!missing.pch) {
                var.names = rep(FALSE, length(cord.X))
            } else {
                var.names = rep(TRUE, length(cord.X))
            }
        } else if (is.vector(var.names, mode = "logical")) {
            if (length(var.names) == 1){
                var.names = rep(var.names,length(cord.X))}
            else if (length(var.names) != length(cord.X))
            stop.message('var.names', sample.X)
            var.names.list = unlist(lapply(seq_len(length(var.names)), function(x){if(var.names[x]){rownames(cord.X[[x]])}
                else {pch[(ind.group[x] + 1) : ind.group[x + 1]]}}))
        } else if (is.list(var.names)) {
            if (length(var.names) != length(cord.X))
            stop.message('var.names', sample.X)
            if (sum(sapply(seq_len(length(var.names)), function(x){if(!lapply(var.names, is.logical)[[x]]){
                if(is.null(ind.var.sel[[x]])){
                    length(var.names[[x]])
                } else {
                    length(var.names[[x]][ind.var.sel[[x]]])
                }
            } else {0}})) !=
            sum(sapply(seq_len(length(var.names)), function(x){if(!lapply(var.names, is.logical)[[x]]){nrow(cord.X[[x]])}else {0}}))){
                stop.message('var.names', sample.X)
            }
            var.names.list = unlist(sapply(seq_len(length(var.names)), function(x){if(lapply(var.names, is.logical)[[x]]){
                if (var.names[[x]]) {
                    row.names(cord.X[[x]])
                } else {
                    pch[(ind.group[x] + 1) : ind.group[x + 1]]
                }
            } else {
                if (is.null(ind.var.sel[[x]])){
                    as.character(var.names[[x]])
                } else {
                    as.character(var.names[[x]])[ind.var.sel[[x]]]
                }
            }
            }))
            var.names = sapply(var.names, function(x){if(is.logical(x)){x}else{TRUE}})
        } else {
            stop.message('var.names', sample.X)
        }
    suppressMessages(require("ellipse"));
        
    #-- Start: Computation ellipse
        circle = list()
        circle[[1]] = ellipse(0, levels = 1, t = 1)
        circle[[2]] = ellipse(0, levels = 1, t = rad.in)
        circle = data.frame(do.call("rbind", circle), "Circle" = c(rep("Main circle", 100), rep("Inner circle", 100)))
        
        #-- End: Computation ellipse
        #-- Start: data set
        df = data.frame(do.call(rbind, cord.X), "Block" = paste0("Block: ", unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))}))))
        if (style=="3d")
        names(df)[seq_len(3)] = c("x", "y","z")
        else
        names(df)[seq_len(2)] = c("x", "y")
        df$names = as.vector(var.names.list)
        df$pch = pch; df$cex = cex; df$col = col; df$font = font
        if(missing.pch)
        df$pch=1
        if (overlap)
        {
            df$Overlap = title
            df$Block = factor(unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))})))
            if(style %in%c("ggplot2","lattice"))
            title=NULL # to avoid double title
        } else {
            df$Overlap = df$Block
            if(style %in%c("ggplot2","lattice"))
            df$Block = factor(unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))})))
        }
        if (cutoff != 0){
            if(style=="3d")
            df = df[abs(df$x) > cutoff | abs(df$y) > cutoff | abs(df$z) > cutoff, ,drop = FALSE]
            else
            df = df[abs(df$x) > cutoff | abs(df$y) > cutoff, ,drop = FALSE]
            ind.group = c(0, cumsum(table(df$Block)[unique(df$Block)])) # add unique to have names of cumsum matching the order of the blocks in df
        }
        if (nrow(df) == 0)
        stop("Cutoff value very high for the components ", comp1, " and ", comp2, ".No variable was selected.")

        #-- End: data set
        #save(list=ls(),file="temp.Rdata")
        #-- Start: ggplot2
        if (style == "ggplot2" &  plot)
        {
            Block = NULL# R check
            # visible variable issues for x, y and Circle
            # according to http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
            # one hack is to set to NULL first.
            x = y = Circle = NULL

            #-- Initialise ggplot2
            p = ggplot(df, aes(x = x, y = y, color = Block)) +
            labs(title = title, x = X.label, y = Y.label) + theme_bw()
            for (i in levels(df$Block))
            p = p + geom_point(data = subset(df, df$Block == i), size = 0, shape = 0)

            #-- Display sample or var.names
            for (i in seq_len(length(var.names))){
                if (var.names[i]) {
                    p = p + geom_text(data = df[c((ind.group[i] + 1) : ind.group[i + 1]), ],
                    label = df[c((ind.group[i] + 1) : ind.group[i + 1]), "names"],
                    size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                    color = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                    fontface = df[c((ind.group[i] + 1) : ind.group[i + 1]), "font"])
                } else {
                    p = p + geom_point(data = df[c((ind.group[i] + 1) : ind.group[i + 1]), ],
                    size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                    shape = df[c((ind.group[i] + 1) : ind.group[i + 1]), "pch"],
                    color = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"])
                }
            }

            #-- Modify scale colour - Change X/Ylabel - split plots into Blocks
            p = p + scale_colour_manual(values = unique(col)[match(levels(factor(as.character(df$Block))), levels(df$Block))], name = "Features", breaks = levels(df$Block),labels=c("Resistome","Microbiome"))
            p = p + scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1))
            p = p + facet_wrap(~ Overlap, ncol = 2, as.table = TRUE)

            #-- Legend
            if (!legend)
            {
                p = p + theme(legend.position="none")
            } else {
                p = p + guides(colour = guide_legend(override.aes = list(shape = 19,
                size = unique(df$cex))))
            }
            #-- abline
            if (abline)
            p = p + geom_vline(aes(xintercept = 0), linetype = 2,
            colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,
            colour = "darkgrey")

            #-- circle correlation
            for (i in c("Main circle", "Inner circle")){
                p = p + geom_path(data = subset(circle, Circle == i),
                aes_string(x = "x", y = "y"), color = "Black")
            }
            p = p + theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,
               axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,
               legend.text=black.bold.10.text,legend.title = black.bold.11.title);    
            print(p)
        }
        #-- End: ggplot2

        #-- Start: 3d
        #-- End: graphics
        if(plot){
            return(invisible(df))}
        else
        return(df)
    
}
###############################################################################
ComputeColorGradientCorr <- function(nd.vec, centered=TRUE){
    library(RColorBrewer)
    if(sum(nd.vec<0, na.rm=TRUE) > 0){
    centered <- T;
    }else{
    centered <- F;
    }
    color <- c(grDevices::colorRampPalette(c("#BD0313", "#E32636", "#FF7783", "#FF9DA6"))(50), grDevices::colorRampPalette(c("#add8e6","#003366"))(50));
    breaks <- generate_breaks(nd.vec, length(color), center = centered);
    return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}
###############################################################################
#Minerva: R package (mictools: Author: albanese et al.)
##############################################
mictoolsc <- function(x, alpha=9, C=5, seed=654331, nperm=nperm, p.adjust.method=p.adjust.method){
    ## Setup
    bins <- seq(0,1, length.out = 10000 + 1)
    varnames <- colnames(x)
    if (is.null(varnames))
    varnames <- paste("Var", 1:ncol(x))

    ## Compute null distribution for TIC
    ticnull <- mictools_null(x, alpha=alpha, C=C, nperm=nperm, seed=seed)

    ## Compute histogram of the null distribution and right tailed area
    histtic <- graphics::hist(ticnull, breaks = bins, plot = FALSE, right=FALSE)
    histcumsum <- rev(cumsum(rev(histtic$counts)))

    ## Compute observed values
    ## ticmicmat <- mictools_pstats(x, alpha=alpha, C=C, est="mic_e")
    ticmicmat <- pstats(x, alpha=alpha, C=C, est="mic_e")

    # ix <- t(utils::combn(1:ncol(x),2))
    ticmicmat <- data.frame(TIC=ticmicmat[, 4], I1=ticmicmat[, 1], I2=ticmicmat[,2], 
    VarNames1=varnames[ticmicmat[, 1]], VarNames2=varnames[ticmicmat[, 2]], stringsAsFactors = FALSE);
    nulldist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))], 
     NullCount=histtic$counts, NullCumSum=histcumsum, stringsAsFactors = FALSE); 
    ## Histogram of observed values
    histtictrue <- graphics::hist(ticmicmat[,1], breaks = bins, plot=FALSE, right=FALSE)
    obscumtic <- rev(cumsum(rev(histtictrue$counts)));
    obsdist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))],
    Count=histtictrue$counts, CountCum=obscumtic, stringsAsFactors = FALSE);  
    ## One-dimensional linear interpolation
    pval <- stats::approx(bins[1:(length(bins)-1)], y = histcumsum, xout=ticmicmat[,1], method="linear", rule=2)$y / (histcumsum[1] + 1);
    pval.df <- data.frame(pval=pval, I1=ticmicmat[,2], I2=ticmicmat[,3], 
    Var1=varnames[ticmicmat[,2]], Var2=varnames[ticmicmat[,3]], 
    stringsAsFactors = FALSE);
    pval.df$adj.P.Val <- stats::p.adjust(pval.df$pval, method=p.adjust.method); 
    ## Return values as list
    retlist <- list(tic=ticnull, nulldist=nulldist, 
    obstic=ticmicmat, obsdist=obsdist, pval=pval.df); 
    return(retlist)
}
###############################################################################
mic_strength1<- function(x, pval, alpha=NULL, C=5, pthr=0.05, pval.col=NULL){
        nbins <- c(1,    25,   50,   250,   500, 1000, 2500, 5000, 10000, 40000);
        alphas <- c(0.85, 0.80, 0.75, 0.70, 0.65, 0.6,  0.55, 0.5,  0.45,  0.4);
        if (is.null(alpha))
        alpha <- alphas[.bincode(nrow(x), nbins, include.lowest=TRUE)]
        ##  Check parameters.....
        pval.col <- check_pvalcol(pval, pval.col)
        ## Select association by pvalue
        pix <- pval[, pval.col[1]] < pthr
        ## Compute mic in res[,2]
        res <- sapply(which(pix), function(y, x.sub, vartouse, alpha, C){
        xtmp <- x.sub[,vartouse[y, 1]];
        ytmp <- x.sub[,vartouse[y, 2]];
        mine_stat(xtmp, ytmp, alpha=alpha, C=C, est="mic_e", measure = "mic")
        }, x.sub=x, alpha=alpha, C=C, vartouse=pval[, pval.col[2:3]])
        ## Add variable index of the MIC computed
        res.df <- data.frame(TicePval=pval[pix, pval.col[1]], MIC=res, 
             I1=pval[pix, pval.col[2]], I2=pval[pix, pval.col[3]], 
             stringsAsFactors = FALSE)
        return(res.df);
}
###############################################################################
check_pvalcol<-function(pval, pval.col){
    if (length(pval.col) > 3)
        pval.col <- pval.col[1:3]
    if (is.character(pval.col)){
        mx <- match(pval.col, names(pval))
        pval.col <- mx[!is.na(mx)]
    }
    if ((length(pval.col) == 0) || (is.null(pval.col)))
        pval.col <- 1:3
    if (length(pval.col) == 1)
        pval.col <- seq(pval.col, pval.col+2)
    if (length(pval.col)==2)
        pval.col <- c(pval.col, pval.col[length(pval.col)] + 1)
    if (max(pval.col) > ncol(pval))
        stop(paste0("Not available column "), max(pval.col), call. = FALSE)
    return(pval.col)
}
#############################
require(ggplot2);
black.bold.10.text <- element_text(face = "bold", color = "black", size = 10);
black.bold.8.text <- element_text(face = "bold", color = "black", size = 8);
black.bold.11.title <- element_text(face = "bold", color = "black", size = 12);
bold.10.text.x <- element_text(face = "bold", size = 10);
#############################

##############Core resistome source code: Microbiome R package#################
core <- function(x, detection, prevalence, include.lowest=FALSE) {
    xorig <- x
    
    # TODO: add optional renormalization such that the core member
    # abundances would
    # sum up to 1 ?
    taxa <- core_members(x, detection, prevalence,
                       include.lowest=include.lowest)
    prune_taxa(taxa, xorig)
}
###############################################
core_members <- function(x, detection=detection, prevalence=prevalence,
                         include.lowest=FALSE) {
  # Pick taxa x samples matrix
    x <- abundances(x)
    if (include.lowest) {
        taxa <- names(which(prevalence(x, detection,
                                   include.lowest=include.lowest) >= prevalence))
    } else {
    taxa <- names(which(prevalence(x, detection,
                                   include.lowest=include.lowest) > prevalence))
    }
    taxa
}
##########################################
abundances <- function(x, transform="identity") {
  # Pick the OTU data
    if (class(x) == "phyloseq") {
    # Pick OTU matrix
        otu <- get_taxa(x)
    # Ensure that taxa are on the rows
        if(!taxa_are_rows(x) && ntaxa(x) > 1 && nsamples(x) > 1) {
            otu <- t(otu)
        }
        if(ntaxa(x) == 1) {
            otu <-matrix(otu, nrow=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }
        if(nsamples(x) == 1) {
            otu <-matrix(otu, ncol=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }
    } else if(is.vector(x)) {
        otu <- as.matrix(x, ncol=1)
    } else{
        # If x is not vector or phyloseq object then let us assume it is a
        # taxa x samples
        # count matrix
        otu <- x
    }
  # Apply the indicated transformation
  if (!transform == "identity") {
        otu <- transform(otu, transform)
  }
  otu
}
#######################################
plot_core <- function(x, prevalences=seq(, 1, 1, 0.1), detections=20,
                      plot.type="lineplot", colours=NULL, # gray(seq(0, 1, length=5)),
                      min.prevalence=NULL, taxa.order=NULL, horizontal=FALSE) {
    if (length(detections) == 1) {
     detections <- 10^seq(log10(0.001), log10(max(abundances(x),
                                                 na.rm=TRUE)), length=detections)
    }
    if (plot.type == "lineplot") {
        # Calculate the core matrix (prevalences x abundance thresholds)
        coremat <- core_matrix(x, prevalences, detections)
        res <- core_lineplot(coremat)
    } else if (plot.type == "heatmap") {
        # Here we use taxon x abundance thresholds table indicating prevalences
        res <- core_heatmap(abundances(x),
                        dets=detections, cols=colours, 
                        min.prev=min.prevalence, taxa.order=taxa.order)
    }
    p <- res$plot;
    if (horizontal) {
        p <- p + coord_flip() + theme(axis.text.x=element_text(angle=90)); 
    }
    p+theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
}
#####################################
core_heatmap <- function(x, dets, cols, min.prev, taxa.order){
    suppressMessages(require(reshape));
    data <- x
    #colours <- gray(seq(0, 1, length=5)), 
    DetectionThreshold <- Taxa <- Prevalence <- NULL
    # Prevalences with varying dets
    prev <- lapply(dets, function(th) {
        prevalence(data, detection=th)
    })
    prev <- do.call("cbind", prev)
    colnames(prev) <- as.character(dets)

    # # Exclude rows and cols that never exceed the given prevalence
    if (!is.null(min.prev)) {
        prev <- prev[rowMeans(prev > min.prev) > 0,
                   colMeans(prev > min.prev) > 0]
    }
    df <- as.data.frame(prev)
    df$ID <- rownames(prev)
    df <- melt(df, "ID")
    names(df) <- c("Features", "DetectionThreshold", "Prevalence")
    df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    if(is.null(taxa.order)) {
        o <- names(sort(rowSums(prev)))
    } else {
        o <- taxa.order
    }
    df$Features <- factor(df$Features, levels=o)
    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x=DetectionThreshold, y=Features, fill=Prevalence))
    p <- p + geom_tile()
    p <- p + xlab("Detection Threshold")
    p <- p + scale_x_log10()
    if (!is.null(cols)) {
        p <- p + scale_fill_gradientn("Prevalence",
                breaks=seq(from=0, to=1, by=0.1), 
                colours=cols,
                limits=c(0, 1))
    }
    return(list(plot=p, data=df))
}
####################################
prevalence <- function(x, detection=0, sort=FALSE, count=FALSE,
    include.lowest=FALSE) {
    if(is.null(detection)) {
        detection <- (-Inf)
    }
    if(is.null(x)) {
        warning("x is NULL - returning NULL")
        return(NULL)
    }
    # Convert phyloseq to matrix
    if(class(x) == "phyloseq") {
        x <- abundances(x)
    }
    if(is.vector(x)) { 
        if(include.lowest) {
            prev <- sum(x >= detection)
        } else {
            prev <- sum(x > detection)
        }   
    } else if(is.matrix(x) || is.data.frame(x)) {
        
        if(include.lowest) {
            prev <- rowSums(x >= detection)
        } else{
            prev <- rowSums(x > detection)
        }
    }
    if(!count) {
        prev <- prev/prevalence_nsamples(x)
    }
    if(sort) {
        prev <- rev(sort(prev))
    }  
    prev 
}
##################################################
# Internal auxiliary function
prevalence_nsamples <- function(x) {
    if(is.vector(x)) {
        n <- length(x)
    } else if(is.matrix(x) || is.data.frame(x)) {
        n <- ncol(x)
    } 
    n 
}
#################################################

################################################################################