##################################################
## R scripts for ResistoXplorer 
## Description: R code for ARG list module (ARG-microbial host network exploration)
## some codes adopted from mirnet tool
###################################################
###############################################################################
SetupGeneListData <- function(mirs, dbType, idType){
    dataSet$listData <- TRUE;
    data.org <<- dataSet$org <- dbType;
    dataSet$idType <- idType;
    current.msg <<- NULL;
    lines <- strsplit(mirs, "\r|\n|\r\n")[[1]];
    mir.lists <- strsplit(lines, "\\s+");
    mir.mat <- do.call(rbind, mir.lists);
    if(dim(mir.mat)[2] == 1){ # add 0
        mir.mat <- cbind(mir.mat, rep(0, nrow(mir.mat)));
    }else if(dim(mir.mat)[2] > 2){
        mir.mat <- mir.mat[,1:2];
        current.msg <<- "More than two columns found in the list. Only first two columns will be used.";
    }
    if(idType == "mir_id"){ 
        mir.vec <- gsub("mir", "miR", mir.vec);
    }
    mir.vec <- mir.mat[,1];
    mir.mat <- mir.mat[,-1, drop=F];
    rownames(mir.mat) <- mir.vec;
    dataSet$mir.orig <- mir.mat;
    dataSet$mir.orig.no <- nrow(mir.mat);
    dataSet <<- dataSet;
    return (nrow(mir.mat));
}
###############################################################################
LoadDbLibrary <-function(fun.type){
    db <- dataSet$org;
    db.rda<-readRDS(paste(lib.path, "db/",db, "/", tolower(db), "_", fun.type, ".rds", sep=""));
    current.setlink <- ""; #not sure about this
    ms.list <- strsplit(as.character(db.rda[,2]),", "); #use ", " as seperator 
    names(ms.list) <- db.rda[,1];
    current.mset <- ms.list;
    set.ids<- names(current.mset);
    current.setlink <<- current.setlink;
    current.setids <<- set.ids;
    current.geneset <<- current.mset;
    current.universe <<- unique(unlist(current.geneset));
}
###############################################################################
CreateResisNets <- function(net.type){
    col.nm <- "Species";
    if (dataSet$org=="AMP" || dataSet$org=="MARDy"){   
        nd.ids <- c(dataSet$mir.res[, "Gene"],dataSet$mir.res[, "Species"]);
    } else {
        nd.ids <- c(dataSet$mir.res[, "Accession"],dataSet$mir.res[, "Species"]);
    }
    dataSet$mirnet <<- net.type;
    dataSet$net.tp <<- col.nm;
    my.nodes <- dataSet$mir.res[, c("Gene",col.nm)];
    nd.nms <- c(dataSet$mir.res[, "Gene"], dataSet$mir.res[, col.nm]);
    names(nd.ids) <- nd.nms;
    dups <- duplicated(nd.ids); #note using unique will lose the names attribute
    dataSet$node.anot <- nd.ids[!dups];
    dataSet$node.anot <<-tapply(dataSet$node.anot, names(dataSet$node.anot), FUN = paste, collapse = ', ');
    library(igraph);
    mir.graph <- simplify(graph.data.frame(my.nodes, directed=FALSE));

    # add node expression value
    match.index <- match(tolower(V(mir.graph)$name), tolower(rownames(dataSet$mir.mapped)));
    
    #matched ARGs should have expression value but all the species should have same expression value
    expr.vals <- dataSet$mir.mapped[match.index, 1];
    nonmatch.index <- which(is.na(expr.vals));
    expr.vals[nonmatch.index] <- 0;
    names(expr.vals)[nonmatch.index] <- V(mir.graph)$name[nonmatch.index];
    mir.graph <- set.vertex.attribute(mir.graph, "abundance", index = V(mir.graph), value = expr.vals);
    mir.graph <<- mir.graph;
    substats <- DecomposeMirGraph(mir.graph, 2, net.type);
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        mir.query <- nrow(dataSet$mir.mapped);
        #mir.query <- nrow(dataSet$mir.orig); #original query
        mir.count <- length(unique(my.nodes[,1]));#matched mir
        tgt.count <- length(unique(my.nodes[,2]));#matched target
        return(c(mir.query, mir.count, tgt.count, ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}
# decompose to individual connected subnetworks, discard very small ones (defined by minNodeNum)
DecomposeMirGraph <- function(gObj, minNodeNum = 2, net.type){
    comps <- decompose.graph(gObj, min.vertices=minNodeNum);
    if(length(comps) == 0){
        current.msg <<- paste("No connected nodes found after this filtering!");
        return(NULL);
    }
    # first get stats
    queries <- tolower(rownames(dataSet$mir.mapped));
    net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
    for(i in 1:length(comps)){
        g <- comps[[i]];
        if(vcount(g) > 0){
            net.stats[i,] <- c(
                               vcount(g),
                               ecount(g),
                               sum(queries %in% tolower(V(g)$name))
                              );
        }
    }
    # now sort graph based on node size and add names
    ord.inx <- order(net.stats[,1], decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    comps <- comps[ord.inx];
    if(net.type == "gene2spec"){
        names(comps) <- rownames(net.stats) <- paste("ARGs-taxa", 1:length(comps), sep="");
    } else {
        names(comps) <- rownames(net.stats) <- paste("ARGs-antibiotics", 1:length(comps), sep="");
    }
    net.stats <- cbind(rownames(net.stats), net.stats);
    colnames(net.stats) <- c("Name", "Node", "Edge", "Query");
    
    # note, we report stats for all nets (at least 2 nodes);
    # but only contruct at least min node
    hit.inx <- net.stats$Node >= minNodeNum;
    comps <- comps[hit.inx];
    sub.stats <- NULL;
    json.res <- rep(list(list()), length(comps));
    i <- 0;
    for(nm in names(comps)){
        sub.stats <- c(sub.stats, vcount(comps[[nm]]));   
    }
    # now save the components
    mir.nets <<- comps;
    net.stats <<- net.stats[, -1];  # remove the first name col 
    
    # update the mir.res edge table 
    # both side of the edge must present in all.nodes
    all.nodes <- V(gObj)$name;
    hit.inx <- (dataSet$mir.res[, "Gene"] %in% all.nodes) & (dataSet$mir.res[, dataSet$net.tp] %in% all.nodes);
    dataSet$mir.filtered <<- dataSet$mir.res[hit.inx, ];
    return(sub.stats);
}
###############################################################################
ReduceEdgeDensity <- function(nd.type="all"){
    all.nms <- V(mir.graph)$name;
    edge.mat <- get.edgelist(mir.graph);
    dgrs <- degree(mir.graph);
    nodes2rm <- NULL;
    set.seed(8574);
    if(length(all.nms) > 50){
        # only get top 50 with highest density (degree)
        inx <- rank(-dgrs) < 50;
        seed.vec <- all.nms[inx];
    }else{
        seed.vec <- all.nms;
    }
    paths.list <-list();
    
    # now calculate the shortest paths only between these densely connected nodes
    for(pos in 1:length(seed.vec)){
        paths.list[[pos]] <- get.shortest.paths(mir.graph, seed.vec[pos], seed.vec[-pos])$vpath;
    }
    nds.inxs <- unique(unlist(paths.list));
    nodes2rm <- all.nms[-nds.inxs];
    
    # keep queries
    if(nd.type == "mir"){ # only apply removing to miRNA nodes
        mir.nms <- unique(edge.mat[,1]);
        nodes2rm <- nodes2rm[nodes2rm %in% mir.nms];
    }else if(nd.type=="other"){
        my.nms <- unique(edge.mat[,2]);
        nodes2rm <- nodes2rm[nodes2rm %in% my.nms];
    }else{
        #nothing to do
    }
    path.list <- NULL; gc();
    nodes2rm <- unique(nodes2rm);
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(mir.graph, 2, "gene2spec");
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,"miRNA"])), length(unique(dataSet$mir.filtered[,"Gene"])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}
###############################################################################
FilterResNet <- function(nd.type, min.dgr, min.btw){
    all.nms <- V(mir.graph)$name;
    edge.mat <- get.edgelist(mir.graph);
    dgrs <- degree(mir.graph);
    nodes2rm.dgr <- nodes2rm.btw <- NULL;
    if(nd.type == "mir"){
        mir.nms <- unique(edge.mat[,1]);
        hit.inx <- all.nms %in% mir.nms;
    }else if(nd.type=="other"){
        my.nms <- unique(edge.mat[,2]);
        hit.inx <- all.nms %in% my.nms;
    }else{ # all
        hit.inx <- rep(TRUE, length(all.nms));
    }
    if(min.dgr > 0){
        rm.inx <- dgrs <= min.dgr & hit.inx;
        nodes2rm.dgr <- V(mir.graph)$name[rm.inx];
    }
    if(min.btw > 0){
        btws <- betweenness(mir.graph);
        rm.inx <- btws <= min.btw & hit.inx;
        nodes2rm.btw <- V(mir.graph)$name[rm.inx];
    }
    nodes2rm <- unique(c(nodes2rm.dgr, nodes2rm.btw));
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(mir.graph, 2, "gene2spec");
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,"Species"])), length(unique(dataSet$mir.filtered[,"Gene"])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}
###############################################################################
FilterResNetByList <- function(ids, id.type){
    lines <- strsplit(ids, "\r|\n|\r\n")[[1]];
    lines<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", lines, perl=TRUE);
    
    # need to first convert to correct id used in the graph 
    hit.inx <- dataSet$mir.res[, id.type] %in% lines;
    nodes2rm <- unique(dataSet$mir.res[!hit.inx, "miRNA"]);
    nodes2rm <- nodes2rm[nodes2rm %in% V(mir.graph)$name];    # make sure they are in the igraph object
    mir.graph <- simplify(delete.vertices(mir.graph, nodes2rm));
    current.msg <<- paste("A total of", length(nodes2rm) , "was reduced.");
    substats <- DecomposeMirGraph(mir.graph, 2,"gene2spec");
    if(!is.null(substats)){
        mir.graph <<- mir.graph;
        return(c(nrow(dataSet$mir.orig), length(unique(dataSet$mir.filtered[,1])), length(unique(dataSet$mir.filtered[,3])), ecount(mir.graph), length(mir.nets), substats));
    }else{
        return(0);
    }
}
#######################################################################################
convertIgraph2JSON <- function(g, filenm){
    nms <- V(g)$name;
    lbls <- as.character(dataSet$node.anot[nms]);
    
    #need to chrck print(lbls)
    # setup shape (species square, gene circle)
    shapes <- rep("square", length(nms));
    
    # get edge data
    edge.mat <- get.edgelist(g);
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], type=rep("arrow", nrow(edge.mat)));
    
    # now get coords
    pos.xy <- PerformLayOut(g, "Default");
    node.btw <- as.numeric(betweenness(g));
    node.dgr <- as.numeric(degree(g));
    node.exp <- as.numeric(get.vertex.attribute(g, name="abundance", index = V(g)));
    if(vcount(g) > 1000){
        minSize = 1;
    }else if(vcount(g) > 300){
        minSize = 2;
    }else{
        minSize = 3;
    }
    # node size to betweenness values
    node.sizes <- as.numeric(rescale2NewRange(node.exp, minSize, 8));
    
    # update gene node size, shape and color
    mir.inx <- nms %in% edge.mat[,2];
    
    # slightly highlight mir node in general
    shapes[mir.inx] <- "circle";
    if(anal.type == "gene"){ #highlight ARGs if they are the query
        node.sizes[mir.inx] <- node.sizes[mir.inx] + 2;
    }else{
        node.sizes[mir.inx] <- node.sizes[mir.inx] + 0.4;
    }
    # update mir node color
    node.cols <- rep("yellow", length(node.dgr));
    node.cols[mir.inx] <- "#77579b"; # dark blue
    topo.colsw <- node.cols;
    node.cols[mir.inx] <- "#77579b"; #changed same for both the background
    topo.colsb <- node.cols;
    #get the node color based on class and mechanism of genes
    #cls.mech.data <- as.data.frame(dataSet$mir.res.orig[indx1,c("Gene","Class","Mechanism")]);
    #cls.var<-factor(cls.mech.data$Class);
    #cls.len<-length(cls.var); 
    #x.colors <- rep(col_vector,length.out=cls.len);
    # color based on expression
    bad.inx <- is.na(node.exp) | node.exp==0;
    if(!all(bad.inx)){
        exp.val <- node.exp;
        #cols <- color_scale("#78ff4d", "#FA8072", "#ebebeb");
        cols1 <- color_scale("#269b00", "#b30000", "#333333");
        node.colsb.exp <- cols1[sapply(node.exp, getIdx, min(exp.val), max(exp.val))];
        node.colsw.exp <- cols1[sapply(node.exp, getIdx, min(exp.val), max(exp.val))];
        #node.colsb.exp <- getExpColors(node.exp, c("#78ff4d", "#FA8072", "#ebebeb")); 
        #node.colsw.exp <- getExpColors(node.exp, c("#269b00", "#b30000", "#333333"));
        node.colsb.exp[bad.inx] <- "#d3d3d3"; 
        node.colsw.exp[bad.inx] <- "#c6c6c6"; 
    }else{
        node.colsb.exp <- rep("#d3d3d3",length(node.exp)); 
        node.colsw.exp <- rep("#c6c6c6",length(node.exp)); 
    }
    # now create the json object
    nodes <- vector(mode="list");
    for(i in 1:length(node.sizes)){
        nodes[[i]] <- list(
                  id=nms[i], 
                  size=node.sizes[i], 
                  type=shapes[i],
                  url=lbls[i],
                  colorb=topo.colsb[i],
                  colorw=topo.colsw[i],
                  x=pos.xy[i,1], 
                  y=pos.xy[i,2],
                  attributes=list(
                    expr = node.exp[i],
                    expcolb=node.colsb.exp[i],
                    expcolw=node.colsw.exp[i],
                    degree=node.dgr[i], # actual degree in orginal network
                    between=node.btw[i])
                );
    }
    # save node table
    nd.tbl <- data.frame(Id=nms, Label=nms, Degree=node.dgr, Betweenness=node.btw);
    fileNm <- paste("node_table_", substring(filenm, 0, nchar(filenm)-5), ".csv", sep="")
    write.csv(nd.tbl, file=fileNm, row.names=FALSE);
    
    # covert to json
    library(RJSONIO);
    netData <- list(resnet=dataSet$mirnet, restarget=dataSet$mirtarget, organism=dataSet$org, nodes=nodes, edges=edge.mat);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    
    # also save to GraphML
    write.graph(g, file="resistonet.graphml", format="graphml");
}

##' @importFrom grDevices colorRampPalette
color_scale <- function(c1="grey", c2="red", c3="green") {
    pal <- colorRampPalette(c(c1, c2, c2))
    colors <- pal(100)
    return(colors)
}
###############################################################################
# col vec is for low high null
getExpColors <- function(nd.vec, col.vec){
    nvec <- rep("", length(nd.vec));
    m.inx <- is.null(nd.vec) | is.na(nd.vec);
    nvec[m.inx] <- col.vec[3];
    pos.inx <- nd.vec > 0;
    nvec[pos.inx] <- col.vec[2];
    nvec[!pos.inx] <- col.vec[1];
    as.character(nvec);
}
getIdx <- function(v, MIN, MAX) {
    if ( MIN == MAX ) {
        return(100)
    }
    intervals <- seq(MIN, MAX, length.out=100)
    max(which(intervals <= v))
}
###############################################################################
PrepareGraphML <- function(net.nm){
    write.graph(mir.nets[[net.nm]], file=paste(net.nm, ".graphml", sep=""), format="graphml");
}
###############################################################################
PrepareMirNet <- function(mir.nm, file.nm){
    my.mirnet <- mir.nets[[mir.nm]];
    current.mirnet <<- my.mirnet;
    convertIgraph2JSON(my.mirnet, file.nm);
    return(1);
}
###############################################################################
PerformLayOut <- function(g, algo){
    vc <- vcount(g);
    if(algo == "Default"){
        if(vc > 3000) {
          pos.xy <- layout.lgl(g, maxiter = 100);
        }else if(vc > 2000) {
          pos.xy <- layout.lgl(g, maxiter = 150);
        }else if(vc > 1000) {
          pos.xy <- layout.lgl(g, maxiter = 200);
        }else if(vc < 150){
          pos.xy <- layout.kamada.kawai(g);
        }else{
          pos.xy <- layout.fruchterman.reingold(g);
        }
    }else if(algo == "FrR"){
        pos.xy <- layout.fruchterman.reingold(g);
    }else if(algo == "random"){
        pos.xy <- layout.random(g);
    }else if(algo == "lgl"){
        if(vc > 3000) {
          pos.xy <- layout.lgl(g, maxiter = 100);
        }else if(vc > 2000) {
          pos.xy <- layout.lgl(g, maxiter = 150);
        }else {
          pos.xy <- layout.lgl(g, maxiter = 200);
        }
    }else if(algo == "gopt"){
        # this is a slow one
        if(vc > 3000) {
          maxiter = 50;
        }else if(vc > 2000) {
          maxiter = 100;
        }else if(vc > 1000) {
          maxiter = 200;
        }else{
          maxiter = 500;
        }
        pos.xy <- layout.graphopt(g, niter=maxiter);
    }
    pos.xy;
}
###############################################################################
UpdateNetworkLayout <- function(algo, filenm){
    pos.xy <- PerformLayOut(current.mirnet, algo);
    nms <- V(current.mirnet)$name;
    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id=nms[i], 
                  x=pos.xy[i,1], 
                  y=pos.xy[i,2]
        );
    }
    # now only save the node pos to json
    library(RJSONIO);
    netData <- list(nodes=nodes);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    return(filenm);
}
###############################################################################
GetNetNames <- function(){
    rownames(net.stats);
}
###############################################################################
GetNetStats <- function(){
    as.matrix(net.stats);
}
###############################################################################
GetNetsNameString <- function(){
    paste(rownames(net.stats), collapse="||");
}
###############################################################################
UpdateSubnetStats <- function(){
    old.nms <- names(mir.nets);
    net.stats <- ComputeSubnetStats(mir.nets);
    ord.inx <- order(net.stats[,1], decreasing=TRUE);
    net.stats <- net.stats[ord.inx,];
    rownames(net.stats) <- old.nms[ord.inx];
    net.stats <<- net.stats;
}
###############################################################################
ComputeSubnetStats <- function(comps){
    net.stats <- as.data.frame(matrix(0, ncol = 3, nrow = length(comps)));
    colnames(net.stats) <- c("Node", "Edge", "Query");
    queries <- rownames(dataSet$mir.mapped);
    for(i in 1:length(comps)){
        g <- comps[[i]];
        net.stats[i,] <- c(vcount(g),ecount(g),sum(queries %in% V(g)$name));
    }
    return(net.stats);
}
###############################################################################
# from to should be valid nodeIDs
GetShortestPaths <- function(from, to){
    paths <- get.all.shortest.paths(current.mirnet, from, to)$res;
    if(length(paths) == 0){
        return (paste("No connection between the two nodes!"));
    }
    path.vec <- vector(mode="character", length=length(paths));
    for(i in 1:length(paths)){
        path.inx <- paths[[i]]; 
        path.ids <- V(current.mirnet)$name[path.inx];
        
        #path.sybls <- V(current.mirnet)$Label[path.inx];
        path.sybls <- path.ids;
        pids <- paste(path.ids, collapse="->");
        psbls <- paste(path.sybls, collapse="->");
        path.vec[i] <- paste(c(pids, psbls), collapse=";")
    }
    if(length(path.vec) > 50){
        path.vec <- path.vec[1:50];
    }
    all.paths <- paste(path.vec, collapse="||");
    return(all.paths);
}
###############################################################################
ExtractMirNetModule <- function(nodeids){
    set.seed(8574);
    nodes <- strsplit(nodeids, ";")[[1]];
    g <- current.mirnet;
    
    # try to see if the nodes themselves are already connected
    hit.inx <- V(g)$name %in% nodes; 
    gObj <- induced.subgraph(g, V(g)$name[hit.inx]);

    # now find connected components
    comps <-decompose.graph(gObj, min.vertices=1);
    if(length(comps) == 1){ # nodes are all connected
        g <- comps[[1]];
    }else{
        # extract modules
        paths.list <-list();
        sd.len <- length(nodes);
        for(pos in 1:sd.len){
            paths.list[[pos]] <- get.shortest.paths(g, nodes[pos], nodes[-(1:pos)])$vpath;
        }
        nds.inxs <- unique(unlist(paths.list));
        nodes2rm <- V(g)$name[-nds.inxs];
        g <- simplify(delete.vertices(g, nodes2rm));
    }
    nodeList <- get.data.frame(g, "vertices");
    if(nrow(nodeList) < 3){
        return ("NA");
    }
    module.count <- module.count + 1;
    module.nm <- paste("module", module.count, sep="");
    colnames(nodeList) <- c("Id", "Label");
    ndFileNm = "mirnet_node_list.csv";
    write.csv(nodeList, file=ndFileNm, row.names=F, quote=F);
    edgeList <- get.data.frame(g, "edges");
    edgeList <- cbind(rownames(edgeList), edgeList);
    colnames(edgeList) <- c("Id", "Source", "Target");
    edgFileNm = "mirnet_edge_list.csv";
    write.csv(edgeList, file=edgFileNm, row.names=F, quote=F);
    filenm <- paste(module.nm, ".json", sep="");
    convertIgraph2JSON(g, filenm);
    
    # record the module 
    mir.nets[[module.nm]] <<- g;
    UpdateSubnetStats();
    module.count <<- module.count
    return (filenm);
}
###############################################################################
# exclude nodes in current.net (networkview)
ExcludeNodes <- function(nodeids, filenm){
    nodes2rm <- strsplit(nodeids, ";")[[1]];
    current.mirnet <- delete.vertices(current.mirnet, nodes2rm);
    
    # need to remove all orphan nodes
    bad.vs<-V(current.mirnet)$name[degree(current.mirnet) == 0];
    current.mirnet <<- delete.vertices(current.mirnet, bad.vs);
    
    # return all those nodes that are removed 
    nds2rm <- paste(c(bad.vs, nodes2rm), collapse="||");
    
    # update topo measures
    node.btw <- as.numeric(betweenness(current.mirnet));
    node.dgr <- as.numeric(degree(current.mirnet));
    node.exp <- as.numeric(get.vertex.attribute(current.mirnet, name="abundance", index = V(current.mirnet)));
    nms <- V(current.mirnet)$name;
    nodes <- vector(mode="list");
    for(i in 1:length(nms)){
        nodes[[i]] <- list(
                  id = nms[i], 
                  expr = node.exp[i],
                  degree = node.dgr[i], 
                  between = node.btw[i],
                  expr = node.exp[i]
                );
    }
    # now only save the node pos to json
    library(RJSONIO);
    netData <- list(deletes=nds2rm,nodes=nodes);
    sink(filenm);
    cat(toJSON(netData));
    sink();
    return(filenm);
}
###############################################################################
PerformMirTargetEnrichAnalysis <- function(fun.type, file.nm, IDs, algo){
    perm.num <- 1000;   
    # prepare lib
    LoadDbLibrary(fun.type);
    mirnet.type <- dataSet$mirnet;
    nodeList <- get.data.frame(current.mirnet, "vertices");
    if(tolower(fun.type) == 'mirfamily'){
        print("Mapping miRNA family from nodes.")
        hit.inx <- dataSet$mir.filtered$miRNA %in% nodeList[,1];
        mir.query <- unique(dataSet$mir.filtered$miRNA[hit.inx]);
        my.data <- unique(dataSet$mir.filtered[hit.inx,c("Accession", "miRNA")]); # The original dataset contains miRNA Accession number, you can consider it as entrez id.
        ora.vec <- my.data[hit.inx,"Accession"];
        sybls <- my.data[hit.inx,"miRNA"];
        names(ora.vec) <- sybls;
    } else{
        #ora.vec<-unique(dataSet$mir.filtered$Gene);
        hit.inx <- dataSet$mir.res$Gene %in% nodeList[,1];
        mir.query <- unique(dataSet$mir.res$Gene[hit.inx]);
        #my.data <- unique(dataSet$mir.filtered[hit.inx,c("Entrez", "Gene")]);
        ora.vec <- mir.query;
        #sybls <- my.data[hit.inx,"Gene"];
        names(ora.vec) <- mir.query;
    }
    #q.vec <-  unlist(strsplit(IDs, "; "));
    #ora.vec <- ora.vec[q.vec];
    ora.vec <- ora.vec[!is.na(ora.vec)];
    ora.nms <- names(ora.vec);

    # prepare for the result table
    set.size<-length(current.geneset);
    res.mat<-matrix(0, nrow=set.size, ncol=4);
    colnames(res.mat)<-c("Total", "Expected", "Hits", "Pval");
    rownames(res.mat)<-names(current.geneset);

    # not all query genes can be used, need to cut query to only the universe covered 
    hits.inx <- tolower(ora.vec) %in% tolower(current.universe);
    ora.vec <- ora.vec[hits.inx];
    ora.nms <- ora.nms[hits.inx];
    q.size<-length(ora.vec);
    # get the matched query for each pathway
    hits.query <- lapply(current.geneset, 
        function(x) {
            ora.nms[ora.vec %in% x];
        }
    );
    hit.num<-unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
    names(hits.query) <- names(current.geneset);

    # total unique gene number
    uniq.count <- length(current.universe);
    
    # unique gene count in each pathway
    set.size <- unlist(lapply(current.geneset, length));
    res.mat[,1] <-set.size;
    res.mat[,2] <-q.size*(set.size/uniq.count);
    res.mat[,3] <-hit.num;
    # standard hypergeometric tests use lower.tail = F for P(X>x)
    raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
    #res.mat[,4]<- raw.pvals;
    fdr.pvals <- p.adjust(raw.pvals, "fdr");
    res.mat[,4] <- fdr.pvals;
    res.mat <- res.mat[hit.num>0,,drop = F];
    hits.query <- hits.query[hit.num>0];
    if(nrow(res.mat)> 1){
        # order by p value
        ord.inx <-order(res.mat[,4]);
        res.mat <- signif(res.mat[ord.inx,],3);
        hits.query <- hits.query[ord.inx];

        # return all the hits disgarding p value cutoff
        #  imp.inx <- res.mat[,4] <= 0.05;        
        #  if(sum(imp.inx) < 10){ # too little left, give the top ones
        #      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
        #      res.mat <- res.mat[1:topn,];
        #      hits.query <- hits.query[1:topn];
        #  }else{
        #      res.mat <- res.mat[imp.inx,];
        #      hits.query <- hits.query[imp.inx];
        #      if(sum(imp.inx) > 120){
        #          # now, clean up result, synchronize with hit.query
        #          res.mat <- res.mat[1:120,];
        #          hits.query <- hits.query[1:120];
        #      }
        #  }
    }
    #get gene symbols
    resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
    if(nrow(resTable) == 0){
        current.msg <<- "No hits found for your query!";
        print(current.msg);
        return(0);
    }
    current.msg <<- "Functional enrichment analysis was completed";
    
    # write json
    fun.anot <- hits.query; 
    hit.num = resTable[,4]; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
    fun.ids <- as.vector(current.setids[names(fun.anot)]); 
    pval <- resTable[,5]; if(length(pval) ==1) { pval <- matrix(pval) };
    if(algo == "emp"){
        hit.inx <- pval == 0;
        pval[hit.inx] <- paste("<", 1/perm.num);
    }
    json.res <- list(
                    fun.anot = fun.anot,
                    fun.ids = fun.ids,
                    pval = pval,
                    hit.num = hit.num
        );
    json.mat <- toJSON(json.res, .na='null');
    json.nm <- paste(file.nm, ".json", sep="");    
    sink(json.nm)
    cat(json.mat);
    sink();

    # write csv
    # csv.nm <- paste(file.nm, ".csv", sep="");
    write.csv(resTable, file="resistonet_enrichment.csv", row.names=F);
    return(1);
}
###############################################################################
###############################################################################
CalculateMirTargetSet <- function(nms, operation){
    nms <- strsplit(nms, ";")[[1]];
    valid.inx <- nms %in% dataSet$mir.filtered$miRNA;
    nms <- nms[valid.inx];
    if(length(nms) == 0){
        print("No valid ARG or microbe (taxa) found!");
        return("error||No valid ARG or microbe (taxa) were selected!");
    }
    hit.inx <- dataSet$mir.filtered$miRNA == nms[1];
    targets <- dataSet$mir.filtered[hit.inx, 3]; # targets at third column
    for(i in 2:length(nms)){
        hit.inx <- dataSet$mir.filtered$miRNA == nms[i];
        if(operation == "intersect"){
            targets <- intersect(targets, dataSet$mir.filtered[hit.inx,3]);
        }else if(operation == "union"){
            targets <- union(targets, dataSet$mir.filtered[hit.inx,3]);
        }
    }
    # include original queries
    return(paste(unique(c(nms, targets)), collapse="||"));
}
###############################################################################
###############################################################################
GetAnotNames <- function(){
    return(rownames(dataSet$data.anot));
}
###############################################################################
GetSeqnm <- function(rowid){
    inx <- which(rownames(dataSet$mir.res) == rowid);
    mirid <- dataSet$mir.res[inx, "miRNA"];
    return(mirid);
}
###############################################################################
GetSeq <- function(mir.id){
    inx <- which(rownames(dataSet$mir.res) == mir.id);
    seq <- dataSet$mir.res[inx, "Accession"];
    return(seq);
}
###############################################################################
RemoveMirEntry <- function(mir.id) {
    inx <- which(rownames(dataSet$mir.res) == mir.id);
    if(inx < 1){
        return(0);
    } else{
        dataSet$mir.res <<- dataSet$mir.res[-inx,];
        return(1);
    }
}
###############################################################################
RemoveMirNode <- function(node.id) {
    # node ID is name of either gene or miRNA
    row.ids <- NULL;
    # first try mir names
    hits <- dataSet$mir.res[, 1] == node.id;
    if(sum(hits) > 0){
        row.ids <- rownames(dataSet$mir.res)[hits];
    }else{
        hits <- dataSet$mir.res[, 3] == node.id;
        if(sum(hits) > 0){
            row.ids <- rownames(dataSet$mir.res)[hits];
        }
    }
    if(is.null(row.ids)){
        return("NA");
    }else{
        dataSet$mir.res <<- dataSet$mir.res[!hits,];
        return(row.ids);
    }
}
###############################################################################
GetUniqueSourceNames <- function(argType){
    db.path <- paste("../../data/browse.json");
    library("RJSONIO");
    browse <- fromJSON(db.path);
    browse <<- browse;
    res <- sort(names(browse[[argType]]));
    return(res);
}
###############################################################################
GetUniqueClassNames <- function(argType, source){
    res <- sort(names(browse[[argType]][[source]]));
    return(res);
}
###############################################################################
GetUniqueSpeciesNames <- function(argType){
    res <- unique(unlist(browse[[argType]][[1]]));
    return(res);
}
###############################################################################
GetUpdateClassNames <- function(argType, source){
    res <- sort(names(browse[[argType]][[source]]));
    return(res);
}
###############################################################################
SetupFilterList <- function(exp, miranda, tarpmir){
   dataSet$parameter <- c(exp, miranda, tarpmir);
   dataSet <<- dataSet;
}
###############################################################################
PerformGeneMapping <- function(){
    dbType <- dataSet$org;
    idType <- dataSet$idType;
    mir.mat <- dataSet$mir.orig;
    idVec <- rownames(mir.mat);
    source.vec <- "";
    #read database from library
    db<-readRDS(paste(lib.path, "db/",dbType, "/", tolower(dbType), "_genelist.rds", sep=""));
    #default name is Name and Organism
    ind.gene <-  which (colnames(db)=="Name");
    ind.spec <-  which (colnames(db)=="Organism");
    colnames(db)[c(ind.gene,ind.spec)] <- c("Gene", "Species");
    GeneSpeciesMapping(idType,idVec,db);
    #no of matches
    hit.num <-length(which(!is.na(dataSet$hit.inx)));
    if(hit.num == 0){
        current.msg <<- "No hits found in the database. Please check your input.";
        return(0);
    } else{
        res <- db[na.omit(dataSet$hit.inx),];
        #the column names are variable between databases
        if(dbType =="AMP"){
            varcol.nm <- "Phenotype";                
        } else if (dbType =="CARD") {
            varcol.nm <- "Family";
        }else {
            varcol.nm <- "Class";
        }
        res <- res[,c("Gene","Species","Accession",varcol.nm,"Mechanism","Reference")];    
        hit.num <- length(unique(res$Gene));
        
        # Species can be NA, which have to be removed
        res <- delete.na(res);
        arg.target.no <- length(res$Species);
        arg.target.no.uniq <- length(unique(res$Species));            
    } 
    dataSet$mir.mapped <- mir.mat;
    write.csv(res, file="arg_species_target.csv", row.names=FALSE);    
    dataSet$mir.res <- res;
    dataSet$mirtarget <- "gene";
    dataSet$arg.target.no <- arg.target.no;
    dataSet$hit.num <- hit.num;
    dataSet$arg.target.no.uniq <- arg.target.no.uniq;
    dataSet <<- dataSet; 
    return(1);
}
###############################################################################
GeneSpeciesMapping <- function(idType,idVec,db){
    qvec <- idVec;
    # variables to record results
    hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
    match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
    match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
    if(idType == "gene"){
        hit.inx <-which(db$Gene %in% qvec==TRUE);
        #hit.inx <- match(qvec,db$Gene);
        match.values <- db$Gene[hit.inx];
        match.state[!is.na(hit.inx)] <- 1;
    }
    # empty memory
    dataSet$hit.inx <- hit.inx;
    dataSet$hit.values <- match.values;
    dataSet$match.state <- match.state;
    dataSet <<- dataSet;
}
###############################################################################
GetInputGeneListStat <- function(){
    return(c(dataSet$mir.orig.no,dataSet$hit.num,dataSet$arg.target.no,dataSet$arg.target.no.uniq));
}
###############################################################################
GetMirResCol <- function(colInx){
    res <- dataSet$mir.res[,colInx];
    return(res); 
}
###############################################################################
GetMirResRowNames <- function(){
    rownames(dataSet$mir.res);
}
###############################################################################
delete.na <- function(DF, n=0) {
    DF[rowSums(is.na(DF)) <= n,]
}
###############################################################################
################################################################################