########################################################################
## R script for ResistoXplorer
## Description: Data/resource management functions
########################################################################
#######Composition plot#######
##############################
PlotStackedAreaComposition <- function(barplotName, viewOpt, genelvl, metadata, imgOpt, colpalopt, calcmeth, feat_cnt, format="png", dpi=72){
    suppressMessages(require(reshape));
    
  #using filtered data
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows =TRUE);
    }
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj),sample_table);
    if(viewOpt=="smpl_grp"){
        data <- as.data.frame(t(otu_table(data1)));
        data <- cbind(data,variable=data1@sam_data[[metadata]]);
        data <- aggregate(. ~variable,data,sum);
        gp_nm <- rownames(data)<-data[,1];
        data <- data[,-1];
        data <- data[order(rownames(data)),]; 
        clsLbl <- sort(unique(factor(data1@sam_data[[metadata]])));
        colvec <- NULL;
    }else {
        smpl_nm <- sample_names(data1); 
        data <- data.frame(otu_table(data1));

        # reorder data based on groups
        sam <- sample_data(data1);
        clsLbl <- factor(sam[[metadata]]);
        ord.inx <- order(clsLbl);
        smpl_nm <- smpl_nm[ord.inx];
        clsLbl <- clsLbl[ord.inx];
        colvec <- as.numeric(clsLbl)+1;
        data <- t(data[,ord.inx]);
    }
    data_tax <- tax_table(data1);       
    
    #reshaping data
    taxa_nm <- as.data.frame(data_tax[,genelvl]);
    taxa_nm <- as.matrix(taxa_nm);
    y <- which(is.na(taxa_nm)==TRUE);

    #converting NA values to unassigned; before order it to last position using ZZZ as its name
    taxa_nm[y] <- "ZZZ";   
    colnames(data) <- taxa_nm[,1];
    nms <- colnames(data); 
    data <- as.matrix(data); 
    data <- data %*% sapply(unique(nms),"==",nms); 
    data <- data.frame(data);
    data <- data[ ,order(names(data))];
    indx <- which(colnames(data)=="ZZZ");
    colnames(data)[indx] <- "NA";
    if(calcmeth=="sum"){
        ind <- which(colSums(data)>feat_cnt);
        ind1 <- which(colSums(data)<feat_cnt);
    } else {
        dt <- apply(data,2,median);
        ind <- which(dt>feat_cnt);
        ind1 <- which(dt<feat_cnt);
    }
    if(length(ind)==0){
        current.msg <<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    }
    if(length(ind1)>0){
        colnames(data)[ind1] <-"Others";
        data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
    }
    feat_no <- ncol(data);
    
    #adjust height according to number of legends
    h <- 540;
    if(feat_no < 10){
        h <- h;
    }else if(feat_no < 20){
        h <- h+100;
    }else if (feat_no < 50){
        h <- h+200;
    }else if (feat_no < 100){
        h <- h+400;
    }else if (feat_no > 100){
        h <- h+500;
    } 
    # width calculation
    a <- nsamples(data1);
    min.w <- 540;
    if(a<50){
        w <- a*38+50;
    }else{
        w <- a*20+50;
    }
    if(w <min.w){
        w <- min.w;
    }
    write.csv(t(data), file="feature_abundance.csv");
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step); 
    data <- data[order(data[,2]),];
    data$variable <- gsub("\\.", " ",data$variable);
    data <- data[,-1];
    if(viewOpt=="smpl_grp"){
        data$step <- rep(1:length(gp_nm),feat_no);
        lbl <- gp_nm;
    }else {
        data$step <- rep(1:a,feat_no);
        lbl <- smpl_nm;
    }
    #color schema
    color_var <- levels(factor(data$variable));
    x <- length(color_var); 
    if(colpalopt=="grad"){
        indx <- which(color_var=="NA");
       
        #color schema for ggplot
        x.colors <- hcl(h=seq(15,375,length=(x+1)),l=65,c=100)[1:x];
        x.colors[indx] <- "#666666";
    }else if(colpalopt=="cont21"){
        x.colors <- rep(custom_col21,length.out=x);
    }else if(colpalopt=="cont28"){
        x.colors <- rep(custom_final28,length.out=x);
    }else if(colpalopt=="cont74"){
        x.colors <- rep(custom_col74,length.out=x);
    }else{
        x.colors<-rep(custom_col42,length.out=x);
    }
    imgSet$stack <- barplotName;
    imgSet <<- imgSet;
    Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
    box <- ggplot(data,aes(x=step,y=value))+theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust =1,vjust= 1))+
        geom_area(aes(fill=variable),position='fill')+
        scale_x_continuous(breaks=seq(1,length(unique(data$step)),1),labels=lbl)+
        labs(y="Relative abundance",x="",fill=genelvl)+
        scale_fill_manual(values=c(x.colors)) + 
        theme(legend.position="bottom",legend.box = "vertical")+
        theme(axis.text.x = element_text(colour=colvec),axis.title.x=element_blank());
    box <- box + theme(axis.text.x = bold.10.text.x,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);    
    print(box);
    dev.off();
    analSet$stack<-data;
    analSet$stack.genelvl<-genelvl;
    analSet$plot <- "Stacked Area";
    analSet <<- analSet;
    return(1); 
}
###############################################################################
PlotStackedBar <- function(barplotName,taxalevel,metadata,imgOpt,colpalopt,calcmeth,feat_cnt,format="png",dpi=72){
    suppressMessages(require(reshape));
    
    #using filtered data
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows =TRUE);
    }
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj),sample_table);
    data <- as.data.frame(otu_table(data1));
    data_tax <- tax_table(data1);
 
    # reorder data based on groups
    sam <- sample_data(data1);
    smpl_nm <- sample_names(data1); 
    clsLbl <- factor(sam[[metadata]]);
    ord.inx <- order(clsLbl);
    smpl_nm <- smpl_nm[ord.inx];

    clsLbl <- clsLbl[ord.inx];
    colvec <- as.numeric(clsLbl)+1;
    data <- t(data[,ord.inx]);
    
    #reshaping data
    taxa_nm <- as.data.frame(data_tax[,taxalevel]);
    taxa_nm <- as.matrix(taxa_nm);
    y <- which(is.na(taxa_nm)==TRUE);
   
    #converting NA values to unassigned; before order it to last position using ZZZ as its name
    taxa_nm[y] <- "ZZZ";  
    colnames(data) <- taxa_nm[,1];
    nms <- colnames(data);
    data <- as.matrix(data);
    data <- data %*% sapply(unique(nms),"==",nms); 
    data <- as.data.frame(data);
    data <- data[ , order(names(data))];
    indx <- which(colnames(data)=="ZZZ");
    colnames(data)[indx] <-"NA";
    if(calcmeth=="sum"){
        ind <- which(colSums(data)>feat_cnt);
        ind1 <- which(colSums(data)<feat_cnt);
    } else {
        dt <- apply(data,2,median);
        ind <- which(dt>feat_cnt);
        ind1 <- which(dt<feat_cnt);
    }
    if(length(ind)==0){
        current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    }
    if(length(ind1)>0){
        colnames(data)[ind1] <- "Others";
        data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
    }
    yLbl <- "Actual Abundance";
    if(imgOpt=="barnorm"){
        data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
        data <- as.data.frame(t(data));
        yLbl <- "Relative Abundance";
    }
    # height according to number of legends
    feat_no <- ncol(data);
    h <- 540;
    if(feat_no < 10){
        h <- h;
    } else if(feat_no<20){
        h <- h+100;
    } else if (feat_no<50){
        h <- h+200;
    } else if(feat_no<100){
        h <- h+400;
    } else if(feat_no>100){
        h <- h+500;
    } 
    # width calculation
    a <- nsamples(data1);
    min.w <- 480;
    if(a<50){
        w <- a*38+50;
    }else{
        w <- a*21.5+50;
    }
    if(w<min.w){
        w <- min.w;
    }
    write.csv(t(data), file="feature_abundance.csv");
    data$step <- factor(rownames(data));
    data$step <- rep(1:length(data$step));
    data <- melt(data,id='step');
    data <- data[order(data[,2]),];
    data$variable <- gsub("\\.", " ",data$variable);
    
    #color schema
    color_var <- levels(factor(data$variable));
    x <- length(color_var); 
    if(colpalopt=="grad"){
        indx <- which(color_var=="NA");
        
        #color schema for ggplot
        x.colors <- hcl(h=seq(15,375, length=(x+1)),l=65,c=100)[1:x];
        x.colors[indx] <-"#666666";
    }else if(colpalopt=="cont21"){
        x.colors <- rep(custom_col21, length.out=x);
    }else if(colpalopt=="cont28"){
        x.colors <- rep(custom_final28, length.out=x);
    }else if(colpalopt=="cont74"){
        x.colors <- rep(custom_col74, length.out=x);
    }else {
        x.colors <- rep(custom_col42,length.out=x);
    }
    imgSet$stack <- barplotName;
    imgSet <<- imgSet;
    Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi); 
    box <- ggplot(data,aes(x = step, y = value, fill = variable))+ geom_bar(stat="identity", position="stack")+theme_bw()+
            theme(legend.position="bottom")+scale_fill_manual(values=c(x.colors))+
            theme(axis.text.x = element_text(angle = 45, hjust =1,vjust= 1))+labs(y=yLbl,x="")+
            scale_x_discrete(limits=c(smpl_nm)) +labs(y=yLbl,x="",fill=taxalevel)+
            theme(axis.text.x = element_text(colour=colvec),axis.title.x=element_blank(),panel.background = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank());
    box <- box + theme(legend.box = "vertical",axis.text.x = bold.10.text.x,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);        
    print(box);
    dev.off();
    analSet$stack.genelvl <- taxalevel;
    analSet$plot <- "Stacked Bar";
    analSet$stack <- data;
    analSet <<- analSet;
    return(1);
}
###############################################################################
PlotStackedBarSamGrp <- function(barplotName, genelvl, metadata, imgOpt, colpalopt, calcmeth, feat_cnt, format="png", dpi=72){
    suppressMessages(require(reshape));
    
    #using filtered data
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows=TRUE);
    }
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj), sample_table);
    yLbl <- "Actual Abundance";
    data <- as.data.frame(t(otu_table(data1)));
    data <- cbind(data,variable=data1@sam_data[[metadata]]);
    data <- aggregate(. ~variable,data,sum);
    gp_nm <- rownames(data)<-data[,1];
    data <- data[,-1];
    data <-data[order(rownames(data)),];
    clsLbl <- sort(unique(factor(data1@sam_data[[metadata]])));
    data_tax <-tax_table(data1);       
    
    #reshaping data
    taxa_nm <- as.data.frame(data_tax[,genelvl]);
    taxa_nm <- as.matrix(taxa_nm);
    y <- which(is.na(taxa_nm)==TRUE);
   
    #converting NA values to unassigned; before order it to last position using ZZZ as its name
    taxa_nm[y] <-"ZZZ";   
    colnames(data) <-taxa_nm[,1];
    nms <- colnames(data);
    data <- as.matrix(data);
    data <- data %*% sapply(unique(nms),"==",nms); 
    data <- as.data.frame(data);
    data <- data[ , order(names(data))];
    indx <- which(colnames(data)=="ZZZ");
    colnames(data)[indx] <- "NA";
    if(calcmeth=="sum"){
        ind <- which(colSums(data)>feat_cnt);
        ind1 <- which(colSums(data)<feat_cnt);
    } else{
        dt <- apply(data,2, median);
        ind <- which(dt>feat_cnt);
        ind1 <- which(dt<feat_cnt);
    }
    if(length(ind)==0){
        current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    }
    if(length(ind1)>0){
        colnames(data)[ind1] <- "Others";
        data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
    }
    if(imgOpt=="barnorm"){
        data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
        data <- as.data.frame(t(data));
        yLbl <- "Relative Abundance";
    }
    feat_no <- ncol(data);
   
    #adjust height according to number of legends
    #h<-300;
    h <- length(clsLbl)*80;
    w <- 960;
    if(feat_no<10){
        h <- h+40;
    } else if (feat_no<20){
        h <- h+60;
    } else if(feat_no<50){
        h <- h+120;
    } else if(feat_no<100){
        h <- h+250;
    } else if(feat_no<200){
        h <- h+500;
    } else if(feat_no>200){
        h <- h+600;
    } 
    write.csv(t(data), file="feature_abundance.csv");
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step); 
    data <- data[order(data[,2]),];
    data <- data[,-1];
    data$variable <- gsub("\\.", " ",data$variable);
    data$step <- rep(gp_nm,feat_no);
    
   #color schema
    color_var <- levels(factor(data$variable));
    x <- length(color_var); 
     if(colpalopt=="grad"){
        indx <- which(color_var=="NA");
        #color schema for ggplot
        x.colors <- hcl(h=seq(15,375,length=(x+1)),l=65,c=100)[1:x];
        x.colors[indx] <- "#666666";
    }else if(colpalopt=="cont21"){
        x.colors <- rep(custom_col21,length.out=x);
    }else if(colpalopt=="cont28"){
        x.colors<- rep(custom_final28,length.out=x);
    }else if(colpalopt=="cont74"){
        x.colors <- rep(custom_col74,length.out=x);
    }else{
        x.colors <- rep(custom_col42,length.out=x);
    }   
    imgSet$stack <- barplotName;
    imgSet <<- imgSet;
    Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
    box <- ggplot(data,aes(x = step, y = value, fill = variable))+ geom_bar(stat="identity", position="stack", width = 0.4)+theme_bw()+ 
            theme(legend.position="bottom")+labs(y=yLbl,x="",fill=genelvl)+scale_fill_manual(values=c(x.colors))+
            theme(axis.text.x = element_text(angle = 0,vjust=0.5))+coord_flip()+
            theme(panel.background = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank());
    box <- box + theme(legend.box = "vertical",axis.text.x = bold.10.text.x,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);        
    print(box);
    dev.off();
    analSet$stack <- data;
    analSet$stack.genelvl <- genelvl;
    analSet$plot <- "Stacked Bar";
    analSet <<- analSet;
    return(1); 
}
###############################################################################
PlotStackedBarSam <- function(barplotName, genelvl, samplnm, imgOpt, feat_cnt, format="png", dpi=72){
    suppressMessages(require(reshape));
    
    #using filtered data
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows=TRUE);
    }
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj), sample_table);
    yLbl <- "Actual Abundance";
    data <- as.matrix(otu_table(data1));
    data <- data[, samplnm, drop=FALSE];
    data <- t(data);
    data_tax <- tax_table(data1);       
    
    #reshaping data
    taxa_nm <- as.data.frame(data_tax[,genelvl]);
    taxa_nm <- as.matrix(taxa_nm);
    y <- which(is.na(taxa_nm)==TRUE);
   
    #converting NA values to unassigned; before order it to last position using ZZZ as its name
    taxa_nm[y] <- "ZZZ"; 
    colnames(data) <- taxa_nm[,1];
    nms <- colnames(data);
    data <- data %*% sapply(unique(nms),"==",nms); 
    data <- as.data.frame(data);
    data <- data[ , order(names(data))];
    indx <- which(colnames(data)=="ZZZ");
    colnames(data)[indx] <-"NA";
    ind <- which(colSums(data)>feat_cnt);
    ind1<- which(colSums(data)<feat_cnt);
    if(length(ind)==0){
        current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    }
    if(length(ind1)>0){
        colnames(data)[ind1] <- "Others";
        data <- as.data.frame(t(rowsum(t(data),group = colnames(data))));
    }
    if(imgOpt=="barnorm"){
        data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
        data <- as.data.frame(t(data));
        yLbl <- "Relative Abundance";
    }
    feat_no <- ncol(data);
    
    #adjust height according to number of legends
    w <- 600; 
    write.csv(t(data), file="feature_abundance.csv");
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step); 
    data <- data[order(data[,2]),];
    data$variable <- gsub("\\.", " ",data$variable);
    data <- data[,-1];
    a <- feat_no;
    if(length(a) < 50){
        h <- feat_no*30;
    }else{
        h <- feat_no*22;
    }
    #sorting by descending
    imgSet$stack <- barplotName;
    imgSet <<- imgSet;
    Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
    box <- ggplot(data, aes(x=reorder(variable,value),y=value))+geom_bar(stat="identity",width=0.6,fill="steelblue")+theme_bw()+
        theme(axis.text.x = element_text(angle = 0,vjust=0.5))+
        labs(y=yLbl,x="",fill=genelvl)+coord_flip()+
        theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none");
    box <- box + theme(axis.text.x = bold.10.text.x,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    analSet$stack <- data;
    analSet$stack.genelvl <- genelvl;
    analSet$plot <- "Stacked Bar";
    analSet <<- dataSet;
    return(1); 
}
#############################################################################

######## D3 Visual Exploration ###################
PrepareSampleD3Plot <- function(viewOpt,smplNm,abunviewOpt,jsonNm){
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows =TRUE);
    }
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj));
    
    #get count data and tax_data (hierarchy) filtered data
    contdata <- as.data.frame(otu_table(data1));
    tax_dt <- as.data.frame(tax_table(data1));
    genelvlnm <- colnames(tax_dt);
    ngenelvl <- length(genelvlnm);
    sankey_data <- merge(tax_dt, contdata, by=0, all=TRUE);
    
    #parse it according the view (overall)
    sankey_data <- sankey_data[,c(colnames(tax_dt),smplNm)];
    
    #tax_dt is in factor type; so convert it to character datatype
    sankey_data[,1:ngenelvl] <- apply(sankey_data[,1:ngenelvl],2,as.character);
    
    #remove all genes with abundance is zero
    ind <- which(sankey_data[,smplNm]==0);     
    if(length(ind)>0){
        sankey_data <- sankey_data[-ind,];
    }
    c1 <- noquote(paste(genelvlnm, collapse='+'));
    my.formula <- as.formula(paste(".~", c1));
    data <- aggregate(my.formula,data=sankey_data,sum);
    if(abunviewOpt=="relative"){
        data <- ShowCountorRelative (data,viewOpt,smplNm);
    }
    write.csv(data, file="Hierarchical_feature_abundance.csv");
    sankeyobj <- list();
    if(viewOpt=="sankey"){
        nsrc <- ngenelvl-1; 
        src <- data[,c(1:nsrc)];
        src <- Reduce(c,src);
        ntrgt <- ngenelvl;
        trgt <- data[,c(2:ntrgt)];
        trgt <- Reduce(c,trgt);
        links=data.frame(source=c(src), target=c(trgt),value=c(data[,smplNm]));
        suppressMessages(require(dplyr));
        nodes=data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique());
        #links$IDsource=match(links$source, nodes$name)-1;
        #links$IDtarget=match(links$target, nodes$name)-1;
        sankeyobj$links <- links;
        sankeyobj$nodes <- nodes;
        require(jsonlite);
        json.obj <- jsonlite::toJSON(sankeyobj);
    } else {
        require(RJSONIO);
        json.obj <- RJSONIO::toJSON(list(name=smplNm,children=makeSunburstData(data,viewOpt)));
    }
    #save every thing in json file 
    sink(jsonNm);
    cat(json.obj);
    sink();
    return(1);
}
####################################
makeSunburstData <- function(x,viewOpt){
  if(ncol(x)>2){
    listSplit <- split(x[-1],x[1],drop=T)
    lapply(names(listSplit),function(y){list(name=y,children=makeSunburstData(listSplit[[y]],viewOpt))})
  }else{
    if(viewOpt=="treemap"){
        lapply(seq(nrow(x[1])),function(y){list(name=x[,1][y],value=x[,2][y])})
    }else{
        lapply(seq(nrow(x[1])),function(y){list(name=x[,1][y],size=x[,2][y])})
    }
  }
}
###########################################
PrepareSampleGrpD3Plot <- function(viewOpt,metadata,clslevel,cntmeth,feat_cnt,feat_filtmeth,abunviewOpt,jsonNm){
    data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows =TRUE);
    }
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj));
    #get sample names according to metadata and clslevel selected by user
    sam_data <- data.frame(sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE));
    smpl_nm <- rownames(subset(sam_data,sam_data[metadata]==clslevel));    
    contdata <- as.data.frame(otu_table(data1));
    contdata <- contdata[,smpl_nm];
    
    #filter genes based on count cutoff
    if(feat_filtmeth=="sum"){
        filt_ind<-which(rowSums(contdata)>feat_cnt);
        filt_ind1<-which(rowSums(contdata)<feat_cnt);
    } else{
        dt <- apply(contdata,1,median);
        filt_ind <- which(dt>feat_cnt);
        filt_ind1 <- which(dt<feat_cnt);
    }
    if(length(filt_ind)==0){
        current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    } else {
        contdata <- contdata[filt_ind,];
    }
    if(cntmeth=="sum"){
        contdata$Total_count <- rowSums(contdata);
    } else{
        #have to round the count in case of mean
        contdata$Total_count <- rowMeans(contdata);
        contdata$Total_count <- round(contdata$Total_count);
    }
    contdata$Total_count <- rowSums(contdata);
    
    #remove all genes with abundance is zero
    ind <- which(contdata[,"Total_count"]==0);    
    if(length(ind)>0){
        contdata<-contdata[-ind,];
    }
    tax_dt <- as.data.frame(tax_table(data1));
    genelvlnm <- colnames(tax_dt);
    ngenelvl <- length(genelvlnm);
    sankey_data <- merge(tax_dt, contdata, by=0, all=TRUE);
    
    #parse it according the view (overall)
    sankey_data <- sankey_data[,c(colnames(tax_dt),"Total_count")];
    
    #tax_dt is in factor type; so convert it to character datatype
    sankey_data[,1:ngenelvl]<-apply(sankey_data[,1:ngenelvl],2,as.character);
    c1 <- noquote(paste(genelvlnm, collapse='+'));
    my.formula <- as.formula(paste(".~", c1));
    data <- aggregate(my.formula,data=sankey_data,sum);
    if(abunviewOpt=="relative"){
        data <- ShowCountorRelative(data,viewOpt,"Total_count");
    }
    write.csv(data, file="Hierarchical_feature_abundance.csv");
    sankeyobj <- list();
    if (viewOpt=="sankey"){
        nsrc <- ngenelvl-1; 
        src <- data[,c(1:nsrc)];
        src <- Reduce(c,src);
        ntrgt <- ngenelvl;
        trgt <- data[,c(2:ntrgt)];
        trgt <- Reduce(c,trgt);
        links <- data.frame(source=c(src), target=c(trgt),value=c(data[,"Total_count"]));
        suppressMessages(require(dplyr));
        nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique());
        #links$IDsource=match(links$source, nodes$name)-1;
        #links$IDtarget=match(links$target, nodes$name)-1;
        sankeyobj$links <- links;
        sankeyobj$nodes <- nodes;
        require(jsonlite);
        json.obj <- jsonlite::toJSON(sankeyobj);
    } else {
        require(RJSONIO);
        json.obj <- RJSONIO::toJSON(list(name=clslevel,children=makeSunburstData(data,viewOpt)));
    }
    #save every thing in json file 
    sink(jsonNm);
    cat(json.obj);
    sink();
    return(1);
}
####################################
PrepareOverviewD3Plot <- function(viewOpt,cntmeth,feat_cnt,feat_filtmeth,abunviewOpt,jsonNm){
   data <- dataSet$filt.data;
    if(class(dataSet$filt.data)=="matrix"){
       data <- otu_table(data,taxa_are_rows =TRUE);
    }
    data1 <- merge_phyloseq(data,tax_table(dataSet$proc.phyobj));
    
    #get count data and tax_data (hierarchy) filtered data
    contdata <- as.data.frame(otu_table(data1));
    
    #filter genes based on count cutoff
    if(feat_filtmeth=="sum"){
        filt_ind <- which(rowSums(contdata)>feat_cnt);
        filt_ind1 <- which(rowSums(contdata)<feat_cnt);
    } else{
        dt <- apply(contdata,1,median);
        filt_ind <- which(dt>feat_cnt);
        filt_ind1 <- which(dt<feat_cnt);
    }
    if(length(filt_ind)==0){
        current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
        return(0);
    } else {
        contdata <- contdata[filt_ind,];
    }
    #calculate gene abundance
    if(cntmeth=="sum"){
        contdata$Total_count <- rowSums(contdata);
    } else {
        contdata$Total_count <- rowMeans(contdata);
        contdata$Total_count <- round(contdata$Total_count);
    }
    #remove all genes with abundance is zero (already done)
    ind <- which(contdata[,"Total_count"]==0);    
    if(length(ind)>0){
        contdata <- contdata[-ind,];
    }
    tax_dt <- as.data.frame(tax_table(data1));
    genelvlnm <- colnames(tax_dt);
    ngenelvl <- length(genelvlnm);
    sankey_data <- merge(tax_dt, contdata, by=0, all=TRUE);
    
    #parse it according the view (overall)
    sankey_data <- sankey_data[,c(colnames(tax_dt),"Total_count")];
    
    #tax_dt is in factor type; so convert it to character datatype
    sankey_data[,1:ngenelvl] <- apply(sankey_data[,1:ngenelvl],2,as.character);
    c1 <- noquote(paste(genelvlnm, collapse='+'));
    my.formula <- as.formula(paste(".~", c1));
    data <- aggregate(my.formula,data=sankey_data,sum);
    if(abunviewOpt=="relative"){
        data <- ShowCountorRelative(data,viewOpt,"Total_count");
    }
    write.csv(data, file="Hierarchical_feature_abundance.csv");
    sankeyobj <- list();
    if (viewOpt=="sankey"){
        nsrc <- ngenelvl-1; 
        src <- data[,c(1:nsrc)];
        src <- Reduce(c,src);
        ntrgt <- ngenelvl;
        trgt <- data[,c(2:ntrgt)];
        trgt <- Reduce(c,trgt);
        links <- data.frame(source=c(src), target=c(trgt),value=c(data[,"Total_count"]));
        suppressMessages(require(dplyr));
        nodes <- data.frame(name=c(as.character(links$source), as.character(links$target)) %>% unique());
        #links$IDsource=match(links$source, nodes$name)-1;
        #links$IDtarget=match(links$target, nodes$name)-1;
        sankeyobj$links <-links;
        sankeyobj$nodes <-nodes;
        require(jsonlite);
        json.obj <- jsonlite::toJSON(sankeyobj);
    } else{
        require(RJSONIO);
        json.obj <- RJSONIO::toJSON(list(name="All samples",children=makeSunburstData(data,viewOpt)));
    }
    #save every thing in json file 
    sink(jsonNm);
    cat(json.obj);
    sink();
    return(1);
}
##################################
round_preserve_sum <- function(x, digits = 0) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
}
#######################################
ShowCountorRelative <- function(data, viewOpt, abuncol) {
    #show relative % abundance or absolute count
    digits <- 0;
    if (viewOpt=="sankey"){
        digits <- 2;
    }
    data[,abuncol] <- round_preserve_sum(data[,abuncol]/ sum(data[,abuncol])*100,digits);
    return(data);
}
########################################

###########Alpha diversity#############
#######################################
PlotAlphaData <- function(bargraphName, distName, metadata, generank, colPal, format="png", dpi=72){
    set.seed(13133);
    if(generank=="Feature"){
        data <- dataSet$proc.phyobj;
    }else{
        #merging at taxonomy levels
        data <- fast_tax_glom_first(dataSet$proc.phyobj,generank);
    }
    #bargraphName = paste(bargraphName, ".", format, sep="");
    imgSet$alpha <- bargraphName;
    imgSet <<- imgSet;
    
    #reordering the sample (getting reordered rownames)
    sam <- sample_data(data);
    sam <- sam[order(sam[[metadata]])];
    smplord <- rownames(sam);
    smpl.num <- length(smplord);
    if(smpl.num<25){
        width <- 600
    }else if(smpl.num >= 25 & smpl.num <=50){
        width <- 750
    }else{
        width <- 950
    }
    Cairo(file=bargraphName, width=width, height=450, type=format, bg="white",dpi=dpi);
    ylab <- paste0("Alpha Diversity measure"," (",generank,"-level)");
    if(distName=="pielou"||distName=="evar"){
        alpha_data <- evenness(data, distName);
        alpha_data <- data.frame(alpha_data, sample_data(data));
        alpha_data <- tidyr::gather(alpha_data, key = "Measure", value = "value", distName);
        alpha_data$sample <- rownames(sample_data(data));
        CLASS <- alpha_data[,metadata];
        box <- ggplot(data = alpha_data, aes(x = sample, y = value, color= CLASS)) + facet_wrap(~Measure, scale = "free") +
                geom_point() +labs(y=ylab,color = metadata) +scale_x_discrete(limits =c(smplord));
    } else{
        box <- plot_richness(data,color=metadata,measures = distName) +scale_x_discrete(limits=c(smplord))+labs(y=ylab);
    }
    analSet$alpha <- box$data;
    write.csv(analSet$alpha, file="alphadiversity.csv");
    box=box+theme_bw()+theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1));
    box$layers <- box$layers[-1];
    box <- box+geom_point(size=3, alpha=0.7)+ theme(axis.text.x = black.bold.8.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title)+labs(x="Samples");
    
    #getting scale for plot (using same for boxplot also)
    ylimits <<- layer_scales(box)$y$range$range;
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal);
    }
    print(box);
    dev.off();
    analSet$alpha.meth <-distName;
    analSet$alpha.metadata <-metadata;
    analSet$alpha.genelvl <-generank;
    analSet <<- analSet;
    return(1);
}
########################################
PlotAlphaBoxData <- function(boxplotName, distName, metadata, colPal, generank, format="png", dpi=72){
    set.seed(1313397);
    data <- analSet$alpha;
    CLASS <- data[,metadata];
    imgSet$alpha.box <- boxplotName;
    imgSet <<- imgSet;
    ylab <- paste0("Alpha Diversity measure"," (",generank,"-level)");
    Cairo(file=boxplotName,width=500, height=400, type=format, bg="white",dpi=dpi);
    box1=ggplot(data,aes(data[,metadata],data$value))+stat_boxplot(geom ='errorbar',width=0.2)+geom_boxplot(aes(fill=CLASS, outlier.shape=1),position = position_dodge(width = 0),width=0.3)+ geom_jitter(alpha=0.6,width=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust =1,vjust= 1))+labs(y=ylab,x=metadata,fill = metadata)+coord_cartesian(ylim = c(ylimits[1],ylimits[2]));
    box1 <- box1 + theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    if(colPal!="default"){
        box1 <- box1+scale_fill_brewer(palette=colPal);
    }
    print(box1);
    dev.off();
    return(1);
}
########################################
PerformAlphaDiversityComp <- function(opt, metadata){
    data <- analSet$alpha;
    cls <- data[,metadata];
    x <- data$value;
    stat.info <- NULL;
    if(length(levels(cls)) > 2){
        if(opt=="tt"){
            res <- anova(aov(x ~ cls));
            stat.info <- paste("P-value: ", signif(res$"Pr(>F)"[1], 5), "; [ANOVA] F-value: ", signif(res$"F value"[1], 5), sep="");   
        }else{
            res <- kruskal.test(x ~ cls);
            stat.info <- paste("P-value: ", signif(res$p.value, 5), "; [Kruskal-Wallis] statistic: ", signif(res$statistic, 5) , sep="");   
        }
    }else{
        inx1 <- which(cls==levels(cls)[1]);
        inx2 <- which(cls==levels(cls)[2]);
        if(opt=="tt"){
            res <- t.test(x[inx1], x[inx2]);
            stat.info <- paste("P-value: ", signif(res$p.value, 5), "; [T-test] statistic: ", signif(res$statistic, 5), sep="");  
        }else{
            res <- wilcox.test(x[inx1], x[inx2]);
            stat.info <- paste("P-value: ", signif(res$p.value, 5), "; [Mann-Whitney] statistic: ", signif(res$statistic, 5), sep="");  
        }
    }
    analSet$alpha.stat.info <- stat.info;
    analSet<<- analSet;
    return(stat.info);
}
###############################################################################

###########Rarefaction analysis###############
PerformRarefactionCurves <- function(imgName,step,linecolor, linetype, format="png",dpi=72){
    # using orignal dataset
    data.src <- "orig";
    data_rare <- readRDS("data.orig"); 
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    data_rare <- merge_phyloseq(otu_table(data_rare,taxa_are_rows =TRUE),sample_table);
    if(min(table(factor(sample_table [[linecolor]]))) < 3 | min(table(factor(sample_table [[linetype]]))) < 3){
      current.msg <<-"Too many groups to be displayed - please select a more meaningful group option with at least 3 samples per group.";
      return(0);
    }
    #get good's coverage index
    #goods_coverage <- ComputeGoods(data_rare)
    #write.csv(goods_coverage, "goods_coverage.csv", row.names = FALSE, quote = FALSE);
    feat_no <- nsamples(data_rare);

    #adjust height according to number of legends
    w <- 500;
    if(feat_no < 10){
      w <- w;
    } else if (feat_no < 20){
      w<-w+100;
    } else if (feat_no < 50){
      w<-w+200;
    } else if (feat_no < 100){
      w<-w+350;
    } else if (feat_no > 100){
      w<-w+475;
    }
    Cairo(file=imgName, width = w, height = 540, type=format, bg="white", dpi=dpi);
    linecolor <- ifelse(linecolor == "none", "NULL", linecolor);
    linetype <- ifelse(linetype == "none", "NULL", linetype);
    box <- ggrare2(data_rare,
                   data.src = data.src,
                   color = linecolor,
                   label = "Sample",
                   linetype = linetype,
                   se = FALSE,  # this is not to meaningful
                   step = step);
    print(box);
    dev.off();
    return(1);
}

###########################################
## Rarefaction curves (ggrare code)
###########################################
ggrare2 <- function(physeq_object, data.src, label = NULL, color = NULL, 
            plot = TRUE, linetype = NULL, se = FALSE, step=5) {
    x <- methods::as(otu_table(physeq_object), "matrix")
    if (taxa_are_rows(physeq_object)) { x <- t(x) }
    ## This script is adapted from vegan `rarecurve` function
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)
    step_new = floor(max(tot)/as.integer(step))
    rarefun <- function(i) {
      #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
      n <- seq(1, tot[i], by = step_new)
      if (n[length(n)] != tot[i]) {
        n <- c(n, tot[i])
      }
      y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
      if (nrow(y) != 1) {
        rownames(y) <- c(".S", ".se")
        return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
      } else {
        return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
      }
    }
    f_n <- paste(data.src, step, "rds", sep = ".");
    if(file.exists(f_n)){
      df <- readRDS(file = f_n);
    } else {
      out <- lapply(seq_len(nr), rarefun)
      df <- do.call(rbind, out);
      saveRDS(df, file = f_n);
    }
    # Get sample data
    if (!is.null(sample_data(physeq_object, FALSE))) {
      sdf <- methods::as(sample_data(physeq_object), "data.frame")
      sdf$Sample <- rownames(sdf)
      data <- merge(df, sdf, by = "Sample")
      labels <- data.frame(x = tot, y = S, Sample = rownames(x))
      labels <- merge(labels, sdf, by = "Sample")
    }
    # Add, any custom-supplied plot-mapped variables
    if (length(color) > 1) {
      data$color <- color
      names(data)[names(data) == "color"] <- deparse(substitute(color))
      color <- deparse(substitute(color))
    }
    if (length(label) > 1) {
      labels$label <- label
      names(labels)[names(labels) == "label"] <- deparse(substitute(label))
      label <- deparse(substitute(label))
    }
    if (length(linetype) > 1) {
      data$linetype <- linetype
      names(data)[names(data) == "linetype"] <- deparse(substitute(linetype))
      linetype <- deparse(substitute(linetype))
    }
    p <- ggplot2::ggplot(data = data,
         ggplot2::aes_string(x = "Size",
         y = ".S",
         group = "Sample",
         color = color,
         linetype = linetype))
    p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Feature Richness")+theme_bw()+theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    if (!is.null(label)) {
      p <- p + ggplot2::geom_text(data = labels,
                                  ggplot2::aes_string(x = "x",
                                                      y = "y",
                                                      label = label,
                                                      color = color),
                                  size = 4, hjust = 0) +
        scale_x_continuous(expand = c(0, 0, 0.3, 0))
    }
    p <- p + ggplot2::geom_line()
    if (se) { ## add standard error if available
      p <- p +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                                 ymax = ".S + .se",
                                                 color = NULL,
                                                 fill = color),
                             alpha = 0.2)
    }
    invisible(p)
}
#############################################################################

###########Heatmap#############
#######################################
PlotHeatmap <- function(plotNm, smplDist, clstDist, palette,metadata,generank,viewOpt,doclust,format="png",showfeatname,isscalefeatures,rowV=F,colV=T,var.inx=NA,border=T,width=NA,dpi=72){
    set.seed(2805614);
    #used for color pallete
    variable <<- metadata; 
    data <- dataSet$norm.phyobj;
    dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(data,dataSet$taxa_table);
    
    #if more than 1500 features will be present;subset to most abundant=>1500 features.
    #OTUs already in unique names;
     if(ntaxa(data)>1500){
        data = prune_taxa(names(sort(taxa_sums(data), TRUE))[1:1500], data);
        viewOpt == "overview";
    }
    if(generank=="Feature"){
        data <- data;
        
        # sometimes the names are super big, need to prune it
        nm <- taxa_names(data);
        nm <- substr(nm, 1,20);
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <-nm;
    }else{
        #merging at taxonomy levels
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        y <- which(is.na(nm)==TRUE);
        
        #converting NA values to unassigned
        nm[y] <- "Unassigned";
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <-nm;
        
        #all NA club together
        data1 <- (t(sapply(by(data1,rownames(data1),colSums),identity)));
        nm <- rownames(data1);
    }
    rownames(data1) <- nm;
    
    # arrange samples on the basis of slected experimental factor and using the same for annotation also
    annotation <- data.frame(sample_data(data));
    ind <- which(colnames(annotation)!=metadata && colnames(annotation)!="sample_id");
    if(length(ind)>0){
        ind1 <- ind[1];
        annotation <- annotation[order(annotation[,metadata],annotation[,ind1]),];
    }else{
        annotation <- annotation[order(annotation[,metadata]),];
    }
    #there is an additional column sample_id which need to be removed first
    annotation <- subset(annotation,select = -sample_id);
    sam.ord <- rownames(annotation);
    data1 <- data1[,sam.ord];

    # set up colors for heatmap
    if(palette=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(palette == "heat"){
        colors <- heat.colors(256);
    }else if(palette == "topo"){
        colors <- topo.colors(256);
    }else if(palette == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
    }else if(palette == "spectral"){
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(256));
    }else if(palette == "viridis") {
        suppressMessages(require(viridis));
        colors <- rev(viridis(15));
    }else if(palette == "plasma") {
        suppressMessages(require(viridis));
        colors <- rev(plasma(15));
    }else if(palette == "inferno") {
        suppressMessages(require(viridis));
        colors <- rev(inferno(15));
    }else if(palette == "magma") {
        suppressMessages(require(viridis));
        colors <- rev(magma(15));
    }else if(palette == "cividis") {
        suppressMessages(require(viridis));
        colors <- rev(cividis(15));
    }else{
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }
     #setting the size of plot
    if(is.na(width)){
        minW <- 1000;
        myW <- ncol(data1)*12 + 300;
        if(myW < minW){
            myW <- minW;
        }   
        w <- round(myW/72,2);
    }
    myH <- nrow(data1)*10 + 150;
    h <- round(myH/72,2);
    if(viewOpt == "overview"){
        if(is.na(width)){
        minW <- 800;
        myW <- ncol(data1)*8 + 300;
        if(myW < minW){
            myW <- minW;
        }   
        w <- round(myW/72,2);
    }
    myH <- nrow(data1)*8 + 150;
    h <- round(myH/72,2);
    }
    if(border){
        border.col <- "grey60";
    }else{
        border.col <- NA;
    }
    imgSet$heatmap <- plotNm;
    imgSet <<- imgSet;
    Cairo(file = plotNm, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    require(pheatmap);
    
    # set up color schema for samples
    if(palette == "gray"){
        cols <- GetColorSchema(T);
        uniq.cols <- unique(cols);
    }else{
        cols <- GetColorSchema();
        uniq.cols <- unique(cols);
    }
    if(doclust == "T"){
        rowV <- T;
    }
    if(showfeatname == "T"){
        showfeatname <- T;
    } else {
        showfeatname <- F; 
    }
    if(isscalefeatures == "T"){
        scaleOpt <-"row";
    } else{
        scaleOpt <-"none"; 
    }
    if (smplDist == "bray"||smplDist== "jaccard"){
        require(vegan);
        rowDist <- vegdist(data1, method = smplDist,diag = T);
        smplDist <- vegdist(t(data1), method = smplDist,diag = T);
    } else {
        rowDist <- smplDist;
    }
    pheatmap(data1, 
        annotation=annotation, 
        fontsize=8, fontsize_row=8, 
        clustering_distance_rows=rowDist,
        clustering_distance_cols=smplDist,
        clustering_method=clstDist,
        show_rownames=showfeatname,
        border_color=border.col,
        cluster_rows=colV, 
        cluster_cols=rowV,
        scale=scaleOpt,
        color=colors,
        fontface="bold"
    );
    dev.off();
    # storing for Report Generation
    analSet$heatmap <- data1;
    analSet$heatmap.dist <- smplDist;
    analSet$heatmap.clust <- clstDist;
    analSet$heat.genelvl <- generank;
    analSet <<- analSet;
    return(1);
}
###############################################################################

###Dendrogram######
###################################
PlotTreeGraph <- function(plotNm,distnm,clstDist,metadata,generank,format="png", dpi=72, width=NA){ 
    set.seed(2805619);
    variable <<- metadata;
    data <- dataSet$norm.phyobj;
    dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(data,dataSet$taxa_table);
    if(generank!="Feature"){
        data <- fast_tax_glom_first(data,generank);
    }
    hc.cls <- as.factor(sample_data(data)[[variable]]);

    # must call distance within the phyloseq package
    dist.mat <- phyloseq::distance(data,distnm,type = "samples");

    # build the tree    
    hc_tree <- hclust(dist.mat, method=clstDist);
    imgSet$tree <- plotNm;
    imgSet <<- imgSet;
    if(is.na(width)){
        w <- minH <- 650;
        myH <- nsamples(data)*16 + 150;
        if(myH < minH){
            myH <-minH;
        }   
        w <- round(w/72,2);
        h <- round(myH/72,2);
    }
    Cairo(file=plotNm, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(mar=c(4,2,2,10));
    clusDendro <- as.dendrogram(hc_tree);
    cols <- GetColorSchema();
    names(cols) <- sample_names(data);
    labelColors <- cols[hc_tree$order];
    colLab <- function(n){
        if(is.leaf(n)) {
            a <- attributes(n);
            labCol <- labelColors[a$label];
            attr(n, "nodePar") <- 
                        if(is.list(a$nodePar)){
                            c(a$nodePar,lab.col = labCol,pch=NA)
                        }else{
                            list(lab.col = labCol,pch=NA)
                        }
        }
        n
    }
    clusDendro <- dendrapply(clusDendro, colLab);
    plot(clusDendro,horiz=T,axes=T);
    par(cex=1);
    legend.nm <- unique(as.character(hc.cls));
    legend.nm <- gsub("\\.", " ",legend.nm)
    legend("topleft", legend = legend.nm, pch=15, col=unique(cols), bty = "n");
    dev.off();
    analSet$tree <- hc_tree;
    analSet$tree.dist <- distnm;
    analSet$tree.clust <- clstDist;
    analSet$tree.genelvl <- generank;
    analSet <<- analSet;
    return(1);
}
###############################################################################

###########Ordination analysis######
#######################################
PlotBetaData <- function(barplotNm,ordmeth,distName,colopt,metadata,showlabel,generank,alphaopt,colPal,format="png",dpi=72){    
    metadata <- metadata;
    #using normalized data
    taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(dataSet$norm.phyobj,taxa_table);          
    #merging at profile annotation levels
    if(generank!="Feature"){
        data <- fast_tax_glom_first(data,generank);
    }
    if(colopt=="alphadiv") {
        data1 <- dataSet$proc.phyobj;
        box <- plot_richness(data1,measures = alphaopt);
        alphaboxdata <- box$data;
        sam_nm <- sample_names(data);
        alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,]; 
        alphaval <- alphaboxdata$value;
        sample_data(data)$alphaopt <-alphaval;
        indx <- which(colnames(sample_data(data))=="alphaopt");
        colnames(sample_data(data))[indx] <- alphaopt;
    }else {
        data <- data;
    }
    #random_tree <- phy_tree(createRandomTree(ntaxa(data),rooted=TRUE,tip.label=taxa_names(data)));
    #data <- merge_phyloseq(data,random_tree);
    ord <- ordinate(data,method = ordmeth,distName);   
    
    #barplotNm = paste(barplotNm, ".", format, sep="");
    imgSet$beta2d <- barplotNm;
    imgSet <<- imgSet;
    Cairo(file=barplotNm, width=720, height=500, type=format, bg="white",dpi=dpi);
    if(colopt=="alphadiv") {
        box = plot_ordination(data,ord,color=alphaopt)+ scale_colour_gradient(low="green", high="red");
    }else{
        box = plot_ordination(data,ord,color=metadata);
    }
    box$layers <- box$layers[-1];
    if(showlabel=="samnm"){
        box = box+geom_text(aes(label=sample_id),hjust=0.5, vjust=2,size=3,fontface="bold");
        box=box+geom_point(size =4,alpha=0.6)+theme_bw();
    }else if(showlabel=="none"){
        box=box+geom_point(size =4,alpha=0.8)+theme_bw();
    }else{
        showlabel <<- showlabel;
        bx_data <<- data.frame(box$data);
        box=box+geom_text(aes(label=bx_data[ ,showlabel]),hjust=0.5, vjust=2,size=3,fontface="bold");
        box=box+geom_point(size =4,alpha=0.6)+theme_bw();
    }
    #used for area color for ellipse
    if(colopt=="expfac"){
        sam_data <- sample_data(data);
        clsLbl <- sam_data[[metadata]];
        box=box+ stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = clsLbl), show.legend=FALSE);
    }
    if(colPal!="default" & colopt!="alphadiv"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    box <- box+ theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    
    #saving info for report generation
    analSet$beta <- data;
    analSet$beta.meth <- ordmeth;
    analSet$beta.dist <- distName;
    analSet$beta.genelvl <- generank;
    analSet <<- analSet;
}

###########Ordination analysis (3D plots)########
#########################################################
# perform PCoA analysis 
PCoA3D.Anal <- function(ordMeth,distName,generank,colopt,variable,alphaopt,jsonNm){
    require(vegan);
    variable <<- variable;
    taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(dataSet$norm.phyobj,taxa_table);       
   
    #merging at taxonomy levels
    if(generank!="Feature"){
        data <- fast_tax_glom_first(data,generank);
    }
    if(colopt=="alphadiv") {
        data1 <- dataSet$proc.phyobj;
        box <- plot_richness(data1,measures = alphaopt);
        alphaboxdata <- box$data;
        sam_nm <- sample_names(data);
        alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,]; 
        alphaval <- alphaboxdata$value;
        sample_data(data)$alphaopt <- alphaval;
        indx <- which(colnames(sample_data(data))=="alphaopt");
        colnames(sample_data(data))[indx] <- alphaopt;
    }else{
        data <- data;
    }
    datacolby <<- data;
    if(ordMeth=="NMDS"){
        if(distName=="wunifrac"){
            GP.ord <- ordinate(data,ordMeth,"unifrac",weighted=TRUE,k=3);
        }else{
            GP.ord <- ordinate(data,ordMeth,distName,k=3);
        }
    }else{
        if(distName=="wunifrac"){
            GP.ord <- ordinate(data,ordMeth,"unifrac",weighted=TRUE);
        }else{
            GP.ord <- ordinate(data,ordMeth,distName);
        }
    }
    # obtain variance explained
    sum.pca <- GP.ord;
    imp.pca <- sum.pca$values;
    #imp.pca<-data.frame(imp.pca);
    std.pca <- imp.pca[1,]; # eigen values
    var.pca <- imp.pca[,2]; # variance explained by each PC
    cum.pca <- imp.pca[5,]; # cummulated variance explained
    sum.pca <- append(sum.pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
    pca3d <- list();
    if(ordMeth=="NMDS"){
        pca3d$score$axis <- paste("NMDS", 1:3 , sep="");
        coord <- sum.pca$points;
        write.csv(signif(coord,5), file="ordination_score.csv");
        list2 <- rep(as.numeric(0),nrow(coord));
        coord <- cbind(coord, list2);
        coords <- data.frame(t(signif(coord[,1:3], 5)));
    }else{
        pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(sum.pca$variance[1:3], 3), "%)", sep="");
        coords <- data.frame(t(signif(sum.pca$vectors[,1:3], 5)));
        write.csv(signif(sum.pca$vectors,5), file="ordination_score.csv");  
    }
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    pca3d$score$name <- sample_names(dataSet$norm.phyobj);
    col.type <- "factor";
    if(colopt=="alphadiv") {
        cls <- sample_data(data)[[alphaopt]];
        col.type <- "gradient";
        cols <- ComputeColorGradient(cls);
    }else {
        cls <- factor(sample_data(dataSet$norm.phyobj)[[variable]]);
        
        # now set color for each group
        cols <- unique(as.numeric(cls)) + 1;
    }
    pca3d$score$type <- col.type;
    pca3d$score$facA <- cls;
    rgbcols <- col2rgb(cols);
    cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    pca3d$score$colors <- cols;
    require(RJSONIO);
    json.obj <- toJSON(pca3d);
    sink(jsonNm);
    cat(json.obj);
    sink();
}
###############################################################################
PreparePCA4Shotgun <- function(imgName,imgName2,colopt,alphaopt,generank, format="json", inx1, inx2, inx3,variable,showlabel,colPal,format2d="png",dpi=72){
    set.seed(13134);
    #imgName2 = paste(imgName2, ".", format2d, sep="");
    imgSet$pca <- imgName2;
    imgSet <<- imgSet;
    #using normalized data
    taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(dataSet$norm.phyobj,taxa_table);         
    #merging at taxonomy levels
    if(generank!="Feature"){
        data <- fast_tax_glom_first(data,generank);
    }
    if(colopt=="alphadiv") {
        data1 <- dataSet$proc.phyobj;
        box <- plot_richness(data1,measures = alphaopt);
        alphaboxdata <- box$data;
        sam_nm <- sample_names(data);
        alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,]; 
        alphaval <- alphaboxdata$value;
        sample_data(data)$alphaopt <- alphaval;
        indx <- which(colnames(sample_data(data))=="alphaopt");
        colnames(sample_data(data))[indx] <- alphaopt;
        cls <- sample_data(data)[[alphaopt]];
        col.type <- "gradient";
        cols <- ComputeColorGradient(cls);
    }else{
        data <- data;
    }
    dat <- as.matrix(otu_table(data));
    pca3d <- list();
    pca <- prcomp(t(dat), center=T, scale=T);
    imp.pca <- summary(pca)$importance;
    write.csv(signif(pca$x,5), file="ordination_score.csv");
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
    
    #3D
    coords <- data.frame(t(signif(pca$x[,1:3], 5)));
    colnames(coords) <- NULL; 
    pca3d$score$type <- "factor";
    pca3d$score$xyz <- coords;
    pca3d$score$name <- sample_names(data);
    sam_data <- data.frame(sample_data(data));
    if(colopt=="alphadiv") {
        rgbcols <- col2rgb(cols);
        cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
        variable <- alphaopt;
        pca3d$score$type <- col.type;
    }else {
        cls <- as.character(sample_data(data)[[variable]]);
        variable <- variable;
        # now set color for each group
        cols <- unique(as.numeric(factor(cls)))+1;
        rgbcols <- col2rgb(cols);
        cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
    }
    clsLbl <- sam_data[[variable]];
    pca3d$score$facA <- cls;
    pca3d$score$colors <- cols;
    require(RJSONIO);
    json.obj <- toJSON(pca3d);
    sink(imgName);
    cat(json.obj);
    sink();

    #2D
    require("ggfortify");
    Cairo(file=imgName2, width=720, height=500, type=format2d, bg="white",dpi=dpi);
    label=FALSE;
    if(showlabel=="samnm"){
        label=TRUE;
        box <- autoplot(pca,data=sam_data, colour=variable, size=4,alpha=0.8, label=label)+theme_bw()+labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
    } else if(showlabel=="none") {
        box <- autoplot(pca,data=sam_data,colour=variable,size=4,alpha =0.8,label = label)+theme_bw()+labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
    } else{
       grplbl <- sam_data[ ,showlabel];
       clsLbl <- clsLbl;
       box <- autoplot(pca, data=sam_data, colour=variable, size=4, alpha =0.8, label = label)+geom_text(aes(label=grplbl, colour=clsLbl))+theme_bw()+labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
    }
    if(colopt=="alphadiv") {
        box <-box+scale_colour_gradient(low="green", high="red");
    }else{
        box <-box+ stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = clsLbl), show.legend=FALSE);
    }
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    box <- box + theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    analSet$pca <- pca;
    analSet <<- analSet;
}
###############################################################################
PerformCategaoryComp <- function(method,distnm,variable){
    require("vegan");
    data <- dataSet$proc.phyobj;
    #random_tree <- phy_tree(createRandomTree(ntaxa(data),rooted=TRUE,tip.label=taxa_names(data)));
    #data <- merge_phyloseq(data,random_tree);  
    data <- transform_sample_counts(data, function(x) x/sum(x));
    data.dist <- phyloseq::distance(data, method=distnm);
    group <- get_variable(data,variable);
    stat.info <- "";
    resTab <- list();
    if(method=="adonis"){
        sampledf <- data.frame(sample_data(data));
        res <- adonis(data.dist~sampledf[ ,variable],data = sampledf);
        resTab <- res$aov.tab[1,];
        stat.info <- paste("[PERMANOVA] F-value: ", signif(resTab$F.Model, 5),  "; R-squared: ", signif(resTab$R2, 5), "; P-value < ", signif(resTab$Pr, 5), sep="");
    }else if(method=="anosim"){
        anosim <- anosim(data.dist,group=group);
        resTab$Rval <- anosim$statistic;
        resTab$pval <- anosim$signif;   
        stat.info <- paste("[ANOSIM] R: ", signif(resTab$Rval, 5), "; P-value < ", signif(resTab$pval, 5), sep="");
    }else if (method=="permdisp") {
        beta <- betadisper(data.dist,group=group);
        resTab <- anova(beta);
        stat.info <- paste("[PERMDISP] F-value: ", signif(resTab$"F value"[1], 5), "; P-value: ", signif(resTab$"Pr(>F)"[1], 5), sep="");   
    }
    analSet$stat.info <- stat.info;
    analSet <<- analSet;
    return(1);
}
###############################################################################
###########Core Resistome analysis########
#######################################
CoreResistoAnalysis <- function(imgName,preval,detection,generank,palette,viewOpt,analOpt, expFact, group,format="png",dpi=72,width=NA){
    require("tidyr");
    data <- dataSet$proc.phyobj;
    expFact <- expFact;
    group <- group;
  
    if(!analOpt == "all_samples"){
        data <- eval(parse(text = paste("phyloseq:::subset_samples(data,", expFact, "==", "\"", group, "\"", ")", sep="")))
    
        # check min 2 reps 
        samples_left <- nsamples(data)
        if(samples_left<2){
          current.msg <<- "More than 2 replicates are required in your group!"
          return(0)
        }
    }
    if(generank=="Feature"){
        data <- otu_table(data,taxa_are_rows=T);
        nm <- taxa_names(data);
        nm <- substr(nm, 1,20);
        if(anyDuplicated(nm)){ 
            nm <- make.unique(as.character(nm), sep = "_");
            taxa_names(data) <- nm;
        }
    }else{
        #merging at taxonomy levels
        data <- fast_tax_glom_first(data,generank);
        tax_table <- tax_table(data);
        nm <- as.character(tax_table(data)[,generank]);
        y <- which(is.na(nm)==TRUE);
       
        #converting NA values to unassigned
        nm[y] <- "Unassigned";
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
        data <- otu_table(data1,taxa_are_rows=T);
    }
    data.compositional <- transform_sample_counts(data,function(x) x / sum(x));
    data.core <- core(data.compositional, detection = detection, prevalence = preval);
    core.nm <- as.data.frame(prevalence(data.compositional, detection = detection, sort = TRUE));
    colnames(core.nm)[1] <- "Prevelance";
    fileName <- "core_resistome.csv";
    write.csv(core.nm,file=fileName);    
    imgSet$core <- imgName;
    imgSet <<- imgSet;
   
    #if more than 1500 features will be present;subset to most abundant=>1500 features.
    #OTUs already in unique names;
     if(ntaxa(data.core)>1500){
        data.core = prune_taxa(names(sort(taxa_sums(data.core), TRUE))[1:1500], data.core);
        viewOpt == "overview";
    }
    #setting the size of plot
    if(is.na(width)){
        minW <- 800;
        myW <- 10*18 + 200;
        if(myW < minW){
            myW <- minW;
        }   
        w <- round(myW/72,2);
    }
    myH <- nrow(data.core)*18 + 150;
    h <- round(myH/72,2);
    if(viewOpt == "overview"){
        if(is.na(width)){
            if(w >9.3){
                w <- 9.3;
            }
        }
        if(h > w){
            h <- w;
        }
    }
    Cairo(file=imgName, unit="in",width=w, height=h, type=format, bg="white",dpi=dpi);
    
    # set up colors for heatmap
    if(palette=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(10);
    }else if(palette == "heat"){
        colors <- heat.colors(10);
    }else if(palette == "topo"){
        colors <- topo.colors(10);
    }else if(palette == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"), space="rgb")(10);
    }else if(palette == "spectral"){
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(10));
    }else if(palette == "viridis") {
        suppressMessages(require(viridis));
        colors <- rev(viridis(15));
    }else if(palette == "plasma") {
        suppressMessages(require(viridis));
        colors <- rev(plasma(15));
    }else if(palette == "inferno") {
        suppressMessages(require(viridis));
        colors <- rev(inferno(15));
    }else if(palette == "magma") {
        suppressMessages(require(viridis));
        colors <- rev(magma(15));
    }else if(palette == "cividis") {
        suppressMessages(require(viridis));
        colors <- rev(cividis(15));
    }else{
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(10));
    }
    p <- plot_core(data.core,plot.type = "heatmap",colours = colors,prevalences = seq(.05, 1, .05),detections = 10^seq(log10(1e-3), log10(0.2), length = 10))+ xlab("Detection Threshold (Relative Abundance (%))") + ylab(generank)+guides(fill = guide_legend(keywidth = 1.5, keyheight = 1));
    analSet$core <- as.matrix(core.nm);
    analSet$core.genelvl <- generank;
    analSet <<- analSet;
    print(p);
    dev.off();
    return(1);
}

###############################################################################
#######Correlation analysis########################
##########################################
PrepareCorrHeatMap <- function(imgName, format="png", width=NA, cor.method, 
                colors_cntrst, viewOpt,generank,fix.col, no.clst, top, topNum,datatype){
    main <- xlab <- ylab <- NULL;
    taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(dataSet$norm.phyobj,taxa_table);
    if(generank=="Feature"){
        data1 <- as.matrix(otu_table(data));  
    } else {
        #merging at taxonomy levels
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    data <- t(data1);
    analSet$abund_data <-data;
    if(ncol(data) > 1500){
        filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
        rk <- rank(-filter.val, ties.method='random');
        data <- as.data.frame(data[,rk <=1500]);
    }
    if(generank=="Feature"){
        colnames(data) <- substr(colnames(data), 1, 20);
    }
    corr.mat <- cor(data, method=cor.method);

    # use total abs(correlation) to select
    if(top){
        cor.sum <- apply(abs(corr.mat), 1, sum);
        cor.rk <- rank(-cor.sum);
        var.sel <- cor.rk <= topNum;
        corr.mat <- corr.mat[var.sel, var.sel];
    }
    # set up parameter for heatmap
    suppressMessages(require(RColorBrewer));
    suppressMessages(require(gplots));
    if(colors_cntrst=="gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colors_cntrst == "heat"){
        colors <- heat.colors(256);
    }else if(colors_cntrst == "topo"){
        colors <- topo.colors(256);
    }else if(colors_cntrst == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"))(256);
    }else if(colors_cntrst == "spectral"){
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(10));
    }else if(colors_cntrst == "viridis") {
        suppressMessages(require(viridis));
        colors <- rev(viridis(15));
    }else if(colors_cntrst == "plasma") {
        suppressMessages(require(viridis));
        colors <- rev(plasma(15));
    }else if(colors_cntrst == "inferno") {
        suppressMessages(require(viridis));
        colors <- rev(inferno(15));
    }else if(colors_cntrst == "magma") {
        suppressMessages(require(viridis));
        colors <- rev(magma(15));
    }else if(colors_cntrst == "cividis") {
        suppressMessages(require(viridis));
        colors <- rev(cividis(15));
    }else{
        colors <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256));
    }
    #used for network edges  
    first_color <- colors[1];
    last_color <- colors[length(colors)];
    if(viewOpt == "overview"){
        if(ncol(corr.mat) > 50){
            myH <- ncol(corr.mat)*8 + 40;
        }else if(ncol(corr.mat) > 20){
            myH <- ncol(corr.mat)*8 + 60;
        }else{
            myH <- 650;
        }
        h <- round(myH/72,2);
        w <- h;
    }else{
        if(ncol(corr.mat) > 50){
            myH <- ncol(corr.mat)*12 + 40;
        }else if(ncol(corr.mat) > 20){
            myH <- ncol(corr.mat)*12 + 60;
        }else{
            myH <- 700;
        }
        h <- round(myH/72,2);

        if(is.na(width)){
            w <- h;
        }else if(width == 0){
            w <- h <- 8.2;
            imgSet$corr.heatmap <-imgName;
        }else{
            w <- h <- 8.2;
        }
        # to prevent too small
        min.w <- 4.8;
        if(w < min.w){
            w <- h <- min.w;
        }
    }
    Cairo(file = imgName, unit="in", dpi=72, width=w, height=h, type=format, bg="white");
    imgSet$cor.heat<-imgName;
    imgSet <<- imgSet;
    if(no.clst){
        rowv=FALSE;
        colv=FALSE;
        dendro= "none";
    }else{
        rowv=TRUE;
        colv=TRUE;
        dendro= "both";
    }
    require(pheatmap);
    if(fix.col){
        breaks <- seq(from = -1, to = 1, length = 257);
        pheatmap(corr.mat, 
            fontsize=8, fontsize_row=8, 
            cluster_rows = colv, 
            cluster_cols = rowv,
            color = colors,
            breaks = breaks
            );
    }else{
        pheatmap(corr.mat, 
            fontsize = 8, fontsize_row = 8, 
            cluster_rows = colv, 
            cluster_cols = rowv,
            color = colors,
            fontface= "bold"
        );
    }
    dev.off();
    analSet$cor.heatmat <- corr.mat;
    analSet$colors_cntrst <- colors_cntrst;
    analSet$edge_col_low <- first_color;
    analSet$edge_col_high <- last_color;
    analSet$colors <- colors;
    analSet$corheat.genelvl <- generank;
    analSet$corheat.meth <- cor.method;
    analSet <<- analSet;
    write.csv(signif(corr.mat,5), file="correlation_table.csv");
    return(1);
}
###############################################################################

###########DE(feature boxplot)#############
#######################################
PlotBoxData <- function(boxplotName, feat, format="png", dpi=72){
    require(ggplot2);
    data <- analSet$boxdata;
    a <- data[,feat];
    ind <- which(a=="0");
    a[ind] <- 0.1;
    data$log_feat <- log(a);
    
    #boxplotName = paste(boxplotName,".",format, sep="");
    Cairo(file=boxplotName,width=720, height=360, type=format, bg="white",dpi=dpi);
    box=ggplot(data,aes(x=data$class, y = data[,feat]))+stat_boxplot(geom ='errorbar')+ 
            geom_boxplot(aes(fill=class, outlier.shape=21))+ geom_jitter(alpha=0.6,width=0.1)+scale_color_brewer(palette="Set1")+scale_fill_brewer(palette="Set2")+theme_bw()+labs(y="Abundance",x= de.var) + ggtitle("Orignal Count") +theme(plot.title = element_text(hjust=0.5),legend.position="none",axis.text.x = element_text(angle = 45,hjust =1,vjust= 1));
    box1=ggplot(data,aes(x=data$class, y = data$log_feat))+stat_boxplot(geom ='errorbar')+ 
            geom_boxplot(aes(fill=class, outlier.shape=21))+ geom_jitter(alpha=0.6,,width=0.1)+scale_color_brewer(palette="Set1")+scale_fill_brewer(palette="Set2")+theme_bw()+labs(y="",x=de.var, fill= de.var) + ggtitle("Log-transformed Count")+ theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(angle = 45,hjust =1,vjust= 1));
    require(grid)
    require(gridExtra)
    grid.arrange(ggplotGrob(box), ggplotGrob(box1),ncol=2,nrow=1,top=feat);  
    dev.off();
}
#######################################
#################Differential analysis##################################
###########EdgeR/DESeq2########
#######################################
PerformRNAseqDE <- function(opts,p.lvl,variable,generank){
    claslbl <- as.factor(sample_data(dataSet$norm.phyobj)[[variable]]);

    # build phyloseq obj in fly
    filt.dataphy <- dataSet$filt.data;
    filt.dataphy < -apply(filt.dataphy,2,as.integer);
    filt.dataphy <- otu_table(filt.dataphy,taxa_are_rows =TRUE);
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    filt.dataphy <- merge_phyloseq(filt.dataphy,sample_table);
    taxa_names(filt.dataphy) <- rownames(dataSet$filt.data);
    data <- filt.dataphy;
    dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(data,dataSet$taxa_table);
    if(generank=="Feature"){
        data <- data;
        nm <- taxa_names(data);
    }else{
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
        data1 <- otu_table(data1,taxa_are_rows=T);
        data <- merge_phyloseq(data1,sample_data(data));
        nm <- taxa_names(data);
    }
    dat3t <- as.data.frame(t(otu_table(data))); 
    colnames(dat3t) <- nm;
    if(opts=="DESeq2"){
        # only for small data set (< 100)
        if(length(claslbl) > 100){
            current.msg <<- "Only EdgeR is supported for sample size over 100."; 
            return(0);
        }else{
            require("DESeq2");      
            my.formula <- as.formula(paste("~", variable));
            #converting from phyloseq object to deseq
            diagdds = phyloseq_to_deseq2(data,my.formula);
            geoMeans = apply(counts(diagdds), 1, gm_mean);
            diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans);
            diagdds = DESeq(diagdds, test="Wald", fitType="parametric");
            res = results(diagdds, independentFiltering = FALSE, cooksCutoff =  Inf);
            sigHits <- which(res$padj < p.lvl);
            de.Num <- length(sigHits);
            if(de.Num == 0){
                current.msg <<- "No significant genes were identified using the given adjusted P-value cutoff."; 
            } else{
                current.msg <<- paste("A total of", de.Num, "significant features were identified!");
            }
            resTable <- res[,c("log2FoldChange" ,"lfcSE","pvalue","padj")];
            resTable <- signif(data.matrix(resTable), digits=5); 
            colnames(resTable) <- c("log2FC","lfcSE","Pvalues","FDR");
            analSet$anal.type <- "deseq";
        }
    } else{
        #using by filtered data ,RLE normalization within it.
        dge <- phyloseq_to_edgeR(data,variable);
        et = exactTest(dge);
        tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue");
        res = tt@.Data[[1]];
        de.Num <- sum(res$FDR < p.lvl);
        if(de.Num == 0){
            current.msg <<- "No significant genes were identified using the given adjusted P-value cutoff."; 
        } else{
            current.msg <<- paste("A total of", de.Num, "significant features were identified!");
        }
        resTable <- res[,c("logFC","logCPM","PValue","FDR")];
        resTable <- signif(data.matrix(resTable), digits = 5); 
        colnames(resTable) <- c("log2FC","logCPM","Pvalues","FDR");
        analSet$anal.type <- "edgr";
    }
    resTable <- data.frame(resTable);
    ord.inx <- order(resTable$Pvalues);
    resTable <- resTable[ord.inx, , drop=FALSE];
    write.csv(resTable, file="rnaseq_de.csv");
    if(nrow(resTable) > 500){
        resTable <- resTable[1:500, ];
    }
    analSet$rnaseq$resTable<-analSet$resTable<-data.frame(resTable);
    
    #only getting the names of DE features
    diff_ft <<- rownames(resTable)[1:de.Num]; 
    generank <<- generank;
    
    #individual boxplot for features
    sigfeat <- rownames(resTable);
    box_data <- as.data.frame(dat3t[ ,sigfeat]);
    colnames(box_data) <- sigfeat;    
    box_data$class <- claslbl;
    analSet$boxdata <- box_data;
    analSet$sig.count <- de.Num;
    analSet$rnaseq.genelvl <- generank;
    analSet$rnaseq.meth <- opts;
    tree_data <<- data;
    de.var <<- variable;
    analSet <<- analSet;
    return(1);
}
###############################################################################
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
################################
edgeRnorm <- function(x,method){
    # Enforce orientation.
    # See if adding a single observation, 1, 
    # everywhere (so not zeros) prevents errors
    # without needing to borrow and modify 
    # calcNormFactors (and its dependent functions)
    # It did. This fixed all problems. 
    # Can the 1 be reduced to something smaller and still work?
    x = x + 1;

    # Now turn into a DGEList
    y = DGEList(counts=x, remove.zeros=TRUE);

    # Perform edgeR-encoded normalization, using the specified method (...)
    z = edgeR::calcNormFactors(y, method=method);

    # A check that we didn't divide by zero inside `calcNormFactors`
    if(!all(is.finite(z$samples$norm.factors)) ){
            current.msg <<- paste("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors.");
            return(0);
          }
    return(z)
}
###################################
####################################
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
    require("edgeR")
    require("phyloseq")
    # Enforce orientation.
    if( !taxa_are_rows(physeq) ){ physeq <-t(physeq) }
    x = as(otu_table(physeq), "matrix")
    # Add one to protect against overflow, log(0) issues.
    x = x + 1
    
    # Check `group` argument
    if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
        
      # Assume that group was a sample variable name (must be categorical)
        group = get_variable(physeq, group)
    }
    # Define gene annotations (`genes`) as tax_table
    taxonomy = tax_table(physeq, errorIfNULL=FALSE)
    if( !is.null(taxonomy) ){
        taxonomy = data.frame(as(taxonomy, "matrix"))
    } 

    # Now turn into a DGEList
    y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
    z = edgeR::calcNormFactors(y, method="RLE");
    # Estimate dispersions
    return(estimateTagwiseDisp(estimateCommonDisp(z)))
}
###############################################################################

###########MetagenomeSeq ########
#######################################
PerformMetagenomeSeqAnal <- function(p.lvl,variable,generank,model){
    require("metagenomeSeq");
    filt.dataphy <- dataSet$filt.data;
    filt.dataphy <- apply(filt.dataphy,2,as.integer);
    filt.dataphy <- otu_table(filt.dataphy,taxa_are_rows =TRUE);
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    filt.dataphy <- merge_phyloseq(filt.dataphy,sample_table);
    taxa_names(filt.dataphy) <- rownames(dataSet$filt.data);
    data <- filt.dataphy;
    cls <- as.factor(sample_data(data)[[variable]]);
    lvl <- length(levels(cls));
    if(generank!="Feature"){          
        dataSet$taxa_table<-tax_table(dataSet$proc.phyobj);
        data <- merge_phyloseq(data,dataSet$taxa_table);
    }else{
        data <- data;
    }
    if(generank=="Feature"){
        data <- data;
        nm <- taxa_names(data);
    }else{
        #merging at taxonomy levels
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
       
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
        data1 <- otu_table(data1,taxa_are_rows=T);
        data <- merge_phyloseq(data1,sample_data(data));
        nm <- taxa_names(data);
    }
    data <- phyloseq_to_metagenomeSeq(data);
    data <- cumNorm(data,p=cumNormStat(data));
    mod <- model.matrix(~phenoData(data)@data[,variable]);
    if (model=="zigfit"){
        fit <- fitZig(data, mod);
    }else {
        if(length(levels(cls)) > 2){
            current.msg <<- paste( "More than two group present in sample variable. This model can only be used with two groups.");
            return(0);
        } else {
            fit <- fitFeatureModel(data, mod);
        }
    }
    x <- MRfulltable(fit, number = nrow(assayData(data)$counts),adjustMethod = "BH");
    x <- x[!is.na(rownames(x)), ];
    rownames(x) <- gsub(":1", "", x = rownames(x), fixed = TRUE);
    x$OTUnames <- as.character(rownames(x))
    if(!is.null(tax_table(data, errorIfNULL = FALSE))) {
        #Attach the bacterial taxonomy to the table, if available
        TAX = data.frame(tax_table(data));
        TAX$OTUnames <- as.character(rownames(TAX));
        res = merge(x, TAX, by = "OTUnames")
    } else {
        res = x;
    }
    na.indx <- which(is.na(res$pvalues)=="TRUE");
    
    #NA introduced so remove those
    if(length(na.indx)>0){
        res <- res[-na.indx,];
    }
    # Sort and return #sighits return TRUE or FALSE
    sigHits <- res$adjPvalues<=p.lvl;
    de.Num <- length(which(sigHits));
    if(de.Num == 0){
        current.msg <<- paste( "No significant genes were identified using the given adjusted P-value cutoff. Please change the cutoff limit."); 
    }else{
        current.msg <<- paste("A total of", de.Num, "significant genes were identified!")
    }
    if(model=="ffm"){
        resTable <- res[,c("pvalues","adjPvalues","logFC")];
        resTable <- signif(resTable[,c(3,1,2)], digits = 5);
        colnames(resTable) <- c("logFC","Pvalues","FDR");
    }else {
        resTable <- res[,c("pvalues","adjPvalues")];
        resTable <- signif(resTable[,1:2], digits = 5);
        colnames(resTable) <- c("Pvalues","FDR");
    }
    ord.inx <- order(resTable$Pvalues);
    resTable <- resTable[ord.inx, , drop=FALSE];
    write.csv(resTable, file="metageno_de_output.csv");   
    if(nrow(resTable) > 500){
        resTable <- resTable[1:500, ];
    }
    analSet$metagenoseq$resTable <- analSet$resTable <- data.frame(resTable);
    
    #only getting the names of DE features
    diff_ft <<- rownames(resTable)[1:de.Num];
    sigfeat <- rownames(resTable);
    
    #prepare individual boxplot
    box_data <- MRcounts(data);   
    #subset only diff. Abundant features
    box_data <- box_data[sigfeat, ];

    #samples in rows
    box_data <- t(box_data);
    box_data <- data.frame(box_data);  
    colnames(box_data) <- sigfeat;
    claslbl <- pData(data)[ ,variable];
    box_data$class <- unlist(claslbl);
    analSet$boxdata <- box_data;
    analSet$sig.count <- de.Num;
    de.var <<- variable;
    analSet$anal.type <- "metagseq";
    analSet$metageno.genelvl <- generank;
    analSet <<- analSet;
    return(1);
}
###############################################################################

###########LEfSe ########
#######################################
PerformLefseAnal <- function(p.lvl,lda.lvl,variable,generank){
    require("MASS");
    filt.data <- dataSet$filt.data;
    # cpm normalization
    data <- cpm.default(filt.data, log=FALSE);
    data <- apply(data,2,as.integer);
    filt.dataphy <- otu_table(data,taxa_are_rows =TRUE);
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    data <- merge_phyloseq(filt.dataphy,sample_table);
    taxa_names(data) <- rownames(filt.data);
    if(generank!="Feature"){          
        dataSet$taxa_table<-tax_table(dataSet$proc.phyobj);
        data <- merge_phyloseq(data,dataSet$taxa_table);
    }else{
        data <- data;
    }
    claslbl <<- as.factor(sample_data(dataSet$proc.phyobj)[[variable]]);
    if(generank=="Feature"){
        data <- data;
        tax_orig <<- nm <- taxa_names(data);
        dat3t <- as.data.frame(t(otu_table(data)));
    }else{
        #merging at taxonomy levels
        data <- fast_tax_glom_first(data,generank);
        tax_orig <<- taxa_names(data);
        nm <- as.character(tax_table(data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        dat3t <- as.data.frame(t(t(sapply(by(data1,rownames(data1),colSums),identity))));
    }
    set.seed(56290);
    #KW rank sum test
    rawpvalues <- apply(dat3t,2,function(x) kruskal.test(x,claslbl)$p.value);
    ord.inx <- order(rawpvalues);
    rawpvalues <- rawpvalues[ord.inx];
    clapvalues <- p.adjust(rawpvalues, method ="BH");
    dat3t <- dat3t[,ord.inx];
    if(length(rawpvalues) > 500){
        rawpvalues <- rawpvalues[1:500];
        clapvalues<-clapvalues[1:500];
        dat3t <- dat3t[,1:500];
    };
    wil_data <- as.data.frame(dat3t);
    #if no subclass within classes then no wilcoxon rank sum test
    wil_datadf <- wil_data;

    #Linear Discriminant analysis (LDA)
    ldares <- lda(claslbl~ .,data = wil_datadf);
    ldamean <- data.frame(t(ldares$means));
    class_no <<- length(unique(claslbl));
    ldamean$max <- apply(ldamean[,1:class_no],1,max);
    ldamean$min <- apply(ldamean[,1:class_no],1,min);
    ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
    ldamean$Pvalues <- signif(rawpvalues,digits=5);
    ldamean$FDR <- signif(clapvalues,digits=5);
    resTable <- ldamean;
    
    # it seems lda add ` around names containing dash "-", need to strip this off
    rawNms <- rownames(resTable); 
    rownames(resTable) <- gsub("`", '', rawNms);
    de.Num <- sum(clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
    if(de.Num == 0){
        current.msg <<- "No significant features were identified with given criteria."; 
    }else{
        current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
    }
    # sort by p value
    ord.inx <- order(resTable$Pvalues, resTable$LDAscore);
    resTable <- resTable[ord.inx, ,drop=FALSE];
    #p-values column to appear first; then FDR and then others
    resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
    resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
    
    #only getting the names of DE features
    diff_ft <<- rownames(resTable)[1:de.Num];
    resTable$max <- resTable$min <- NULL;
    
    #space in column names have replaced with (.), so revert back to normal;
    colnames(resTable) <- gsub("\\."," ",colnames(resTable));
    
    #if only two groups are present in sample variable
    cls.lbls <- levels(claslbl);
    grp.dat <- resTable[, cls.lbls];
    res.cls <- as.character(apply(grp.dat, 1, function(x){cls.lbls[which.max(x)]}));
    if(class_no==2){
        indx<-which(res.cls==unique(claslbl)[1]);
        resTable$LDAscore[indx]<--abs(resTable$LDAscore[indx]);
    }
    write.csv(resTable, file="lefse_de_output.csv");
    analSet$lefse$resTable <- analSet$resTable <-resTable;
    
    #subset dataset for bar plot visualization (LDA Score)
    ldabar <- as.data.frame(rownames(resTable));
    ldabar[,2] <- resTable$LDAscore;
    ldabar[,3] <- res.cls;
    
    #visualizing top 25 features based on LDA score
    if (nrow(ldabar)>25) {
        ldabar <- ldabar[1:25,];
    }
    ldabar <<- ldabar;
    
    #preparing data for indvidual box plot
    sigfeat <<- rownames(resTable);
    generank <<- generank;
    box_data <- as.data.frame(wil_datadf[, sigfeat]);
    colnames(box_data) <- sigfeat;
    box_data$class <- claslbl;
    analSet$boxdata <- box_data;
    analSet$sig.count <- de.Num;
    analSet$anal.type <- "lefse";
    de.var <<- variable;
    analSet$lefse.genelvl<-generank;
    analSet <<- analSet;
    return(1);
}
################################
PlotLEfSeSummary <- function(imgName,format="png", dpi=72) {
    set.seed(280561493);
    ldabar <- ldabar;
    Cairo(file=imgName, width=600, height=560, type=format, bg="white",dpi=dpi);
    box = ggplot(ldabar, aes(x=reorder(ldabar[,1],ldabar[,2]), y=ldabar[,2], fill=ldabar[,3]))+ geom_bar(stat="identity",width=0.8) + coord_flip() + labs(y="LDA score",x="Features",fill="Class")+theme_bw()+scale_color_brewer(palette="Set1");
    print(box);
    dev.off();
}
###############################################################################

##ANCOM/ANCOR
####source: users.ugent.be/~shawinke/ABrokenPromise/03_diffAbundDetect.html
##########################################
PerformANCOM <- function(p.lvl,generank,variable){
    claslbl <- as.factor(sample_data(dataSet$norm.phyobj)[[variable]]);
    filt.dataphy <- dataSet$filt.data;
    filt.dataphy <- apply(filt.dataphy,2,as.integer);
    filt.dataphy <- otu_table(filt.dataphy,taxa_are_rows =TRUE);
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    filt.dataphy <- merge_phyloseq(filt.dataphy,sample_table);
    taxa_names(filt.dataphy) <- rownames(dataSet$filt.data);
    data <- filt.dataphy;
    dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(data,dataSet$taxa_table);  
    if(generank=="Feature"){
        data <- data;
    }else{
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1<-as.matrix(otu_table(data));
        rownames(data1)<- nm;
        
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
        data1 <- otu_table(data1,taxa_are_rows=T);
        data <- merge_phyloseq(data1,sample_data(data));
    }
    data <- t(data);
    data <- otu_table(data)@.Data;
    data1 = as.data.frame(cbind(data, group= sample_data(dataSet$norm.phyobj)[[variable]]));
    
    #add pseudo count before performing log ratio
    data1[data1 == 0] <- 0.001;
    
    #multcorr=2 corresponds to BH correction
    res <- ANCOM(data1,sig = p.lvl, multcorr = 2, ncore = 1);
    saveRDS(res,"ancom.rds");
    taxaDet = res$detected;
    
    #taxaDet = gsub(".TP","-TP",gsub("X","",res$detected));
    de.Num <- length(res$detected);
    if(taxaDet == "No significant features detected" & de.Num == 1){
        de.Num <- 0;
        current.msg <<- "No significant genes were identified using the given adjusted P value cutoff."; 
    } else{
         current.msg <<- paste("A total of", de.Num, "significant features were identified!");
    }
    resTable <- matrix(1, ncol=3, nrow=ncol(data));
    colnames(resTable) = c("rawP","adjP", "W-stats");
    rownames(resTable) = colnames(data);
    if(de.Num == 0){
        resTable<-resTable;
    }else {
        resTable[taxaDet,"adjP"] <- 1e-4;     
    }
    resTable[,"W-stats"] <- res$W;
    resTable <- as.data.frame(resTable);
    ord.inx <- order(resTable$adjP,-resTable[,"W-stats"]);
    resTable <- resTable[ord.inx, , drop=FALSE];
    feat.nm <- rownames(resTable); 
    resTable <-data.frame(W=resTable[,-c(1,2)]);
    rownames(resTable)<- feat.nm;
    write.csv(resTable, file="ANCOM_de.csv");

    if(nrow(resTable) > 500){
        resTable <- resTable[1:500, ];
    }
    analSet$ancom$resTable<-analSet$resTable<-resTable;
    
    #only getting the names of DE features
    diff_ft <<- rownames(resTable)[1:de.Num]; 
    generank <<- generank; 
    
    #individual boxplot for features
    sigfeat <- rownames(resTable);
    box_data <- as.data.frame(data1[ ,sigfeat]);
    colnames(box_data) <- sigfeat;    
    box_data$class <- claslbl;
    analSet$boxdata <- box_data;
    analSet$sig.count <- de.Num;
    analSet$ancom.genelvl<- generank;
    analSet$anal.type <- "ancom";
    de.var <<- variable;
    analSet <<- analSet;
    return(1);
}
###############################################

##ALDEx2 test (Fernandes et al., 2014)
####Source: users.ugent.be/~shawinke/ABrokenPromise/03_diffAbundDetect.html
##########################################
PerfromaldexTest = function(p.lvl,generank,variable,diffOpt,mc.samples){
    set.seed(5623009);    
    claslbl <- as.character(sample_data(dataSet$norm.phyobj)[[variable]]);
    filt.dataphy <- dataSet$filt.data;
    filt.dataphy <- apply(filt.dataphy,2,as.integer);
    filt.dataphy <- otu_table(filt.dataphy,taxa_are_rows =TRUE);
    sample_table <- sample_data(dataSet$proc.phyobj, errorIfNULL=TRUE);
    filt.dataphy <- merge_phyloseq(filt.dataphy,sample_table);
    taxa_names(filt.dataphy) <- rownames(dataSet$filt.data);
    data <- filt.dataphy;
    dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
    data <- merge_phyloseq(data,dataSet$taxa_table); 
    if(generank=="Feature"){
        data <- data;
    }else{
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        
        #all NA club together
        data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
        data1 <- otu_table(data1,taxa_are_rows=T);
        data <- merge_phyloseq(data1,sample_data(data));
    }
    if (!taxa_are_rows(data)){
        data <- t(data);
    }
    data <- data.frame(otu_table(data)@.Data);
    suppressMessages(require("ALDEx2"));
    test <-"t";
    if(length(unique(claslbl)) > 2){
        test <-"kw";
    }
    res <- aldex(data, conditions = claslbl,mc.samples=mc.samples, test = test);
    if (test=="t"){
        if (diffOpt=="both"){
            sigHits <- res$we.eBH<=p.lvl & res$wi.eBH<=p.lvl;
            de.Num <- length(which(sigHits));
            resTable = cbind(res$we.ep, res$we.eBH,res$wi.ep, res$wi.eBH);
            colnames(resTable) = c("we_Pvalues","we_adjPvalues","wi_Pvalues","wi_adjPvalues");
        } else {
            if (diffOpt=="t"){
                sigHits <- res$we.eBH<=p.lvl;
                de.Num <- length(which(sigHits));
                resTable = cbind(res$we.ep, res$we.eBH);

            }else {
                sigHits <- res$wi.eBH<=p.lvl;
                de.Num <- length(which(sigHits));
                resTable = cbind(res$wi.ep, res$wi.eBH);
            } 
            colnames(resTable) = c("Pvalues","adjPvalues");
        }
    } else {
        if (diffOpt=="both"){
            sigHits <- res$kw.eBH<=p.lvl & res$glm.eBH<=p.lvl;
            de.Num <- length(which(sigHits));
            resTable = cbind(res$kw.ep, res$kw.eBH,res$glm.ep, res$glm.eBH);
            colnames(resTable) = c("kw_Pvalues","kw_adjPvalues","glm_Pvalues","glm_adjPvalues");
        } else {
            if (diffOpt=="kw"){
                sigHits <- res$kw.eBH<=p.lvl;
                de.Num <- length(which(sigHits));
                resTable = cbind(res$kw.ep, res$kw.eBH);

            }else {
                sigHits <- res$glm.eBH<=p.lvl;
                de.Num <- length(which(sigHits));
                resTable = cbind(res$glm.ep, res$glm.eBH);
            } 
            colnames(resTable) = c("Pvalues","adjPvalues");
        }
    }
    if(de.Num == 0){
        current.msg <<- "No significant features were identified with given criteria."; 
    }else{
        current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
    }
    resTable <- as.data.frame(resTable);
    rownames(resTable) = rownames(res);
    if (diffOpt=="both"){
        if (test=="t"){
            ord.inx <- order(resTable$we_Pvalues);
        } else {
             ord.inx <- order(resTable$kw_Pvalues);
        }
    } else {
        ord.inx <- order(resTable$Pvalues);
    }
    resTable <- resTable[ord.inx, , drop=FALSE];
    write.csv(resTable, file="aldex2_de_output.csv");   
    if(nrow(resTable) > 500){
        resTable <- resTable[1:500, ];
    }
    analSet$aldex$resTable <- analSet$resTable <-resTable;
    
    #only getting the names of DE features
    diff_ft <<- rownames(resTable)[1:de.Num];
    sigfeat <- rownames(resTable);
    generank <<- generank; 
    
    #prepare individual boxplot
    box_data <- data[sigfeat, ];
    
    #samples in rows
    box_data <- t(box_data);
    box_data <- data.frame(box_data);  
    colnames(box_data) <- sigfeat;
    box_data$class <- claslbl;
    analSet$boxdata <- box_data;
    analSet$sig.count <- de.Num;
    de.var <<- variable;
    analSet$anal.type <- "aldex";
    analSet$aldex.genelvl <- generank;
    analSet <<- analSet;
    return(1);
}
###########################################
GetALDEXClassInfo <- function(variable){
    clslbl <- as.factor(sample_data(dataSet$norm.phyobj)[[variable]]);
    return(length(levels(clslbl)));
}
########################################################################

########### Random Forest #############
#######################################
# random forests
RF.Anal <- function(treeNum, tryNum, randomOn, variable, generank){
    suppressMessages(require(randomForest));

    # set up random numbers
    if(is.null(analSet$random.seeds)){
        analSet$random.seeds <- GetRandomNumbers();
        analSet$cur.inx <- 0;
        analSet$rn.seed <- analSet$random.seeds[1];
    }
    if(randomOn == -1){
        rn.sd <- 123456;
    }else if(randomOn == 0){ # keep current
        rn.sd <- analSet$rn.seed;
    }else{ # random on
        cur.inx <- analSet$cur.inx + 1;
        rn.sd <- analSet$random.seeds[cur.inx];        
        analSet$cur.inx <- cur.inx;
    }
    set.seed(rn.sd);
    # save the seed
    analSet$rn.seed <- rn.sd;
    analSet <<- analSet;       

    #merging at taxonomy levels
    if(generank =="Feature"){
        data <- dataSet$norm.phyobj
        data1 <- as.matrix(t(otu_table(data)));      
    } else{
        dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
        data <- merge_phyloseq(dataSet$norm.phyobj,tax_table(dataSet$proc.phyobj)); 
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        y <- which(is.na(nm)==TRUE);
        #converting NA values to unassigned
        nm[y] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        #all NA club together
        data1<-sapply(by(data1,rownames(data1),colSums),identity);
    }
    data.impfeat <<- data1;
    cls <- sample_data(dataSet$norm.phyobj)[[variable]];
    variable <<- variable;
    rf_out <- randomForest(data1,cls, ntree = treeNum, mtry = tryNum, importance = TRUE, proximity = TRUE);
    
    # set up named sig table for display
    impmat <- rf_out$importance;
    impmat <- impmat[rev(order(impmat[,"MeanDecreaseAccuracy"])),]
    sigmat <- impmat[,"MeanDecreaseAccuracy", drop=F];
    sigmat <- signif(sigmat, 5);
    write.csv(sigmat,file="randomforests_sigfeatures.csv");
    analSet$cls <- cls;
    analSet$rf <- rf_out;
    analSet$rf.sigmat <- sigmat;
    analSet <<- analSet;
}
###############################################################################
# plot variable importance ranked by MeanDecreaseAccuracy
PlotRF.Classify <- function(imgName,format="png", dpi=72, width=NA){
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 9;
    }else{
        w <- width;
    }
    h <- w*6/9;
    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    par(mar=c(4,4,3,2));
    cols <- rainbow(length(levels(analSet$cls))+1);
    plot(analSet$rf, main="Random Forest classification", col=cols);
    legend("topright", legend = c("Overall", levels(analSet$cls)), lty=2, lwd=1, col=cols);
    dev.off();
}
###############################################################################
# plot variable importance ranked by MeanDecreaseAccuracy
PlotRF.VIP <- function(imgName,format="png", dpi=72, width=NA){
    vip.score <- rev(sort(analSet$rf$importance[,"MeanDecreaseAccuracy"]));
    if(is.na(width)){
        w <- 9;
    }else if(width == 0){
        w <- 8;
    }else{
        w <- width;    
    }
    h <- w*7/8;
    Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    PlotImpVar(vip.score,"MeanDecreaseAccuracy");
    dev.off();
}
###############################################################################
# get the OOB error for the last signif
GetRFOOB <- function(){
    errors = analSet$rf$err.rate;
    nrow = dim(errors)[1];
    signif(errors[nrow, 1],3);
    print(signif(errors[nrow, 1],3));
}
###############################################################################
GetSigTable.RF <- function(){
    GetSigTable(analSet$rf.sigmat, "Random Forest");
}
###############################################################################
# significance measure, double[][]
GetRFSigMat <- function(){
    return(CleanNumber(analSet$rf.sigmat))
}
###############################################################################
GetRFSigRowNames <- function(){
    rownames(analSet$rf.sigmat);
}
###############################################################################
GetRFSigColNames <- function(){
    colnames(analSet$rf.sigmat);
}
###############################################################################
# return double[][] confusion matrix
GetRFConfMat <- function(){
	signif(analSet$rf$confusion,3);
}
###############################################################################
GetRFConfRowNames <- function(){
	rownames(analSet$rf$confusion);
}
###############################################################################
GetRFConfColNames <- function(){
	colnames(analSet$rf$confusion);
}
###############################################################################
PlotImpVar <- function(imp.vec, xlbl, feat.num=15, color.BW=FALSE){
    cls.len <- length(levels(analSet$cls));
    if(cls.len == 2){
        rt.mrg <- 5;
    }else if(cls.len == 3){
        rt.mrg <- 6;
    }else if(cls.len == 4){
        rt.mrg <- 7;
    }else if(cls.len == 5){
        rt.mrg <- 8;
    }else if(cls.len == 6){
        rt.mrg <- 9;
    }else{
        rt.mrg <- 11;
    }
    op <- par(mar=c(5,9,2,rt.mrg)); # set right side margin with the number of class
    if(feat.num <= 0){
        feat.num = 15;
    }
    if(feat.num > length(imp.vec)){
        feat.num <- length(imp.vec);
    }
    # first get the top subset
    imp.vec <- rev(sort(imp.vec))[1:feat.num];

    # reverser the order for display
    imp.vec <- sort(imp.vec);
    
    # as data should already be normalized, use mean/median should be the same
    # mns is a list contains means of all vars at each level
    # conver the list into a matrix with each row contains var averages across different lvls
    data1 <- as.data.frame(data.impfeat);
    mns <- by(data1[, names(imp.vec)], analSet$cls, 
                    function(x){ # inner function note, by send a subset of dataframe
                        apply(x, 2, mean, trim=0.1)
                    });
    mns <- t(matrix(unlist(mns), ncol=feat.num, byrow=TRUE));
    vip.nms <- names(imp.vec);
    
    # sometimes the name is superlong; so have to strip it
    vip.nms <- substr(vip.nms, 1, 20);
    names(imp.vec) <- NULL;

    # modified for B/W color
    dotcolor <- ifelse(color.BW, "darkgrey", "blue");
    dotchart(imp.vec, bg=dotcolor, xlab= xlbl, cex=1.3);
    mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1)
    axis.lims <- par("usr"); # x1, x2, y1 ,y2

    # get character width
    shift <- 2*par("cxy")[1];
    lgd.x <- axis.lims[2] + shift;
    x <- rep(lgd.x, feat.num);
    y <- 1:feat.num;
    par(xpd=T);
    suppressMessages(require(RColorBrewer));
    nc <- ncol(mns);

    # modified for B/W color
    colorpalette <- ifelse(color.BW, "Greys", "RdYlGn");
    col <- colorRampPalette(brewer.pal(10, colorpalette))(nc); # set colors for each class
    if(color.BW) col <- rev(col);

    # calculate background
    bg <- matrix("", nrow(mns), nc);
    for (m in 1:nrow(mns)){
        bg[m,] <- (col[nc:1])[rank(mns[m,])];
    }
    cls.lbl <- levels(analSet$cls);
    for (n in 1:ncol(mns)){
        points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
        # now add label
        text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
        # shift x, note, this is good for current size
        x <- x + shift/1.25;
    }
    # now add color key, padding with more intermediate colors for contiuous band
    col <- colorRampPalette(brewer.pal(25, colorpalette))(50)
    if(color.BW) col <- rev(col);
    nc <- length(col);
    x <- rep(x[1] + shift, nc);
    shifty <- (axis.lims[4]-axis.lims[3])/3;
    starty <- axis.lims[3] + shifty;
    endy <- axis.lims[3] + 2*shifty;
    y <- seq(from = starty, to = endy, length = nc);
    points(x,y, bty="n", pch=15, col=rev(col), cex=2);
    text(x[1], endy+shifty/8, "High");
    text(x[1], starty-shifty/8, "Low");
    par(op);
}
###############################################################################

###SVM (adopted from MetaboAnalyst R)####
##########################################################
RSVM.Anal <- function(variable,generank,cvType){
   set.seed(7575677);
   #merging at taxonomy levels
   if(generank =="Feature"){
        data <- dataSet$norm.phyobj;
        data1 <- as.data.frame(t(otu_table(data)));      
    } else {
        dataSet$taxa_table <- tax_table(dataSet$proc.phyobj);
        data <- merge_phyloseq(dataSet$norm.phyobj,tax_table(dataSet$proc.phyobj)); 
        data <- fast_tax_glom_first(data,generank);
        nm <- as.character(tax_table(data)[,generank]);
        y <- which(is.na(nm)==TRUE);
        
        #converting NA values to unassigned
        nm[y] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(data));
        rownames(data1) <- nm;
        
        #all NA club together
        data1 <- sapply(by(data1,rownames(data1),colSums),identity);
        data1 <- as.data.frame(data1);
    }
    data.impfeat <<- data1;
    cls <- sample_data(dataSet$norm.phyobj)[[variable]];
    variable <<- variable;  
    ladder = CreateLadder(ncol(data1));
    svm.out <- RSVM(data1, cls, ladder, CVtype=cvType);
    
    # calculate important features
    ERInd <- max(which(svm.out$Error == min(svm.out$Error)))
    MinLevel <- svm.out$ladder[ERInd];
    FreqVec <- svm.out$SelFreq[, ERInd];
    SelInd <- which(rank(FreqVec) >= (svm.out$ladder[1]-MinLevel));
    FreqInd <- svm.out$SelFreq[SelInd, ERInd]
    names(FreqInd) <- names(data1)[SelInd];
    #create a sig table for display
    sig.var <- rev(sort(FreqInd));
    sig.var <- as.matrix(sig.var); # 1-column matrix
    colnames(sig.var) <- "Freqency";
    write.csv(sig.var, file="svm_sigfeatures.csv");
    # add sorted features frequencies as importance indicator
    svm.out <- append(svm.out, list(sig.mat=sig.var, best.inx=ERInd));
    analSet$cls <- cls;
    analSet$svm <- svm.out;
    analSet <<- analSet;
}
###############################################################################
PlotRSVM.Classify <- function(imgName, format="png", dpi=72, width=NA){
  res <- analSet$svm$Error;
  edge <- (max(res)-min(res))/100; # expand y uplimit for text
  if(is.na(width)){
    w <- 8;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w*6/8;  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  plot(res,type='l',xlab='Number of variables (levels)',ylab='Error Rate',
       ylim = c(min(res)-5*edge, max(res)+18*edge), axes=F,
       main="Recursive SVM classification")
  text(res,labels =paste(100*round(res,3),'%'), adj=c(-0.3, -0.5), srt=45, xpd=T)
  points(res, col=ifelse(1:length(res)==analSet$svm$best.inx,"red","blue"));
  axis(2);
  axis(1, 1:length(res), names(res));
  dev.off();
}
###############################################################################
PlotRSVM.Cmpd <- function(imgName, format="png", dpi=72, width=NA){  
  sigs <- analSet$svm$sig.mat;
  data <- sigs[,1];
  if(is.na(width)){
    w <- 8;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w*7/8;  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  PlotImpVar(data, "Frequency");
  dev.off();
}
###############################################################################
CreateLadder <- function(Ntotal, Nmin=5){
  x <- vector()
  x[1] <- Ntotal
  # note SVM is very computationally intensive, large step first 
  # first descend with 0.5 -> 50 var left
  # then descend with 0.6 -> 25 var left
  # then desend with 0.75 -> 5 var
  for( i in 1:100 ){
    if(x[i]>200){
      pRatio = 0.4
    }else if(x[i]>50){
      pRatio = 0.5
    }else if(x[i]>25){
      pRatio = 0.6
    }else{
      pRatio = 0.75
    }
    pp <- round(x[i] * pRatio)
    if( pp == x[i] ){
      pp <- pp-1
    }
    if( pp >= Nmin ) {
      x[i+1] <- pp
    } else{
      break
    }
  }
  x
}
###############################################################################
RSVM <- function(x, y, ladder, CVtype, CVnum=0){
    ## check if y is binary response
    Ytype <- names(table(y))
    if(length(Ytype) != 2){
    print("ERROR!! RSVM can only deal with 2-class problem")
    return(0)
    }
    ## class mean
    m1 <- apply(x[ which(y==Ytype[1]), ], 2, mean)
    m2 <- apply(x[ which(y==Ytype[2]), ], 2, mean)
    md <- m1-m2

    yy <- vector(length=length(y))
    yy[which(y==Ytype[1])] <- 1
    yy[which(y==Ytype[2])] <- -1
    y <- yy

    ## check ladder
    if(min(diff(ladder)) >= 0){
        print("ERROR!! ladder must be monotonously decreasing")
        return(0);
    }
    if(ladder[1] != ncol(x) ){
        ladder <- c(ncol(x), ladder)
    }
    nSample <- nrow(x)
    nGene <- ncol(x)
    SampInd <- seq(1, nSample)
    if(CVtype == "LOO"){
        CVnum <- nSample
    }else{
    if(CVnum == 0 ){
      CVnum <- nSample
    }
    }
    ## vector for test error and number of tests
    ErrVec <- vector(length=length(ladder))
    names(ErrVec) <- as.character(ladder);
    nTests <- 0
    SelFreq <- matrix(0, nrow=nGene, ncol=length(ladder))
    colnames(SelFreq) <- paste("Level", ladder);

    ## for each CV
    for(i in 1:CVnum){
    ## split data
    if(CVtype == "LOO"){
      TestInd <- i
      TrainInd <- SampInd[ -TestInd]
    }else{
      if(CVtype == "bootstrape"){
        TrainInd <- sample(SampInd, nSample, replace=T);
        TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))];
      }else{
        ## Nfold
        CVtype <- as.integer(CVtype);
        TrainInd <- sample(SampInd, nSample*(CVtype-1)/CVtype);
        TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))];
      }
    }
    nTests <- nTests + length(TestInd)

    ## in each level, train a SVM model and record test error
    xTrain <- x[TrainInd, ]
    yTrain <- y[TrainInd]
    xTest  <- x[TestInd,]
    yTest  <- y[TestInd]

    ## index of the genes used in the
    SelInd <- seq(1, nGene)
    for(gLevel in 1:length(ladder))
    {
      ## record the genes selected in this ladder
      SelFreq[SelInd, gLevel] <- SelFreq[SelInd, gLevel] +1

      ## train SVM model and test error
      ###################################################################################
      ## note the scale is changed to T or it never returns sometime for unscaled data ###
      ## note: the classification performance is idenpendent of about scale is T/F  #####
      ## for "LOO", the test data should be as.data.frame, matrxi will trigger error #####
      ###################################################################################
      svmres <- e1071::svm(xTrain[, SelInd], yTrain, scale=T, type="C-classification", kernel="linear" )
      if( CVtype == "LOO" ){
        svmpred <- predict(svmres, as.data.frame(xTest[SelInd], nrow=1) )
      }else{
        svmpred <- predict(svmres, xTest[, SelInd] )
      }
      ErrVec[gLevel] <- ErrVec[gLevel] + sum(svmpred != yTest )

      ## weight vector
      W <- t(svmres$coefs*yTrain[svmres$index]) %*% svmres$SV * md[SelInd]
      rkW <- rank(W)

      if( gLevel < length(ladder) ){
        SelInd <- SelInd[which(rkW > (ladder[gLevel] - ladder[gLevel+1]))]
      }
    }
    }
    ret <- list(ladder=ladder, Error=ErrVec/nTests, SelFreq=SelFreq);
    ret;
}
###########################################################################

###############################################################################