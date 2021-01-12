########################################################################
## R script for ResistoXplorer
## Description: Data/resource management functions
########################################################################

#######Procrustes Analysis########################
##########################################
PerformProcrustes <- function(imgNm,taxrank,generank,ordmeth,distName,colOpt,colPal,variable,showlabel,format="png",dpi=72){
    set.seed(123455);
    suppressMessages(require(vegan));
    suppressMessages(require(ade4));
    
  #normalized data
    m.data <- dataSet$m.norm.phyobj;
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    m.data<-merge_phyloseq(m.data,dataSet$m.taxa_table);
    if(taxrank=="Feature"){
        m.data <- m.data;
        m.data.norm <- t(as.matrix(otu_table(m.data)));
    }else{
        #merging at taxonomy levels
        m.data <- fast_tax_glom_first(m.data,taxrank);
        nm <- as.character(tax_table(m.data)[,taxrank]);
        
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1<-as.matrix(otu_table(m.data));
        rownames(data1) <- nm;
        m.data.norm <- as.matrix(t(t(sapply(by(data1,rownames(data1),colSums),identity))));
    }
    r.data <- dataSet$norm.phyobj;
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.data <- merge_phyloseq(r.data,dataSet$taxa_table);
    if(generank=="Feature"){
        r.data <- r.data;
        r.data.norm <- t(as.matrix(otu_table(r.data)));
    }else{
        #merging at taxonomy levels
        r.data <- fast_tax_glom_first(r.data,generank);
        nm <- as.character(tax_table(r.data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(r.data));
        rownames(data1) <- nm;
        r.data.norm <- as.matrix(t(t(sapply(by(data1,rownames(data1),colSums),identity))));
    }
    #2 {perform ordination analysis vegan based on ordmeth and distName}
    if(ordmeth=="PCA"){
      #need to transform the microbiome data in order to apply PCA
        #m.hel <- decostand(m.data.norm, method = "hellinger");
        m.ord <- rda(m.data.norm);
        #r.hel <- decostand(r.data.norm, method = "hellinger");
        r.ord <- rda(r.data.norm);
    } else if(ordmeth=="PCoA"){
        m.dis <- vegdist(m.data.norm, method = distName);
        add <- !(is.euclid(m.dis));
        m.ord <- cmdscale(m.dis, k = nrow(m.data.norm)-1, eig = TRUE, add = add);
        r.dis <- vegdist(r.data.norm, method = distName);
        add <-  !(is.euclid(r.dis));
        r.ord <- cmdscale(r.dis, k = nrow(r.data.norm)-1, eig = TRUE, add = add);    
    }else{
        m.ord <- metaMDS(m.data.norm,k=3,binomial=TRUE,distance= distName);
        r.ord <- metaMDS(r.data.norm,k=3,binomial=TRUE,distance= distName);
    }
    #only first two axis are used to avoid warning messages
    proc.res <- procrustes(r.ord,m.ord,symmetric = FALSE, scores = "sites",choices = c(1,2,3));
    proc.res.data <- proc.res;
    Cairo(file=imgNm, width=720, height=500, type=format, bg="white",dpi=dpi);
    
    #make sure the sample names are in same order
    m_sam_nm <- sample_names(dataSet$m.norm.phyobj);
    suppressMessages(require(RColorBrewer));
    colnames(proc.res.data$Yrot) <- colnames(proc.res.data$X);
    proc.res.data$Yrot <- proc.res.data$Yrot[m_sam_nm,];
    proc.res.data$X <- proc.res.data$X[m_sam_nm,];
    proc.res.data$X <- as.data.frame(proc.res.data$X);
    proc.res.data$Yrot <- as.data.frame(proc.res.data$Yrot);
    proc.res.data$Yrot$Data <- rep("microbiome",nrow(proc.res.data$Yrot));
    proc.res.data$X$Data <- rep("resistome",nrow(proc.res.data$X));
    proc.res.data$X$Sample <- rownames(proc.res.data$X);
    proc.res.data$Yrot$Sample <- rownames(proc.res.data$Yrot);
    if (colOpt!="omics"){        
        proc.res.data$X[,variable] <- sample_data(dataSet$norm.phyobj)[m_sam_nm][[variable]]; 
        proc.res.data$Yrot[,variable] <- sample_data(dataSet$m.norm.phyobj)[m_sam_nm][[variable]];
    }
    proc.res.data <- rbind(proc.res.data$X,proc.res.data$Yrot);
    proc.score <- signif(proc.res.data[,1:3],5);
    proc.score$Data <- proc.res.data$Data;
    write.csv(proc.score,row.names=proc.res.data[,"Sample"], file="procrustes_ordination_score.csv");
    proc.score <- NULL;
    suppressMessages(require(dplyr));
    if (colOpt=="omics"){
      box <- proc.res.data %>% ggplot(aes(proc.res.data[,1],proc.res.data[,2], color=Data,shape=Data))+
            geom_point(aes(fill=Data),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2");
    } else {
        grp <- proc.res.data[,variable];
        box <- proc.res.data %>% ggplot(aes(proc.res.data[,1],proc.res.data[,2],color=grp,shape=Data))+
            geom_point(aes(fill=grp),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2",fill=variable,color=variable);
    }
    box <- box+geom_line(aes(group = Sample),color="grey")+theme_bw()+scale_shape_manual(values=c(22,21));
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    if(showlabel=="samnm"){
       box <- box+geom_text(aes(label=proc.res.data[,"Sample"]),hjust=0.5, vjust=2,size=3,fontface="bold");
    } else if(showlabel!="none") {
       smplnm <- proc.res.data[,"Sample"]
       clslbl <- sample_data(dataSet$m.norm.phyobj)[smplnm][[showlabel]];
       box<-box+geom_text(aes(label=clslbl),hjust=0.5, vjust=2,size=3,fontface="bold");
    }
    box <- box+ theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);  
    m.ord <<-m.ord;
    r.ord <<- r.ord;
    ordmeth <<- ordmeth;
    m_sam_nm <<- m_sam_nm;
    dataSet$merge.data <- proc.res.data;
    variable <<-variable;
    #m.coords <<- proc.res$Yrot[,1:3];
    #r.coords <<- proc.res$X[,1:3];
    dataSet <<- dataSet;
    dev.off();
}
###############################################################################
Plot3DIDAScore <- function(jsonNm,ordmeth,colOpt,variable,format="json"){
    pca3d <- list();
    pca3d$score$axis <- paste(ordmeth, 1:3 , sep="");
    # "m_" have been just added internal for rbind the coordinates
    data <- dataSet$merge.data;
    #if (ncol(coords)==2){ #case in rCCA
        #coords<-cbind(coords,Z=rep(0,nrow(coords)));
    #}
    coords <- data.frame(t(signif(data[,c(1,2,3)], 5)));
    #sample_names<-colnames(coords);
    sample_names <- data$Sample;
    colnames(coords) <- NULL; 
    pca3d$score$xyz <- coords;
    if(colOpt=="omics"){
        cls <- as.character(data[,"Data"]);
        clsb <- as.character(data[,"Data"]);
    }else {
        cls <- as.character(data[,variable]);
        clsb <- as.character(data[["Data"]]);
        # now set color for each group   
    }
    # connecting samples based on thier names
    clsc <- as.character(data[,"Sample"]);
    pca3d$score$name <- clsc;
    cols <- unique(as.numeric(cls)) + 1;
    pca3d$score$facA <- cls;
    pca3d$score$facB <- clsb;
    pca3d$score$facC <- clsc;
    variable<<-variable;
    
    # now set color for each group
    grp.num <- length(unique(cls)) + 1;
    cols <- 1:grp.num + 1;
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
PerformprocrustesProtest <- function(nperm){
    #input: 2 matrices (ordination scores)
    #perform protest
    protest.res <- protest(m.ord,r.ord,permutations = nperm,scores = "sites");
    pro.stat.info <- paste("[Procrustes] Sum of Squares = ", signif(protest.res$ss, 4), "; Correlation coefficient (squared m12) = ", signif(protest.res$scale, 4), "; P-value < ", signif(protest.res$signif, 5), sep="");
    dataSet$pro.stat.info <- pro.stat.info;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################
#######Coinertia Analysis########################
##########################################
PerformCIA <- function(imgNm,taxrank,generank,ordmeth,distName,colOpt,colPal,variable,showlabel,format="png",dpi=72){
    set.seed(123458);
    suppressMessages(require (ade4));
    #normalized data  
    m.data <- dataSet$m.norm.phyobj;
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    m.data <- merge_phyloseq(m.data,dataSet$m.taxa_table);
    if(taxrank=="Feature"){
        m.data <- m.data;
        m.data.norm <- t(as.matrix(otu_table(m.data)));
    }else{
        #merging at taxonomy levels
        m.data <- fast_tax_glom_first(m.data,taxrank);
        nm <- as.character(tax_table(m.data)[,taxrank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(m.data));
        rownames(data1) <- nm;
        m.data.norm <- as.matrix(t(t(sapply(by(data1,rownames(data1),colSums),identity))));
    }
    r.data <- dataSet$norm.phyobj;
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.data <- merge_phyloseq(r.data,dataSet$taxa_table);
    if(generank=="Feature"){
        r.data<-r.data;
        r.data.norm <- t(as.matrix(otu_table(r.data)));
    }else{
        #merging at taxonomy levels
        r.data <- fast_tax_glom_first(r.data,generank);
        nm <- as.character(tax_table(r.data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned"; 
        data1 <- as.matrix(otu_table(r.data));
        rownames(data1) <- nm;
        r.data.norm <- as.matrix(t(t(sapply(by(data1,rownames(data1),colSums),identity))));
    }
    #2 {perform ordination analysis vegan based on ordmeth and distName}
    if(ordmeth=="PCA"){
      #need to transform the microbiome data in order to apply PCA
        mcoa <- dudi.pca(m.data.norm, scannf = FALSE, nf = 3);
        rcoa <- dudi.pca(r.data.norm, scannf = FALSE, nf = 3);
    } else{
        suppressMessages(require(vegan));
        #have to provide distance matrix as input to dudi.pca
        m.dis <- vegdist(m.data.norm, method = distName);
        mcoa <- dudi.pco(m.dis, scannf = FALSE, nf = 2);
        r.dis <- vegdist(r.data.norm, method = distName);
        rcoa <- dudi.pco(r.dis, scannf = FALSE, nf = 2);    
    }
    #only first two axis are used to avoid warning messages
    cia.res <- coinertia(rcoa,mcoa, scannf = F, nf = 3);
    Cairo(file=imgNm, width=720, height=500, type=format, bg="white",dpi=dpi);
    cia.res.data <- cia.res;
    suppressMessages(require(RColorBrewer));      
    m_sam_nm<-sample_names(dataSet$m.norm.phyobj);
    cia.res.data$mY <- cia.res.data$mY[m_sam_nm,];
    cia.res.data$mX <- cia.res.data$mX[m_sam_nm,];
    cia.res.data$mX <- as.data.frame(cia.res.data$mX);
    cia.res.data$mY <- as.data.frame(cia.res.data$mY);
    cia.res.data$mY$Data <- rep("microbiome",nrow(cia.res.data$mY));
    cia.res.data$mX$Data <- rep("resistome",nrow(cia.res.data$mX));
    cia.res.data$mX$Sample <- rownames(cia.res.data$mX);
    cia.res.data$mY$Sample <- rownames(cia.res.data$mY);
    if(colOpt!="omics"){        
        cia.res.data$mX[,variable] <- sample_data(dataSet$norm.phyobj)[m_sam_nm][[variable]]; 
        cia.res.data$mY[,variable] <- sample_data(dataSet$m.norm.phyobj)[m_sam_nm][[variable]];
    }
    cia.res.data <- rbind(cia.res.data$mX,cia.res.data$mY);
    cia.score <- signif(cia.res.data[,1:3],5);
    cia.score$Data <- cia.res.data$Data;
    write.csv(cia.score,row.names=cia.res.data[,"Sample"], file="coinertia_ordination_score.csv");
    cia.score <- NULL;
    suppressMessages(require(dplyr));
    if(colOpt=="omics"){
      box <- cia.res.data %>% ggplot(aes(cia.res.data[,1],cia.res.data[,2], color=Data,shape=Data))+
            geom_point(aes(fill=Data),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2");
    } else {
        grp <- cia.res.data[,variable];
        box <- cia.res.data %>% ggplot(aes(cia.res.data[,1],cia.res.data[,2],color=grp,shape=Data))+
            geom_point(aes(fill=grp),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2",fill=variable,color=variable);
    }
    box <- box+geom_line(aes(group = Sample),color="grey")+theme_bw()+scale_shape_manual(values=c(22,21));
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    if(showlabel=="samnm"){
       box<-box+geom_text(aes(label=cia.res.data[,"Sample"]),hjust=0.5, vjust=2,size=3,fontface="bold");
    } else if(showlabel!="none") {
       smplnm <- cia.res.data[,"Sample"]
       clslbl <- sample_data(dataSet$m.norm.phyobj)[smplnm][[showlabel]];
       box<-box+geom_text(aes(label=clslbl),hjust=0.5, vjust=2,size=3,fontface="bold");
    }
    box <- box+ theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    dataSet$mcoa <- mcoa;
    dataSet$rcoa <- rcoa;
    ordmeth <<- ordmeth;
    m_sam_nm <<- m_sam_nm;
    dataSet$merge.data<- cia.res.data;
    variable <<-variable;
    #m.coords <<- cia.res$mY[,1:3];
    #r.coords <<- cia.res$mX[,1:3];
    dataSet <<- dataSet;
}
###############################################################################
PerformCIARVtest <- function(nperm){
    dataSet$rcoa$tab <- dataSet$rcoa$tab[rownames(dataSet$mcoa$tab),];
    rv.sig.res <-RV.rtest(dataSet$mcoa$tab,dataSet$rcoa$tab,nperm);
    cia.stat.info <- paste("[Coinertia] RV coefficient = ", signif(rv.sig.res$obs, 5), "; P-value = ", signif(rv.sig.res$pvalue, 5), sep="");
    dataSet$cia.stat.info <- cia.stat.info;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################
######Pairwise correlation##
################################
PairwiseIntegrativeCorrelation <- function(taxrank,generank,cor.meth,cor.threshold,p.lvl,p.adj.method){
    #using only filtered data 
    m.data <- otu_table(dataSet$m.data.prefilt,taxa_are_rows =TRUE);
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    m.data <- merge_phyloseq(m.data,dataSet$m.taxa_table,sample_table);
    
    #perform correlation on relative abundance
    m.data <- transform_sample_counts(m.data, function(OTU) OTU/sum(OTU));
    if(taxrank=="Feature"){
        m.data <- m.data;
        m.data <- as.matrix(otu_table(m.data));
        nm <- substr(rownames(m.data), 1,20);
        rownames(m.data) <- nm;
    }else{
        #merging at taxonomy levels
        m.data<-fast_tax_glom_first(m.data,taxrank);
        nm<-as.character(tax_table(m.data)[,taxrank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_taxa"; 
        data1 <- as.matrix(otu_table(m.data));
        rownames(data1) <- nm;
        m.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    r.data <- otu_table(dataSet$data.prefilt,taxa_are_rows =TRUE);
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.data <- merge_phyloseq(r.data,dataSet$taxa_table,sample_table);
    #perform correlation of relative abundance
    r.data <- transform_sample_counts(r.data, function(OTU) OTU/sum(OTU));
    if(generank=="Feature"){
        r.data <- r.data;
        r.data <- as.matrix(otu_table(r.data));
        nm <- substr(rownames(r.data), 1,20);
        rownames(r.data) <- nm;
    }else{
        #merging at taxonomy levels
        r.data <- fast_tax_glom_first(r.data,generank);
        nm <- as.character(tax_table(r.data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_ARG"; 
        data1 <- as.matrix(otu_table(r.data));
        rownames(data1) <- nm;
        r.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    #ARGs and species that were not found in more than half of samples, to alleviate the bias from potential joint-ranking of zero values by Spearman’s rank
    #(artificial association bias (Li et al., 2015)
    r.dataFilt <- (r.data[rowSums(r.data != 0) > ncol(r.data)/2,]);
    m.dataFilt <- (m.data[rowSums(m.data != 0) > ncol(m.data)/2,]); 
    r.featnm <- rownames(r.dataFilt);
    #additional filtering can be done
    #if(ncol(data) > 1000){
      #filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
      #rk <- rank(-filter.val, ties.method='random');
      #data <- as.data.frame(data[,rk <=1000]);
    #}
    #the sample names are same in both the datsets
    data <- t(rbind(r.dataFilt, m.dataFilt));
    res <- associate(data,NULL,cor.meth,p.lvl,p.adj.method,mode="table",
            filter.self.correlations=TRUE,cth=NULL, order=FALSE, n.signif=0); 
    colnames(res) <- c("Features1", "Features2", "Correlation", "Adj.P.value");        
    #simplify the results and remove duplicates
    uid <- unique(unlist(res[c("Features1", "Features2")], use.names=FALSE));
    swap <- match(res[["Features1"]], uid) > match(res[["Features2"]], uid);
    idx = !duplicated(data.frame(
            V1 = ifelse(swap, res[["Features2"]], res[["Features1"]]),
            V2 = ifelse(swap, res[["Features1"]], res[["Features2"]])));
    res <- res[idx, , drop=FALSE]
    sigHits <-res$Adj.P.value<=p.lvl;
    de.Num <- length(which(sigHits)); 
    if(de.Num == 0){
        current.msg <<- "No significant correlation were identified with given adjusted P-value cutoff.";
        return(0);
    } else {
        sig_df <- res[res$Adj.P.value <= p.lvl,];
        sigcor <- abs(sig_df$Correlation) >= cor.threshold; 
        de.Numc <- length(which(sigcor)); 
        if(de.Numc == 0){
            current.msg <<- "No significant correlation were identified with given correlation coeffecient cutoff.";
            return(0);
        }else {
            corr.res <- sig_df[abs(sig_df$Correlation) >= cor.threshold,];
            corr.res[,3] <- round(corr.res[,3], digits=5);
            #corr.res[,4] <- round(corr.res[,4], digits=5);
            if (cor.meth=="spearman"){
                write.csv(corr.res, file="spearman_correlation.csv"); 
                saveRDS(corr.res,"spearman_correlation_net.rds");
            } else {
                write.csv(corr.res, file="pearson_correlation.csv"); 
                saveRDS(corr.res,"pearson_correlation_net.rds");
            }    
        }
    }
    dataSet$r.featnm <- r.featnm;
    dataSet$assoc.msg <- current.msg;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################
PerformCClassoCorrelation  <- function(taxrank,generank,cor.threshold,p.lvl,k_max,n_boot){
    #using only filtered data 
    m.data <- otu_table(dataSet$m.data.prefilt,taxa_are_rows =TRUE);
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    m.data <- merge_phyloseq(m.data,dataSet$m.taxa_table,sample_table);
    #perform correlation on relative abundance
    if(taxrank == "Feature"){
        m.data <- m.data;
        m.data <- as.matrix(otu_table(m.data));
        nm <- substr(rownames(m.data), 1,20);
        rownames(m.data) <- nm;
    }else{
        #merging at taxonomy levels
        m.data <- fast_tax_glom_first(m.data,taxrank);
        nm <- as.character(tax_table(m.data)[,taxrank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_taxa"; 
        data1 <- as.matrix(otu_table(m.data));
        rownames(data1) <- nm;
        m.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    #add pseudo count
    m.data[m.data == 0] <- 0.05;
    m.data <- sweep(m.data,2, colSums(m.data), FUN="/");
    r.data <- otu_table(dataSet$data.prefilt,taxa_are_rows =TRUE);
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.data <- merge_phyloseq(r.data,dataSet$taxa_table,sample_table);
    if(generank=="Feature"){
        r.data <- r.data;
        r.data <- as.matrix(otu_table(r.data));
        nm <- substr(rownames(r.data), 1,20);
        rownames(r.data) <- nm;
    }else{
        #merging at taxonomy levels
        r.data <- fast_tax_glom_first(r.data,generank);
        nm <- as.character(tax_table(r.data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_ARG"; 
        data1 <- as.matrix(otu_table(r.data));
        rownames(data1) <- nm;
        r.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    #add pseudo count
    r.data[r.data == 0] <- 0.05;
    r.data <- sweep(r.data,2, colSums(r.data), FUN="/");
    r.featnm <- rownames(r.data);
    #the sample names are same in both the datasets
    data <- t(rbind(r.data,m.data));
    featnm <- colnames(data);
    #CCLasso
    res <- cclasso(data,k_max = k_max, n_boot = n_boot);
    p.mat <-res$p_vals; #p-values
    corr.mat <- res$cor_w; #correlation value
    rownames(corr.mat)<-colnames(corr.mat)<-rownames(p.mat)<-colnames(p.mat) <-featnm;
    
    #there can be sometimes NA's in cor.mat and p.mat
    na.inx <- which(is.na(corr.mat));
    if (length(na.inx)> 0){
        p.mat[is.na(p.mat)]<-1;#replace NA in  p.mat to 1
        corr.mat[is.na(corr.mat)]<-0; #replace NA in corr.mat to 0
    }
    res <- flattenCorrMatrix(corr.mat,p.mat);
    sigHits <- res$P.value<=p.lvl;
    de.Num <- length(which(sigHits)); 
    if(de.Num == 0){
        current.msg <<- "No significant correlation were identified with given P-value cutoff.";
        return(0);
    } else {
        sig_df <- res[res$P.value <= p.lvl,];
        sigcor <- abs(sig_df$Correlation) >= cor.threshold; 
        de.Numc <- length(which(sigcor)); 
        if(de.Numc == 0){
            current.msg <<- "No significant correlation were identified with given correlation coeffecient cutoff.";
            return(0);
        }else {
            cor.res_filt <- sig_df[abs(sig_df$Correlation) >= cor.threshold,];
            cor.res_filt[,3] <- round(cor.res_filt[,3], digits=5);
            #cor.res_filt[,4] <- round(cor.res_filt[,4], digits=5);
            saveRDS(cor.res_filt,"cclasso_correlation_net.rds");
            write.csv(cor.res_filt, file="cclasso_correlation.csv"); 
        }
    }
    dataSet$r.featnm <- r.featnm;
    dataSet$cclsso.msg <- current.msg;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    Features1 = rownames(cormat)[row(cormat)[ut]],
    Features2 = rownames(cormat)[col(cormat)[ut]],
    Correlation  =(cormat)[ut],
    P.value = pmat[ut]
    );
}
###############################################################################
PerformMICorr <- function(taxrank,generank,cor.threshold,p.lvl,padj.meth,nperm){
    suppressMessages(require("minerva"));
    #using only filtered data 
    m.data <- otu_table(dataSet$m.data.prefilt,taxa_are_rows =TRUE);
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    m.data <- merge_phyloseq(m.data,dataSet$m.taxa_table,sample_table);
    #perform correlation on relative abundance
    if(taxrank=="Feature"){
        m.data <- m.data;
        m.data <- as.matrix(otu_table(m.data));
        nm <- substr(rownames(m.data), 1,20);
        rownames(m.data) <- nm;
    }else{
        #merging at taxonomy levels
        m.data <- fast_tax_glom_first(m.data,taxrank);
        nm <- as.character(tax_table(m.data)[,taxrank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_taxa"; 
        data1 <- as.matrix(otu_table(m.data));
        rownames(data1) <- nm;
        m.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    #add pseudo count
    #m.data[m.data == 0] <- 0.05;
    m.data <- sweep(m.data,2, colSums(m.data), FUN="/");
    r.data <- otu_table(dataSet$data.prefilt,taxa_are_rows =TRUE);
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.data <- merge_phyloseq(r.data,dataSet$taxa_table,sample_table);
    if(generank=="Feature"){
        r.data <- r.data;
        r.data <- as.matrix(otu_table(r.data));
        nm <- substr(rownames(r.data), 1,20);
        rownames(r.data)<-nm;
    }else{
        #merging at taxonomy levels
        r.data <- fast_tax_glom_first(r.data,generank);
        nm <- as.character(tax_table(r.data)[,generank]);
        #converting NA values to unassigned
        nm[is.na(nm)] <- "Unassigned_ARG"; 
        data1<-as.matrix(otu_table(r.data));
        rownames(data1) <- nm;
        r.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
    #add pseudo count
    #r.data[r.data == 0] <- 0.05;
    r.data <- sweep(r.data,2, colSums(r.data), FUN="/");
    r.featnm <- rownames(r.data);
    #the sample names are same in both the datasets
    data <- t(rbind(r.data,m.data));
    #tic measure and p-value
    tic.res <- mictoolsc(data, nperm=nperm,p.adjust.method = padj.meth);
    #calculate MIC(correlation coeff.); use adj-p value cutoff
    res <- mic_strength1(data, tic.res$pval, pval.col=c(6,4, 5),pthr=p.lvl); #MIC coeffecient for only significant pairs
    res <- res [,c(3,4,2,1)];
    colnames(res) <- c("Features1", "Features2", "Correlation", "Adj.P.value");        
    de.Num<-length(nrow(res)); 
    if(de.Num == 0){
        current.msg <<- "No significant correlation were identified with given adjusted P-value cutoff.";
        return(0);
    } else {
        #sig_df <- res[res$Adj.P.value <= p.lvl,];#already filtered
        sigcor <- abs(res$Correlation) >= cor.threshold; 
        de.Numc <- length(which(sigcor)); 
        if(de.Numc == 0){
            current.msg <<- "No significant correlation were identified with given correlation coeffecient cutoff.";
            return(0);
        }else {
            corr.res <- res[abs(res$Correlation) >= cor.threshold,];
            corr.res[,3] <- round(corr.res[,3], digits=5);
            #corr.res[,4] <- round(corr.res[,4], digits=5);
            write.csv(corr.res, file="mic_correlation.csv"); 
            saveRDS(corr.res,"mic_correlation_net.rds");
        }    
    }
    dataSet$r.featnm <- r.featnm;
    dataSet$mic.msg <- current.msg;
    dataSet <<- dataSet;
    return(1);    
}
###############################################################################
PlotCorrelationNet <- function(layoutOpt,nodeszOpt,cor.meth){
    if (cor.meth == "spearman"){
        cor.res <- readRDS("spearman_correlation_net.rds");
    } else if (cor.meth=="pearson") {
        cor.res <- readRDS("pearson_correlation_net.rds");
    } else if (cor.meth=="mic") {
         cor.res <- readRDS("mic_correlation_net.rds");
    }else {
         cor.res <- readRDS("cclasso_correlation_net.rds");
    }
    edge.list <- cor.res;
    edge.list = edge.list[,c("Features1", "Features2"), drop=FALSE];
    suppressMessages(require("igraph"));
    g = graph_from_data_frame(edge.list, directed = FALSE, vertices = NULL);
    E(g)$weight = abs(cor.res[, "Correlation"]);
    E(g)$correlation = cor.res[, "Correlation"];
    corr_matrix = get.adjacency(g,sparse=FALSE);
    colnames(edge.list) = c("Source", "Target")
    nodes = unique(c(edge.list[,1], edge.list[,2]));
    node.list = data.frame(Id=nodes, Name=nodes);    
    node.dgr <- as.numeric(degree(g));
    node.btw <- as.numeric(betweenness(g));
    
    # node size to degree
    if(vcount(g)>500){
        min.size = 1;
    }else if(vcount(g)>200){
        min.size = 2;
    }else{
        min.size = 2;
    }
    if (nodeszOpt=="dgr"){
        node.sizes <- node.dgr;
    } else {
        node.sizes <- node.btw;
    }
    node.sizes <- as.numeric(rescale2NewRange(node.sizes, min.size, 8));
    nms <- V(g)$name;
    #shapes and color would be based on wether the name is a microbe or resistome
    shapes <- rep("circle", length(nms));
    mat.inx <- nms %in% dataSet$r.featnm;
    shapes[mat.inx] <- "square";
    colors <- rep("#77579b", length(nms));
    colors[mat.inx]<-"#ffcc00";
    edge.width <- as.numeric(rescale2NewRange(E(g)$weight,1, 10));
    #edges info
    edge.mat <- get.edgelist(g);
    edge.mat1 = data.frame(edge.mat);
    edge.mat1$color = ComputeColorGradientCorr(E(g)$correlation);
    edge.mat1 = as.matrix(edge.mat1);
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3], size=edge.width, correlation = E(g)$correlation);
    #different layouts
    pos.xy <- PerformLayOut(g,layoutOpt);
    nodes <- vector(mode="list");
    for(i in 1:length(node.sizes)){
        nodes[[i]] <- list(
        id = nms[i],
        label=nms[i],
        type=shapes[i],
        size=node.sizes[i],
        color=colors[i],
        x = pos.xy[i,1],
        y = pos.xy[i,2]
    );
    }
  #nodes measures
    nd.tbl <- data.frame(Id=nms, Degree=node.dgr, Betweenness=round(node.btw,2));
    # order 
    ord.inx <- order(nd.tbl[,2], nd.tbl[,3], decreasing = TRUE)
    nd.tbl <- nd.tbl[ord.inx, ];
    netData <- list(nodes=nodes, edges=edge.mat);
    if (cor.meth=="spearman"){
        write.csv(nd.tbl, file="spearman_node_table.csv", row.names=FALSE);
        sink("spear_corr_network.json");
    } else if (cor.meth=="pearson") {
        write.csv(nd.tbl, file="pearson_node_table.csv", row.names=FALSE);
        sink("pear_corr_network.json");
    } else if (cor.meth=="mic") {
        write.csv(nd.tbl, file="mic_node_table.csv", row.names=FALSE);
        sink("mic_corr_network.json");
    } else {
        write.csv(nd.tbl, file="cclasso_node_table.csv", row.names=FALSE);
        sink("cclasso_corr_network.json");
    }
    suppressMessages(require(RJSONIO));
    cat(toJSON(netData));
    sink();
    return(1);
}
###############################################################################
###########################
##### rCCA (MixOmics R package)
###########################
InteromicsrCCA <- function(taxrank,generank,ncomp,analmethod,cvmethod,colOpt,colPal,variable,showlabel,imgNm1,format="png",dpi=72){
    set.seed(11313);
    #using only filtered data (microbiome)
    m.data <- otu_table(dataSet$m.filt.data,taxa_are_rows =TRUE);
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    m.dataphy <- merge_phyloseq(m.data,dataSet$m.taxa_table,sample_table);
    #merging at taxonomy levels
    m.dataphy <- fast_tax_glom_first(m.dataphy,taxrank);
    nm <- as.character(tax_table(m.dataphy)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Unassigned"; 
    data1 <- as.matrix(otu_table(m.dataphy));
    rownames(data1) <- nm;
    m.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    
    #using only filtered data (resistome)
    r.data <- otu_table(dataSet$filt.data,taxa_are_rows =TRUE);
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.dataphy <- merge_phyloseq(r.data,dataSet$taxa_table,sample_table);    
    #merging at functional levels
    r.dataphy <- fast_tax_glom_first(r.dataphy,generank);
    nm <- as.character(tax_table(r.dataphy)[,generank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Unassigned"; 
    data1 <- as.matrix(otu_table(r.dataphy));
    rownames(data1) <- nm;
    r.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    
    #add pseudo count (offset)
    suppressMessages(require("mixOmics"));
    r.data[r.data == 0] <- 0.05;
    m.data[m.data == 0] <- 0.05;
    #CLR normalization
    m.data.norm <- logratio.transfo(as.matrix(t(m.data)), logratio = 'CLR', offset = 0);    
    r.data.norm <- logratio.transfo(as.matrix(t(r.data)), logratio = 'CLR', offset = 0); 
    
    #functions require Sample by feature matrix
    if (analmethod == "shrinkage"){
        res <- rcc(r.data.norm, m.data.norm, ncomp = ncomp, method = 'shrinkage');
    } else {
        cv <- tune.rcc(r.data.norm,m.data.norm, grid1 = seq(0.001, 1, length = 5), grid2 = seq(0.001, 1, length = 5), 
            validation = cvmethod,folds = 5, plot=FALSE);
        res <- rcc(r.data.norm, m.data.norm, ncomp = ncomp, method = analmethod,lambda1 = cv$opt.lambda1,lambda2 = cv$opt.lambda2);       
    }
    Cairo(file=imgNm1, width=720, height=500, type=format, bg="white",dpi=dpi);
    rcc.res.data <- res;
    suppressMessages(require(RColorBrewer));      
    m_sam_nm <- sample_names(m.dataphy);
    rcc.res.data$variates$Y <- rcc.res.data$variates$Y[m_sam_nm,];
    rcc.res.data$variates$X <- rcc.res.data$variates$X[m_sam_nm,];
    rcc.res.data$variates$X <- as.data.frame(rcc.res.data$variates$X);
    rcc.res.data$variates$Y <- as.data.frame(rcc.res.data$variates$Y);
    rcc.res.data$variates$Y$Data <- rep("microbiome",nrow(rcc.res.data$variates$Y));
    rcc.res.data$variates$X$Data <- rep("resistome",nrow(rcc.res.data$variates$X));
    rcc.res.data$variates$X$Sample <- rownames(rcc.res.data$variates$X);
    rcc.res.data$variates$Y$Sample <- rownames(rcc.res.data$variates$Y);
    if (colOpt!="omics"){        
        rcc.res.data$variates$X[,variable] <- sample_data(r.dataphy)[m_sam_nm][[variable]]; 
        rcc.res.data$variates$Y[,variable] <- sample_data(m.dataphy)[m_sam_nm][[variable]];
    }
    rcc.res.data <- rbind(rcc.res.data$variates$X,rcc.res.data$variates$Y);
    rcc.score <- signif(rcc.res.data[,1:3],5);
    rcc.score$Data <- rcc.res.data$Data;
    write.csv(rcc.score,row.names=rcc.res.data[,"Sample"], file="rCCA_ordination_score.csv");
    rca.score <- NULL;
    suppressMessages(require(dplyr));
    if (colOpt=="omics"){
      box <- rcc.res.data %>% ggplot(aes(rcc.res.data[,1],rcc.res.data[,2], color=Data,shape=Data))+
            geom_point(aes(fill=Data),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2");
    } else {
        grp <- rcc.res.data[,variable];
        box <- rcc.res.data %>% ggplot(aes(rcc.res.data[,1],rcc.res.data[,2],color=grp,shape=Data))+
            geom_point(aes(fill=grp),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2",fill=variable,color=variable);
    }
    box <- box+geom_line(aes(group = Sample),color="grey")+theme_bw()+scale_shape_manual(values=c(22,21));
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    if(showlabel=="samnm"){
       box <- box+geom_text(aes(label=rcc.res.data[,"Sample"]),hjust=0.5, vjust=2,size=3,fontface="bold");
    } else if(showlabel!="none") {
       smplnm <- rcc.res.data[,"Sample"]
       clslbl <- sample_data(m.dataphy)[smplnm][[showlabel]];
       box <- box+geom_text(aes(label=clslbl),hjust=0.5, vjust=2,size=3,fontface="bold");
    }
    box <- box+ theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    dataSet$cca.res <- res;
    dataSet$merge.data <- rcc.res.data;
    dataSet$omics.meth <-"rcca";
    m_sam_nm <<- m_sam_nm;
    #m.coords <<- res$variates$Y[,1:3];
    #r.coords <<- res$variates$X[,1:3];
    dataSet$cca.X <- r.data.norm;
    dataSet$cca.Y <- m.data.norm;
    variable <<-variable;
    dataSet <<- dataSet;
    return(1);
}
###############################################
####sPLS ###################
###############################################
PerformintgromsPLS <- function(taxrank,generank,ncomp,analmethod,colOpt,colPal,variable,showlabel,imgNm1,format="png",dpi=72){
    set.seed(1131345);
    #using only filtered data (microbiome)
    m.data <- otu_table(dataSet$m.filt.data,taxa_are_rows =TRUE);
    dataSet$m.taxa_table <- tax_table(dataSet$m.taxa_table);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    m.dataphy <- merge_phyloseq(m.data,dataSet$m.taxa_table,sample_table);
    #merging at taxonomy levels
    m.dataphy <- fast_tax_glom_first(m.dataphy,taxrank);
    nm <- as.character(tax_table(m.dataphy)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Unassigned"; 
    data1 <- as.matrix(otu_table(m.dataphy));
    rownames(data1) <- nm;
    m.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    
    #using only filtered data (resistome)
    r.data <- otu_table(dataSet$filt.data,taxa_are_rows =TRUE);
    dataSet$taxa_table <- tax_table(dataSet$taxa_table);
    r.dataphy <- merge_phyloseq(r.data,dataSet$taxa_table,sample_table);    
    #merging at functional levels
    r.dataphy <- fast_tax_glom_first(r.dataphy,generank);
    nm <- as.character(tax_table(r.dataphy)[,generank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Unassigned"; 
    data1 <- as.matrix(otu_table(r.dataphy));
    rownames(data1) <- nm;
    r.data <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    
    #add pseudo count (offset)
    suppressMessages(require("mixOmics"));
    r.data[r.data == 0] <- 0.05;
    m.data[m.data == 0] <- 0.05;
    #CLR normalization
    m.data.norm <- logratio.transfo(as.matrix(t(m.data)), logratio = 'CLR', offset = 0);    
    r.data.norm <- logratio.transfo(as.matrix(t(r.data)), logratio = 'CLR', offset = 0); 
    res <- spls(r.data.norm,m.data.norm, ncomp = ncomp, mode = analmethod);
    spls.res.data <- res;
    Cairo(file=imgNm1, width=720, height=500, type=format, bg="white",dpi=dpi);    
    suppressMessages(require(RColorBrewer));      
    m_sam_nm <- sample_names(m.dataphy);
    spls.res.data$variates$Y <- spls.res.data$variates$Y[m_sam_nm,];
    spls.res.data$variates$X <- spls.res.data$variates$X[m_sam_nm,];
    spls.res.data$variates$X <- as.data.frame(spls.res.data$variates$X);
    spls.res.data$variates$Y <- as.data.frame(spls.res.data$variates$Y);
    spls.res.data$variates$Y$Data <- rep("microbiome",nrow(spls.res.data$variates$Y));
    spls.res.data$variates$X$Data <- rep("resistome",nrow(spls.res.data$variates$X));
    spls.res.data$variates$X$Sample <- rownames(spls.res.data$variates$X);
    spls.res.data$variates$Y$Sample <- rownames(spls.res.data$variates$Y);
    if (colOpt!="omics"){        
        spls.res.data$variates$X[,variable] <- sample_data(r.dataphy)[m_sam_nm][[variable]]; 
        spls.res.data$variates$Y[,variable] <- sample_data(m.dataphy)[m_sam_nm][[variable]];
    }
    spls.res.data <- rbind(spls.res.data$variates$X,spls.res.data$variates$Y);
    spls.score <- signif(spls.res.data[,1:3],5);
    spls.score$Data <- spls.res.data$Data;
    write.csv(spls.score,row.names=spls.res.data[,"Sample"], file="sPLS_ordination_score.csv");
    spls.score <- NULL;
    suppressMessages(require(dplyr));
    if (colOpt=="omics"){
      box <- spls.res.data %>% ggplot(aes(spls.res.data[,1],spls.res.data[,2], color=Data,shape=Data))+
            geom_point(aes(fill=Data),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2");
    } else {
        grp <- spls.res.data[,variable];
        box <- spls.res.data %>% ggplot(aes(spls.res.data[,1],spls.res.data[,2],color=grp,shape=Data))+
            geom_point(aes(fill=grp),size=4,alpha =0.8)+labs(x="Component 1",y="Component 2",fill=variable,color=variable);
    }
    box <- box+geom_line(aes(group = Sample),color="grey")+theme_bw()+scale_shape_manual(values=c(22,21));
    if(colPal!="default"){
        box <- box+scale_color_brewer(palette=colPal)+scale_fill_brewer(palette=colPal);
    }
    if(showlabel=="samnm"){
       box<-box+geom_text(aes(label=spls.res.data[,"Sample"]),hjust=0.5, vjust=2,size=3,fontface="bold");
    } else if(showlabel!="none") {
       smplnm <- spls.res.data[,"Sample"]
       clslbl <- sample_data(m.dataphy)[smplnm][[showlabel]];
       box<-box+geom_text(aes(label=clslbl),hjust=0.5, vjust=2,size=3,fontface="bold");
    }
    box <- box+ theme(axis.text.x = black.bold.10.text,axis.text.y = black.bold.10.text,axis.title.y=black.bold.11.title,axis.title.x=black.bold.11.title,axis.title=black.bold.11.title,legend.text=black.bold.10.text,legend.title = black.bold.11.title);
    print(box);
    dev.off();
    dataSet$spls.res <- res;
    dataSet$merge.data <- spls.res.data;
    dataSet$omics.meth <- "spls";
    m_sam_nm <<- m_sam_nm;
    #m.coords <<- res$variates$Y[,1:3];
    #r.coords <<- res$variates$X[,1:3];
    variable <<-variable;
    dataSet$cca.X <- r.data.norm;
    dataSet$cca.Y <- m.data.norm;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################
PlotCIM <- function(imgNm,format="png",distNm,distalgo,colors_cntrst,viewOpt,width=NA,dpi=72,border=T){
    if(dataSet$omics.meth=="spls"){
        res <- dataSet$spls.res;
    } else{
        res <- dataSet$cca.res;
    }
    #use the cim function from MixOmics package to get plot info and redo in pheatmap
    #the issue is unnecessary plotting of image
    Cairo(file = "cim_default.png", unit="in", dpi=22, width=200, height=200, type="png", bg="white");
    cim.res <- cim(res,dist.method = c(distNm,distNm),clust.method=c(distalgo,distalgo));
    dev.off();
    #get dendrogram ; already computed by cim
    row.dend <- cim.res$ddr;
    col.dend <- cim.res$ddc;
    #covert into hclust object
    row.clust <- as.hclust(row.dend);
    col.clust <- as.hclust(col.dend);
    
    #reorder the correlation matrix based on orignal data
    ordimp.x <- colnames(dataSet$cca.X);
    ordimp.y <- colnames(dataSet$cca.Y);
    cor.mat.orig <- cim.res$mat.cor;
    cim.res <- NULL;
    cor.mat.ord <- cor.mat.orig[ordimp.x,ordimp.y]; 
    
    # transpose the data (rows are microbes and column is resistome)
    cor.mat.ord <- t(cor.mat.ord);
    write.csv(cor.mat.ord,"check.csv");
    
    # set up parameter for heatmap
    suppressMessages(require(RColorBrewer));
    suppressMessages(require(gplots));
    if(colors_cntrst == "gbr"){
        colors <- colorRampPalette(c("green", "black", "red"), space="rgb")(256);
    }else if(colors_cntrst == "heat"){
        colors <- heat.colors(256);
    }else if(colors_cntrst == "topo"){
        colors <- topo.colors(256);
    }else if(colors_cntrst == "gray"){
        colors <- colorRampPalette(c("grey90", "grey10"))(256);
    }else if(colors_cntrst == "spectral"){
        suppressMessages(require(RColorBrewer));
        colors <- rev(colorRampPalette(brewer.pal(10, "Spectral"))(256));
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
    first_color <- colors[1];
    last_color <- colors[length(colors)];
    if(viewOpt == "overview"){
        if(nrow(cor.mat.ord) > 3000){
            h <- 50;
        }else if(nrow(cor.mat.ord) > 2000 & nrow(cor.mat.ord) < 3000){
            h <- 30;
        }else if(nrow(cor.mat.ord) > 1000 & nrow(cor.mat.ord) < 2000){
            h <- 20;
        }else if(nrow(cor.mat.ord) > 500 & nrow(cor.mat.ord) < 1000){
            h <- 15;
        }else if(nrow(cor.mat.ord) > 50){
            myH <- nrow(cor.mat.ord)*12 + 40;
            h <-0;
        }else if(nrow(cor.mat.ord) > 20){
            myH <- nrow(cor.mat.ord)*10 + 60;
            h <- 0;
        }else{
            h <- 9;
        }
        if (h==0){
            h <- round(myH/72,2);
        }

        if(ncol(cor.mat.ord) > 3000){
            w <- 50;
        }else if(ncol(cor.mat.ord) > 2000 & ncol(cor.mat.ord) < 3000){
            w <- 30;
        }else if(ncol(cor.mat.ord) > 1000 & ncol(cor.mat.ord) < 2000){
            w <- 20;
        }else if(ncol(cor.mat.ord) > 500 & ncol(cor.mat.ord) < 1000){
            w <- 15;
        }else if (ncol(cor.mat.ord) > 50){
            myW <- ncol(cor.mat.ord)*12 + 40;
            w <- 0;
        }else if(ncol(cor.mat.ord) > 20){
            myW <- ncol(cor.mat.ord)*10 + 60;
            w <- 0;
        }else{
            w <- 9;
        }
        if (w==0){
            w <- round(myW/72,2);
        }
    }else{
        if(nrow(cor.mat.ord) > 3000){
            myH <- nrow(cor.mat.ord)*5;
        }else if(nrow(cor.mat.ord) > 2000 & nrow(cor.mat.ord) < 3000){
            myH <- nrow(cor.mat.ord)*4;
        }else if(nrow(cor.mat.ord) > 1000 & nrow(cor.mat.ord) < 2000){
            myH <- nrow(cor.mat.ord)*3;
        }else if(nrow(cor.mat.ord) > 500 & nrow(cor.mat.ord) < 1000){
            myH <- nrow(cor.mat.ord)*2;
        }else if(nrow(cor.mat.ord) > 50){
            myH <- nrow(cor.mat.ord)*15 + 40;
        }else if(nrow(cor.mat.ord) > 20){
            myH <- nrow(cor.mat.ord)*12 + 60;
        }else{
            myH <- 792;
        }
        h <- round(myH/72,2);
        if(ncol(cor.mat.ord) > 3000){
            myH <- ncol(cor.mat.ord)*5;
        }else if(ncol(cor.mat.ord) > 2000 & ncol(cor.mat.ord) < 3000){
            myH <- ncol(cor.mat.ord)*4;
        }else if(ncol(cor.mat.ord) > 1000 & ncol(cor.mat.ord) < 2000){
            myH <- ncol(cor.mat.ord)*3;
        }else if(ncol(cor.mat.ord) > 500 & ncol(cor.mat.ord) < 1000){
            myH <- ncol(cor.mat.ord)*2;
        }else if(ncol(cor.mat.ord) > 50){
            myW <- ncol(cor.mat.ord)*15 + 40;
        }else if(ncol(cor.mat.ord) > 20){
            myW <- ncol(cor.mat.ord)*12 + 60;
        }else{
            myW <- 792;
        }
        w <- round(myW/72,2);
        # to prevent too small
        min.w <- 4.8;
        if(w < min.w){
            w <- h <- min.w;
        }
    }
    if(border){
        border.col <- "grey60";
    }else{
        border.col <- NA;
    }
    Cairo(file = imgNm, unit="in", dpi=72, width=w, height=h, type=format, bg="white");
    require(pheatmap);
    pheatmap(cor.mat.ord, 
            fontsize=8, fontsize_row=8, 
            cluster_rows = col.clust, 
            cluster_cols = row.clust,
            color = colors,
            border_color = border.col,
            fontface="bold"
        );
    dev.off();
}
###############################################################################
PlotScreePlot <- function(imgNm,format="png",dpi=72){
    if(dataSet$omics.meth=="spls"){
        res <- dataSet$spls.res;
    } else{
        res <- dataSet$cca.res;
    }
    Cairo(file=imgNm, width=680, height=600, type=format, bg="white",dpi=dpi);
    plot(res, scree.type = "barplot");
    dev.off();
}
###############################################################################
PlotCorrCircleplot <- function(imgNm,compx,compy,colPal,cor_thres,format="png",dpi=72){
    if(dataSet$omics.meth=="spls"){
        res <- dataSet$spls.res;
    } else{
        res <- dataSet$cca.res;
    }
    #colors
    if(colPal=="default"){
        col1 <- "#F8766D"
        col2 <- "#00BFC4"
    } else {
        dis.cols <- brewer.pal(2, colPal);
        col1 <- dis.cols[1];
        col2 <- dis.cols[2];
    }
    Cairo(file=imgNm, width=720, height=500, type=format, bg="white",dpi=dpi);
    ccr.img <- plotVar2(res,comp = c(compx,compy),cex=c(3.2,3.2),col=c(col2,col1),
               cutoff= cor_thres, title="",legend=TRUE);
    dev.off();           
    if (!exists("ccr.img")){
        current.msg<<-"Cutoff value very high for the respective components.No variable was selected.";
        return (0);
    }else {
        return (1);
    }
}
###############################################################################