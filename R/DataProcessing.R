##################################################
## R script for ResistoXplorer
## Description: Data/resource management functions
## Author: Achal Dhariwal, UiO
###################################################

# init resources for analysis
Init.Data <-function(){
    dataSet <<- analSet <<- list(); 
    imgSet <<- list();
    current.msg <<- current.org <<- "";
    msg.vec <<- vector(mode="character");
    module.count <<- 0;
    lib.path <<- "../../data/lib/";
    # preload some general package
    library('Cairo');
    CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
    print("init ResistoXplorer!");
}
###############################################################################
SetAnalType <- function(analType){
    anal.type <<- analType;
}
###############################################################################
GetClassInfo <- function(){
    return(levels(dataSet$cls));
}
###############################################################################

### Read Metadata file###########
#################################
ReadSampleTable <- function(dataName) { 
    msg <- NULL;
    mydata <- .readDataTable(dataName);
    if(is.null(mydata) || class(mydata) == "try-error"){
        current.msg <<- "Failed to read in the metadata! Please make sure that your metadata file is in the right format and do not have empty cells or NA.";
        return(0);
    }
    # look for #NAME, store in a list
    sam.nm <- substr(colnames(mydata[1]),1,5); 
    sam.nm <- tolower(sam.nm);
    sam.inx <- grep("#name",sam.nm);
    if(length(sam.inx) > 0){ 
        smpl_nm <- mydata[,1];
        smpl_var <- colnames(mydata[-1]);
    }else{
        current.msg <<- "Please make sure you have label #NAME in your sample data file!";
        return(0);
    }
    # converting to character matrix as duplicate row names not allowed in data frame.
    mydata <- as.matrix(mydata[,-1]);
    rownames(mydata) <- smpl_nm;
    colnames(mydata) <- smpl_var;

    # empty cell or NA cannot be tolerated in metadata
    na.inx <- is.na(mydata);
    na.msg <- NULL;
    if(sum(na.inx) > 0){
        mydata[na.inx] <- "Unknown";
        na.msg <- paste("A total of", sum(na.inx), "empty or NA values were replaced by 'Unknown'.");
    }
    dataSet$sample_data <- data.frame(mydata);
    dataSet <<- dataSet;
    current.msg <<- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");
    dataSet$smpl.msg <- current.msg;
    dataSet <<- dataSet;
    return (1);
}
###############################################################################

### Read Abundance (count data) file###########
#################################
ReadAbundanceData <- function(dataName,db,analtype,datatype) {
    mydata <- .readDataTable(dataName);
    if(is.null(mydata) || class(mydata) == "try-error"){
        current.msg <<- "Failed to read in the abundance data! Please make sure the abundance table is in the right format and do not have empty cells or NA.";
        return("F");
    }
    # look for #NAME, store in a list
    #just to check the name of labels.
    sam.nm <-substr(colnames(mydata[1]),1,5); 
    sam.nm <-tolower(sam.nm);
    sam.inx <- grep("^#name",sam.nm);
    if(length(sam.inx) > 0){ 
        smpl_nm <- colnames(mydata[-1]);
    }else{
        current.msg <<- "Please make sure you have label #NAME in your abundance data file!";
        return(0);
    }
    dat.nms <- mydata[,1];
    mydata <- as.matrix(mydata[,-1]);
    if(mode(mydata)=="character"){
        current.msg <<- paste("Errors in parsing your data as numerics - possible reason: comma as decimal separator?");
        return(0);
    }
    rownames(mydata) <- dat.nms;
    dataSet$name <- basename(dataName);
    if(analtype == "omics"){
        if(datatype == "microme"){
            dataSet$m.data.orig <- mydata;            
        }else{
            dataSet$data.orig <- mydata;
        }
    }else{
        dataSet$data.orig <- mydata;
    }
    dataSet$smpl_nm <- smpl_nm;
    current.msg <<- paste("A total of ", ncol(mydata) , " samples and ", nrow(mydata), " features are present in abundance data file.");
    dataSet$db.type <-db;
    dataSet <<- dataSet;
    return (1);
}
###############################################################################

### Read Annotation file###########
#################################
ReadPofileTable <- function(dataName,analtype,datatype) {
    msg <- NULL;
    mydata <- .readDataTable(dataName);
    if(is.null(mydata) || class(mydata) == "try-error"){
        current.msg <<- "Failed to read in the annotation data! Please make sure the data is in the right format!";
        return(0);
    }
    # look for #ANNOTATION, store in a list
    sam.nm <- substr(colnames(mydata[1]),1,11); 
    sam.nm <- tolower(sam.nm);
    sam.inx <- grep("^#annotation",sam.nm);
    if(length(sam.inx)>0){ 
        tax_nm <- mydata[,1];
        tax_rank <- colnames(mydata[-1]);
    }else{
        current.msg <<- "No labels #ANNOTATION found in your annotation file!";
        return(0);
    }
    # converting to character matrix as duplicate row names not allowed in data frame.
    mydata <- as.matrix(mydata[,-1]);
    rownames(mydata) <- tax_nm;
    colnames(mydata) <- tax_rank;
    current.msg <<- paste("Annotation file has total of ", nrow(mydata), "features and", ncol(mydata), "annotation levels. Features which are not present in abundance table will be removed from further anaylsis.");
    if(analtype == "omics"){
        if(datatype == "microme"){
            dataSet$m.taxa_table <- mydata;
        }else{
            dataSet$taxa_table <- mydata;
        }
    }else{
        dataSet$taxa_table <- mydata;
    }
    dataSet <<- dataSet;
    return (1);
}
###############################################################################

### Get Annotation info from precompiled DB###########
#################################
MapDbAnnotation <- function (db,data.proc) {          
        if(db=="ResFinder"||db=="ARDB"||db=="AMP"||db=="CARD"||db=="BacMet"){
            gene_list.db <- readRDS(paste(lib.path, "db/", db, "/", tolower(db), "_annot.rds", sep=""));           
        }else if(db=="MEGARes1"){ #V:1.0
            gene_list.db <- readRDS(paste(lib.path, "db/MEGARes/megares1_","annot.rds", sep=""));
        }else if(db=="MEGARes2_drug"){V:2.0
            gene_list.db <- readRDS(paste(lib.path, "db/MEGARes/", tolower(db), "_annot.rds", sep=""));
        }else if(db=="MEGARes2_full"){ #V:2.0
            gene_list.db <- readRDS(paste(lib.path, "db/MEGARes/", tolower(db), "_annot.rds", sep=""));
        }else{
            gene_list.db <- readRDS(paste(lib.path, "db/Others/", tolower(db), "_annot.rds", sep=""));
            #other databases: AMRfinder,ARG-ANNOT, deepARGDb etc.
        }
        rowInd <- which(rownames(data.proc) %in% rownames(gene_list.db)==TRUE);    
        if(length(rowInd)==0){ 
            return(NULL);    
        } else if(length(rowInd)< round((length(rownames(data.proc))*0.1),0)) {
             return(0); 
        } else{
            data.proc <- data.proc[rowInd,];
            #create a dummy taxa_table
            taxa_tableInd <- which(rownames(gene_list.db) %in% rownames(data.proc)==TRUE);
            taxa_table <- gene_list.db[taxa_tableInd,]; 
            taxa_table <- taxa_table[row.names(data.proc), ]; 
            data.map <- list(data.proc=data.proc,taxa_table=taxa_table);
            return(data.map);
        }
}
###############################################################################

###########Data Summary#############
#######################################
SanityChecking <- function(db,profinfo.type){
   #check for compositionality and sparsity based on initial data   
   iscomp <- all(abs(colSums(dataSet$data.orig) - 1) < 1e-06);
   sparsity  <- signif((length(which(dataSet$data.orig==0))/length(dataSet$data.orig))*100,4);
   # remove features with all the zeros
   feat.sums <- apply(dataSet$data.orig, 1, function(x){sum(x>0, na.rm=T)}); 
   gd.inx <- feat.sums > 0; # occur in at least 1 samples
   if(length(which(gd.inx=="TRUE"))==0){
       current.msg <<-"No feature read counts are present in more than one samples.";
       return(0);
   }
   data.proc <- dataSet$data.orig[gd.inx, ];
   gd_feat <- nrow(data.proc);
   
   #zero variance filter (not applied)
   # filtering the constant features here (check for columns with all constant (var=0))
   # not sure though (should we apply it or not)
   #varCol <- apply(dataSet$data.proc, 1, var, na.rm=T);
   #constCol <- varCol == 0 | is.na(varCol);
   # making copy of data.proc and proc.phyobj(phyloseq)
   #data.proc<-dataSet$data.proc[!constCol, ];
   #if(length(data.proc)==0){
       #current.msg <<-"All features are found to be constant and have been removed from data. No data left after such processing.";
       #return(0);
   #}
   #saveRDS(data.proc, file="data.proc.orig"); # save an copy
   #write.csv(dataSet$data.prefilt,"phylo.txt") 
   # now get stats
   gene_no.orig <- nrow(dataSet$data.orig);

   # no of samples in abundance table
   sample_no <- ncol(dataSet$data.orig);
   
   #metadata sanity checking
   #initial no of experimental factors
   init_exp_fac <-ncol(dataSet$sample_data);
   dataSet$sample_data <- dataSet$sample_data[sapply(dataSet$sample_data, function(col) length(unique(col))) > 1];
   if(ncol(dataSet$sample_data)==0){
        current.msg <<- "No experimental factor have more than one group. Please provide variables with at least two groups.";
        return(0);
   }
   #converting all sample variables to factor type
   character_vars <- lapply(dataSet$sample_data, class) == "character";
   dataSet$sample_data[, character_vars] <- lapply(dataSet$sample_data[, character_vars], as.factor);
   
   #check whether sample variables are continuous and removing them
   dataSet$sample_data <- dataSet$sample_data[sapply(dataSet$sample_data, function(x) length(x)/length(levels(x))) > 2]; 
   if(ncol(dataSet$sample_data)==0){
        current.msg <<-"No experimental factor have discrete or categorical values. Please provide discrete variables.";
        return(0);
    }
    disc_exp_fac <- ncol(dataSet$sample_data);
    cont_exp_fac <- init_exp_fac - disc_exp_fac;
    
    #abundance and metadata file cross checking
    init_meta_smpl_no <- length(rownames(dataSet$sample_data));
    if (length(colnames(data.proc))> init_meta_smpl_no){
        current.msg <<-"The number of samples present in your metadata file is less than abundance file! 
                        Please make sure that you provide at least same or more number of samples in metadata file.";
        return(0);
    }
    indx <- match(colnames(data.proc), rownames(dataSet$sample_data));
    if(all(is.na(indx))){
        current.msg <<- "No sample name matches between abundance and metadata file !";
        return(0);
    }else if (all(!is.na(indx))){
        dataSet$sample_data <- dataSet$sample_data[indx,];
        abund_meta_sample_no <- nrow(dataSet$sample_data);
    }else{
        # remove unmatched sample names and then update both metadata and count table
        indx1 <- indx[!is.na(indx)];
        dataSet$sample_data <- dataSet$sample_data[indx1, ,drop=FALSE];
        data.proc <- data.proc[,rownames(dataSet$sample_data)];
        abund_meta_sample_no <- nrow(dataSet$sample_data);
    }
   # now store data.orig to RDS 
   saveRDS(dataSet$data.orig, file="data.orig");
   
    #annotation file 
   if(profinfo.type == "db"){
       init_annot_feat_no <- 0; 
       data.map<-MapDbAnnotation(db,data.proc);
        if (is.null(data.map)){
            current.msg <<- "No hits found in the selected database! Please make sure you have selected the correct annotation database and your feature names are in correct required format.";
            return(0);
        }else if(!is.list(data.map) & length(data.map) == 1){
            current.msg <<- "Less than 10% of features found in the selected annotation database! Please either make sure that your all feature names are in correct required format or upload your own annotation file";
            return(0);
        }else {
           data.proc <- data.map$data.proc;
           taxa_table <-data.map$taxa_table;
           dataSet$taxa_table <-taxa_table;
        }
   } else{
        init_annot_feat_no <- nrow(dataSet$taxa_table);
        rowInd <- which(rownames(data.proc) %in% rownames(dataSet$taxa_table)==TRUE);  
        if(length(rowInd)==0){
            current.msg <<- "No feature names matches between abundance and annotation file !";
            return(0);
        }else if (length(rowInd)< round(nrow(data.proc)*0.1,0)){
            current.msg <<- "Less than 10% of feature names matched between annotation and abundance file!";
            return(0);
        }else{
            data.proc <- data.proc[rowInd,];
            #create a dummy taxa_table
            taxa_tableInd <- which(rownames(dataSet$taxa_table) %in% rownames(data.proc)==TRUE);
            taxa_table <- dataSet$taxa_table[taxa_tableInd,]; 
            taxa_table <- taxa_table[row.names(data.proc), ]; 
            dataSet$taxa_table <- taxa_table;
        }
   }
   #stats
   rem.gene_no <- nrow(data.proc);
   smpl.sums <- apply(data.proc, 2, sum);
   tot_size <- sum(smpl.sums);
   smin <- min(smpl.sums);
   smean <- mean(smpl.sums);
   smax <- max(smpl.sums);
   dataSet$data.prefilt <- dataSet$data.proc<-data.proc;
   dataSet$anot.type <- profinfo.type;
   dataSet$anot.db <- db;
   annot.lvl <- ncol(dataSet$taxa_table);
   #write.csv(dataSet$data.prefilt,"filt_data.txt")
   #write.csv(dataSet$taxa_table,"taxa_data.txt")
   #dataSet$data.proc2<-data.proc2;
   dataSet <<- dataSet;
   return(c(1,gene_no.orig,gd_feat,rem.gene_no,sample_no,init_exp_fac,tot_size,smean,smax,smin,sparsity,iscomp,disc_exp_fac,cont_exp_fac,init_meta_smpl_no,abund_meta_sample_no,init_annot_feat_no, annot.lvl));   
}
###############################################################################

### Data Summary (Integration)###########
#################################
SanityOmicsChecking <- function(){
   #first do microbiome data 
   m.feat.sums <- apply(dataSet$m.data.orig, 1, function(x){sum(x>0, na.rm=T)});
   m.comp <- length(which(colSums(dataSet$m.data.orig) > 1));
   if (m.comp == 0){
       miscomp <- TRUE;
    }else{
       miscomp <- FALSE;
    }
   m.sparsity  <- signif((length(which(dataSet$m.data.orig == 0))/length(dataSet$m.data.orig))*100,4);
   m.gd.inx <- m.feat.sums > 1; # occur in at least 1 samples
   if(length(which(m.gd.inx=="TRUE"))==0){
       current.msg <<-"Reads occur in only one sample in taxonomic profile.  All these are considered as artifacts and have been removed from data. No data left after such processing.";
       return(0);
   }
   m.data.proc <- dataSet$m.data.orig[m.gd.inx, ];
   m.gd_feat <- nrow(m.data.proc);
   
   # now get stats
   m.gene_no.orig <- nrow(dataSet$m.data.orig);

   # from abundance table
   m.sample_no <- ncol(dataSet$m.data.orig);
   
   # now store data.orig to RDS 
   saveRDS(dataSet$m.data.orig, file="m.data.orig");
   
   m.smpl.sums <- apply(m.data.proc, 2, sum);
   m.tot_size <- sum(m.smpl.sums);
   m.smin <- min(m.smpl.sums)
   m.smean <- mean(m.smpl.sums)
   m.smax <- max(m.smpl.sums); 
   #####resistome data#########
   r.feat.sums <- apply(dataSet$data.orig, 1, function(x){sum(x>0, na.rm=T)});
   r.comp <- length(which(colSums(dataSet$data.orig) > 1));
   if(r.comp == 0){
       riscomp <- TRUE;
    }else{
       riscomp <- FALSE;
    }
   r.sparsity <- signif((length(which(dataSet$data.orig==0))/length(dataSet$data.orig))*100,4);
   r.gd.inx <- r.feat.sums > 1; # occur in at least 1 samples
   if(length(which(r.gd.inx=="TRUE"))==0){
       current.msg <<- "Reads occur in only one sample in resistome profile.  All these are considered as artifacts and have been removed from data. No data left after such processing.";
       return(0);
   }
   data.proc <- dataSet$data.orig[r.gd.inx, ];
   r.gd_feat <- nrow(data.proc);
   
   # now get stats
   r.gene_no.orig <- nrow(dataSet$data.orig);
   r.sample_no <- ncol(dataSet$data.orig);
   
   # now store data.orig to RDS 
   saveRDS(dataSet$data.orig, file="data.orig");
   r.smpl.sums <- apply(data.proc, 2, sum);
   r.tot_size <-sum(r.smpl.sums);
   r.smin <- min(r.smpl.sums)
   r.smean <- mean(r.smpl.sums)
   r.smax <- max(r.smpl.sums); 
   
   #additional info
   m.samnames <- colnames(m.data.proc);
   r.samnames <- colnames(data.proc);
   indx <- match(m.samnames,r.samnames);
    if(all(is.na(indx))){
        current.msg <<- "Please make sure that sample names match between taxonomic and resistomic profile.";
        return(0);
    }
    match_sam_no <- length(which(!is.na(indx)));
    
    # keep only samples that matches in both the data
    indx2 <- match(r.samnames,m.samnames);
    if(length(which(is.na(indx)))=="0" && length(which(is.na(indx)))=="0"){
        all_sam_match <- TRUE;
    } else{
        all_sam_match <- FALSE;
    }
    r.sam.keep <- which(!is.na(indx2));
    data.proc <- data.proc[,r.sam.keep];
    m.sam.keep <- which(!is.na(indx));
    m.data.proc <- m.data.proc[,m.sam.keep];
   #############
   exp.fact <- 0;
   dis.exp.fact <- 0;
    
   #*need to check whether both metadata file and abundance files have same no of samples 
   sample_data <- dataSet$sample_data;
   if (exists("sample_data")){
     
       #converting all sample variables to factor type
       character_vars <- lapply(dataSet$sample_data, class) == "character";
       dataSet$sample_data[, character_vars] <- lapply(dataSet$sample_data[, character_vars], as.factor);
       
       #first total experimental factor
       exp.fact <- ncol(dataSet$sample_data); 
      
       #check whether sample variables are continuous and removing them
       dataSet$sample_data <- dataSet$sample_data[sapply(dataSet$sample_data, function(x) length(x)/length(levels(x))) > 2];
       if(ncol(dataSet$sample_data)==0){
            current.msg <<- "No sample variable have discrete or categorical values. Please provide discrete variables.";
            return(0);
        }
       dis.exp.fact <- ncol(dataSet$sample_data);
   }
   sample_data <-NULL;
   dataSet$m.data.prefilt <- dataSet$m.data.proc <-m.data.proc;
   dataSet$data.prefilt <- dataSet$data.proc <-data.proc;
   dataSet <<- dataSet;  
   return(c(1,m.gene_no.orig,m.gd_feat,m.sparsity,m.tot_size,m.smean,m.smax,m.smin,miscomp,
          riscomp,r.gene_no.orig,r.gd_feat,r.sparsity,r.tot_size,r.smean,r.smax,r.smin,m.sample_no,
          r.sample_no,exp.fact,dis.exp.fact,all_sam_match,match_sam_no,match_sam_no));   
}
###############################################################################

####Library Size Overview######
#################################
PlotLibSizeView <- function(dataType, imgName, format="png", dpi=72){
    if(dataType == "microme"){
        data_bef <- data.matrix(dataSet$m.data.proc);
    }else{
        data_bef <- data.matrix(dataSet$data.proc);
    }
    smpl.sums <- colSums(data_bef);
    names(smpl.sums) <- colnames(data_bef);
    smpl.sums <- sort(smpl.sums);
    smpl.sums <- rev(smpl.sums);
    vip.nms <- names(smpl.sums);
    names(smpl.sums) <- NULL;
    vip.nms <- substr(vip.nms, 1, 20);
    myH <- ncol(data_bef)*25 + 50;
    Cairo(file=imgName, width=600, height=myH, type=format, bg="white",dpi=dpi);
    xlim.ext <- GetExtendRange(smpl.sums, 10);
    par(mar=c(4,10,4,2));
    dotchart(smpl.sums, col="forestgreen", xlim=xlim.ext, pch=19, xlab="Read Counts", main="Library Size Overview");
    mtext(side=2, at=1:length(vip.nms), vip.nms, las=2, line=1)
    text(x=smpl.sums, y=1:length(smpl.sums), labels= round(smpl.sums), col="blue", pos=4, xpd=T)
    dev.off();   
}
###############################################################################

### Data Filtration ###########
#################################
# filter data based on low counts in high percentage samples
# note, first is abundance, followed by variance
ApplyAbundanceFilter <- function(filt.opt, count,oneminCount, isonefeatfilt, smpl.perc, analtype, datatype){
    msg <- NULL;
    if(analtype == "omics"){
        if(datatype == "microme"){
            data <- dataSet$m.data.prefilt;  
        }else{
            data <- dataSet$data.prefilt;
        }
    }else{
        data <- dataSet$data.prefilt;
        if(isonefeatfilt=="T"){
            gd.inx <- apply(data, 1, function(x){sum(x>0, na.rm=T)})==1; #features having count only in one samples
            if(length(which(gd.inx=="TRUE"))!=0){#if there are features having count in only one samples
                one.count <- 2;
                rsum.inx <- rowSums(data)<= one.count; #
                rm.inx <- gd.inx&rsum.inx==TRUE;
                one.featno <- length(which(rm.inx)==TRUE);
                data <- data[!rm.inx,];
                msg <- paste("A total of ", one.featno, " features having counts less than or equals to ", one.count , " in just one sample were removed", ".", sep=""); 
            } else {
            one.featno <-0;
            msg <- paste("None of the features are present in just one sample", ".", sep=""); 
            }
        }
    }
    #this data is used for sample categorial comparision further
    rmn_feat <- nrow(data);
    if(count==0){# no low-count filtering
        kept.inx <- rep(TRUE, rmn_feat);
    }else{
        if(filt.opt=="prevalence"){
            rmn_feat <- nrow(data);
            minLen <- smpl.perc*ncol(data);
            kept.inx <- apply(data, MARGIN = 1,function(x) {sum(x > count) >=minLen});  
        }else if(filt.opt == "mean"){
            filter.val <- apply(data, 1, mean, na.rm=T);
            kept.inx <- filter.val > count;
        }else if(filt.opt == "median"){
            filter.val <- apply(data, 1, median, na.rm=T);
            kept.inx <- filter.val > count;
        }
    }
    data <- data[kept.inx, ];
    if(datatype == "microme"){
        dataSet$m.filt.data <-data;
    } else{
        dataSet$filt.data <-data;
    }
    dataSet <<- dataSet;
    current.msg <<- paste(c(msg,"A total of ", sum(!kept.inx), " low abundance features were removed based on ", filt.opt, ".", sep="")); 
    return(1);
}
# filter data based on low abundace or variance
# note this is applied after abundance filter
ApplyVarianceFilter <- function(filtopt, filtPerct,analtype,datatype){
    if(analtype=="omics"){
        if(datatype=="microme"){
            data <- dataSet$m.filt.data; 
        }else{
            data <- dataSet$filt.data;
        }
    }else{
        data <- dataSet$filt.data;
    }
    msg <- NULL;
    
    rmn_feat <- nrow(data);

    filter.val <- nm <- NULL;

    if(filtPerct==0){# no low-count filtering
        remain <- rep(TRUE, rmn_feat);
    }else{
        if(filtopt == "iqr"){
            filter.val <- apply(data, 1, IQR, na.rm=T);
            nm <- "IQR";
        }else if(filtopt == "sd"){
            filter.val <- apply(data, 1, sd, na.rm=T);
            nm <- "standard deviation";
        }else if (filtopt == "cov"){
            sds <- apply(data, 1, sd, na.rm=T);
            mns <- apply(data, 1, mean, na.rm=T);
            filter.val <- abs(sds/mns);
            nm <- "Coeffecient of variation";
        }
        # get the rank
        rk <- rank(-filter.val, ties.method='random');
        var.num <- nrow(data);
        remain <- rk < var.num*(1-filtPerct);
    }
    data <- data[remain,];
    if(datatype=="microme"){
        dataSet$m.filt.data <-data;
    }else{
        dataSet$filt.data<-data;
    }
    rm.msg1 <- paste("A total of ", sum(!remain), " low variance features were removed based on ", filtopt, ".", sep="");
    rm.msg2 <- paste("The number of features remains after the data filtering step:", nrow(data));   
    current.msg <<- paste(c(current.msg,rm.msg1,rm.msg2), collapse=" ");
    dataSet$filt.msg <-current.msg;
    dataSet <<- dataSet;
    return(1);
}
###############################################################################

### Data Normalization###########
#################################
# note, here also update data type array/count
PerformNormalization <- function(rare.opt,scale.opt,transform.opt,analtype,datatype){
    if(analtype=="omics"){
        if(datatype=="microme"){
            data <- data.matrix(dataSet$m.filt.data); 
        }else{
            data <- data.matrix(dataSet$filt.data);
        }
    }else{
        data <- data.matrix(dataSet$filt.data);
    }
    tax_nm <- rownames(data); 
    msg <- NULL;
    library(phyloseq);
    if(rare.opt!="none"){
        otu.tab <- PerformRarefaction(data, rare.opt);
        msg <- c(msg, paste("Performed data rarefaction."));
    }else{
        msg <- c(msg, paste("No data rarefaction was performed."));
    }
    if(scale.opt!="none"){
        if(scale.opt=="cpm"){    
            data <- cpm.default(data, log=FALSE);    
            msg <- c(msg, paste("Performed count per million (CPM) normalization."));
        }else if(scale.opt=="rel"){ #repetitive so will be removed
            data <- apply(data, 2, function(x) {x/max(sum(x), 1e-32)});
            msg <- c(msg, paste("Performed proportion (relative) normalization"));
        }else if(scale.opt=="logcpm"){
            data <- cpm.default(data,log=TRUE); 
            msg <- c(msg, paste("Performed log count per million (log CPM) normalization"));
        }else if(scale.opt=="upperquartile"){
            suppressMessages(require("edgeR"));
            otuUQ <- edgeRnorm(data,method="upperquartile");
            data <- as.matrix(otuUQ$counts);
            msg <- c(msg, paste("Performed upperquartile normalization"));
        }else if(scale.opt=="CSS"){
            suppressMessages(require("metagenomeSeq"));
            
          #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
            data1 <- as(data, "matrix");
            dataMR <- newMRexperiment(data1);
            data <- cumNorm(dataMR, p=cumNormStat(dataMR));
            data <- MRcounts(data, norm = T);
            msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
        }else{
            print(paste("Unknown scaling parameter:", scale.opt));
        }
        otu.tab <- otu_table(data, taxa_are_rows =TRUE);
        taxa_names(otu.tab) <- tax_nm;
    }else{
        msg <- c(msg, paste("No data scaling was performed."));
    }
    if(transform.opt!="none"){
        if(transform.opt=="rle"){
            suppressMessages(require("edgeR"));
            otuRLE <- edgeRnorm(data,method="RLE");
            data <- as.matrix(otuRLE$counts);
            msg <- c(msg, paste("Performed RLE Normalization"));
        }else if(transform.opt=="TMM"){
            suppressMessages(require("edgeR"));
            otuTMM <- edgeRnorm(data, method="TMM");
            data <- as.matrix(otuTMM$counts);
            msg <- c(msg, paste("Performed TMM Normalization"));
        }else if(transform.opt=="clr"){ 
            pseudocount <- min(which(data!=0))*0.01;
            data[data==0] <- pseudocount;
            data <- cenLR_rcom(t(data));
            msg <- c(msg, paste("Performed centered log-ratio normalization."));
        }else if(transform.opt=="alr"){
            pseudocount <- min(which(data!=0))*0.01;
            data[data==0]<- pseudocount;
            data <- addLR_rcom(t(data));
            tax_nm <- rownames(data);
            msg <- c(msg, paste("Performed additive log-ratio normalization."));
        }else if(transform.opt=="ilr"){
            suppressMessages(require("compositions"));
            pseudocount <- min(which(data != 0))*0.01;
            data[data==0]<- pseudocount;
            colnm <- colnames(data);
            data <- as.data.frame(compositions::ilr(t(data)));
            tax_nm <- rownames(data);
            colnames(data) <- colnm;
            
            # the features name get distorted so no idea
            msg <- c(msg, paste("Performed isometric log-ratio normalization."));
        }else if(transform.opt=="hell"){ 
            suppressMessages(require("vegan"));
            data <- decostand(data, method="hellinger", MARGIN=2);
            msg <- c(msg, paste("Performed Hellinger transformation."));
        }else{
            print(paste("Unknown scaling parameter:", transform.opt));
        }
        otu.tab <- otu_table(data,taxa_are_rows =TRUE);
        taxa_names(otu.tab) <- tax_nm;
    }else{
        msg <- c(msg, paste("No data transformation was performed."));
    }
    # create phyloseq obj
    #after rarefaction the OTU sequences changes automatically
    dataSet$sample_data$sample_id <- rownames(dataSet$sample_data);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    phy.obj <- merge_phyloseq(otu.tab,sample_table);
    if(datatype=="microme"){
        dataSet$m.norm.phyobj<-phy.obj; 
    } else{
        dataSet$norm.phyobj<-phy.obj; 
    }
    current.msg <<- paste(msg, collapse=" ");
    dataSet$norm.msg <- current.msg;
    dataSet <<- dataSet;
    return(1);
}
######################################
PerformRarefaction <- function(data, rare.opt){
    data <- data.matrix(data);
    tax_nm <- rownames(data);

    # data must be count data (not contain fractions)        
    data <-round(data);
    
    # create phyloseq obj
    otu.tab <- otu_table(data,taxa_are_rows =TRUE);
    taxa_names(otu.tab)<-tax_nm;
    
    #random_tree<-phy_tree(createRandomTree(ntaxa(otu.tab),rooted=TRUE,tip.label=tax_nm));
    dataSet$sample_data$sample_id <- rownames(dataSet$sample_data);
    sample_table <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    phy.obj <- merge_phyloseq(otu.tab,sample_table);  
    msg <- NULL;
    
    #first doing rarefaction
    if(rare.opt=="rarewo"){
        phy.obj <- rarefy_even_depth(phy.obj, replace=FALSE,rngseed = T);
        msg <- c(msg, paste("Rarefaction without replacement to minimum library depth."));
    }
    data<- otu_table(phy.obj);
    return(data);
}
###############################################################################


###########Phyloseq object creation#############
#######################################
CreatePhyloseqObj <- function(){
    require(phyloseq);
    data.proc <- dataSet$data.proc;
    
    #standard name to be used
    #prepare data for phyloseq visualization.
    data.proc<-apply(data.proc,2,as.integer);
    
    #constructing phyloseq object for aplha diversity.
    data.proc<- otu_table(data.proc,taxa_are_rows = TRUE);
    taxa_names(data.proc)<- rownames(dataSet$data.proc);
    
    #sanity check: sample names of both abundance and metadata should match.
    smpl_var <- colnames(dataSet$sample_data);
    indx <- match(colnames(dataSet$data.proc), rownames(dataSet$sample_data));
    if(all(is.na(indx))){
        current.msg <<- "Please make sure that sample names matches in both the files.";
        return(0);
    }
    dataSet$sample_data <- sample_data(dataSet$sample_data, errorIfNULL=TRUE);
    taxa_table <- tax_table(as.matrix(dataSet$taxa_table));
    taxa_names(taxa_table) <- rownames(dataSet$taxa_table);
    dataSet$taxa_table <- taxa_table;
    dataSet$proc.phyobj <- phyloseq(data.proc,dataSet$sample_data,dataSet$taxa_table);
    
    #saveRDS(dataSet$proc.phyobj, file="proc.phyobj.orig");
    dataSet <<- dataSet;
    return(1);
}

###############################################################################

################################################################################


