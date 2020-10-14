# R script corresponding to publication: Felipe-Lucia et al. 2020. PNAS. Land-use intensity alters networks between biodiversity, ecosystem functions and services. 
# Author: Maria Felipe-Lucia, with contributions from Angel de Frutos and Caterina Penone.
# Contact: maria.felipe.lucia@gmail.com.

## SECTIONS: Set Rstudio for an overview of the script sections: Edit > Folding > Collapse all. #### 


## DATA REQUIRED ####
# Data to reproduce the analyses is available at https://www.bexis.uni-jena.de/PublicData/SearchPublicData.aspx 
# Use the ID number for the search. 
# Assembled datasets available: 
  # Grassland biodiversity (diversity data ID 21766; species information ID 21726), 
  # Grassland ecosystem functions and services (ID 27087), 
  # Forest biodiversity (diversity data ID 24607; species information ID 24608), 
  # Forest ecosystem functions and services (ID 24367).
# Individual datasets: see IDs provided in parentheses in the Supplementary Methods of the paper (Felipe-Lucia et al. 2020 PNAS).
# Please note that some datasets will only be available after an embargo period of around 3 years after data collection.

# Datasets and versions used in the analyses:
# ID 21766: here I used the 190411_EP_species_diversity_GRL.txt # check last update!!
# ID 21726: here I used the 190411_EP_species_info_GRL.txt      # check last update!!
# ID 21688: 21688_Assembled functions grassland EPs (2008-2011) for MF synth_1.5.21
# ID 27087 and 27088: I used May 2019 version      # (more updated than BExis (ID 21688) # check last update!!
# ID xxxx: butterflies_charismatic         # (from Gossner et al paper under review)
# ID 25086: LUI data from Bexis            # https://www.bexis.uni-jena.de/LuiTool/LuiTool.aspx?DatasetId=25086
  #input data: new datasets > standardized > regional > all Explos > years 2008-2015 > LUI overall > plots: EPs >>> (original file name is LUI_newSet_reg_ov_16.05.2019+171744)
# ID 20826: convertPlots                   # 20826_Basic Information of all Experimental Plots (EPs)_1.7.6
# ID xxxx: env_var_Grass                   # used by Soliveres et al. 2016 Nature.
# ID 24607: I used 180607_forestDiv.RAW.txt            # check last update!!
# ID 24608: I used 170313_forestSP_info.txt            # check last update!!
# ID 24426: 24426_Cercozoa and Endomyxa (Rhizaria, protists), Illumina Sequences, all EPs, grassland and forest) 2011_1.2.10
# ID 24367: forest_ES.raw                   # from Felipe-Lucia en al. 2018 Nat. Comm. 
# ID 17166: 17166_MinSoil_2011_Mineral_Soil_Enzymes_1.1.4
# ID 16466: 16466_ForMI - Forest Management Intensity Index_1.3.2
# ID 11603: 11603_GPslopeAspectElevation_6.3.1\\11603_forestEPVIP
# ID 20826: 20826_Basic Information of all Experimental Plots (EPs)_1.7.6
# ID 14686: 14686_MinSoil_2011_Mineral_Soil_Texture_1.9.6
# ID 14447: 14447_MinSoil_2011_Mineral_Soil_pH_1.10.12
# ID xxxxx: Forest_soil_depth              # assembled by Pete Manning/Ingo Schöning

# ID 21707 : multidiv: an R function to calculate multidiversity or multifunctionality


##  R FUNCTIONS used that are not implemented in R packages ####
 # Multifunctionality (multidiv.R) ####

# from Eric Allan Github, available at https://github.com/eric-allan/multidiversity.

#### a function to calculate multifunctionality
#### threshold can be a proportion of the maximum e.g. 0.5 or can be the median or mean
#### it can also be a vector of the same length as the number of diversities to
#### allow different thresholds for different diversities
#### threshold = FALSE calculates the mean of the scaled diversities or functions
#### scaling can be by any function specified in "sc", i.e. max, mean, sd etc., max is the default
#### if the maximum should be a mean of the top n values specify sc = "maxx", mx = n
#### centering by the mean is possible with cent = TRUE, to take z-scores of the diversities/processes, use sc="sd" and cent = TRUE 
#### "by" specifies a vector of the same length as nrow(x) to be used to split the data and use different
#### thresholds for the groups in "by"
#### "weights" allows different weightings for the different groups and should be a vector of the same length as the number of
#### groups (shorter vectors will be recycled)
#### to weight a diversity as 0 and drop it from the calculation, code the weight as "NA"
#### the default is weights = 1: all diversities weighted the same

### a function to calculate max based on several values (e.g. a mean of top 5) the number of values is specified by mx
maxx <- function(x, n, ...){     
  return(mean(rev(sort(x))[1:n], ...))}

### the multidiversity function
multidiv <- function(x, threshold=FALSE, sc = "max", mx = 1, cent = FALSE, by =FALSE, weights = 1){
  
  result <- matrix(nrow=nrow(x), ncol=2)
  weights <- rep_len(weights, length.out=ncol(x))   ### expand weights
  
  ### split the data on "by"
  if(any(by!=FALSE)){
    
    xs <- split(x, by)
    xst <- list()
    
    for(i in 1:length(unique(by))){
      if(mx > 1){
        xst[[i]] <- scale(xs[[i]], scale = apply(xs[[i]], 2, maxx, mx, na.rm=T), center = cent)}
      else{
        xst[[i]] <- scale(xs[[i]], scale = apply(xs[[i]], 2, match.fun(sc), na.rm=T), center = cent) 
      }}
    x.stand <- do.call("rbind", xst)
  }
  
  ### otherwise standardise globally
  else{
    if(mx > 1){
      x.stand <- scale(x, scale = apply(x, 2, maxx, mx, na.rm=T), center = cent)}
    else{
      x.stand <- scale(x, scale = apply(x, 2, match.fun(sc), na.rm=T), center = cent)
    }}
  
  #### sum diversities measured on each plot
  #### NA weighted diversities are dropped from the calculation
  x2 <- sweep(x, 2, weights, "*")  ## remove diversities with NA threshold from calc. of how many measured
  gm <- apply(x2, 1, function(x)(sum(complete.cases(x))))
  
  
  ### no threshold: average standardised values
  if(FALSE %in% threshold){  ### prevent error message if threshold is vector
    x.stand2 <- sweep(x.stand, 2, weights, "*")
    m <- apply(x.stand2, 1, mean, na.rm=T)
  }
  
  else{
    
    if (any(c("median","mean") %in% threshold)){  ### mean or median threshold
      
      tf <- match.fun(threshold)
      
      x.thresh <- apply(x.stand, 2, function(x)(1*(x > tf(x, na.rm=T))))
      gg <- apply(x.thresh, 1, sum, na.rm=T)
      m <- gg/gm
      
    }
    
    else{ 
      
      x.thresh <- 1*sweep(x.stand, 2, threshold, ">")  ### does each variable pass threshold?
      
      weights2 <- matrix(rep(weights, nrow(x.thresh)), nrow=nrow(x.thresh), byrow=TRUE)
      x.thresh <- x.thresh*weights2   ### multiply by weights
      
      gg <- apply(x.thresh, 1, sum,na.rm=T)
      m <- gg/gm
    }}
  
  result[, 1] <- m
  result[, 2] <- gm
  
  colnames(result)<-c("m", "groups measured")
  
  return(result)
}


 # Partial correlations (pcor.R) ####

# Contrary to the package pcor, this functions allows for NAs
# from http://www.yilab.gatech.edu/pcor.R or http://www.yilab.gatech.edu/pcor.html


#pcor test
pcor.test <- function(x,y,z,use="mat",method="p",na.rm=T){
  # The partial correlation coefficient between x and y given z
  #
  # pcor.test is free and comes with ABSOLUTELY NO WARRANTY.
  #
  # x and y should be vectors
  #
  # z can be either a vector or a matrix
  #
  # use: There are two methods to calculate the partial correlation coefficient.
  #	 One is by using variance-covariance matrix ("mat") and the other is by using recursive formula ("rec").
  #	 Default is "mat".
  #
  # method: There are three ways to calculate the correlation coefficient, 
  #	    which are Pearson's ("p"), Spearman's ("s"), and Kendall's ("k") methods.
  # 	    The last two methods which are Spearman's and Kendall's coefficient are based on the non-parametric analysis.
  #	    Default is "p".
  #
  # na.rm: If na.rm is T, then all the missing samples are deleted from the whole dataset, which is (x,y,z).
  #        If not, the missing samples will be removed just when the correlation coefficient is calculated.
  #	   However, the number of samples for the p-value is the number of samples after removing 
  #	   all the missing samples from the whole dataset.
  #	   Default is "T".
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(use == "mat"){
    p.use <- "Var-Cov matrix"
    pcor = pcor.mat(x,y,z,method=method,na.rm=na.rm)
  }else if(use == "rec"){
    p.use <- "Recursive formula"
    pcor = pcor.rec(x,y,z,method=method,na.rm=na.rm)
  }else{
    stop("\'use\' should be either \"rec\" or \"mat\"!\n")
  }
  
  # print the method
  if(gregexpr("p",method)[[1]][1] == 1){
    p.method <- "Pearson"
  }else if(gregexpr("s",method)[[1]][1] == 1){
    p.method <- "Spearman"
  }else if(gregexpr("k",method)[[1]][1] == 1){
    p.method <- "Kendall"
  }else{
    stop("\'method\' should be \"pearson\" or \"spearman\" or \"kendall\"!\n")
  }
  
  # sample number
  n <- dim(na.omit(data.frame(x,y,z)))[1]
  
  # given variables' number
  gn <- dim(z)[2]
  
  # p-value
  if(p.method == "Kendall"){
    statistic <- pcor/sqrt(2*(2*(n-gn)+5)/(9*(n-gn)*(n-1-gn)))
    p.value <- 2*pnorm(-abs(statistic))
    
  }else{
    statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
    p.value <- 2*pnorm(-abs(statistic))
  }
  
  data.frame(estimate=pcor,p.value=p.value,statistic=statistic,n=n,gn=gn,Method=p.method,Use=p.use)
}			

# By using var-cov matrix
pcor.mat <- function(x,y,z,method="p",na.rm=T){
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }
  
  data <- data.frame(x,y,z)
  
  if(na.rm == T){
    data = na.omit(data)
  }
  
  xdata <- na.omit(data.frame(data[,c(1,2)]))
  Sxx <- cov(xdata,xdata,m=method)
  
  xzdata <- na.omit(data)
  xdata <- data.frame(xzdata[,c(1,2)])
  zdata <- data.frame(xzdata[,-c(1,2)])
  Sxz <- cov(xdata,zdata,m=method)
  
  zdata <- na.omit(data.frame(data[,-c(1,2)]))
  Szz <- cov(zdata,zdata,m=method)
  
  # is Szz positive definite?
  zz.ev <- eigen(Szz)$values
  if(min(zz.ev)[1]<0){
    stop("\'Szz\' is not positive definite!\n")
  }
  
  # partial correlation
  Sxx.z <- Sxx - Sxz %*% solve(Szz) %*% t(Sxz)
  
  rxx.z <- cov2cor(Sxx.z)[1,2]
  
  rxx.z
}

# By using recursive formula
pcor.rec <- function(x,y,z,method="p",na.rm=T){
  # 
  
  x <- c(x)
  y <- c(y)
  z <- as.data.frame(z)
  
  if(dim(z)[2] == 0){
    stop("There should be given data\n")
  }
  
  data <- data.frame(x,y,z)
  
  if(na.rm == T){
    data = na.omit(data)
  }
  
  # recursive formula
  if(dim(z)[2] == 1){
    tdata <- na.omit(data.frame(data[,1],data[,2]))
    rxy <- cor(tdata[,1],tdata[,2],m=method)
    
    tdata <- na.omit(data.frame(data[,1],data[,-c(1,2)]))
    rxz <- cor(tdata[,1],tdata[,2],m=method)
    
    tdata <- na.omit(data.frame(data[,2],data[,-c(1,2)]))
    ryz <- cor(tdata[,1],tdata[,2],m=method)
    
    rxy.z <- (rxy - rxz*ryz)/( sqrt(1-rxz^2)*sqrt(1-ryz^2) )
    
    return(rxy.z)
  }else{
    x <- c(data[,1])
    y <- c(data[,2])
    z0 <- c(data[,3])
    zc <- as.data.frame(data[,-c(1,2,3)])
    
    rxy.zc <- pcor.rec(x,y,zc,method=method,na.rm=na.rm)
    rxz0.zc <- pcor.rec(x,z0,zc,method=method,na.rm=na.rm)
    ryz0.zc <- pcor.rec(y,z0,zc,method=method,na.rm=na.rm)
    
    rxy.z <- (rxy.zc - rxz0.zc*ryz0.zc)/( sqrt(1-rxz0.zc^2)*sqrt(1-ryz0.zc^2) )
    return(rxy.z)
  }			
}	


 # Raw correlations (flattenCorrMatrix) ####
  # Function to transform correlations matrices to "list" format (i.e. from/to/corr) 
  # http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software#correlation-matrix-with-significance-levels-p-value 
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
library(Hmisc)
  
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
   )
  }

 # Network metrics (developed specifically for the analyses of this paper) ####

    # Function to calculate network metrics; windows move 1 plot at a time until a window reaches the plot with highest LUI to avoid rolling over the low LUI plots.
    # if a particular block of networks give problems can be excluded modifying the function where indicated
    # usage:  metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=vertex_info_grass, estimate="POS",correlation="RAW")
            # xxxxxxxxxxx <- metrics$Select.Net.metrics_final
            # xxxxxxxxxxx <- metrics$Select.Node.metrics_final
            # xxxxxxxxxxx <- metrics$Select.Net.matrix_final
            # xxxxxxxxxxx <- metrics$Select.Net
            # xxxxxxxxxxx <- metrics$edges.list

    # np:  number of plots
    # mws: moving window size
    # nn = c(1:(np-(mws-1)))  # network number
    # data: your dataset
    # vertices: vertex info --> can be subset as Vertices=subset (vertex_info_grass,Node_type =="BD" | Node_type =="EF")
    # estimate: "POS" or "NEG" # select correlation type
    # correlation = "RAW"or "SIGN" or "PARTIAL"
       # if correlation = PARTIAL, requires function "pcor.R"
       # if correlation = RAW or SIGN, requires funtion "flattenCorrMatrix" + load "Himsc"

Network_metrics <- function(np,mws,data,Vertices,estimate,correlation) { ####
  nn = c(1:(np-(mws-1))) 
  Subset <- list () 
  for (i in nn){ 
    order_1 <- c(0:(mws-1))+i
    Subset[length(Subset)+1] <- list(order_1)
  }
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame( row = rownames(cormat)[row(cormat)[ut]],
                column = rownames(cormat)[col(cormat)[ut]],
                cor  =(cormat)[ut],
                p = pmat[ut])}
  
  edges.list<-list()
  pcor.test.temp <- NULL
  for (i in 1:length(Subset)){   
    if (correlation == "RAW") {
      temp_corr <- rcorr(as.matrix(data[Subset[[i]],],type="spearman"))
      temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
      ifelse (estimate=="POS",
              temp_corr <-temp_corr[(temp_corr$cor>0),],
              temp_corr <-temp_corr[(temp_corr$cor<0),])
      corr_list <- merge(temp_corr,Vertices,by.x="row",by.y="Nodes",all=FALSE)
      corr_list <- merge(corr_list,Vertices,by.x="column",by.y="Nodes",all=FALSE)
      
    } else if  (correlation == "SIGN") {  #significant correlations
      temp_corr <- rcorr(as.matrix(data[Subset[[i]],],type="spearman"))
      temp_corr<-flattenCorrMatrix(temp_corr$r, temp_corr$P)
      if (estimate=="POS") { 
        temp_corr <-temp_corr[(temp_corr$cor>0),]
        temp_corr <-temp_corr[(temp_corr$p<=0.05),]
      } else {    
        temp_corr <-temp_corr[(temp_corr$cor<0),]
        temp_corr <-temp_corr[(temp_corr$p<=0.05),]}
      corr_list <- merge(temp_corr,Vertices,by.x="row",by.y="Nodes",all=FALSE)
      corr_list <- merge(corr_list,Vertices,by.x="column",by.y="Nodes",all=FALSE)
      
    } else { # i.e. PARTIAL correlations 
      pcor.test.tabla <-NULL  
      for (z in 1:ncol(data)){ 
        for (j in 1:ncol(data)){
          if (z>j) {
            pcor.test.temp<-pcor.test(data[Subset[[i]],z],data[Subset[[i]],j],data[Subset[[i]],-c(z,j)], method="spearman", na.rm = F)
            pcor.test.temp$from<-names(data)[z]
            pcor.test.temp$to<-names(data)[j]
            pcor.test.tabla[length(pcor.test.tabla)+1] <- list(pcor.test.temp) 
          }
        }
      }
      pcor.test.tabla_temp <- Reduce(rbind, pcor.test.tabla) 
      pcor.test.tabla_temp <-pcor.test.tabla_temp[,c("from","to","estimate","p.value")]
      ifelse(estimate=="POS", 
             pcor.test.tabla_temp <-droplevels(subset(pcor.test.tabla_temp,estimate>0)),
             pcor.test.tabla_temp <-droplevels(subset(pcor.test.tabla_temp,estimate<0)))
      corr_list <- merge(pcor.test.tabla_temp,Vertices,by.x="from",by.y="Nodes",all=FALSE)
      corr_list <- merge(corr_list,Vertices,by.x="to",by.y="Nodes",all=FALSE)
    }
    colnames(corr_list)<- c("to","from","estimate","p.value", "Node_type.from", "Node_type.to")
    
    # drop internal links 
    corr_list_filtered <- droplevels(corr_list[!(corr_list$Node_type.from == "BD" & corr_list$Node_type.to == "BD" ),])
    corr_list_filtered <- droplevels(corr_list_filtered[!(corr_list_filtered$Node_type.from == "EF" & corr_list_filtered$Node_type.to == "EF" ),])
    corr_list_filtered <- droplevels(corr_list_filtered[!(corr_list_filtered$Node_type.from == "ES" & corr_list_filtered$Node_type.to == "ES" ),])
    
    # convert my matrix into an object readble by igraph
    edges<-corr_list_filtered[,c("to","from","estimate")]
    colnames(edges)<- c("from","to","weight")
    edges.list[length (edges.list)+1] <- list(edges)
    print (paste(i, "out of",length(Subset)))
  } # if is finished    
  
  # add window block
  for( w in seq_along(edges.list)){                
    edges.list[[w]]$window_block <- rep(w,nrow(edges.list[[w]]))
  }
  
  edges.list <- Reduce(rbind, edges.list)  
  
  Select.Net.metrics <-list()
  Select.Net.metrics_final<-list() 
  Select.Node.metrics <-list() 
  Select.Net.matrix<-list()
  Select.Net.matrix_final <-list()
  Select.Net <-list()
  
  for (q in nn){ #if a block gives problems, can be excluded (e.g. exclude block 8): nn[-8]
    Net <- graph_from_data_frame(d=subset(edges.list,window_block == q),vertices=Vertices, directed=F)
    
    # save the net object as matrix for checking
    Select.Net.matrix<- as_adjacency_matrix(Net, type="both", attr="weight",sparse=FALSE) 
    Select.Net.matrix_final[length (Select.Net.matrix_final)+1] <- list(Select.Net.matrix)
    
    # save the net object for plotting
    Select.Net[length (Select.Net)+1] <- list(Net) 
    
    ## METRICS: use absolute weights !
    degree<-  degree(Net, mode="all")  # degree: number of links of a node 
    density <-edge_density(Net, loops=F)   # density = connectance: proportion of present edges from all possible edges in the network
    weighted.density <- sum(abs(E(Net)$weight)) / ((vcount(Net)*(vcount(Net)-1))/2) # own adaptation inspired by StackOverflow
    modularity<-modularity(Net,membership(cluster_walktrap(Net)), weights=abs(E(Net)$weight))#(values>0.4 suggest that the network has a modular structure; Newman,2006).
    # evenness
    m<-Select.Net.matrix_final[[q]]
    df<-data.frame(row=rownames(m)[row(m)[upper.tri(m)]], 
                   col=colnames(m)[col(m)[upper.tri(m)]], 
                   corr=m[upper.tri(m)])
    df2 <- droplevels (df [!(df$corr==0),])
    
    H <- vegan::diversity(abs(df2 [,"corr"]))# with absolute values
    S <- dim (df2)[1]
    evenness <- H/log(S)
    
    # create the data frame for NETWORK LEVEL metrics: 
    Net.metrics <- c(weighted.density,modularity, evenness )
    indicators <- c("connectance","modularity", "evenness")
    names (Net.metrics)<- indicators
    Select.Net.metrics[length(Select.Net.metrics)+1] <- list(Net.metrics)
    
    # create the data frame for NODE LEVEL metrics: 
    temp.Node.metrics <- list(as.vector(degree))
    vertex <-  names (degree)                                           
    Node.metrics<- lapply(temp.Node.metrics, setNames, vertex)         
    nodes <- names (degree)                                                                          
    Node.metrics[length(Node.metrics)+1] <- list(nodes)                                                            
    Node.metrics[length(Node.metrics)+1] <- list(rep (q, lengths(Node.metrics)[1])) #info window_block                                      
    names (Node.metrics) <- c("degree","nodes","window_block") 
    Select.Node.metrics[length(Select.Node.metrics)+1] <- list(Node.metrics)
    
    edges<-NULL
    Net<-NULL
    Net.metrics<- NULL
    Node.metrics<- NULL
    Select.Net.matrix<-NULL
  }
  Select.Net.metrics_final <- Reduce(rbind, Select.Net.metrics)  
  Select.Net.metrics_final <- as.data.frame(Select.Net.metrics_final)
  Select.Net.metrics_final$window_block <- rep (nn) #if a block was excluded, adapt this (e.g. block 8 excluded): nn[-8]
  rownames (Select.Net.metrics_final) <-NULL  
  
  Select.Node.metrics_final <- rbindlist(Select.Node.metrics, fill=TRUE,use.names=TRUE)            
  Select.Node.metrics_final <- as.data.frame(Select.Node.metrics_final)
  #reorder columns
  Select.Node.metrics_final<-Select.Node.metrics_final[,c("window_block","nodes","degree")]
  Select.Node.metrics_final$nodes <- as.factor(Select.Node.metrics_final$nodes)
  
  output <- list (Select.Net.metrics_final,
                  Select.Node.metrics_final,
                  Select.Net.matrix_final,
                  Select.Net,
                  edges.list)
  names(output) <- c("Select.Net.metrics_final",
                     "Select.Node.metrics_final",
                     "Select.Net.matrix_final",
                     "Select.Net",
                     "edges.list")
  return (output)
} 


 # GAM Summary ####
k_val      <- function(x) (x$smooth [[1]]$bs.dim)
EDF_val    <- function(x) round (summary(x)$s.table [1],3)
F_val      <- function(x) round (summary(x)$s.table [3],3)
p.value_val<- function(x) round (summary(x)$s.table [4],3)
Adj.R2_val <- function(x) round (summary(x)$r.sq,3)
DE_val     <- function(x) round (summary(x)$dev.expl,3)
N_val      <- function(x) round (summary(x)$n,3)


 # AlignPlots ####

# Source: https://gist.github.com/scbrown86/3af317eb4ab692342e0ae46237dbace2
# Copied from https://stackoverflow.com/a/30414008/1710632 on 2018-02-08
# check also "align_plots()" from library(cowplot)

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}






## Libraries are provided in the respective sections ####

# 1. Assembling the DATASET ####

library(data.table) 
library(reshape)
library(vegan)
library(plyr)
library(dplyr)
library("scales")
source("multidiv.R")# see above and check updated code by Noelle Schenk.

  # 1.1 GRASSLANDS ====

    # 1.1.1 Biodiversity data ----
    #------------------------#

### Raw diversity
allsp_Grass <- fread("...\\Grasslands\\200205_EP_species_diversity_GRL.txt") # check last update!!
summary(allsp_Grass)
allsp_Grass$Plot_bexis<-as.factor(allsp_Grass$Plot_bexis)
allsp_Grass$Plot<-as.factor(allsp_Grass$Plot)
allsp_Grass$Species<-as.factor(allsp_Grass$Species)
allsp_Grass$type<-as.factor(allsp_Grass$type)
allsp_Grass$DataID<-as.factor(allsp_Grass$DataID)
allsp_Grass$Year<-as.factor(allsp_Grass$Year)
allsp_Grass$Dataversion<-as.factor(allsp_Grass$Dataversion)
summary(allsp_Grass)
levels(allsp_Grass$Year)

### Species information
fgs_Grass <- fread("...\\Grasslands\\200205_EP_species_info_GRL.txt")
summary(fgs_Grass)
fgs_Grass[] <- lapply(fgs_Grass, factor)
summary(fgs_Grass)

### merge functional groups with species file
names(allsp_Grass)
names(fgs_Grass)
allspfun_Grass <- merge(allsp_Grass, fgs_Grass, by ="Species", all=TRUE) 
summary(allspfun_Grass) 
summary(allspfun_Grass$Year) # year with more records is 2011

# check which groups have data for multiple years
unique(allspfun_Grass[,c("Year","Group_broad"),with=F])

### keep data for a single year: Remove all data not from 2011 or closest.
allspfun_Grass_clean <- droplevels(allspfun_Grass[!(Group_broad %in% "Plant" & !Year %in% 2011)])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "Birds" & !Year %in% 2011)])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "bacteria.RNA" & !Year %in% 2011)])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "Bats" & !Year %in% 2010)]) # there is not 2011 data, so 2010 is the closest year
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "Protists" & !Year %in% 2011)])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "soilfungi" & !Year %in% 2011)])

#final selection
unique(allspfun_Grass_clean[,c("Year","Group_broad"),with=F])

# check which groups are in each trophic level and there are not repeated groups
unique(allspfun_Grass_clean[,c("Trophic_level","Group_broad"),with=F])

# check if any sp is in two different datasets: 
summary(count(allspfun_Grass_clean,c("DataID","Species"))) 
# remove species with no info on Trophic level or Functional group (all sp should have info)
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!is.na(Trophic_level)]) 

#drop trophic groups with no ecological meaning: soilfungi.other + protist.unknown +  protist.parasite.nonplant (we don't know what they do)
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "soilfungi.other")])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "protist.unknown")])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "protist.parasite.nonplant")])

# remove some groups that will be used as ES
# birds (we consider SR of all birds as ES)
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "Birds")])

#drop groups with less than 10 sp OR less than 130 plots
# tertiary.consumer: bats (8 sp) because birds were dropped before
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Group_broad %in% "Bats")])
# belowground herbivores (5) and predators (6). 
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "belowground.herbivore")])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "belowground.predator")])
#omnivores (3 sp max), snail omnivores (5), and ant omnivores (10 and only in 110 plots) 
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "ant.omnivore")])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "omnivore")])
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "omnivore.snail")])

# merge herbivore (57 sp) + herbivore snail (19) --> rename the trophic group info
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Trophic_level=="herbivore.snail"]<- "herbivore"

# merge decomposers (5 sp) and detritivores (3)  --> rename the trophic group info
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Trophic_level=="detritivore"]<- "decomposer"

# merge myriapods to their corresponding trophic group
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Trophic_level=="myriapod.decomposer"]<- "decomposer"
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Trophic_level=="myriapod.secondary.consumer"]<- "secondary.consumer"

# drop empty levels
summary(allspfun_Grass_clean$Trophic_level)
allspfun_Grass_clean <- droplevels(allspfun_Grass_clean[!(Trophic_level %in% "herbivore.snail")])

#check correlations between lichen, moss and plants
# Calculate sum of values >1 (i.e. presence) by plot and trophic level
subset.autotrophs <- aggregate( (value>0) ~ Plot+Group_broad, data = subset(allspfun_Grass_clean,Trophic_level=="autotroph"), sum)
subset.autotrophs.cast<- cast(subset.autotrophs, Plot~ Group_broad)
cor(subset.autotrophs.cast)
# Separate lichen and mosses from plants --> rename the trophic group info
allspfun_Grass_clean$Group_broad <-as.character(allspfun_Grass_clean$Group_broad) 
allspfun_Grass_clean$Trophic_level <-as.character(allspfun_Grass_clean$Trophic_level)

allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Group_broad=="Lichen"]<- "Lichen"
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Group_broad=="Moss"] <- "Moss"
allspfun_Grass_clean$Trophic_level [allspfun_Grass_clean$Group_broad=="Plant"] <- "Plant"

#back-transform to factor
allspfun_Grass_clean$Trophic_level <-as.factor(allspfun_Grass_clean$Trophic_level)
summary(allspfun_Grass_clean)

# records per Trophic level, not species
summary(allspfun_Grass_clean$Trophic_level)
# cast which groups have only cover or presence/absence data: 
cast(Trophic_level ~ type , data=allspfun_Grass_clean) 

# Select only wanted columns and transform into a plot x trophic level matrix (with SR as value in each column)
allspfun_Grass_cleanSR <- allspfun_Grass_clean

# Calculate sum of values >1 (i.e. presence) by plot and trophic level
allspfun_Grass_cleanSR <- aggregate( (value>0) ~ Plot+Trophic_level, data = allspfun_Grass_cleanSR, sum)
allspfun_Grass_cleanSR<- cast(allspfun_Grass_cleanSR, Plot~ Trophic_level)
head(allspfun_Grass_cleanSR)
summary(allspfun_Grass_cleanSR)

# check normality of the data: some variables are not normal --> use Spearman correlations (alternative: transform data to meet normalization and use Pearson corr).
windows()
par(mfrow= c (4,4))     
for (i in 1:17) {
  hist(allspfun_Grass_cleanSR[,-1][,i],main=names(allspfun_Grass_cleanSR [,-1])[i],xlab=names(allspfun_Grass_cleanSR [,-1])[i])
}
par(op)

    # 1.1.2 EF+ES data ----
    #------------------#

#Total flower cover + Forage Quality + SR_birds from from Bexis: check updates!
EFESBexis <- read.table("...\\21688_Assembled functions grassland EPs (2008-2011) for MF synth_1.5.21\\21688.txt",header=TRUE, sep="\t", dec=".") 
summary(EFESBexis)

# EFES: Synthesis may 2019 version (check last updates in BExis) # note the .csv format
EFES_Grass <- read.table("...\\Grasslands\\may2019_grassland_functions.csv",header=TRUE, sep=";", dec=".") 
summary(EFES_Grass) 
names(EFES_Grass) 

# select columns 
EFES_Grass_clean <-subset(EFES_Grass, select = c(Plot,Plotn,Explo,
                                                 Biomass, # this is the average value
                                                 soilCflxs, Soil.C.stock, # drop SoilOrganicC
                                                 NRI,soilNflxs, #N_leaching_risk2015 should be the inverse of NRI
                                                 PRI,Phosphatase,# P_leaching_risk2015 should be the inverse of PRI
                                                 Groundwater.recharge,  # average value across years
                                                 herbivory.2013,  # transform to percentage log 
                                                 Root.biomass,Root.decomposition,Litter.decomposition,dung.decomposition,
                                                 Parasitoid.traps,pathogen.infection, # transform
                                                 Total.pollinators,
                                                 Aggregation,mAMFhyphae
))
summary(EFES_Grass_clean)

# Transform some data so larger values indicate more EF/ES                                                
##herbivory data : calculate herbivory control
EFES_Grass_clean$herbivory_control_log <-log(1/EFES_Grass_clean$herbivory.2013) # as pathogen regulation in Allan's Ecol.lett. paper
##pathogen.infection: calculate pathogen (pest) control
EFES_Grass_clean$pathogen_control_log <-log(1/EFES_Grass_clean$pathogen.infection) # as pathogen regulation in Allan's Ecol.lett. paper


#Charismatic butterflies as cultural ES (based on Gossner paper under review)
butterflies_charismatic <- read.table("...\\butterflies_charismatic.txt",header=TRUE, sep="\t", dec="." ) # check ID in BExis
summary(butterflies_charismatic)

#get all butterfly data from EPs
Lepidoptera <- droplevels(allspfun_Grass[(Group_broad %in% "Lepidoptera")])
summary(Lepidoptera)
#replace "_" by " " in the species names 
Lepidoptera$Species <- gsub('_', ' ', Lepidoptera$Species)
Lepidoptera$Species<-as.factor(Lepidoptera$Species)
summary(Lepidoptera$Species)
#unique list of sp
Lepidoptera_list<-unique (Lepidoptera[,"Species",with=F])
summary(Lepidoptera_list)

#merge
butterflies<- merge(butterflies_charismatic, Lepidoptera[,c(1:5)], by ="Species", all=F) # merge only the ones matching!
summary(butterflies)
dim(butterflies) #13426
#select the charismatic
butterflies_selected <-droplevels(subset(butterflies,Charismatic>0))

# ES indicator 1: abundance of charismatic butterfly sp per plot
library(reshape)
butterflies_selected <- cast(butterflies_selected,Plot~Species,sum)
butterflies_selected$Abund<-rowSums(butterflies_selected)
dim(butterflies_selected) # only 137 rows!

# ES indicator 2: richness of charismatic butterfly sp per plot
library(vegan)
names(butterflies_selected)
butterflies_selected$SR<-specnumber(butterflies_selected[,c(2:46)]) #drop first and last(Plot and Abund) 
butterflies_selected <-butterflies_selected[,c("Plot","Abund","SR")]
summary(butterflies_selected) # richness values per plot up to 12; abundance, max value =195...

# Merge the datasets
EFES_Grass_clean <- merge (EFES_Grass_clean,EFESBexis[,c("Plot","Total_flower_cover", "forage_quality", "SR_birds")],by.x="Plotn",by.y="Plot", all=TRUE) 
EFES_Grass_clean <- merge (EFES_Grass_clean,butterflies_selected[,c("Plot","Abund")],by.x="Plotn",by.y="Plot", all=TRUE) 
dim(EFES_Grass_clean)

names(EFES_Grass_clean)
# drop unwanted columns
EFES_Grass_clean <-subset(EFES_Grass_clean, select = -c(herbivory.2013, pathogen.infection))
names(EFES_Grass_clean)

# rename EFES
colnames(EFES_Grass_clean)<- c("Plotn","Plot","Exploratory","Biomass" ,"soilCflxs","Soil.C.stock" ,"NRI","soilNflxs" ,"PRI","Phosphatase","Groundwater.recharge", "Root.biomass","Root.decomp" ,"Litter.decomp", "Dung.decomp", "Parasitoid.traps","Total.pollinators" , "Soil.agg" ,"mAMFhyphae", "Herbivory.control","Pathogen.control","Flower.cover",    "Forage.quality", "Pot.bird_watching", "Butterfl.abund") # use "." instead of "-" to avoid problems with "mutate" below.
EFES_Grass_final <-EFES_Grass_clean

# reorder to sort first EF and then ES
EFES_Grass_final<-EFES_Grass_final [,c("Plot","Plotn","Exploratory",
                                       "Soil.agg" ,"mAMFhyphae",  #"Hyphal length"
                                       "Root.biomass","Root.decomp" ,"Litter.decomp", "Dung.decomp",
                                       "soilCflxs","NRI","soilNflxs" ,"PRI","Phosphatase",
                                       "Groundwater.recharge",
                                       "Herbivory.control","Pathogen.control","Parasitoid.traps","Total.pollinators" , 
                                       "Soil.C.stock" , "Biomass" ,   "Forage.quality",
                                       "Flower.cover",  "Pot.bird_watching",  "Butterfl.abund" )] 
summary(EFES_Grass_final)

# check normality of the data : no quite normal but we'll use Spearman correlations
windows() # maximize the window otherwise it doesn't fit
par(mfrow= c(5,5))     
for (i in 1:22) {
  hist(EFES_Grass_final[,4:25][,i],main=names(EFES_Grass_final [,4:25])[i],xlab=names(EFES_Grass_final [,4:25])[i])
}
par(op)
dev.off()


    # 1.1.3 CREATE JOIN DATASET ----
    #-----------------------------#

### merge sp+fg+EF+ES
names(allspfun_Grass_cleanSR)
spfgEFES_Grass <- merge(x=allspfun_Grass_cleanSR, y=EFES_Grass_final, by.x ="Plot", by.y ="Plotn",all=TRUE)
dim(spfgEFES_Grass)
head(spfgEFES_Grass)
summary(spfgEFES_Grass)# consider removing ES with lots of NAs (will depend on the correlations)
# mAMFhyphae: 57
# Soil.aggr: 57
# Parasitoid.traps: 67 
# Total.pollinators: 31
# Flower.cover: 80

# rename columns
names(spfgEFES_Grass)
colnames(spfgEFES_Grass) <- c("Plotn","Bacteria", "AG.decomp","AG.herb", "Lichen","Moss","Plant",  
                              "Plant.pathog", "pollinator","Protist.bacteriv","Protist.eukaryv","Protist.omniv", "Protist.parasite.plant",
                              "Snd.cons" ,"Soilfungi.decomp", "Soilfungi.pathot", "Soilfungi.symb",
                              "Plot","Exploratory",
                              "Soil.agg" ,"mAMFhyphae",
                              "Root.biomass","Root.decomp" ,"Litter.decomp", "Dung.decomp",
                              "soilCflxs","NRI","soilNflxs" ,"PRI","Phosphatase",
                              "Groundwater.recharge",
                              "Herbivory.control","Pathogen.control","Parasitoid.traps","Total.pollinators" , 
                              "Soil.C.stock" , "Biomass" ,   "Forage.quality",
                              "Flower.cover",  "Pot.bird_watching",  "Butterfl.abund")

# drop: pollinator + reorder BD to the logic order
spfgEFES_Grass_order <- spfgEFES_Grass[,c("Plot","Plotn","Exploratory",
                                          "Bacteria",
                                          "Protist.bacteriv","Protist.eukaryv","Protist.omniv","Protist.parasite.plant", 
                                          "Soilfungi.decomp",  "Soilfungi.pathot", "Soilfungi.symb",
                                          "Lichen","Moss","Plant", "Plant.pathog", 
                                          "AG.decomp","AG.herb", 
                                          "Snd.cons" ,
                                          "Soil.agg" ,"mAMFhyphae",
                                          "Root.biomass","Root.decomp" ,"Litter.decomp", "Dung.decomp",
                                          "soilCflxs","NRI","soilNflxs" ,"PRI","Phosphatase",
                                          "Groundwater.recharge",
                                          "Herbivory.control","Pathogen.control","Parasitoid.traps","Total.pollinators" , 
                                          "Soil.C.stock" , "Biomass" ,   "Forage.quality",
                                          "Flower.cover",  "Pot.bird_watching",  "Butterfl.abund"
)] 


    # 1.1.4 add LUI ----
    #---------------#
# get LUI data from Bexis: https://www.bexis.uni-jena.de/LuiTool/LuiTool.aspx?DatasetId=25086
# Bexis> Data>LUI Tool>input data: new datasets > standardized > regional > all Explos > years 2008-2015 > LUI overall > plots: EPs >>> (original file name is LUI_newSet_reg_ov_16.05.2019+171744)

LUI2008_2015 <- read.table("...\\Grasslands\\LUI_2008-2015.txt",header=TRUE, sep="\t", dec=",")
convert_plots<- read.table("...\\Grasslands\\convertPlots.txt",header=TRUE, sep="\t") # check ID
meanLUI <- merge(x=LUI2008_2015[,c("EP.Plotid","LUI")], y=convert_plots, by.x ="EP.Plotid", by.y ="Plot",all=TRUE)
head(meanLUI)
spfgEFES_Grass_order <- merge(x=meanLUI, y=spfgEFES_Grass_order, by.x="Plot0",by.y="Plotn",all=TRUE) #
# order plots by LUI 
spfgEFESLUI_Grass_order <-spfgEFES_Grass_order[order(spfgEFES_Grass_order$LUI),] 


    # 1.1.5 Take Residuals of Environmental variables + SCALE DATA ----
    #----------------------------------------------------------------#
# Filter by environmental variables: from Soliveres et al.2016 Nature paper: region, (!LUI), soil type and depth, pH, TWI (topographic wetness index), elevation

env_var_Grass<- read.table("...\\env_var_Grass.txt",header=TRUE, sep="\t", dec=".") # check ID
head(env_var_Grass)
env_var_Grass <- env_var_Grass[,c( "Plot", "Soil.type","pH", "Soil.depth", "elevation","TWI")]

# merge to the dataset
spfgEFESLUIenv_Grass <- merge(spfgEFESLUI_Grass_order,env_var_Grass, by.x= "Plot0",by.y = "Plot",all = TRUE )

# scale data between 0 and 1 (as in Soliveres et al. 2016 Nature paper)
library("scales") # only for individual vectors, need to combine with lapply
names(spfgEFESLUIenv_Grass)
#Note that the 'scale' function doesn't like to have a number as column name!!
data_Grass_scaled<-lapply(spfgEFESLUIenv_Grass[,-c(1:5,43)],rescale)# -c("Plot","Exploratory","Plot0","EP.Plotid","LUI","Soil.type"): Only for numeric values (not categorical)  
data_Grass_scaled <- as.data.frame (data_Grass_scaled)
data_Grass_scaled <- cbind(spfgEFESLUIenv_Grass[,c(1,4,5,3,43)],data_Grass_scaled) # bring categorical variables back and sorted
summary(data_Grass_scaled) 

# filter ALL data (BD-EF-ES) by the residuals of a model with only Explo+EnvVar as predictors
names(data_Grass_scaled)

var_names <- c("Bacteria", "Protist.bacteriv","Protist.eukaryv","Protist.omniv","Protist.parasite.plant", 
               "Soilfungi.decomp",  "Soilfungi.pathot", "Soilfungi.symb",
               "Lichen","Moss","Plant", "Plant.pathog", 
               "AG.decomp","AG.herb", 
               "Snd.cons" ,
               "Soil.agg" ,"mAMFhyphae",
               "Root.biomass","Root.decomp" ,"Litter.decomp", "Dung.decomp",
               "soilCflxs","NRI","soilNflxs" ,"PRI","Phosphatase",
               "Groundwater.recharge",
               "Herbivory.control","Pathogen.control","Parasitoid.traps","Total.pollinators" , 
               "Soil.C.stock" , "Biomass" ,   "Forage.quality",
               "Flower.cover",  "Pot.bird_watching",  "Butterfl.abund")

model_list_allRes <- list()
residual_list_allRes <- list()

for (i in 1:length (var_names)) { 
  lm(as.formula(paste(var_names [i],"~ Exploratory+Soil.type+pH+Soil.depth+elevation+TWI")),data=data_Grass_scaled) -> model_EnvVar_allRes
  resid(model_EnvVar_allRes)-> residual_EnvVar_allRes
  
  model_list_allRes[length(model_list_allRes)+1] <- list(model_EnvVar_allRes)     
  residual_list_allRes[length(residual_list_allRes)+1] <- list(residual_EnvVar_allRes)   
}

names(model_list_allRes) <- var_names
names(residual_list_allRes) <- var_names

# join to the initial table 
zz <- NULL
zz <- data.frame (c(1:150))
colnames(zz) <- "id"
for (i in 1:length(residual_list_allRes)) { 
  a <- as.data.frame (residual_list_allRes [i])
  a$id <- rownames (a)
  zz <- merge(zz,a, all = TRUE)
}

data_Grass_res <- data.frame (data_Grass_scaled[,c(1:4)],zz [,-1] ) # remove the env_var as we don't need them anymore , keep "Plot0","Plot","Exploratory","LUI" 
dim (data_Grass_res)
head (data_Grass_res)



  # 1.2 FORESTS ====

    # 1.2.1 Biodiversity data ----
    #------------------------#

# Raw diversity
allsp_Forest <- fread("...\\Forests\\180607_forestDiv.RAW.txt") # Check updates
summary(allsp_Forest)
allsp_Forest$plotID<-as.factor(allsp_Forest$plotID)
allsp_Forest$species<-as.factor(allsp_Forest$species)
allsp_Forest$type<-as.factor(allsp_Forest$type)
allsp_Forest$dataID<-as.factor(allsp_Forest$dataID)
allsp_Forest$year<-as.factor(allsp_Forest$year)
summary(allsp_Forest)

#replace "_" by "." in the species names only
allsp_Forest$species <- gsub('_', '.', allsp_Forest$species)
head(allsp_Forest)
levels(allsp_Forest$year)

### Species information
fgs_Forest <- fread("...\\Forests\\170313_forestSP_info.txt") # old version: check last update!!
summary(fgs_Forest)
fgs_Forest[] <- lapply(fgs_Forest, factor)
summary(fgs_Forest)
head(fgs_Forest)
#replace "_" by "." in the species names of the OTUs only
fgs_Forest$species <- gsub('_', '.', fgs_Forest$species)

### merge functional groups with species file
head(allsp_Forest)
head(fgs_Forest)
allspfun_Forest <- merge(allsp_Forest, fgs_Forest[,-"dataID"], by ="species", all=TRUE) #

# check which groups have data for multiple years
unique(allspfun_Forest[,c("year","Trophic.level"),with=F])
unique(allspfun_Forest[,c("year","group"),with=F]) # Note that lichens and bryophytes were recorded half each year, so both years are needed!!

cast(group ~ year , data=allspfun_Forest)
cast(Trophic.level ~ year , data=allspfun_Forest) 

# check which groups are in each trophic level
unique(allspfun_Forest[,c("Trophic.level","group"),with=F])

#Remove duplicate, chose 2011 as a baseline (Lichens & Mosses were recorded half in each year -Need data from both years)
allspfun_Forest_clean <- droplevels(allspfun_Forest[!(group %in% "arthropod" & !year %in% 2008)])
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "plant" & !year %in% 2011)])
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "fungi.deadw" & !year %in% 2010)]) # more data in 2010
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "fungi.root.soil" & !year %in% 2011)]) 
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "bat" & !year %in% 2010)])
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "bird" & !year %in% 2011)])
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "micromammmal" & !year %in% 2008)])# more data in 2008 # but there are only 10 sp

# drop vert.herb (micromammals) (only 6 sp!!)
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(Trophic.level %in% "vert.herb")])

#remove bats because there are too few sp and are the only tertiary.consumer # there are max 11 sp so better drop it 
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(Trophic.level %in% "tertiary.consumer")])

# remove some groups that will be used as ES (and are already in the other dataset)
# birds
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!(group %in% "bird")])
# pollinators from this dataset remain but will be used as EF 

#check correlations between lichen, moss and plants
# Calculate sum of values >1 (i.e. presence) by plot and trophic level
subset.autotrophsForest <- aggregate( (value>0) ~ plotID+group, data = subset(allspfun_Forest,Trophic.level=="autotroph"), sum)
subset.autotrophsForest.cast<- cast(subset.autotrophsForest, plotID~ group)
cor(subset.autotrophsForest.cast) 

# Separate lichen and mosses from plants --> rename the trophic group info
allspfun_Forest_clean$group <-as.character(allspfun_Forest_clean$group) #it gives error if a factor
allspfun_Forest_clean$Trophic.level <-as.character(allspfun_Forest_clean$Trophic.level) #it gives error if a factor

allspfun_Forest_clean$Trophic.level [allspfun_Forest_clean$group=="lichen"]<- "Lichen"
allspfun_Forest_clean$Trophic.level [allspfun_Forest_clean$group=="bryophyte"] <- "Moss"
allspfun_Forest_clean$Trophic.level [allspfun_Forest_clean$group=="plant"] <- "Plant"

#back transform this to factor
allspfun_Forest_clean$group <-as.factor(allspfun_Forest_clean$group) #it gives error if a factor
allspfun_Forest_clean$Trophic.level <-as.factor(allspfun_Forest_clean$Trophic.level)
summary(allspfun_Forest_clean) 

### remove species with no info on Throphic level (not really needed because all sp have info)
allspfun_Forest_clean <- droplevels(allspfun_Forest_clean[!is.na(Trophic.level)])

# Select only wanted columns and transform into a plot x trophic level matrix (with SR as value in each column)
# Calculate sum of values >1 (i.e. presence) by plot and trophic level
allspfun_Forest_cleanSR <- aggregate( (value>0) ~ plotID+Trophic.level, data = allspfun_Forest_clean, sum)
allspfun_Forest_cleanSR<- cast(allspfun_Forest_cleanSR, plotID~ Trophic.level)

# Adding protist data for 2011
all_protist <- fread("...\\24426_Cercozoa and Endomyxa (Rhizaria, protists), Illumina Sequences, all EPs, grassland and forest) 2011_1.2.10\\24426.txt") # 2011 data is twice!!
# remove duplicates
all_protist<-unique(all_protist)
dim(all_protist) #615593

# get sp info from the grasslands datase
protists_info<-droplevels(subset(fgs_Grass,Group_broad=="Protists"))
head(protists_info)
# match OTUs names: replace "_protist24426" by "" in the species names only
protists_info$Species <- gsub('_protist24426', '', protists_info$Species)
summary(protists_info)
protists_info$Species<-as.factor(protists_info$Species)
#drop the rest of oTUs (for 2017) - not needed if using all=FALSE in the merge
protists_info<- protists_info[!grepl("_protist24466", protists_info$Species),]

# merge sp info data
protist_2011 <- merge(all_protist, protists_info, by.x ="OTUs",by.y="Species", all=TRUE) # by =c("species","dataID"): gives many more NAs!!!!

#drop unwanted protist
summary(protist_2011$Trophic_level)
protist_2011 <- droplevels(protist_2011[!(Trophic_level %in% "protist.unknown")])
protist_2011 <- droplevels(protist_2011[!(Trophic_level %in% "protist.parasite.nonplant")])
protist_2011 <- droplevels(protist_2011[!(Trophic_level %in% "protist.mainly.bacterivore")])

# Calculate sum of values ("raw_abund") >1 (i.e. presence) by plot and trophic level
protist_2011agg <- aggregate( (raw_abund>0) ~ MyPlotID+Trophic_level, data = protist_2011, sum)
protist_2011agg<- cast(protist_2011agg, MyPlotID~ Trophic_level)
head(protist_2011agg)

# match plot names: replace "AEF" by "AEW" in all MyPlotID 
protist_2011agg$MyPlotID <- gsub('F', 'W', protist_2011agg$MyPlotID)
summary(protist_2011agg)
dim(protist_2011agg)

# drop grasslands datasets
protist_2011aggF<- protist_2011agg[!grepl("G", protist_2011agg$MyPlotID),]
dim(protist_2011aggF)

## merge all data
allspfun_Forest_cleanSR2 <- merge(allspfun_Forest_cleanSR, protist_2011aggF, by.x ="plotID",by.y="MyPlotID", all=TRUE) # Don't use by =c("species","dataID"), it brings more NAs!
summary(allspfun_Forest_cleanSR2)
dim(allspfun_Forest_cleanSR2) #150

##rename to match the Grasslands dataset
colnames(allspfun_Forest_cleanSR2)<- c( "Plot0","Bacteria","AG.decomp","AG.herb","Lichen","Moss", "AG.omniv",
                                  "Soilfungi.pathot", "Plant" , "Pollinator","Soilfungi.decomp", "Snd.cons","Soilfungi.symb",
                                  "Protist.bacteriv","Protist.eukaryv","Protist.omniv", "Protist.parasite.plant")

#reorder + KEEP pollinators (no other data)
allspfun_Forest_df <-allspfun_Forest_cleanSR2 [,c( "Plot0","Bacteria","Protist.bacteriv","Protist.eukaryv","Protist.omniv", "Protist.parasite.plant",
                                             "Soilfungi.pathot", "Soilfungi.decomp", "Soilfungi.symb",
                                             "Lichen","Moss","Plant","AG.decomp","AG.herb","AG.omniv","Pollinator",  
                                             "Snd.cons")]

# check normality of the data
dim(allspfun_Forest_df)
names(allspfun_Forest_df)
windows()
par(mfrow= c (4,4))     
for (i in 1:16) {
  hist(allspfun_Forest_df[,-1][,i],main=names(allspfun_Forest_df [,-1])[i],xlab=names(allspfun_Forest_df [,-1])[i])
}
par(op)
dev.off()


    # 1.2.2 add EF-ES data ----
    #--------------------------#
# from Felipe-Lucia et al. 2018 Nat.Comm. paper on Forest ES 
EFES_Forest <- read.table("...\\Forests\\forest_ES.raw.txt",header=TRUE, sep="\t", dec=".") ## update with dataset ID

enzymatic_act <- read.table("...\\17166_MinSoil_2011_Mineral_Soil_Enzymes_1.1.4\\17166.txt",header=TRUE, sep="\t", dec="." )

# calculate mini multifunctionality measures for the Carbon cycle enzymatic activities as with the grasslands data
# run "multidiv.R" function (see above and check updated code by Noelle Schenk) # check ID
enzymatic_act$soilCflxs <- multidiv(enzymatic_act[, c("Glu", "N_Ac", "Xyl")], sc="sd", cent=TRUE)[,1] # Note changes from default "threshold=FALSE,sc=maxx, mx =5"

## merge all data
EFES_Forest1 <- merge(x=EFES_Forest, y=enzymatic_act[,c("EP_Plotid","Pho","soilCflxs")], by.x ="Plot", by.y ="EP_Plotid",all=TRUE) 
dim(EFES_Forest1)

EFES_Forest_final <-EFES_Forest1 

colnames(EFES_Forest_final)<-c("Plot","Plot0","Exploratory","Root.decomp" , "Dung.decomp","Mycorrhiza", "N.aval","P.aval",
                               "Pest.control","Soil.C.stock" ,"Temp.reg" ,"Trees.C.stock" ,"Timber"  ,
                               "Edible.fungi" , "Edible.plants" ,"Cultural.plants","Pot.bird_watching","Phosphatase","soilCflxs") 

# check normality of the data: not quite normal -> use Spearman correlations
dim(EFES_Forest_final)
windows()
par(mfrow= c(4,4))     
for (i in 1:16) {
  hist(EFES_Forest_final[,-c(1:3)][,i],main=names(EFES_Forest_final [,-c(1:3)])[i],xlab=names(EFES_Forest [,-c(1:3)])[i])
}
par(op)
dev.off()


    # 1.2.3 CREATE JOIN DATASET ----
    #---------------------------#

### merge sp+fg+EF+ES
head(allspfun_Forest_df)
head(EFES_Forest_final)
spfgEFES_Forest <- merge(x=allspfun_Forest_df, y=EFES_Forest_final, by="Plot0",all=TRUE)
dim(spfgEFES_Forest)

summary(spfgEFES_Forest)# consider removing ES with lots of NAs
# Pollinator: 68
# Root_decomp: 15 
# Mycorrhiza: 39

# reorder columns to the logic order 
names(spfgEFES_Forest)
spfgEFES_Forest <- spfgEFES_Forest[,c( "Plot0","Plot","Exploratory","Bacteria","Protist.bacteriv","Protist.eukaryv","Protist.omniv", "Protist.parasite.plant",
                                       "Soilfungi.pathot","Soilfungi.decomp","Soilfungi.symb", "Lichen","Moss","Plant",
                                       "AG.decomp","AG.herb","AG.omniv","Pollinator","Snd.cons",
                                       "Root.decomp","Dung.decomp","Mycorrhiza","Phosphatase","soilCflxs", "N.aval", "P.aval" , 
                                       "Pest.control","Temp.reg","Soil.C.stock","Trees.C.stock","Timber",
                                       "Edible.fungi","Edible.plants","Cultural.plants","Pot.bird_watching")]


    # 1.2.4 add LUI ----
    #---------------#
ForMI <- read.table("...\\16466_ForMI - Forest Management Intensity Index_1.3.2\\16466.txt",header=TRUE, sep="\t", dec=".")
head(ForMI)

spfgEFESLUI_Forest <- merge(x=ForMI [,c("EP_Plotid", "ForMI")], y=spfgEFES_Forest, by.x ="EP_Plotid", by.y ="Plot",all=TRUE) #
dim(spfgEFESLUI_Forest)

# order plots by LUI (ForMI)
spfgEFESLUI_Forest_order <-spfgEFESLUI_Forest[order(spfgEFESLUI_Forest$ForMI),] 
head (spfgEFESLUI_Forest_order)
# check data is well distributed across the LUI gradient
windows()
plot (spfgEFESLUI_Forest_order$ForMI)


    # 1.2.5 take Residuals of Environmental variables ----
    #------------------------------------------------#
# Filter by environmental variables, as in Felipe-Lucia et al. 2018 Nat. Comm. paper: pH, %clay or soil type, soil depth, elevation, slope
# load environmental variables
Plot_info <- read.table("...\\11603_GPslopeAspectElevation_6.3.1\\11603_forestEPVIP.txt", sep="\t", header=TRUE) 
names(Plot_info)
Plot_info2 <- read.table("...\\20826_Basic Information of all Experimental Plots (EPs)_1.7.6\\20826.txt", sep="\t", header=TRUE) 
names(Plot_info2)
Soil_texture<- read.table("...\\14686_MinSoil_2011_Mineral_Soil_Texture_1.9.6\\14686.txt", sep="\t", header=TRUE) 
names(Soil_texture)
pH<- read.table("...\\14447_MinSoil_2011_Mineral_Soil_pH_1.10.12\\14447.txt", sep="\t", header=TRUE) 
head(pH)
pH$pH <- (pH$pH_1+pH$pH_2)/2
Soil_depth<- read.table("...\\Forest_soil_depth.txt", sep="\t", header=TRUE) # assembled by Pete Manning, no ID yet
summary(Soil_depth)

env_var1<-merge(x=Plot_info[,c("id","Exploratorium","EP_Plotid","rw","hw","elevation","elevationSd","slope","slopeSd","aspect","slopePercent")], y=Plot_info2[,c("EP_PlotID","SoilTypeGerman","SoilTypeWRB")], by.x="EP_Plotid", by.y="EP_PlotID", all = FALSE)
env_var2<-merge(x=env_var1, y=Soil_texture[,c("EP_Plotid","Clay","Fine_Silt","Medium_Silt","Coarse_Silt","Fine_Sand","Medium_Sand","Coarse_Sand")], by.x="EP_Plotid", by.y="EP_Plotid", all = FALSE)
env_var3<-merge(x=env_var2, y=pH[,c("EP_Plotid","pH")], by.x="EP_Plotid", by.y="EP_Plotid", all = FALSE)
env_var4<-merge(x=env_var3, y=Soil_depth, by.x="id", by.y="PlotID", all = FALSE)
dim(env_var4)
summary(env_var4)
env_var_Forest<-env_var4

# merge to the dataset
spfgEFESLUIenv_Forest <- merge(spfgEFESLUI_Forest_order,env_var_Forest [,c("EP_Plotid","elevation","SoilTypeGerman","pH","Soil_depth")],
                               by.x="EP_Plotid", by.y="EP_Plotid",all = TRUE )
dim(spfgEFESLUIenv_Forest)

# scale data --> scale data between 0 and 1 (as in Soliveres et al. 2016 Nature paper)
library("scales") # only for individual vectors, need to combine with lapply
names(spfgEFESLUIenv_Forest)
data_Forest_scaled<-lapply(spfgEFESLUIenv_Forest[,-c(1:4,38)],rescale) # exclude categorical variables ("EP_Plotid","ForMI","Plot0","Exploratory","SoilTypeGerman")
data_Forest_scaled <- as.data.frame (data_Forest_scaled)
data_Forest_scaled <- cbind(spfgEFESLUIenv_Forest[,c(1:4,38)],data_Forest_scaled)# add them back
summary(data_Forest_scaled) 

# filter ALL data (BD-EF-ES) by the residuals of a model with only Explo+EnvVar as predictors
names(data_Forest_scaled)

var_namesF <- c("Bacteria","Protist.bacteriv","Protist.eukaryv","Protist.omniv", "Protist.parasite.plant",
               "Soilfungi.pathot","Soilfungi.decomp","Soilfungi.symb", "Lichen","Moss","Plant",
               "AG.decomp","AG.herb","AG.omniv","Pollinator","Snd.cons",
               "Root.decomp","Dung.decomp","Mycorrhiza","Phosphatase","soilCflxs", "N.aval", "P.aval" , 
               "Pest.control","Temp.reg","Soil.C.stock","Trees.C.stock","Timber",
               "Edible.fungi","Edible.plants","Cultural.plants","Pot.bird_watching")

model_list_allRes <- list()
residual_list_allRes <- list()

for (i in 1:length (var_namesF)) { # EF + ES
  lm(as.formula(paste(var_namesF [i],"~ Exploratory+SoilTypeGerman+pH+Soil_depth+elevation")),data=data_Forest_scaled) -> model_EnvVar_allRes
  resid(model_EnvVar_allRes)-> residual_EnvVar_allRes
  
  model_list_allRes[length(model_list_allRes)+1] <- list(model_EnvVar_allRes)     
  residual_list_allRes[length(residual_list_allRes)+1] <- list(residual_EnvVar_allRes)   
}

names(model_list_allRes) <- var_namesF
names(residual_list_allRes) <- var_namesF

# join to the initial table 
zz <- NULL
zz <- data.frame (c(1:150))
colnames(zz) <- "id"
for (i in 1:length(residual_list_allRes)) {
  a <- as.data.frame (residual_list_allRes [i])
  a$id <- rownames (a)
  zz <- merge(zz,a, all = TRUE)
}

data_Forest_res <- data.frame (data_Forest_scaled[,c(1:4)],zz [,-1] ) # remove the env_var as we don't need them anymore 
dim (data_Forest_res)
head (data_Forest_res)

save.image(".../Mydata.RData")


  # 1.3 Create common "variables_names" for figures and tables ####

# Note that: 
  # "var_names" == Node_name for Grasslands only
  # "var_namesF" == Node_name for Forests only

Variable_name <- c("Bacteria","Protists bacterivore","Protists eukaryvore","Protists omnivore","Protists parasite plant","Soil fungi decomposer","Soil fungi pathotroph","Soil fungi symbiont","Lichens","Bryophytes","Vascular plants","Plant pathogens","Aboveground decomposers","Aboveground herbivores","Aboveground omnivores","Secondary consumers","Root biomass","Root decomposition","Soil C cycling","Nitrogen retention","Soil N cycling","Phosphorus retention","N availability","P availability","Phosphatase","Dung decomposition","Infiltration","Soil C stock","Trees C stock","Temperature regulation","Herbivory control","Pest control","Forage quality","Forage biomass","Timber","Edible fungi","Edible plants","Cultural value plants","Charismatic butterflies","Potential bird-watching")

Node_name <- c("Bacteria","Protist.bacteriv","Protist.eukaryv","Protist.omniv","Protist.parasite.plant","Soilfungi.decomp","Soilfungi.pathot","Soilfungi.symb","Lichen","Moss","Plant","Plant.pathog","AG.decomp","AG.herb","AG.omniv","Snd.cons","Root.biomass","Root.decomp","soilCflxs","NRI","soilNflxs","PRI","N.aval","P.aval","Phosphatase","Dung.decomp","Groundwater.recharge","Soil.C.stock","Trees.C.stock","Temp.reg","Herbivory.control","Pest.control","Forage.quality","Biomass","Timber","Edible.fungi","Edible.plants","Cultural.plants","Butterfl.abund","Pot.bird_watching")

variables_names <-data.frame(Variable_name,Node_name)

## clean RData and save only the final objects ----

rm(list=setdiff(ls(), c("data_Grass_res","data_Grass_res_order","var_names", "meanLUI",
                        "data_Forest_res","data_Forest_res_order","var_namesF", "ForMI","variables_names")))

save.image(".../Mydata.RData")



# 2.NETWORK level analyses ####

# Calculate Network metrics (connectance, modularity and evenness, plus node weighted degree) 
# using partial SPEARMAN correlations (untransformed data)
# and a moving window of 60 plots for grasslands and 50 plots for forests. 
  ## Warnings/Errors with smaller windows:
    #1: In sqrt(1/diag(V)): NaNs produced
    #2: In cov2cor(Sxx.z): diag(.) had 0 or NA entries; non-finite result is doubtful
    #3: In sqrt((n - 2 - gn)/(1 - pcor^2)): NaNs produced
    #4: Error in solve.default(Szz): system is computationally singular: reciprocal condition number = 5.45831e-19 
    #5: Error in pcor.mat(x, y, z, method = method, na.rm = na.rm): 'Szz' is not positive definite!

    # 2.1 GRASSLANDS data ====

      # 2.1.1 Prepare data ----
# Sort plots from low to high LUI
data_Grass_res_order <-data_Grass_res[order(data_Grass_res$LUI),] 
names (data_Grass_res_order)

# Drop data with more than 25 NAs
data_Grass_select<-subset(data_Grass_res_order, select = -c(Plot0,Plot, LUI,Exploratory,
                                                            Soil.agg,mAMFhyphae,Pathogen.control,Parasitoid.traps,
                                                            Flower.cover,Total.pollinators,Litter.decomp)) 

      ## create "vertex_info_grass" ONLY with used variables!! ----
vertex_info_grass <-data.frame(Nodes = names(data_Grass_select), Node_type = c(rep ("BD",15), rep ("EF",8), rep ("ES",7)))

      # 2.1.2 Calculate networks ----
library(igraph)
library(data.table)
library(vegan)
source("pcor.R") # PARTIAL CORRELATIONS that allows for NAs: http://www.yilab.gatech.edu/pcor.html 
# requires function "Network_metrics" (see above)

# np:  number of plots
# mws: moving window size
# nn = c(1:(np-(mws-1)))  # network number
# data: your dataset
# vertices: vertex info --> can be subset as Vertices=subset (vertex_info_grass,Node_type =="BD" | Node_type =="EF")
# estimate: "POS" or "NEG" # select correlation type
# correlation = "PARTIAL" or "RAW" --> if PARTIAL, run "pcor.R" function
# if a particular block of networks give problems can be excluded modifying the function where indicated

# Synergies networks (positive correlations)
metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=vertex_info_grass, estimate="POS",correlation="PARTIAL")
Select.Net.metrics_Grass_pos <- metrics$Select.Net.metrics_final
Select.Nodes_Grass_pos <- metrics$Select.Node.metrics_final
Select.Grass.Net.matrix_pos <- metrics$Select.Net.matrix_final
Select.Grass_pos.NET <- metrics$Select.Net
Select.Grass_pos.edges <- metrics$edges.list

# Trade-offs networks (negative correlations)
metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=vertex_info_grass, estimate="NEG",correlation="PARTIAL")
Select.Net.metrics_Grass_NEG <- metrics$Select.Net.metrics_final
Select.Nodes_Grass_NEG <- metrics$Select.Node.metrics_final
Select.Grass.Net.matrix_NEG <- metrics$Select.Net.matrix_final
Select.Grass.Net_NEG <- metrics$Select.Net
Select.Grass.edges_NEG <- metrics$edges.list


      # 2.1.3 Create variable "LUI_table" ----
#define parameters
np=150
mws=60
nn = c(1:(np-(mws-1))) 

# create the subsets
Subset <- list () 
LUI_values <- list ()
meanLUI <- NULL

for (i in 1:max(nn)){    
  order_1 <- c(0:(mws-1))+i #mws of 60 plots
  Subset[length(Subset)+1] <- list(order_1)
}
# calculate mean LUI for each subset
for (i in 1:max(nn)){    
  meanLUI <- mean(data_Grass_res_order[Subset[[i]],"LUI"]) 
  LUI_values[length(LUI_values)+1] <- list(meanLUI)
}
LUI_table <- Reduce(rbind, LUI_values)
LUI_table<-as.data.frame(LUI_table)
LUI_table$windows_block<-c(1:length(nn)) 
colnames(LUI_table) <-c("mean_LUI","window_block")
rownames(LUI_table)<-NULL

      # 2.1.4 Calculate D: node weighted degree ----
      #-----------------------------------------#
# weighting calculated multiplying by the mean correlation

mean_corr <- NULL
mean_corr_all<- list()
mean_corr_all2<- list()

for (i in 1:length(nn)){ 
  for (j in 1: ncol(Select.Grass.Net.matrix_pos[[i]])){ 
    mean_corr <- sum(Select.Grass.Net.matrix_pos [[i]][,j])/(ncol(Select.Grass.Net.matrix_pos[[i]])-1)
    mean_corr <- as.data.frame(mean_corr)
    mean_corr$node_name <- names(Select.Grass.Net.matrix_pos[[i]][,j])[j]
    mean_corr$window_block<- i 
    colnames(mean_corr)<-c("mean_corr","node_name","window_block")
    mean_corr_all[length(mean_corr_all)+1] <- list(mean_corr)
    
  }
  mean_corr_all <- Reduce(rbind, mean_corr_all)
  mean_corr_all <- as.data.frame(mean_corr_all)
  rownames (mean_corr_all) <- NULL
  mean_corr_all2[length(mean_corr_all2)+1] <- list(mean_corr_all)
  mean_corr_all<- list()
}
mean_corr_all2<- Reduce(rbind, mean_corr_all2)
mean_corr_all2 <- as.data.frame(mean_corr_all2)
head(mean_corr_all2)

# add mean_corr as a weight
Select.Nodes_Grass_pos <- merge(Select.Nodes_Grass_pos,mean_corr_all2,by.x=c("nodes","window_block"),by.y=c("node_name","window_block"))
# D is degree_2w
Select.Nodes_Grass_pos$degree_2w <-Select.Nodes_Grass_pos$degree*Select.Nodes_Grass_pos$mean_corr
Select.Nodes_Grass_pos$habitat<-"Grass"
Select.Nodes_Grass_pos$corr_type<-"pos"

#add LUI
Select.Nodes_Grass_pos<- merge(Select.Nodes_Grass_pos,LUI_table,by="window_block")
# add nodes info (BD-EF-ES)
Select.Nodes_Grass_pos <- merge(Select.Nodes_Grass_pos,vertex_info_grass,by.x="nodes",by.y="Nodes")



# save workspace
save.image(".../MyRData.RData")

    # 2.2 FORESTS ====
      # 2.2.1 Prepare data ----

# order plots by LUI (ForMI)
data_Forest_res_order <-data_Forest_res[order(data_Forest_res$ForMI),] 
head (data_Forest_res_order)
summary(data_Forest_res_order)
#drop columns more than 30 NAs
data_Forest_select<-subset(data_Forest_res_order, select = -c(EP_Plotid,ForMI,Plot0,Exploratory,
                                                              Pollinator, Mycorrhiza))

      ## create "vertex_info_Forest" ONLY with used variables!! ----
names(data_Forest_select)
vertex_info_Forest <-data.frame(Nodes = names(data_Forest_select), Node_type = c(rep ("BD",15), rep ("EF",6), rep ("ES",9)))

      # 2.2.2 Calculate networks  ----

library(igraph)
library(data.table)
library(vegan)
source("pcor.R") # PARTIAL CORRELATIONS that allows for NAs: http://www.yilab.gatech.edu/pcor.html 
# requires function "Network_metrics" (see above)

#Synergies networks (positive correlations)
metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=vertex_info_Forest, estimate="POS",correlation="PARTIAL")
Select.Net.metrics_Forest_pos <- metrics$Select.Net.metrics_final
Select.Node_Forest_pos <- metrics$Select.Node.metrics_final
Select.Forest.Net.matrix_pos <- metrics$Select.Net.matrix_final
Select.Forest_pos.NET <- metrics$Select.Net
Select.Forest_pos.edges <- metrics$edges.list

#Trade-offs networks (negative correlations)
metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=vertex_info_Forest, estimate="NEG",correlation="PARTIAL")
Select.Net.metrics_Forest_NEG <- metrics$Select.Net.metrics_final
Select.Node_Forest_NEG <- metrics$Select.Node.metrics_final
Select.Forest.Net.matrix_NEG <- metrics$Select.Net.matrix_final
Select.Forest.Net_NEG <- metrics$Select.Net
Select.Forest.edges_NEG <- metrics$edges.list

      # 2.2.3 Create variable "LUI_table_Forest" ----
#define parameters
np=150
mws=50
nn = c(1:(np-(mws-1))) 

# create the subsets
Subset <- list () 
LUI_values <- list ()
meanLUI <- NULL

for (i in 1:max(nn)){    
  order_1 <- c(0:(mws-1))+i #mws of 50 plots
  Subset[length(Subset)+1] <- list(order_1)
}
for (i in 1:max(nn)){    
  meanLUI <- mean(data_Forest_res_order[Subset[[i]],"ForMI"], na.rm=TRUE) # TAKING subset of the moving window!!
  LUI_values[length(LUI_values)+1] <- list(meanLUI)
}
LUI_table_Forest <- Reduce(rbind, LUI_values)
LUI_table_Forest<-as.data.frame(LUI_table_Forest)
LUI_table_Forest$windows_block<-c(1:length(nn)) 
colnames(LUI_table_Forest) <-c("mean_LUI","window_block")
rownames(LUI_table_Forest)<-NULL


      # 2.2.4 Calculate D: node weighted degree ----
      #-----------------------------------------#
mean_corr <- NULL
mean_corr_all<- list()
mean_corr_all2<- list()

for (i in 1:length(nn)){ 
  for (j in 1: ncol(Select.Forest.Net.matrix_pos[[i]])){ 
    mean_corr <- sum(Select.Forest.Net.matrix_pos [[i]][,j])/(ncol(Select.Forest.Net.matrix_pos[[i]])-1)
    mean_corr <- as.data.frame(mean_corr)
    mean_corr$node_name <- names(Select.Forest.Net.matrix_pos[[i]][,j])[j]
    mean_corr$window_block<- i 
    colnames(mean_corr)<-c("mean_corr","node_name","window_block")
    mean_corr_all[length(mean_corr_all)+1] <- list(mean_corr)
    
  }
  mean_corr_all <- Reduce(rbind, mean_corr_all)
  mean_corr_all <- as.data.frame(mean_corr_all)
  rownames (mean_corr_all) <- NULL
  mean_corr_all2[length(mean_corr_all2)+1] <- list(mean_corr_all)
  mean_corr_all<- list()
}
mean_corr_all2<- Reduce(rbind, mean_corr_all2)
mean_corr_all2 <- as.data.frame(mean_corr_all2)

# add mean corr as a weight
Selected.Nodes_Forest_pos <- merge(Select.Node_Forest_pos,mean_corr_all2,by.x=c("nodes","window_block"),by.y=c("node_name","window_block"))
# D is degree_2w
Selected.Nodes_Forest_pos$degree_2w <-Selected.Nodes_Forest_pos$degree*Selected.Nodes_Forest_pos$mean_corr
Selected.Nodes_Forest_pos$habitat<-"Forest"
Selected.Nodes_Forest_pos$corr_type<-"pos"

#add LUI
Selected.Nodes_Forest_pos<- merge(Selected.Nodes_Forest_pos,LUI_table_Forest,by="window_block")
# add nodes info (BD-EF-ES)
Selected.Nodes_Forest_pos <- merge(Selected.Nodes_Forest_pos,vertex_info_Forest,by.x="nodes",by.y="Nodes")


# save workspace
save.image(".../MyRData.RData")

    # 2.3 MODELS ====
      # 2.3.1 For synergy networks ----
#add LUI
Select.Net.metrics_Grass_pos<- merge(Select.Net.metrics_Grass_pos,LUI_table,by="window_block")
Select.Net.metrics_Forest_pos<- merge(Select.Net.metrics_Forest_pos,LUI_table_Forest,by="window_block")

## SELECTED GAMs MODELS: best fit with reasonable ecology behind - ie. only one up/down curve (see example of how to check the models below)
library(mgcv)
wdGg <- gam(connectance  ~ s(mean_LUI,k=6),data=Select.Net.metrics_Grass_pos) 
evGg2 <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_pos)
moGg2 <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_pos)

wdFg2 <- gam(connectance  ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_pos) 
evFg2 <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_pos)  #similar to Polynomial
moFg2 <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_pos)  #similar to Polynomial

# Example for how to check the models 
summary(wdGg1)
plot.gam(wdGg1)
par(mfrow=c(2,2))
gam.check(wdGg1)
AIC(wdGg1,wdGg)  

## including example of comparing GAMs with lm, and Polynomial
wdGPol1 <- lm(weighted.density ~ mean_LUI,data=subset(Select.Net.metrics_pos,habitat=="Grass"))
summary(wdGPol1)
wdGPol2 <- lm(weighted.density ~ poly(mean_LUI,2),data=subset(Select.Net.metrics_pos,habitat=="Grass"))
summary(wdGPol2) 
wdGPol3 <- lm(weighted.density ~ poly(mean_LUI,3),data=subset(Select.Net.metrics_pos,habitat=="Grass"))
summary(wdGPol3) 

anova(wdGPol1,wdGPol2,wdGPol3) 
AIC(wdGPol1,wdGPol2,wdGPol3) #m1

anova(wdGPol1,wdGg)
AIC(wdGPol1,wdGg) #gam 


      # 2.3.2 for trade-off networks----
#add LUI
Select.Net.metrics_Grass_NEG<- merge(Select.Net.metrics_Grass_NEG,LUI_table,by="window_block")
Select.Net.metrics_Forest_NEG<- merge(Select.Net.metrics_Forest_NEG,LUI_table_Forest,by="window_block")

## SELECTED GAMs MODELS: best fit with reasonable ecology behind - ie. only one up/down curve (see example of how to check the models below)
library(mgcv)
wdGg_NEG <- gam(connectance  ~ s(mean_LUI,k=6),data=Select.Net.metrics_Grass_NEG) 
evGg_NEG <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_NEG)
moGg_NEG <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_NEG)

wdFg_NEG <- gam(connectance  ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_NEG) 
evFg_NEG <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_NEG)  
moFg_NEG <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_NEG) 

# to check the models 
summary(moFg_NEG)
plot.gam(moFg_NEG)
par(mfrow=c(2,2))
gam.check(moFg_NEG)

    # 2.4 TABLE S2 ==== 

# requires functions from "Summary GAMs tables for SI"
# list final models from 2.3
models <- list(wdGg,moGg2,evGg2,wdFg2,moFg2,evFg2,wdGg_NEG,moGg_NEG,evGg_NEG,wdFg_NEG,moFg_NEG,evFg_NEG) 
Gams_PosNeg <- data.frame(Corr_type = rep(c("Positive", "Negative"), each=6),
                          Habitat = rep(c("Grassland", "Forest"), each=3),
                          Metric =rep(c("Connectance","Modularity", "Evenness"),4) ,
                          k =  unlist (lapply(models, k_val)),
                          EDF =   unlist (lapply(models, EDF_val)),
                          F =     unlist (lapply(models, F_val)),
                          p.value=unlist (lapply(models, p.value_val)),
                          Adj.R2 =unlist (lapply(models, Adj.R2_val)),
                          DE =    unlist (lapply(models, DE_val)),
                          N=      unlist (lapply(models, N_val)))
# save table for SI
write.table(Gams_PosNeg,"...\\Gams_PosNeg.txt",sep="\t", dec=".") 

    # 2.5 FIGURES ====
      # 2.5.1 Fig. 2 ----

library (gratia)
library(ggplot2)
library(mgcViz)

# prepare data and check significance
  # summary(wdGg)
wdGg_conf <-confint(wdGg, parm = "mean_LUI", type = "confidence")
wdGg_NEG_conf <-confint(wdGg_NEG, parm = "mean_LUI", type = "confidence")
moGg_conf <-confint(moGg2, parm = "mean_LUI", type = "confidence")
moGg_NEG_conf <-confint(moGg_NEG, parm = "mean_LUI", type = "confidence") # N.S.
evGg_conf <-confint(evGg2, parm = "mean_LUI", type = "confidence")
evGg_NEG_conf <-confint(evGg_NEG, parm = "mean_LUI", type = "confidence")
wdFg_conf <-confint(wdFg2, parm = "mean_LUI", type = "confidence")
wdFg_NEG_conf <-confint(wdFg_NEG, parm = "mean_LUI", type = "confidence")
moFg_conf <-confint(moFg2, parm = "mean_LUI", type = "confidence") # **
moFg_NEG_conf <-confint(moFg_NEG, parm = "mean_LUI", type = "confidence") # N.S.
evFg_conf <-confint(evFg2, parm = "mean_LUI", type = "confidence")
evFg_NEG_conf <-confint(evFg_NEG, parm = "mean_LUI", type = "confidence")

# PLOT THE FIGURE
plot1 <-ggplot( ) +
  labs(y ="",x="") + 
  geom_line(data=wdGg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = wdGg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdGg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = wdGg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  ggtitle("Grasslands")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01)) 

plot2 <-ggplot( ) +
  labs(y ="",x="") + 
  geom_line(data=moGg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = moGg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  theme_classic()  +
  geom_line(data=moGg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = moGg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))  

plot3 <-ggplot( ) +
  labs(y ="",x="Land-use intensity") + 
  geom_line(data=evGg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = evGg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evGg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = evGg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01)) 

plot4 <-ggplot( ) +
  labs(y ="Connectance (scaled)",x="") +
  geom_line(data=wdFg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = wdFg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdFg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = wdFg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  ggtitle("Forests")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot5 <-ggplot( ) +
  labs(y ="Modularity (scaled)",x="") + 
  geom_line(data=moFg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = moFg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  theme_classic()  +
  geom_line(data=moFg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = moFg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))  

plot6 <-ggplot( ) +
  labs(y ="Evenness (scaled)",x="Land-use intensity") + 
  geom_line(data=evFg_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = evFg_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evFg_NEG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "red") +
  geom_ribbon(data = evFg_NEG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="red",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01)) 

# Multiple graphs on the same page
windows ()
gridPrint(grobs = list(plot4,plot1,plot5,plot2,plot6,plot3), ncol = 2) # to have forest on the left and grasslands on the right columns

      # 2.5.2 Fig. S2----

library(mgcViz)
library(grid)
library(ggplot2)
library(gridExtra)

# A (synergy networks)
plot1 <- plot( sm(getViz(wdFg2), 1) ) + labs(y ="Connectance",x="")+ #ylim (-0.015,0.015)+
  l_ciLine( colour = "black", linetype = 2)+ #l_ciPoly(alpha=0.5) : colour the CI in grey
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+ 
  ggtitle("Forest")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot2 <- plot( sm(getViz(wdGg), 1) ) + labs(y = "",x="")+ #+ labs(y ="weighted.density") ylim (-0.015,0.015)+
  l_ciLine( colour = "black", linetype = 2)+#  l_ciPoly(alpha=0.5) +
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+
  ggtitle("Grassland")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot3 <- plot( sm(getViz(moFg2), 1) ) + labs(y ="Modularity",x="")+ #ylim (-0.15,0.15)+ 
  l_ciLine( colour = "black", linetype = 2)+ # l_ciPoly(alpha=0.5) +
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))

plot4 <- plot( sm(getViz(moGg2), 1) ) +labs(y ="",x="")+  #labs(y ="modularity")+ ylim (-0.15,0.15)+
  l_ciLine( colour = "black", linetype = 2)+ #l_ciPoly(alpha=0.5) +
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))

plot5 <- plot( sm(getViz(evFg2), 1) ) + labs(y ="Evenness")+ # ylim (-0.015,0.015)+ 
  l_ciLine( colour = "black", linetype = 2)+ #l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+labs(x = "Land-use intensity")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot6 <- plot( sm(getViz(evGg2), 1) ) + 
  l_ciLine( colour = "black", linetype = 2) + # ylim (-0.015,0.015)+ l_ciPoly(alpha=0.5)+ 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+labs(y = "")+labs(x = "Land-use intensity")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

# Multiple graphs on the same page
windows ()
gridPrint(grobs = list(plot1,plot2,plot3,plot4,plot5,plot6), ncol = 2)
      
# B (trade-offs networks)
plot1 <- plot( sm(getViz(wdFg_NEG), 1) ) + labs(y ="Connectance",x="")+
  l_ciLine( colour = "black", linetype = 2)+ #l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+ 
  ggtitle("Forest")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot2 <- plot( sm(getViz(wdGg_NEG), 1) ) + labs(y = "",x="")+ #ylim (-0.015,0.015)+
  l_ciLine( colour = "black", linetype = 2)+#l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+
  ggtitle("Grassland")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot3 <- plot( sm(getViz(moFg_NEG), 1) ) + labs(y ="Modularity",x="")+ 
  l_ciLine( colour = "black", linetype = 2)+#l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))  ## No sign

plot4 <- plot( sm(getViz(moGg_NEG), 1) ) +labs(y ="",x="")+ #ylim (-0.15,0.15)+ 
  l_ciLine( colour = "black", linetype = 2)+#l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))## No sign

plot5 <- plot( sm(getViz(evFg_NEG), 1) ) + labs(y ="Evenness")+
  l_ciLine( colour = "black", linetype = 2)+#l_ciPoly(alpha=0.5) + 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+labs(x = "Land-use intensity") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot6 <- plot( sm(getViz(evGg_NEG), 1) ) + 
  l_ciLine( colour = "black", linetype = 2) +#l_ciPoly(alpha=0.5)+ 
  l_points(shape = 19, size = 1, alpha = 0.5) + theme_classic()+l_fitLine(colour = "blue")+labs(y = "")+labs(x = "Land-use intensity")+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

# Multiple graphs on the same page
windows ()
gridPrint(grobs = list(plot1,plot2,plot3,plot4,plot5,plot6), ncol = 2)


# 3.NODE level analyses ####
    # 3.1 Models: LUI effects on D (weighted node degree) - Table S6 ====

# required objects from Section 2: Selected.Nodes_Grass_pos,Selected.Nodes_Forest_pos$nodes
# "weighted node degree" == "degree_2w"

#set contrasts
options(contrasts=c("contr.treatment", "contr.helmert"))

# Check interactions
baseG<-lm( degree_2w ~ mean_LUI, data= Select.Nodes_Grass_pos)
summary(baseG)
GNodes_allmodel1 <-glm( degree_2w ~ mean_LUI*nodes, data= Select.Nodes_Grass_pos) # sign
summary(GNodes_allmodel1)
GNodes_allmodel2 <-glm( degree_2w ~ mean_LUI+nodes, data= Select.Nodes_Grass_pos) # sign
summary(GNodes_allmodel2)
AIC(GNodes_allmodel1,GNodes_allmodel2) #interaction best!

baseF<-lm( degree_2w ~ mean_LUI, data= Selected.Nodes_Forest_pos)
summary(baseF)
FNodes_allmodel1 <-glm( degree_2w ~ mean_LUI*nodes, data= Selected.Nodes_Forest_pos) # sign
summary(FNodes_allmodel1) 
FNodes_allmodel2 <-glm( degree_2w ~ mean_LUI+nodes, data= Selected.Nodes_Forest_pos) # sign
summary(FNodes_allmodel2)
AIC(FNodes_allmodel2,FNodes_allmodel1) #interaction best!

# extract the coeff table
D_Grass <- as.data.frame(summary(GNodes_allmodel1)$coeff)
D_Grass$Habitat <- "Grassland"
D_Forest <-as.data.frame(summary(FNodes_allmodel1)$coeff) 
D_Forest$Habitat <- "Forest"
D_models <-rbind(D_Grass,D_Forest)

# save table for SI - Table S6
write.table(D_models,"...\\D_models.txt",sep="\t", dec=".") 

  
    # 3.2 Modules ====
      # 3.2.1 Prepare data ----
library(igraph)
library(reshape2)
      # 3.2.1.1 Grasslands----
# define the parameters
np = 150 # number of plots
mws = 60 # moving window size
nn = c(1:(np-(mws-1)))  # network number
N = length (V(Select.Grass_pos.NET[[min(nn)]])) # number of nodes - igraph object from Section 2
# "variables_names" created in Section 1
# "vertex_info_grass" created in Section 2

# get modules using the "cluster_louvain" algorithm (this gives a more balanced number of modules than "walktrap")
communities_low<-cluster_louvain(Select.Grass_pos.NET[[min(nn)]], weights = E(Select.Grass_pos.NET[[min(nn)]])$weight)
modules_low <- communities_low$membership 
communities_high<-cluster_louvain(Select.Grass_pos.NET[[max(nn)]], weights = E(Select.Grass_pos.NET[[max(nn)]])$weight)
modules_high <- communities_high$membership 

# create a dataframe 
a <- data.frame (communities_low$names,modules_low)
a <- a[order(a$communities_low),]
b <- data.frame (communities_high$names,modules_high)
b <- b[order(b$communities_high),]
Grass_modules<- cbind(a,b$modules_high)
colnames (Grass_modules) <- c("nodes","modules_low","modules_high")

Grass_modules_long <- reshape2::melt(Grass_modules, id="nodes")

#rename variables for the figure
Grass_modules_long$variable <-with (Grass_modules_long, ifelse (variable  =="modules_low","Low land-use intensity", "High land-use intensity"))
Grass_modules_long$variable <- factor (Grass_modules_long$variable, levels=c ("Low land-use intensity","High land-use intensity"))

# add nodes info
Grass_modules_long <- merge (Grass_modules_long,variables_names, by.x="nodes",by.y= "Node_name")
Grass_modules_long <- merge (Grass_modules_long,vertex_info_grass, by.x="nodes",by.y= "Nodes")
# colour nodes names
Grass_modules_long$color <- ifelse (Grass_modules_long$Node_type=="BD","grey15",ifelse(Grass_modules_long$Node_type=="EF","royalblue","indianred"))

# order the labels by LUI to use the annotate function
Grass_modules_long_low <- subset(Grass_modules_long,variable=="Low land-use intensity")
Grass_modules_long_high <- subset(Grass_modules_long,variable=="High land-use intensity")
Grass_modules_long_high <- cbind(Grass_modules_long_low,Grass_modules_long_high)

#order 1º by value + 2º by variable name at low LUI
Grass_modules_long_high<-Grass_modules_long_high[order(Grass_modules_long_high$value,Grass_modules_long_high$Variable_name),]
#order High Lui by Value
names(Grass_modules_long_high)
Grass_modules_long_high<- Grass_modules_long_high[,c(9:16)]
Grass_modules_long_high<-Grass_modules_long_high[order(Grass_modules_long_high$value),]
Grass_modules_long_high$Variable_name <- factor (Grass_modules_long_high$Variable_name, levels=Grass_modules_long_high$Variable_name)  

# Now, reorder LOW levels for the figure
Grass_modules_long_low <- Grass_modules_long_low[order(Grass_modules_long_low$value),]
Grass_modules_long_low$Variable_name <- factor (Grass_modules_long_low$Variable_name, levels=Grass_modules_long_low$Variable_name)

      # 3.2.1.2 Forests ----

# define the parameters (these were define in network metrics already)
np = 150 # number of plots
mws = 50 # moving window size
nn = c(1:(np-(mws-1)))  # network number
N = length (V(Select.Forest_pos.NET[[min(nn)]])) # number of nodes- igraph object from Section 2
# "variables_names" created in Section 1
# "vertex_info_grass" created in Section 2

# get modules using the "cluster_louvain" algorith
communities_low<-cluster_louvain(Select.Forest_pos.NET[[min(nn)]], weights = E(Select.Forest_pos.NET[[min(nn)]])$weight)
modules_low <- communities_low$membership 
communities_high<-cluster_louvain(Select.Forest_pos.NET[[max(nn)]], weights = E(Select.Forest_pos.NET[[max(nn)]])$weight)
modules_high <- communities_high$membership 

# create a dataframe 
a <- data.frame (communities_low$names,modules_low)
a <- a[order(a$communities_low),]
b <- data.frame (communities_high$names,modules_high)
b <- b[order(b$communities_high),]
Forest_modules<- cbind(a,b$modules_high)
colnames (Forest_modules) <- c("nodes","modules_low","modules_high")

Forest_modules_long <- reshape2::melt(Forest_modules, id="nodes")

#rename variables for the figure
Forest_modules_long$variable <-with (Forest_modules_long, ifelse (variable  =="modules_low","Low land-use intensity", "High land-use intensity"))
Forest_modules_long$variable <- factor (Forest_modules_long$variable, levels=c ("Low land-use intensity","High land-use intensity"))

# add nodes info
Forest_modules_long <- merge (Forest_modules_long,variables_names, by.x="nodes",by.y= "Node_name")
Forest_modules_long <- merge (Forest_modules_long,vertex_info_Forest, by.x="nodes",by.y= "Nodes")
# colour nodes names
Forest_modules_long$color <- ifelse (Forest_modules_long$Node_type=="BD","grey15",ifelse(Forest_modules_long$Node_type=="EF","royalblue","indianred"))

# order the labels by LUI to use the annotate function
Forest_modules_long_low <- subset(Forest_modules_long,variable=="Low land-use intensity")
Forest_modules_long_high <- subset(Forest_modules_long,variable=="High land-use intensity")
Forest_modules_long_high <- cbind(Forest_modules_long_low,Forest_modules_long_high)

#order 1º by value + 2º by variable name at low LUI
Forest_modules_long_high<-Forest_modules_long_high[order(Forest_modules_long_high$value,Forest_modules_long_high$Variable_name),]
#order High Lui by Value
names(Forest_modules_long_high)
Forest_modules_long_high<- Forest_modules_long_high[,c(9:16)]
Forest_modules_long_high<-Forest_modules_long_high[order(Forest_modules_long_high$value),] 
Forest_modules_long_high$Variable_name <- factor (Forest_modules_long_high$Variable_name, levels=Forest_modules_long_high$Variable_name)  

# Now, reorder LOW levels for the figure
Forest_modules_long_low <- Forest_modules_long_low[order(Forest_modules_long_low$value),]
Forest_modules_long_low$Variable_name <- factor (Forest_modules_long_low$Variable_name, levels=Forest_modules_long_low$Variable_name)


      # 3.2.2 Figure 3 ----
library("ggforce")
library (gtable)
library (mgcViz)                        
# add parameters required to plot the figure
Grass_modules_long <- Grass_modules_long [order (Grass_modules_long$variable,Grass_modules_long$nodes),]
Grass_modules_long$id <- c(rep(c(1:N),2)) # different number for each node
Grass_modules_long$idd <- 1 # fix number = 1
Grass_modules_long$value <- as.factor (Grass_modules_long$value ) #it needs to be a factor!
# calculate the parameters to annotate the text in the figure
 # default: ((1, [A],  ([A]+2.5), ([A]+2.5+[B]-0.5),   ([A]+2.5+[B]-0.5+2), ([A]+2.5+[B]-0.5+2+[C]-0.5)) - 0.5   
 rev (table(Grass_modules_long$variable,Grass_modules_long$value)[1,]) # LOW LUI (0 9 6 15)
 # c(1,9,     11.5,17,   19,33.5) -0.5 
 # c(0.5:8.5, 11.0:16.5, 18.5:33) 
 rev (table(Grass_modules_long$variable,Grass_modules_long$value)[2,])  # high LUI: 8 9 5 8  
# c(1,8,     10.5,19,   21,25.5,  27.5,35) -0.5 # info                              
# c(0.5:7.5, 10.0:18.5, 20.5:25.0, 27.0:34.5)  # this one to annotate      
 
 
Forest_modules_long <- Forest_modules_long [order (Forest_modules_long$variable,Forest_modules_long$nodes),]
Forest_modules_long$id <- c(rep(c(1:N),2)) # different number for each node
Forest_modules_long$idd <- 1 # fix number = 1
Forest_modules_long$value <- as.factor (Forest_modules_long$value ) #it needs to be a factor!
# calculate the parameters to annotate the text in the figure
# default: (((1, [A])-0.5,  (+2.5, +[B])-0.5,   (+2, +[C])-0.5) - 0.5  
rev (table(Forest_modules_long$variable,Forest_modules_long$value)[1,]) #(9 5 10 6)
# c(1,9,     11.5,16,   18,27.5,   29.5,35) -0.5 
# c(0.5:8.5, 11.0:15.5, 17.5:27.0, 29.0:34.5)  
rev (table(Forest_modules_long1$variable,Forest_modules_long1$value)[2,])  #high LUI: 14  9  7  
#c(1,14,     16.5,25,   27,33.5) -0.5 # info                                   
#c(0.5:13.5, 16.0:24.5, 26.5:33.0)    # to annotate

##
plotG_mods<- ggplot(Grass_modules_long, aes(variable, id = id, value = idd, split = value)) +
  geom_parallel_sets(aes(fill =Variable_name),alpha = 0.3, axis.width = 0.1) + 
  geom_parallel_sets_axes(axis.width = 0.1 ,fill="white",linetype=1,color="black") + 
  geom_parallel_sets_labels(angle = 0, size = 3,axis.width = 0.1,colour = "white")+ 
  theme_classic () +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
        axis.ticks.x=element_blank(),axis.line.x=element_blank())+ 
  annotate("text", x=0.9 ,y=c(0.5:8.5, 11.0:16.5, 18.5:33) , 
           label=rev(Grass_modules_long_low$Variable_name) ,hjust = 1,size=3, 
           color =rev( Grass_modules_long_low$color), fontface= 2) + #low lUI
  annotate("text", x=2.1 ,y=c(0.5:7.5, 10.0:18.5, 20.5:25.0, 27.0:34.5), 
           label=rev(Grass_modules_long_high$Variable_name) ,hjust = 0,size=3,
           color =rev( Grass_modules_long_high$color), fontface= 2) + # high LUI
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+ 
  ggtitle("Grasslands")+  
  scale_fill_manual(values=c(rep("grey30",15),rep("grey50",6),rep("grey80",9)))+ # number of items in each low LUI module
  scale_y_continuous(expand = c(0, 0.3))+ 
  scale_x_discrete(expand = expand_scale(add = c(0.7, 0.7))) 

plotF_mods<- ggplot(Forest_modules_long, aes(variable, id = id, value = idd, split = value)) +
  geom_parallel_sets(aes(fill =Variable_name),alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1 ,fill="white",linetype=1,color="black") + #width and colour of modules  
  geom_parallel_sets_labels(angle = 0, size = 3,axis.width = 0.1,colour = "white")+ #text modules stratum  
  theme_classic () +
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
        axis.ticks.x=element_blank(),axis.line.x=element_blank())+ #remove x and y axes
  annotate("text", x=0.9 ,y=c(0.5:8.5, 11.0:15.5, 17.5:27.0, 29.0:34.5) , 
           label=rev(Forest_modules_long_low$Variable_name) ,hjust = 1,size=3, 
           color =rev(Forest_modules_long_low$color), fontface= 2) + 
  annotate("text", x=2.1 ,y=c(0.5:13.5, 16.0:24.5, 26.5:33.0), 
           label=rev(Forest_modules_long_high$Variable_name) ,hjust = 0,size=3,
           color =rev( Forest_modules_long_high$color), fontface= 2) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(),# remove axis
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), # remove legend
        plot.title = element_text(hjust = 0.5))+ #center title
  ggtitle("Forests")+  # add title
  scale_fill_manual(values=c(rep("grey50",6),rep("grey30",10),rep("grey80",5),rep("grey1",9)))+ # to color alluviums
  scale_y_continuous(expand = c(0, 0.3))+ # expansion of Y axis
  scale_x_discrete(expand = expand_scale(add = c(0.7, 0.7)))  # distance from the modules to Y axis


#Both plots together
windows () 
gridPrint(grobs = list(plotG_mods,plotF_mods))


source (AlignPlots.r)
plots_Fig3<- AlignPlots(plotG_mods,plotF_mods)

tiff("...\\Fig3.modules.tiff", width = 7, height = 10, units = 'in', res = 300)  
#windows()
do.call(grid.arrange, plots_Fig3)
dev.off()


      # 3.2.3 Modules Table S5 ----
modules_grass <-Grass_modules_long [,c("variable", "value","Variable_name")]
modules_grass_low <- droplevels(subset(modules_grass,variable=="Low land-use intensity"))
modules_grass_low <- modules_grass_low[order(modules_grass_low$value),]
modules_grass_high <- droplevels(subset(modules_grass,variable=="High land-use intensity"))
modules_grass_high <- modules_grass_high[order(modules_grass_high$value),]
modules_grass <-cbind(modules_grass_low,modules_grass_high)
modules_grass$habitat <- "Grassland"

modules_forests <-Forest_modules_long [,c("variable", "value","Variable_name")]
modules_forests_low <- droplevels(subset(modules_forests,variable=="Low land-use intensity"))
modules_forests_low <- modules_forests_low[order(modules_forests_low$value),]
modules_forests_high <- droplevels(subset(modules_forests,variable=="High land-use intensity"))
modules_forests_high <- modules_forests_high[order(modules_forests_high$value),]
modules_forests <-cbind(modules_forests_low,modules_forests_high)
modules_forests$habitat <- "Forest"

modules_table<- rbind(modules_grass,modules_forests ) 
write.table(modules_table,"...\\modules_table.txt",sep="\t", dec=".") 


    # 3.3 Hubs ====
      # 3.3.1 Prepare data ----
library(plyr)
library(dplyr)
library (data.table)

# requires objects from network metrics: Selected.Nodes_Grass_pos, Selected.Nodes_Forest_pos

# for Grass_pos:
#--------------------#

##Check moving windows size and number of networks!! 
nn= 91 

# extract the node with maximum weighted degree ("degree_2w") by node type

## BD hubs
Hub_BD <-NULL
Hub_BD_selected <- list()

for (i in 1:nn){  
  Node_BD = Select.Nodes_Grass_pos  %>% filter(Node_type == "BD")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_BD = data.frame(Node_BD)
  Node_BD <- merge(Node_BD,subset(Select.Nodes_Grass_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_BD = Node_BD %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  Hub_BD_selected[length(Hub_BD_selected)+1] <- list(Hub_BD)
}

Hub_BD_selected <- Reduce(rbind, Hub_BD_selected)
Hub_BD_selected <- as.data.frame(Hub_BD_selected)
# drop nodes that are not hubs 
Hub_BD_selected<-droplevels (subset(Hub_BD_selected,degree_2w>0))

## EF hubs
Hub_EF <-NULL
Hub_EF_selected <- list()

for (i in 1:nn){  
  Node_EF = Select.Nodes_Grass_pos  %>% filter(Node_type == "EF")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_EF = data.frame(Node_EF)
  Node_EF <- merge(Node_EF,subset(Select.Nodes_Grass_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_EF = Node_EF %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  Hub_EF_selected[length(Hub_EF_selected)+1] <- list(Hub_EF)
}

Hub_EF_selected <- Reduce(rbind, Hub_EF_selected)
Hub_EF_selected <- as.data.frame(Hub_EF_selected)
# drop nodes that are not hubs
Hub_EF_selected<-droplevels (subset(Hub_EF_selected,degree_2w>0))

## ES hubs
Hub_ES <-NULL
Hub_ES_selected <- list()

for (i in 1:nn){  
  Node_ES = Select.Nodes_Grass_pos  %>% filter(Node_type == "ES")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_ES = data.frame(Node_ES)
  Node_ES <- merge(Node_ES,subset(Select.Nodes_Grass_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_ES = Node_ES %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  
  Hub_ES_selected[length(Hub_ES_selected)+1] <- list(Hub_ES)
}

Hub_ES_selected <- Reduce(rbind, Hub_ES_selected)
Hub_ES_selected <- as.data.frame(Hub_ES_selected)
# drop nodes that are not hubs 
Hub_ES_selected<-droplevels (subset(Hub_ES_selected,degree_2w>0))

## merge the three hubs
Selected.Hubs_Grass_pos<- rbind(Hub_BD_selected,Hub_EF_selected,Hub_ES_selected)
Selected.Hubs_Grass_pos$habitat<-"Grass"
Selected.Hubs_Grass_pos$corr_type<-"pos"
# add LUI
Selected.Hubs_Grass_posLUI <- merge (Selected.Hubs_Grass_pos,Select.Nodes_Grass_pos[,c("mean_LUI","nodes", "window_block")], by=c("nodes","window_block"), all =FALSE)


# for Forest_pos:
#--------------------#
##Check moving windows size!! 
nn= 101 

# extract the node with maximum weighted degree ("degree_2w") by node type
## BD hubs
Hub_BD <-NULL
Hub_BD_selected <- list()

for (i in 1:nn){  
  Node_BD = Selected.Nodes_Forest_pos  %>% filter(Node_type == "BD")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_BD = data.frame(Node_BD)
  Node_BD <- merge(Node_BD,subset(Selected.Nodes_Forest_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_BD = Node_BD %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  
  Hub_BD_selected[length(Hub_BD_selected)+1] <- list(Hub_BD)
}

Hub_BD_selected <- Reduce(rbind, Hub_BD_selected)
Hub_BD_selected <- as.data.frame(Hub_BD_selected)
# drop nodes that are not hubs
Hub_BD_selected<-droplevels (subset(Hub_BD_selected,degree_2w>0))

## EF hubs
Hub_EF <-NULL
Hub_EF_selected <- list()

for (i in 1:nn){  
  Node_EF = Selected.Nodes_Forest_pos  %>% filter(Node_type == "EF")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_EF = data.frame(Node_EF)
  Node_EF <- merge(Node_EF,subset(Selected.Nodes_Forest_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_EF = Node_EF %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  
  Hub_EF_selected[length(Hub_EF_selected)+1] <- list(Hub_EF)
}

Hub_EF_selected <- Reduce(rbind, Hub_EF_selected)
Hub_EF_selected <- as.data.frame(Hub_EF_selected)
# drop nodes that are not hubs 
Hub_EF_selected<-droplevels (subset(Hub_EF_selected,degree_2w>0))

## ES hubs
Hub_ES <-NULL
Hub_ES_selected <- list()

for (i in 1:nn){  
  Node_ES = Selected.Nodes_Forest_pos  %>% filter(Node_type == "ES")%>% filter(window_block == i)%>% count(Node_type, nodes)
  Node_ES = data.frame(Node_ES)
  Node_ES <- merge(Node_ES,subset(Selected.Nodes_Forest_pos,window_block==i)[,c("nodes","degree_2w","window_block")], by = "nodes")
  Hub_ES = Node_ES %>% filter(degree_2w, degree_2w == max(degree_2w), degree_2w == max(degree_2w))
  
  Hub_ES_selected[length(Hub_ES_selected)+1] <- list(Hub_ES)
}

Hub_ES_selected <- Reduce(rbind, Hub_ES_selected)
Hub_ES_selected <- as.data.frame(Hub_ES_selected)
# drop nodes that are not hubs
Hub_ES_selected<-droplevels (subset(Hub_ES_selected,degree_2w>0))


## merge the three hubs
Selected.Hubs_Forest_pos<- rbind(Hub_BD_selected,Hub_EF_selected,Hub_ES_selected)
Selected.Hubs_Forest_pos$habitat<-"Forest"
Selected.Hubs_Forest_pos$corr_type<-"pos"
# add LUI
Selected.Hubs_Forest_posLUI <- merge (Selected.Hubs_Forest_pos,Selected.Nodes_Forest_pos[,c("mean_LUI","nodes", "window_block")], by=c("nodes","window_block"), all=FALSE)


## merge F+G hubs
Selected.HubsLUI<-rbind(Selected.Hubs_Grass_posLUI,Selected.Hubs_Forest_posLUI)




      # 3.3.2. Figure 4 ====

library(ggplot2)
library(gridExtra)
library(plyr)
library(dplyr)
library (grid)
library (gtable)

# Run Source (00_AlignPlots.r) 

# requires object variable_names (created in 1.Dataset)
Selected.NodesLUI <- rbind(Select.Nodes_Grass_pos,Selected.Nodes_Forest_pos)
Selected.NodesLUI<- merge (Selected.NodesLUI,variables_names, by.x="nodes",by.y= "Node_name") 
Selected.HubsLUI <- merge (Selected.HubsLUI,variables_names, by.x="nodes",by.y= "Node_name") 

Selected.NodesLUI$habitat<-as.factor(Selected.NodesLUI$habitat)
Selected.HubsLUI$habitat<-as.factor(Selected.HubsLUI$habitat)

#rename habitats for the figure
head(Selected.NodesLUI)
Selected.NodesLUI$habitat <-with (Selected.NodesLUI, ifelse (habitat  =="Grass","Grasslands", "Forests"))
Selected.HubsLUI$habitat <-with (Selected.HubsLUI, ifelse (habitat  =="Grass","Grasslands", "Forests"))

#plot the fig.
plot_BD<-(subset(Selected.NodesLUI,Node_type=="BD")) %>% 
  ggplot(aes(x=mean_LUI,y=degree_2w,colour=factor(Variable_name)))+
  facet_grid(Node_type~habitat, scales="free") + #, scales="free_x": common y axis
  guides(color=guide_legend(title="Trophic groups"))+
  labs(x = NULL, y="")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=9,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(0.7, 'lines'),
        axis.text.x.bottom= element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.y = element_blank())+
  geom_smooth( method = "loess",span = 0.3, se= FALSE) #alt, span=0.8 or 0.3

plot_EF<-(subset(Selected.NodesLUI,Node_type=="EF")) %>% 
  ggplot(aes(x=mean_LUI,y=degree_2w,colour=factor(Variable_name)))+
  facet_grid(Node_type~habitat, scales="free") + #, scales="free_x": common y axis
  guides(color=guide_legend(title="Ecosystem functions"))+
  labs(y="Node degree (D)")+ 
  labs(x="")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=9,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(0.7, 'lines'),
        axis.text.x.bottom= element_blank(),
        axis.ticks.x= element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+
  geom_smooth( method = "loess",span = 0.3, se= FALSE) #method="gam", formula = y ~s(x);span=0.3 or 0.8

plot_ES<-(subset(Selected.NodesLUI,Node_type=="ES")) %>% 
  ggplot(aes(x=mean_LUI,y=degree_2w,colour=factor(Variable_name)))+
  facet_grid(Node_type~habitat, scales="free") + #, scales="free_x": common y axis
  guides(color=guide_legend(title="Ecosystem services"))+
  labs(y="")+
  labs(x="Land-use intensity")+
  theme_bw() +  #theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=9,face="bold"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(0.7, 'lines'),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())+
  geom_smooth( method = "loess",span = 0.3, se= FALSE) #method="gam", formula = y ~s(x); span=0.3 or 0.8

windows()
plots_Fig4 <- AlignPlots(plot_BD,plot_EF,plot_ES) 
do.call(grid.arrange, plots_Fig4)


tiff("...\\Fig4.hubs.tiff", width = 8, height = 10, units = 'in', res = 500)  

dev.off()




# 4.SENSITIVITY ANALYSES ####
###############################.


    # 4.1 Strength of the correlations ==== 
      # 4.1.1 Prepare data ---- 

# add LUI data to edge list from Section 2
corr_GrassNEG<- merge(Select.Grass.edges_NEG,LUI_table,by="window_block")
corr_GrassPOS<- merge(Select.Grass_pos.edges,LUI_table,by="window_block")
corr_ForestNEG<- merge(Select.Forest.edges_NEG,LUI_table_Forest,by="window_block")
corr_ForestPOS<- merge(Select.Forest_pos.edges,LUI_table_Forest,by="window_block")

# merge pos+neg
corr_GrassNEG$corr_type<-"NEG"
corr_GrassPOS$corr_type<-"POS"
corr_Grass<- rbind(corr_GrassPOS,corr_GrassNEG)
corr_Grass$corr_type<-as.factor(corr_Grass$corr_type)

corr_ForestNEG$corr_type<-"NEG"
corr_ForestPOS$corr_type<-"POS"
corr_Forest<- rbind(corr_ForestPOS,corr_ForestNEG)
corr_Forest$corr_type<-as.factor(corr_Forest$corr_type)

# add link_type (BD-EF / BD-EF / EF-ES): match nodes to vertex_info
corr_GrassALL<-corr_Grass 
pair_Grass <- corr_GrassALL [,c("from","to")]
pair_Grass [] <- vertex_info_grass$Node_type[match(unlist(pair_Grass), vertex_info_grass$Nodes)]
pair_Grass$pair <- paste0(pair_Grass$from, "-", pair_Grass$to)
pair_Grass[pair_Grass=="ES-EF"]<-"EF-ES"
pair_Grass[pair_Grass=="ES-BD"]<-"BD-ES"
corr_GrassALL <- cbind (corr_GrassALL,pair=pair_Grass[,"pair"])
head (corr_GrassALL)

corr_ForestALL<-corr_Forest 
pair_Forest <- corr_ForestALL [,c("from","to")]
pair_Forest [] <- vertex_info_Forest$Node_type[match(unlist(pair_Forest), vertex_info_Forest$Nodes)]
pair_Forest$pair <- paste0(pair_Forest$from, "-", pair_Forest$to)
pair_Forest[pair_Forest=="ES-EF"]<-"EF-ES"
pair_Forest[pair_Forest=="ES-BD"]<-"BD-ES"
corr_ForestALL <- cbind (corr_ForestALL,pair=pair_Forest[,"pair"])
head (corr_ForestALL)

      # 4.1.2 Models (GAMs) ----
library(mgcv)
# see in section 2 how to check GAMs
summary(gam_Forest_NEG_EF.ES) # Deviance explained < 10% in all cases
gam_GRASS_POS<-gam(weight ~ s(mean_LUI,k=6),data=subset(corr_GrassALL,corr_type=="POS")) 
gam_GRASS_NEG<-gam(weight ~ s(mean_LUI,k=4),data=subset(corr_GrassALL,corr_type=="NEG")) 

gam_GRASS_POS_BD.EF<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_GrassALL, corr_type=="POS" & pair=="BD-EF"))
gam_GRASS_POS_BD.ES<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_GrassALL, corr_type=="POS" & pair=="BD-ES")) 
gam_GRASS_POS_EF.ES<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_GrassALL, corr_type=="POS" & pair=="EF-ES")) 

gam_GRASS_NEG_BD.EF<-gam(weight ~ s(mean_LUI,k=4), data=subset (corr_GrassALL, corr_type=="NEG" & pair=="BD-EF"))
gam_GRASS_NEG_BD.ES<-gam(weight ~ s(mean_LUI,k=4), data=subset (corr_GrassALL, corr_type=="NEG" & pair=="BD-ES")) #n.s.
gam_GRASS_NEG_EF.ES<-gam(weight ~ s(mean_LUI,k=4), data=subset (corr_GrassALL, corr_type=="NEG" & pair=="EF-ES"))

gam_Forest_POS<-gam(weight ~ s(mean_LUI,k=6),data=subset(corr_ForestALL,corr_type=="POS")) 
gam_Forest_NEG<-gam(weight ~ s(mean_LUI,k=8),data=subset(corr_ForestALL,corr_type=="NEG")) 

gam_Forest_POS_BD.EF<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_ForestALL, corr_type=="POS" & pair=="BD-EF"))
gam_Forest_POS_BD.ES<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_ForestALL, corr_type=="POS" & pair=="BD-ES")) 
gam_Forest_POS_EF.ES<-gam(weight ~ s(mean_LUI,k=6), data=subset (corr_ForestALL, corr_type=="POS" & pair=="EF-ES")) 

gam_Forest_NEG_BD.EF<-gam(weight ~ s(mean_LUI,k=8), data=subset (corr_ForestALL, corr_type=="NEG" & pair=="BD-EF"))
gam_Forest_NEG_BD.ES<-gam(weight ~ s(mean_LUI,k=8), data=subset (corr_ForestALL, corr_type=="NEG" & pair=="BD-ES"))
gam_Forest_NEG_EF.ES<-gam(weight ~ s(mean_LUI,k=8), data=subset (corr_ForestALL, corr_type=="NEG" & pair=="EF-ES"))


      # 4.1.3 Figure (Fig. S2)----

#Grassland
windows ()                             
Grass_plot<- plot(weight~mean_LUI, data=corr_GrassALL,col="transparent", ylab="Partial Spearman correlation (rho)",xlab="Land-use intensity",ylim=c(-0.3,0.3)) 
lines(subset (corr_GrassALL, corr_type=="POS")$mean_LUI, fitted(gam_GRASS_POS), col="blue",lwd=2) 
lines(subset (corr_GrassALL, corr_type=="NEG")$mean_LUI, fitted(gam_GRASS_NEG), col="red",lwd=2) 
abline(h=0, col="gray50")

lines(subset (corr_GrassALL, corr_type=="POS" & pair=="BD-EF")$mean_LUI, fitted(gam_GRASS_POS_BD.EF), col="royalblue",lty=2) 
lines(subset (corr_GrassALL, corr_type=="POS" & pair=="BD-ES")$mean_LUI, fitted(gam_GRASS_POS_BD.ES), col="royalblue",lty=3) 
lines(subset (corr_GrassALL, corr_type=="POS" & pair=="EF-ES")$mean_LUI, fitted(gam_GRASS_POS_EF.ES), col="royalblue",lty=4) 

lines(subset (corr_GrassALL, corr_type=="NEG" & pair=="BD-EF")$mean_LUI, fitted(gam_GRASS_NEG_BD.EF), col="indianred",lty=2) 
lines(subset (corr_GrassALL, corr_type=="NEG" & pair=="BD-ES")$mean_LUI, fitted(gam_GRASS_NEG_BD.ES), col="indianred",lty=3) 
lines(subset (corr_GrassALL, corr_type=="NEG" & pair=="EF-ES")$mean_LUI, fitted(gam_GRASS_NEG_EF.ES), col="indianred",lty=4) 

legend("topright", legend=c("Positive", "Negative", "BD-EF","BD-ES","EF-ES"),    ## inset=c(-0.24,0),)
       lty=c(1,1,2,3,4), cex=0.8,col=c("blue","red", "black", "black", "black"),
       box.lty=0)

# Forest
windows ()                             
Forest_plot<- plot(weight~mean_LUI, data=corr_ForestALL,col="transparent", ylab="Partial Spearman correlation (rho)",xlab="Land-use intensity",ylim=c(-0.3,0.3)) 
lines(subset (corr_ForestALL, corr_type=="POS")$mean_LUI, fitted(gam_Forest_POS), col="blue",lwd=2) 
lines(subset (corr_ForestALL, corr_type=="NEG")$mean_LUI, fitted(gam_Forest_NEG), col="red",lwd=2) 
abline(h=0, col="gray50")

lines(subset (corr_ForestALL, corr_type=="POS" & pair=="BD-EF")$mean_LUI, fitted(gam_Forest_POS_BD.EF), col="royalblue",lty=2) 
lines(subset (corr_ForestALL, corr_type=="POS" & pair=="BD-ES")$mean_LUI, fitted(gam_Forest_POS_BD.ES), col="royalblue",lty=3) 
lines(subset (corr_ForestALL, corr_type=="POS" & pair=="EF-ES")$mean_LUI, fitted(gam_Forest_POS_EF.ES), col="royalblue",lty=4) 

lines(subset (corr_ForestALL, corr_type=="NEG" & pair=="BD-EF")$mean_LUI, fitted(gam_Forest_NEG_BD.EF), col="indianred",lty=2) 
lines(subset (corr_ForestALL, corr_type=="NEG" & pair=="BD-ES")$mean_LUI, fitted(gam_Forest_NEG_BD.ES), col="indianred",lty=3) 
lines(subset (corr_ForestALL, corr_type=="NEG" & pair=="EF-ES")$mean_LUI, fitted(gam_Forest_NEG_EF.ES), col="indianred",lty=4) 

legend("topright", legend=c("Positive", "Negative", "BD-EF","BD-ES","EF-ES"),    # inset=c(-0.24,0),)
       lty=c(1,1,2,3,4), cex=0.8,col=c("blue","red", "black", "black", "black"),
       box.lty=0)

# add multiplot
windows () 
gridPrint(grobs = list(Grass_plot,Forest_plot))


    # 4.2 Number of correlations ====
      # 4.2.1 Prepare data ----

corr_GrassALL$group2 <- paste0(corr_GrassALL$corr_type,"_",corr_GrassALL$pair)
corr_GrassALL$group2 <-as.factor(corr_GrassALL$group2)
count(corr_GrassALL, mean_LUI, by =group2) -> corr_countGrass

corr_ForestALL$group2 <- paste0(corr_ForestALL$corr_type,"_",corr_ForestALL$pair)
corr_ForestALL$group2 <-as.factor(corr_ForestALL$group2)
count(corr_ForestALL, mean_LUI, by =group2) -> corr_countForest

#merge both for ggplot
corr_countGrass$habitat<-"Grass"
corr_countForest$habitat<-"Forest"
corr_countALL<- rbind(corr_countGrass,corr_countForest)

# add corr_type and Link_type
corr_countALL$corr_type <-gsub("_.*$","", corr_countALL$by) # keep all characters before "_"
corr_countALL$link_type <-gsub(".*_","", corr_countALL$by) # keep all characters after "_"

#rename columns
colnames(corr_countALL)<- c("mean_LUI", "group2","n", "habitat","corr_type", "link_type")

# calculate the relative number of correlations 
##n_corrected: divided by the max number of correlations 
corr_countALL$n_corrected <-ifelse(corr_countALL$habitat =="Forest" & corr_countALL$link_type=="BD-EF", corr_countALL$n/90,# OR replace 90 by (nrow(subset(vertex_info_Forest, Node_type=="BD")))*(nrow(subset(vertex_info_Forest, Node_type=="EF"))) #i.e. number of BD nodes * number of EF nodes
                                   ifelse(corr_countALL$habitat =="Forest" & corr_countALL$link_type=="BD-ES", corr_countALL$n/135, 
                                          ifelse(corr_countALL$habitat =="Forest" & corr_countALL$link_type=="EF-ES", corr_countALL$n/54,
                                                 ifelse(corr_countALL$habitat =="Grass" & corr_countALL$link_type=="BD-EF", corr_countALL$n/120,
                                                        ifelse(corr_countALL$habitat =="Grass" & corr_countALL$link_type=="BD-ES", corr_countALL$n/105,
                                                               ifelse(corr_countALL$habitat =="Grass" & corr_countALL$link_type=="EF-ES", corr_countALL$n/56,"ERROR"
                                                               ))))))

corr_countALL$n_corrected<-as.numeric(corr_countALL$n_corrected)

      # 4.2.2 Figure (Fig. S3)----
library(grid)
library(ggplot2)
library(gridExtra)
library(mgcViz)

#rename habitats for the figure
head(corr_countALL)
corr_countALL$habitat <-with (corr_countALL, ifelse (habitat  =="Grass","Grassland", "Forest"))

plot_n<-ggplot(subset(corr_countALL, corr_type=="POS"), aes(x=mean_LUI, y=n, colour = factor(link_type))) + # only positive ones are needed (Negative are the inverse)
  geom_point( size = 1) + 
  theme_bw() + 
  labs(x = "")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none" ,
        strip.text.y = element_blank(),
        axis.text.x.bottom= element_blank(),
        axis.ticks.x= element_blank())+
  # xlab("Land-use intensity") + 
  ylab("Total number of positive correlations") + 
  facet_grid(corr_type~habitat, scales="free_x") +
  geom_smooth(method = 'lm', formula = y ~ x) # method = 'loess',span=1

plot_ncorrected<-ggplot(subset(corr_countALL, corr_type=="POS"), aes(x=mean_LUI, y=n_corrected, colour = factor(link_type))) + # only positive ones are needed (Negative are the inverse)
  geom_point( size = 1) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.position=c(0.9,0.22),
        strip.text.y = element_blank(),
        strip.text.x = element_blank()) +
  xlab("Land-use intensity") + 
  ylab("Relative number of positive correlations") + 
  facet_grid(corr_type~habitat, scales="free_x") +
  geom_smooth(method = 'lm', formula = y ~ x) # method = 'loess',span=1

windows () 
gridPrint(grobs = list(plot_n,plot_ncorrected))


tiff("...\\FigS3.tiff", width = 7, height = 10, units = 'in', res = 300)  
dev.off()

# check significance

summary(lm_GRASS_POS_BD.ES)
lm_GRASS_POS<-lm(n_corrected ~ mean_LUI,data=subset(corr_countALL,corr_type=="POS" & habitat=="Grass")) #***
lm_GRASS_POS_BD.EF<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="BD-EF" & habitat=="Grass"))#***
lm_GRASS_POS_BD.ES<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="BD-ES"& habitat=="Grass")) # NS Y 0% 
gam_GRASS_POS_BD.ES<-gam(n_corrected ~ s(mean_LUI), data=subset (corr_countALL, corr_type=="POS" & link_type=="BD-ES"& habitat=="Grass")) #***
plot(gam_GRASS_POS_BD.ES)
gam.check(gam_GRASS_POS_BD.ES)
lm_GRASS_POS_EF.ES<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="EF-ES"& habitat=="Grass")) #***

lm_Forest_POS<-lm(n_corrected ~ mean_LUI,data=subset(corr_countALL,corr_type=="POS"& habitat=="Forest")) #***
lm_Forest_POS_BD.EF<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="BD-EF"& habitat=="Forest")) # *
lm_Forest_POS_BD.ES<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="BD-ES"& habitat=="Forest")) #***
lm_Forest_POS_EF.ES<-lm(n_corrected ~ mean_LUI, data=subset (corr_countALL, corr_type=="POS" & link_type=="EF-ES"& habitat=="Forest")) #***

    # 4.3 Comparing metrics by link types (BD-EF, BD-ES, EF-ES) only for partial positive correlations####

      # 4.3.1 Compute network metrics====
# Grassland 
library(igraph)
library(data.table)
library(vegan)
source("pcor.R") # PARTIAL CORRELATIONS that allows for NAs: http://www.yilab.gatech.edu/pcor.html 
# requires function "Network_metrics" (see above)

metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=subset (vertex_info_grass,Node_type =="BD" | Node_type =="EF"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_GrassPOS_BDEF <- metrics$Select.Net.metrics_final
Node_metrics_GrassPOS_BDEF <- metrics$Select.Node.metrics_final
Network_matrix_GrassPOS_BDEF <- metrics$Select.Net.matrix_final
Net.GrassPOS_BDEF <- metrics$Select.Net
edges.list_GrassPOS_BDEF <- metrics$edges.list

metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=subset (vertex_info_grass,Node_type =="BD" | Node_type =="ES"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_GrassPOS_BDES <- metrics$Select.Net.metrics_final
Node_metrics_GrassPOS_BDES <- metrics$Select.Node.metrics_final
Network_matrix_GrassPOS_BDES <- metrics$Select.Net.matrix_final
Net.GrassPOS_BDES <- metrics$Select.Net
edges.list_GrassPOS_BDES <- metrics$edges.list

metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=subset (vertex_info_grass,Node_type =="EF" | Node_type =="ES"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_GrassPOS_EFES <- metrics$Select.Net.metrics_final
Node_metrics_GrassPOS_EFES <- metrics$Select.Node.metrics_final
Network_matrix_GrassPOS_EFES <- metrics$Select.Net.matrix_final
Net.GrassPOS_EFES <- metrics$Select.Net
edges.list_GrassPOS_EFES <- metrics$edges.list

  #Forest 
metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=subset (vertex_info_Forest,Node_type =="BD" | Node_type =="EF"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_ForestPOS_BDEF <- metrics$Select.Net.metrics_final
Node_metrics_ForestPOS_BDEF <- metrics$Select.Node.metrics_final
Network_matrix_ForestPOS_BDEF <- metrics$Select.Net.matrix_final
Net.ForestPOS_BDEF <- metrics$Select.Net
edges.list_ForestPOS_BDEF <- metrics$edges.list

metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=subset (vertex_info_Forest,Node_type =="BD" | Node_type =="ES"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_ForestPOS_BDES <- metrics$Select.Net.metrics_final
Node_metrics_ForestPOS_BDES <- metrics$Select.Node.metrics_final
Network_matrix_ForestPOS_BDES <- metrics$Select.Net.matrix_final
Net.ForestPOS_BDES <- metrics$Select.Net
edges.list_ForestPOS_BDES <- metrics$edges.list

metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=subset (vertex_info_Forest,Node_type =="EF" | Node_type =="ES"), 
                            estimate="POS",correlation="PARTIAL")
Network_metrics_ForestPOS_EFES <- metrics$Select.Net.metrics_final
Node_metrics_ForestPOS_EFES <- metrics$Select.Node.metrics_final
Network_matrix_ForestPOS_EFES <- metrics$Select.Net.matrix_final
Net.ForestPOS_EFES <- metrics$Select.Net
edges.list_ForestPOS_EFES <- metrics$edges.list  

save.image()

      # 4.3.2 Models ----

#add LUI
Select.Net.metrics_Grass_BDEF<- merge(Network_metrics_GrassPOS_BDEF,LUI_table,by="window_block")
Select.Net.metrics_Grass_BDES<- merge(Network_metrics_GrassPOS_BDES,LUI_table,by="window_block")
Select.Net.metrics_Grass_EFES<- merge(Network_metrics_GrassPOS_EFES,LUI_table,by="window_block")

Select.Net.metrics_Forest_BDEF<- merge(Network_metrics_ForestPOS_BDEF,LUI_table_Forest,by="window_block")
Select.Net.metrics_Forest_BDES<- merge(Network_metrics_ForestPOS_BDES,LUI_table_Forest,by="window_block")
Select.Net.metrics_Forest_EFES<- merge(Network_metrics_ForestPOS_EFES,LUI_table_Forest,by="window_block")

# GAMs 
library(mgcv)
wdGg_BDEF <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Grass_BDEF) 
evGg_BDEF <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_BDEF) #
moGg_BDEF <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_BDEF)

wdGg_BDES <- gam(connectance ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_BDES)  
evGg_BDES <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_BDES)
moGg_BDES <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_BDES)

wdGg_EFES <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Grass_EFES) #k=5
evGg_EFES <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_EFES)
moGg_EFES <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_EFES) #k=5

wdFg_BDEF <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_BDEF) 
evFg_BDEF <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_BDEF)  
moFg_BDEF <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Forest_BDEF) #n.s. #k=4

wdFg_BDES <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_BDES) 
evFg_BDES <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_BDES)  
moFg_BDES <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_BDES)

wdFg_EFES <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_EFES) ##  
evFg_EFES <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_EFES)  
moFg_EFES <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_EFES)

# check models e.g.
summary(moFg_EFES)
plot.gam(wdFg_EFES)
par(mfrow=c(2,2))
gam.check(evGg_BDEF)
AIC(evGg_BDEF1,evGg_BDEF2)

      # 4.3.3 Table S3 ----

# requires functions from "Summary GAMs tables for SI"
models <- list(wdGg_BDEF,moGg_BDEF,evGg_BDEF,wdGg_BDES,moGg_BDES,evGg_BDES,wdGg_EFES,moGg_EFES,evGg_EFES,
               wdFg_BDEF,moFg_BDEF,evFg_BDEF,wdFg_BDES,moFg_BDES,evFg_BDES,wdFg_EFES,moFg_EFES,evFg_EFES) 
Gams_link_type <- data.frame(Habitat = rep(c("Grassland", "Forest"), each=3),
                             Link_type = rep(c("BD-EF", "BD-ES", "EF-ES"),each=6),
                             Metric =rep(c("Connectance","Modularity", "Evenness"),6) ,
                          k =  unlist (lapply(models, k_val)),
                          EDF =   unlist (lapply(models, EDF_val)),
                          F =     unlist (lapply(models, F_val)),
                          p.value=unlist (lapply(models, p.value_val)),
                          Adj.R2 =unlist (lapply(models, Adj.R2_val)),
                          DE =    unlist (lapply(models, DE_val)),
                          N=      unlist (lapply(models, N_val)))
# save table for SI
write.table(Gams_link_type,"...\\Gams_link_type.txt",sep="\t", dec=".") 

      # 4.3.4 Figure (Fig.S4) ----

library (gratia)
library(ggplot2)
library(mgcViz)

#prepare data
wdGg_BDEF_conf <-confint(wdGg_BDEF, parm = "mean_LUI", type = "confidence")
wdGg_BDES_conf <-confint(wdGg_BDES, parm = "mean_LUI", type = "confidence")
wdGg_EFES_conf <-confint(wdGg_EFES, parm = "mean_LUI", type = "confidence")

wdFg_BDEF_conf <-confint(wdFg_BDEF, parm = "mean_LUI", type = "confidence")
wdFg_BDES_conf <-confint(wdFg_BDES, parm = "mean_LUI", type = "confidence")
wdFg_EFES_conf <-confint(wdFg_EFES, parm = "mean_LUI", type = "confidence")

moGg_BDEF_conf <-confint(moGg_BDEF, parm = "mean_LUI", type = "confidence")
moGg_BDES_conf <-confint(moGg_BDES, parm = "mean_LUI", type = "confidence")
moGg_EFES_conf <-confint(moGg_EFES, parm = "mean_LUI", type = "confidence")

moFg_BDEF_conf <-confint(moFg_BDEF, parm = "mean_LUI", type = "confidence")
moFg_BDES_conf <-confint(moFg_BDES, parm = "mean_LUI", type = "confidence")
moFg_EFES_conf <-confint(moFg_EFES, parm = "mean_LUI", type = "confidence")

evGg_BDEF_conf <-confint(evGg_BDEF, parm = "mean_LUI", type = "confidence")
evGg_BDES_conf <-confint(evGg_BDES, parm = "mean_LUI", type = "confidence")
evGg_EFES_conf <-confint(evGg_EFES, parm = "mean_LUI", type = "confidence")

evFg_BDEF_conf <-confint(evFg_BDEF, parm = "mean_LUI", type = "confidence")
evFg_BDES_conf <-confint(evFg_BDES, parm = "mean_LUI", type = "confidence")
evFg_EFES_conf <-confint(evFg_EFES, parm = "mean_LUI", type = "confidence")

# plots (Fig. S4)

plot1 <-ggplot( ) +
  labs(y =" ",x=" ") + 
  theme_classic()  +
  geom_line(data=wdGg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = wdGg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=wdGg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = wdGg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=wdGg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = wdGg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  ggtitle("Grasslands")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot2 <-ggplot( ) +
  labs(y =" ",x=" ") +  
  theme_classic()  +
  geom_line(data=moGg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = moGg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=moGg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = moGg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=moGg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = moGg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))

plot3 <-ggplot( ) +
  labs(y =" ",x="Land-use intensity") + 
  theme_classic()  +
  geom_line(data=evGg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = evGg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=evGg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = evGg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=evGg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = evGg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.025,0.025), breaks = seq(-0.02,0.02, by = 0.02))

plot4 <-ggplot( ) +
  labs(y ="Connectance",x=" ") +
  theme_classic()  +
  geom_line(data=wdFg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = wdFg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=wdFg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = wdFg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=wdFg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = wdFg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  ggtitle("Forest")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.025,0.025), breaks = seq(-0.02,0.02, by = 0.01)) 

plot5 <-ggplot( ) +
  labs(y ="Modularity",x=" ") + 
  theme_classic()  +
  geom_line(data=moFg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = moFg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=moFg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = moFg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=moFg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = moFg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.15,0.15), breaks = seq(-0.1,0.1, by = 0.1))

plot6 <-ggplot( ) +
  labs(y ="Evenness",x="Land-use intensity") + 
  theme_classic()  +
  geom_line(data=evFg_BDEF_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkolivegreen3") +
  geom_ribbon(data = evFg_BDEF_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkolivegreen3",alpha = .15)+ 
  geom_line(data=evFg_BDES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorange2") +
  geom_ribbon(data = evFg_BDES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorange2",alpha = .15)+
  geom_line(data=evFg_EFES_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkorchid") +
  geom_ribbon(data = evFg_EFES_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkorchid",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.025,0.025), breaks = seq(-0.02,0.02, by = 0.02))

windows ()
gridPrint(grobs = list(plot4,plot1,plot5,plot2,plot6,plot3), ncol = 2) # to have forest on the left and grasslands on the right columns



    # 4.4 COMPARING ABOVE-BELOWGROUND #####
      # 4.4.1 Prepare data ----
        # 4.4.1.1 For ABOVEGROUND SPECIES/EF/ES----
          ## create "vertex_info_grass_AG" ONLY with used variables!! ----
data_Grass_select_AG<-subset(data_Grass_select, select = c(Lichen, Moss, Plant,Plant.pathog,
                                                           AG.decomp,AG.herb,Snd.cons, 
                                                           Dung.decomp,Herbivory.control,
                                                           Biomass,Forage.quality,Pot.bird_watching,Butterfl.abund))

vertex_info_grass_AG <-data.frame(Nodes = names(data_Grass_select_AG), Node_type = c(rep ("BD",7), rep ("EF",1), rep ("ES",5)))

          ## create "vertex_info_Forest_AG" ONLY with used variables!! ----
data_Forest_select_AG<-subset(data_Forest_res_order, select = c(Lichen, Moss, Plant,AG.decomp,AG.herb,AG.omniv,Pollinator,Snd.cons,Dung.decomp,
                                                                Pest.control,Temp.reg,Trees.C.stock,Timber,Edible.fungi,Edible.plants,Cultural.plants,Pot.bird_watching))
#drop columns with more than 30 NAs
data_Forest_select_AG<-subset(data_Forest_select_AG, select = -c(Pollinator)) 
vertex_info_Forest_AG <-data.frame(Nodes = names(data_Forest_select_AG), Node_type = c(rep ("BD",7), "EF", rep ("ES",8))) 

       # 4.4.1.2 For BELOWGROUND SPECIES/EF/ES (basically take the rest from AG)----
          ## create "vertex_info_grass_BG" ONLY with used variables!! ----
data_Grass_select_BG<-subset(data_Grass_select, select = -c(Lichen, Moss, Plant,Plant.pathog,AG.decomp,AG.herb,Snd.cons,
                                                            Dung.decomp,Herbivory.control,
                                                            Biomass,Forage.quality,Pot.bird_watching,Butterfl.abund))

vertex_info_grass_BG <-data.frame(Nodes = names(data_Grass_select_BG), Node_type = c(rep ("BD",8), rep ("EF",7), rep ("ES",2))) 

          ## create "vertex_info_Forest_AG" ONLY with used variables!! ----
data_Forest_select_BG<-subset(data_Forest_res_order, select = -c(EP_Plotid,ForMI,Plot0,Exploratory,
                                                                 Lichen, Moss, Plant,AG.decomp,AG.herb,AG.omniv,Pollinator,Snd.cons,
                                                                 Dung.decomp,Pest.control,Temp.reg,Trees.C.stock,
                                                                 Timber,Edible.fungi,Edible.plants,Cultural.plants,Pot.bird_watching))
#drop columns with more than 30 NAs
data_Forest_select_BG<-subset(data_Forest_select_BG, select = -c(Mycorrhiza))
vertex_info_Forest_BG <-data.frame(Nodes = names(data_Forest_select_BG), Node_type = c(rep ("BD",8), rep ("EF",5), rep ("ES",1)))


       # 4.4.2 Networks ----
  library(igraph)
  library(data.table)
  library(vegan)
  source("pcor.R") # PARTIAL CORRELATIONS that allows for NAs: http://www.yilab.gatech.edu/pcor.html 
  # requires function "Network_metrics" (see above)
  
  Network_metrics_GrassAG <- Network_metrics (np=150, mws=60,data=data_Grass_select_AG,Vertices=vertex_info_grass_AG, estimate="POS",correlation="PARTIAL")
  Select.Net.metrics_Grass_AG <- Network_metrics_GrassAG$Select.Net.metrics_final
  Select.Node.metrics_Grass_AG <- Network_metrics_GrassAG$Select.Node.metrics_final
  Select.Grass.Net.matrix_AG <- Network_metrics_GrassAG$Select.Net.matrix_final
  Select.Grass_AG.NET <- Network_metrics_GrassAG$Select.Net
  Select.Grass_AG.edges <- Network_metrics_GrassAG$edges.list
  
  Network_metrics_ForestAG <- Network_metrics (np=150, mws=50,data=data_Forest_select_AG,Vertices=vertex_info_Forest_AG, estimate="POS",correlation="PARTIAL")
  Select.Net.metrics_Forest_AG <- Network_metrics_ForestAG$Select.Net.metrics_final
  Select.Node.metrics_Forest_AG <- Network_metrics_ForestAG$Select.Node.metrics_final
  Select.Forest.Net.matrix_AG <- Network_metrics_ForestAG$Select.Net.matrix_final
  Select.Forest_AG.NET <- Network_metrics_ForestAG$Select.Net
  Select.Forest_AG.edges <- Network_metrics_ForestAG$edges.list
  
  Network_metrics_GrassBG <- Network_metrics (np=150, mws=60,data=data_Grass_select_BG,Vertices=vertex_info_grass_BG, estimate="POS",correlation="PARTIAL")
  Select.Net.metrics_Grass_BG <- Network_metrics_GrassBG$Select.Net.metrics_final
  Select.Node.metrics_Grass_BG <- Network_metrics_GrassBG$Select.Node.metrics_final
  Select.Grass.Net.matrix_BG <- Network_metrics_GrassBG$Select.Net.matrix_final
  Select.Grass_BG.NET <- Network_metrics_GrassBG$Select.Net
  Select.Grass_BG.edges <- Network_metrics_GrassBG$edges.list
  
  Network_metrics_ForestBG <- Network_metrics (np=150, mws=50,data=data_Forest_select_BG,Vertices=vertex_info_Forest_BG, estimate="POS",correlation="PARTIAL")
  Select.Net.metrics_Forest_BG <- Network_metrics_ForestBG$Select.Net.metrics_final
  Select.Node.metrics_Forest_BG <- Network_metrics_ForestBG$Select.Node.metrics_final
  Select.Forest.Net.matrix_BG <- Network_metrics_ForestBG$Select.Net.matrix_final
  Select.Forest_BG.NET <- Network_metrics_ForestBG$Select.Net
  Select.Forest_BG.edges <- Network_metrics_ForestBG$edges.list
  
  save.image()
  
       # 4.4.3 Models ----
  
  #add LUI
  Select.Net.metrics_Grass_AG<- merge(Select.Net.metrics_Grass_AG,LUI_table,by="window_block")
  Select.Net.metrics_Forest_AG<- merge(Select.Net.metrics_Forest_AG,LUI_table_Forest,by="window_block")
  Select.Net.metrics_Grass_BG<- merge(Select.Net.metrics_Grass_BG,LUI_table,by="window_block")
  Select.Net.metrics_Forest_BG<- merge(Select.Net.metrics_Forest_BG,LUI_table_Forest,by="window_block")
  
  # see how to do model selection in the Main models section (2.3)
  library(mgcv)
  
  wdG_AG <- gam(connectance ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_AG)
  wdG_BG <- gam(connectance ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_BG)
  wdF_AG <- gam(connectance ~ s(mean_LUI,k=7),data=Select.Net.metrics_Forest_AG)
  wdF_BG <- gam(connectance ~ s(mean_LUI,k=7),data=Select.Net.metrics_Forest_BG)
  
  evG_AG <- gam(evenness ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_AG) #
  evG_BG <- gam(evenness ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_BG)
  evF_AG <- gam(evenness ~ s(mean_LUI,k=5),data=Select.Net.metrics_Forest_AG)
  evF_BG <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_BG)
  
  moG_AG <- gam(modularity ~ s(mean_LUI,k=3),data=Select.Net.metrics_Grass_AG) 
  moG_BG <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_BG)#
  moF_AG <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_AG)
  moF_BG <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_BG)
  
  # check models e.g.
  summary(evG_AG1)
  plot.gam(evG_AG1)
  par(mfrow=c(2,2))
  gam.check(evG_AG)
  AIC(evG_AG1,evG_AG)
  
       # 4.4.4 Table S4 ----
models <- list(wdG_AG,moG_AG,evG_AG,wdG_BG,moG_BG,evG_BG,wdF_AG,moF_AG,evF_AG,wdF_BG,moF_BG,evF_BG) 
Gams_AGBG <- data.frame(Habitat = rep(c("Grassland", "Forest"), each=6),
                          Compartment = rep(c("Aboveground", "Belowground"), each=3),
                          Metric =rep(c("Connectance","Modularity", "Evenness"),4) ,
                          k =  unlist (lapply(models, k_val)),
                          EDF =   unlist (lapply(models, EDF_val)),
                          F =     unlist (lapply(models, F_val)),
                          p.value=unlist (lapply(models, p.value_val)),
                          Adj.R2 =unlist (lapply(models, Adj.R2_val)),
                          DE =    unlist (lapply(models, DE_val)),
                          N=      unlist (lapply(models, N_val)))
# save table for SI
write.table(Gams_AGBG,"...\\Gams_AGBG.txt",sep="\t", dec=".") 


       # 4.4.5 Figure (Fig. S5)----

library (gratia)
library(ggplot2)
library(mgcViz)

summary(evF_BG)

wdG_AG_conf <-confint(wdG_AG, parm = "mean_LUI", type = "confidence")
wdG_BG_conf <-confint(wdG_BG, parm = "mean_LUI", type = "confidence")

moG_AG_conf <-confint(moG_AG, parm = "mean_LUI", type = "confidence") 
moG_BG_conf <-confint(moG_BG, parm = "mean_LUI", type = "confidence") 

evG_AG_conf <-confint(evG_AG, parm = "mean_LUI", type = "confidence") 
evG_BG_conf <-confint(evG_BG, parm = "mean_LUI", type = "confidence")

wdF_AG_conf <-confint(wdF_AG, parm = "mean_LUI", type = "confidence")
wdF_BG_conf <-confint(wdF_BG, parm = "mean_LUI", type = "confidence")

moF_AG_conf <-confint(moF_AG, parm = "mean_LUI", type = "confidence")
moF_BG_conf <-confint(moF_BG, parm = "mean_LUI", type = "confidence") 

evF_AG_conf <-confint(evF_AG, parm = "mean_LUI", type = "confidence")
evF_BG_conf <-confint(evF_BG, parm = "mean_LUI", type = "confidence") #*

plot1 <-ggplot( ) +
  labs(y ="Connectance",x="Land-use intensity") + # change y label and see limits:
  geom_line(data=wdG_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = wdG_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdG_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = wdG_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  ggtitle("Grassland")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.02,0.02), breaks = seq(-0.02,0.02, by = 0.02))

plot2 <-ggplot( ) +
  labs(y ="Modularity",x="Land-use intensity") + # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=moG_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = moG_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+
  theme_classic()  +
  geom_line(data=moG_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = moG_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.2,0.2), breaks = seq(-0.2,0.2, by = 0.2))

plot3 <-ggplot( ) +
  labs(y ="Evenness",x="Land-use intensity") +  # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=evG_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = evG_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evG_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = evG_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, by = 0.04))

plot4 <-ggplot( ) +
  labs(y ="Connectance",x="Land-use intensity") + # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=wdF_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = wdF_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdF_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = wdF_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  ggtitle("Forest")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.02,0.02), breaks = seq(-0.02,0.02, by = 0.02))

plot5 <-ggplot( ) +
  labs(y ="Modularity",x="Land-use intensity") +  # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=moF_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = moF_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+
  theme_classic()  +
  geom_line(data=moF_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = moF_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.2,0.2), breaks = seq(-0.2,0.2, by = 0.2))

plot6 <-ggplot( ) +
  labs(y ="Evenness",x="Land-use intensity") +# change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=evF_AG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgreen") +
  geom_ribbon(data = evF_AG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgreen",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evF_BG_conf, aes(x = mean_LUI,y = est), size = 1,colour = "brown") +
  geom_ribbon(data = evF_BG_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="brown",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.04,0.04), breaks = seq(-0.04,0.04, by = 0.04))


windows ()
gridPrint(grobs = list(plot4,plot1,plot5,plot2,plot6,plot3), ncol = 2) # to have forest on the left and grasslands on the right columns








    # 4.5 RANDOm networks #####
      # 4.5.1 Compute the networks ====

# create 100 randomizations of the data set (with each random data I run the nn networks)

library(igraph)
library(data.table)
library(vegan)
source("pcor.R") # PARTIAL CORRELATIONS that allows for NAs: http://www.yilab.gatech.edu/pcor.html 
# requires function "Network_metrics" (see above)

       # Grasslands ----
Network_metrics_Grass_random <- list()

for (w in 1:100){ 
  # set seed to fix random sample per network. If computation stops, add the number of loops to continue with the sequence and avoid repetitions of the same sample 
  set.seed(123+w) #set.seed(123+w+75) #get stuck in the loop 75, so I jump one and continue from there
  # Randomize the order of the data frame
  random_data <- data_Grass_select[sample(1:nrow(data_Grass_select)), ]
  #adjust function parameters as needed
  Network_metrics_random <- Network_metrics (np=150, mws=60,data=random_data,Vertices=vertex_info_grass, estimate="POS",correlation="PARTIAL")
  Network_metrics_random <- Network_metrics_random$Select.Net.metrics_final

  Network_metrics_Grass_random[length(Network_metrics_Grass_random)+1] <- list(Network_metrics_random)

  print(paste ("Random Network",w))
  }

nn= 91

## part 1
Random.Net.metrics_Grass_p1 <- Reduce(rbind, Network_metrics_Grass_random)
Random.Net.metrics_Grass_p1 <- as.data.frame(Random.Net.metrics_Grass_p1)
rownames (Random.Net.metrics_Grass_p1) <-NULL
#add additional variables
Random.Net.metrics_GrassLUI_p1<-Random.Net.metrics_Grass_p1
          #--> up to 50 of loop 74. delete the last and resume after 73 --------------------------
Random.Net.metrics_GrassLUI_p1$mean_LUI<- rep(LUI_table[,"mean_LUI"],74) #1º chunk of data has 6734 (i.e. 74)
Random.Net.metrics_GrassLUI_p1$window_block<- rep(LUI_table[,"window_block"],74)
Random.Net.metrics_GrassLUI_p1$sim_run<- rep(c(1:74),each=nn)
Random.Net.metrics_GrassLUI_p1$habitat <- "Grass"
dim(Random.Net.metrics_GrassLUI_p1)

save.image("...\\random_p1.RData")

## part 2
Random.Net.metrics_Grass_p2 <- Reduce(rbind, Network_metrics_Grass_random)
Random.Net.metrics_Grass_p2 <- as.data.frame(Random.Net.metrics_Grass_p2)
rownames (Random.Net.metrics_Grass_p2) <-NULL
#add additional variables
Random.Net.metrics_GrassLUI_p2<-Random.Net.metrics_Grass_p2
          #--> up to 38 of loop 91.------------- 
Random.Net.metrics_GrassLUI_p2$mean_LUI<- rep(LUI_table[,"mean_LUI"],93) # 2º part has 8463 (i.e. 93)
Random.Net.metrics_GrassLUI_p2$window_block<- rep(LUI_table[,"window_block"],93)
Random.Net.metrics_GrassLUI_p2$sim_run<- rep(c(1:93),each=nn)
Random.Net.metrics_GrassLUI_p2$habitat <- "Grass"
summary(Random.Net.metrics_GrassLUI_p2)

save.image("...\\random_p2.RData")

#load RData from the first chunk of loops

#calculate how many loops I need from each chunk:
# discard the last one from the first, ie. 1-73
100-73 #=27
Random.Net.metrics_GrassLUI_total<- rbind (subset(Random.Net.metrics_GrassLUI_p1,sim_run<74),subset(Random.Net.metrics_GrassLUI_p2,sim_run<28))
dim(Random.Net.metrics_GrassLUI_total) #9100

save.image("...random_all.RData")


       #  Forest ----
Network_metrics_Forest_random <- list()

for (w in 1:100){ 
  # set seed to fix random sample per network. If computation stops, add the number of loops to continue with the sequence and avoid repetitions of the same sample 
  set.seed(123+w) 
  # Randomize the order of the data frame
  random_data <- data_Forest_select[sample(1:nrow(data_Forest_select)), ]
  #adjust function parameters as needed
  Network_metrics_random <- Network_metrics (np=150, mws=50,data=random_data,Vertices=vertex_info_Forest, estimate="POS",correlation="PARTIAL")
  Network_metrics_random <- Network_metrics_random$Select.Net.metrics_final
  
  Network_metrics_Forest_random[length(Network_metrics_Forest_random)+1] <- list(Network_metrics_random)
  
  print(paste ("Random Network",w))
}

Random.Net.metrics_Forest <- Reduce(rbind, Network_metrics_Forest_random)
Random.Net.metrics_Forest <- as.data.frame(Random.Net.metrics_Forest)
rownames (Random.Net.metrics_Forest_pos) <-NULL
#add additional variables
Random.Net.metrics_ForestLUI<-Random.Net.metrics_Forest
Random.Net.metrics_ForestLUI$mean_LUI<- rep(LUI_table_Forest[,"mean_LUI"],100)
Random.Net.metrics_ForestLUI$window_block<- rep(LUI_table_Forest[,"window_block"],100)
Random.Net.metrics_ForestLUI$sim_run<- rep(c(1:100),each=nn)
Random.Net.metrics_ForestLUI$habitat <- "Forest"


save.image("...\\random_10062020.RData")


      # 4.5.2 Plot the figure (Fig. S6)----
Random.Net.metrics <- rbind(Random.Net.metrics_GrassLUI_total,Random.Net.metrics_ForestLUI)

library(reshape)
Random.Net.metrics_long <- melt(Random.Net.metrics, id=c("habitat","mean_LUI","sim_run")) 
Random.Net.metrics_long$sim_run<-as.factor(Random.Net.metrics_long$sim_run)

head(Select.Net.metrics_Forest_pos)
Select.Net.metrics_Grass_pos$habitat <- "Grass"
Select.Net.metrics_Forest_pos$habitat <- "Forest"
Select.Net.metrics_pos<- rbind(Select.Net.metrics_Grass_pos,Select.Net.metrics_Forest_pos)
Select.Net.metrics_pos_long <- melt(Select.Net.metrics_pos, id=c("habitat","mean_LUI")) 

library(ggplot2)
## mean and CI
p= ggplot(subset(Random.Net.metrics_long,variable==c("connectance","modularity","evenness")), aes(x=mean_LUI, y=value)) + 
  geom_smooth(color="limegreen",method = 'loess') +  # geom_line(aes(group = sim_run), size = 0.5,alpha=0.1,color="red") 
  theme_minimal() +
  xlab("Land-use intensity") + 
  ylab("") +
  facet_wrap(variable~habitat, scales="free", ncol=2) 
p2 = p + geom_smooth(data= subset(Select.Net.metrics_pos_long,variable==c("connectance","modularity","evenness")),aes(x=mean_LUI, y=value)) #,method = 'gam', formula = y ~ s(x, bs = "tp")# "cs"
windows()
print(p2) 


    # 4.6 Comparing significant vs. non-significant positive correlations (Fig. S7) ----
#--------------------------------------------------------------------#
# only possible using "raw" = "non-partial correlations" 
# (partial correlations are already very conservative and subsetting the significant ones would result in very few correlations left)

library(mgcv)

       # 4.6.1 Compute the networks ====

         # ALL (raw non-partial positive Spearman correlations) ----

library(igraph)
library(data.table)
library(vegan)
# requires function "Network_metrics"
# requires function "Raw correlations"
library("Hmisc")

metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=vertex_info_grass, estimate="POS",correlation="RAW")
Select.Net.metrics_Grass_RAWpos <- metrics$Select.Net.metrics_final

metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=vertex_info_Forest, estimate="POS",correlation="RAW")
Select.Net.metrics_Forest_RAWpos <- metrics$Select.Net.metrics_final

         # SIGNIFICANT raw non-partial positive Spearman correlations ----
metrics <- Network_metrics (np=150, mws=60,data=data_Grass_select,Vertices=vertex_info_grass, estimate="POS",correlation="SIGN")
Select.Net.metrics_Grass_SIGNpos <- metrics$Select.Net.metrics_final

metrics <- Network_metrics (np=150, mws=50,data=data_Forest_select,Vertices=vertex_info_Forest, estimate="POS",correlation="SIGN")
Select.Net.metrics_Forest_SIGNpos <- metrics$Select.Net.metrics_final

       # 4.6.2 MODELS ====

#add LUI
Select.Net.metrics_Grass_RAWpos<- merge(Select.Net.metrics_Grass_RAWpos,LUI_table,by="window_block")
Select.Net.metrics_Forest_RAWpos<- merge(Select.Net.metrics_Forest_RAWpos,LUI_table_Forest,by="window_block")
Select.Net.metrics_Grass_SIGNpos<- merge(Select.Net.metrics_Grass_SIGNpos,LUI_table,by="window_block")
Select.Net.metrics_Forest_SIGNpos<- merge(Select.Net.metrics_Forest_SIGNpos,LUI_table_Forest,by="window_block")

# For ALL correlations (all signif+non-sign)
summary(moFg_SIGNpos)
wdGg_RAWpos <- gam(connectance ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_RAWpos) 
evGg_RAWpos <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_RAWpos)
moGg_RAWpos <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_RAWpos)

wdFg_RAWpos <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_RAWpos) 
evFg_RAWpos <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_RAWpos)  
moFg_RAWpos <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_RAWpos) #n.s.

# For SIGNIFICANT correlations only 
wdGg_SIGNpos <- gam(connectance ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_SIGNpos) 
evGg_SIGNpos <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Grass_SIGNpos)
moGg_SIGNpos <- gam(modularity ~ s(mean_LUI,k=5),data=Select.Net.metrics_Grass_SIGNpos)

wdFg_SIGNpos <- gam(connectance ~ s(mean_LUI,k=6),data=Select.Net.metrics_Forest_SIGNpos) 
evFg_SIGNpos <- gam(evenness ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_SIGNpos)  
moFg_SIGNpos <- gam(modularity ~ s(mean_LUI,k=4),data=Select.Net.metrics_Forest_SIGNpos) 

# for checking the models use e.g.
AIC(evGg_SIGNpos1,evGg_SIGNpos2)
summary(evGg_SIGNpos1)
plot.gam(evGg_SIGNpos1)
par(mfrow=c(2,2))
gam.check(evGg_SIGNpos1)



       # 4.6.3 Table to SI (Table S7)====

models <- list(wdGg_RAWpos,moGg_RAWpos,evGg_RAWpos,wdFg_RAWpos,moFg_RAWpos,evFg_RAWpos,wdGg_SIGNpos,moGg_SIGNpos,evGg_SIGNpos,wdFg_SIGNpos,moFg_SIGNpos,evFg_SIGNpos) 
Gams_Raw <- data.frame(Habitat = rep(c("Grassland", "Forest"), each=3),
                       Corr_sign = rep(c("Raw", "Sign"), each=6),
                       Metric =rep(c("Connectance","Modularity", "Evenness"),4) ,
                        k =  unlist (lapply(models, k_val)),
                        EDF =   unlist (lapply(models, EDF_val)),
                        F =     unlist (lapply(models, F_val)),
                        p.value=unlist (lapply(models, p.value_val)),
                        Adj.R2 =unlist (lapply(models, Adj.R2_val)),
                        DE =    unlist (lapply(models, DE_val)),
                        N=      unlist (lapply(models, N_val)))
# save table for SI
write.table(Gams_Raw,"...\\Gams_Raw.txt",sep="\t", dec=".") 


       # 4.6.4 FIGURE S7 ====
library (gratia)
library(ggplot2)
library(mgcViz)

summary(evFg_SIGNpos)
  
wdGg_RAWpos_conf <-confint(wdGg_RAWpos, parm = "mean_LUI", type = "confidence")
wdGg_SIGNpos_conf <-confint(wdGg_SIGNpos, parm = "mean_LUI", type = "confidence")

moGg_RAWpos_conf <-confint(moGg_RAWpos, parm = "mean_LUI", type = "confidence")
moGg_SIGNpos_conf <-confint(moGg_SIGNpos, parm = "mean_LUI", type = "confidence") 

evGg_RAWpos_conf <-confint(evGg_RAWpos, parm = "mean_LUI", type = "confidence")
evGg_SIGNpos_conf <-confint(evGg_SIGNpos, parm = "mean_LUI", type = "confidence")

wdFg_RAWpos_conf <-confint(wdFg_RAWpos, parm = "mean_LUI", type = "confidence")
wdFg_SIGNpos_conf <-confint(wdFg_SIGNpos, parm = "mean_LUI", type = "confidence")

moFg_RAWpos_conf <-confint(moFg_RAWpos, parm = "mean_LUI", type = "confidence")   # non significant
moFg_SIGNpos_conf <-confint(moFg_SIGNpos, parm = "mean_LUI", type = "confidence")

evFg_RAWpos_conf <-confint(evFg_RAWpos, parm = "mean_LUI", type = "confidence")
evFg_SIGNpos_conf <-confint(evFg_SIGNpos, parm = "mean_LUI", type = "confidence")

# fig
plot1 <-ggplot( ) +
  labs(y =" ",x=" ") +  # change y label and see limits:
  geom_line(data=wdGg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = wdGg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdGg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = wdGg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  ggtitle("Grassland")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot2 <-ggplot( ) +
  labs(y =" ",x=" ") +# change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=moGg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = moGg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+
  theme_classic()  +
  geom_line(data=moGg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = moGg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.25,0.25), breaks = seq(-0.2,0.2, by = 0.2))

plot3 <-ggplot( ) +
  labs(y =" ",x="Land-use intensity") + # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=evGg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = evGg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evGg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = evGg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot4 <-ggplot( ) +
  labs(y ="Connectance (scaled)",x=" ") + # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=wdFg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = wdFg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=wdFg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = wdFg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  ggtitle("Forest")+ theme(plot.title = element_text(size=13, face="bold",hjust=0.5))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))

plot5 <-ggplot( ) +
  labs(y ="Modularity (scaled)",x=" ") +# change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=moFg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = moFg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+
  theme_classic()  +
  geom_line(data=moFg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = moFg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.25,0.25), breaks = seq(-0.2,0.2, by = 0.2))

plot6 <-ggplot( ) +
  labs(y ="Evenness (scaled)",x="Land-use intensity") +  # change y label and see limits: + ylim (-0.1,0.1)
  geom_line(data=evFg_RAWpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "darkgrey") +
  geom_ribbon(data = evFg_RAWpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="darkgrey",alpha = .15)+ 
  theme_classic()  +
  geom_line(data=evFg_SIGNpos_conf, aes(x = mean_LUI,y = est), size = 1,colour = "blue") +
  geom_ribbon(data = evFg_SIGNpos_conf, aes(x = mean_LUI,y = NULL, ymin = lower, ymax = upper),fill="blue",alpha = .15)+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = c(-0.015,0.015), breaks = seq(-0.01,0.01, by = 0.01))


windows ()
gridPrint(grobs = list(plot4,plot1,plot5,plot2,plot6,plot3), ncol = 2) # to have forest on the left and grasslands on the right columns



# 5. FIGURE 1 ####
library(igraph)
library (ggplot2)
# function to align the vertex labels nicely (from https://kieranhealy.org/blog/archives/2011/02/18/aligning-labels-in-circular-igraph-layouts/)
radian.rescale <- function(x, start=0, direction=1) { #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

#create dataset
vertex_info3 <- data.frame (Nodes = c("BD1","BD2","BD3","EF1","EF2","EF3","ES1","ES2","ES3"), Node_type = c(rep("BD",3), rep("EF",3),rep("ES",3)))    


##CONNECTANCE LOW LUI
edges_lowLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                     rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                            to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                   rep(c("ES1","ES2","ES3"),3)), 
                            weight = c(3.5,1.5,0,1,0,2,
                                       0,2,3,0,2,0,
                                       1,3,0,0,4,0.5,
                                       1,1,0,
                                       0,2,0.5,
                                       3,1,0.5))
#exclude 0s so it doesn't plot them!
edges_lowLUI3 <- edges_lowLUI3 [(edges_lowLUI3$weight!=0),]

#plot the network
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
Net <-lowLUI.Net3
#V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
lab.locs <- radian.rescale(x=1:length(V(lowLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

lowLUI.Net3<- plot(lowLUI.Net3, # only change this line
                   edge.width=E(Net)$weight*5,
                   vertex.color=V(Net)$vertex_colour, 
                   vertex.label.dist=2.5,
                   vertex.label.cex=1.1,
                   vertex.label.font=2, #bold
                   vertex.label.degree=lab.locs,
                   vertex.label.color="black", 
                   #vertex.size=V(Net)$vertex.size*8,
                   layout=layout_in_circle(Net))



# network metrics
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
lowLUI_density <-edge_density(lowLUI.Net3, loops=F) #0.5


##MODULARITY + modules LOW LUI
edges_lowLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                     rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                            to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                   rep(c("ES1","ES2","ES3"),3)), 
                            weight = c(
                              0,0,0,0,4,0,
                              0,0,0,0,4,1.5,
                              3.5,3,0,0,0,0,
                              0,0,0,
                              0,0,0,
                              3,0,0))
#exclude 0s so it doesn't plot them!
edges_lowLUI3 <- edges_lowLUI3 [(edges_lowLUI3$weight!=0),]

#plot the network
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
Net <-lowLUI.Net3
#V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
lab.locs <- radian.rescale(x=1:length(V(lowLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

lowLUI.Net3<- plot(lowLUI.Net3, # only change this line
                   edge.width=E(Net)$weight*5,
                   vertex.color=V(Net)$vertex_colour, 
                   vertex.label.dist=2.5,
                   vertex.label.cex=1.1,
                   vertex.label.font=2, #bold
                   vertex.label.degree=lab.locs,
                   vertex.label.color="black", 
                   #ertex.size=V(Net)$vertex.size*8,
                   layout=layout_in_circle(Net))


# network metrics
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
lowLUI_modularity<-modularity(lowLUI.Net3,membership(cluster_walktrap(lowLUI.Net3)), weights=abs(E(lowLUI.Net3)$weight)) #0.61


### EVENNESS LOW LUI
edges_lowLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                     rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                            to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                   rep(c("ES1","ES2","ES3"),3)), 
                            weight = c(1,5,0,0,0,0,
                                       0,0,0.1,2,0,1,
                                       0,0,0,0,4,1.5,
                                       0,0.7,0,
                                       0,0,0,
                                       1,0,0))
#exclude 0s so it doesn't plot them!
edges_lowLUI3 <- edges_lowLUI3 [(edges_lowLUI3$weight!=0),]

#plot the network
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
Net <-lowLUI.Net3
#V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
lab.locs <- radian.rescale(x=1:length(V(lowLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

lowLUI.Net3<- plot(lowLUI.Net3, # only change this line
                   edge.width=E(Net)$weight*5,
                   vertex.color=V(Net)$vertex_colour, 
                   vertex.label.dist=2.5,
                   vertex.label.cex=1.1,
                   vertex.label.font=2, #bold
                   vertex.label.degree=lab.locs,
                   vertex.label.color="black", 
                   #ertex.size=V(Net)$vertex.size*8,
                   layout=layout_in_circle(Net))

#evenness value
diversity(edges_lowLUI3 [,"weight"])/log(length(V(lowLUI.Net3))) #0.85



## HIGH LUI: CONNECTANCE, MODULARITY, MODULES
edges_highLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                      rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                             to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                    rep(c("ES1","ES2","ES3"),3)), 
                             weight = c(0,0,2,0,0,0,
                                        0,1,0,0,1,0,
                                        1,1,0,0,0,3,
                                        1,0,0,
                                        0,0,0,
                                        1,1,0))
#exclude 0s so it doesn't plot them!
edges_highLUI3 <- edges_highLUI3 [(edges_highLUI3$weight!=0),]

#plot the network
highLUI.Net3 <- graph_from_data_frame(d=edges_highLUI3,vertices=vertex_info3, directed=F)
Net <-highLUI.Net3
#V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
lab.locs <- radian.rescale(x=1:length(V(highLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

highLUI.Net3<- plot(highLUI.Net3, # only change this line
                    edge.width=E(Net)$weight*5,
                    vertex.color=V(Net)$vertex_colour, 
                    vertex.label.dist=2.5,
                    vertex.label.cex=1.1,
                    vertex.label.font=2, #bold
                    vertex.label.degree=lab.locs,
                    vertex.label.color="black",
                    #vertex.size=V(Net)$vertex.size*8,
                    layout=layout_in_circle(Net))




## EVENNESS HIGH LUI
edges_highLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                      rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                             to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                    rep(c("ES1","ES2","ES3"),3)), 
                             weight = c(1,1,1,0,1,1,
                                        0,1,0,0,1,0,
                                        1,1,0,1,0,2,
                                        1,0,1,
                                        0,1,0,
                                        1,1,0))
#exclude 0s so it doesn't plot them!
edges_highLUI3 <- edges_highLUI3 [(edges_highLUI3$weight!=0),]

#plot the network
highLUI.Net3 <- graph_from_data_frame(d=edges_highLUI3,vertices=vertex_info3, directed=F)
Net <-highLUI.Net3
#V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
lab.locs <- radian.rescale(x=1:length(V(highLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

highLUI.Net3<- plot(highLUI.Net3, # only change this line
                    edge.width=E(Net)$weight*5,
                    vertex.color=V(Net)$vertex_colour, 
                    vertex.label.dist=2.5,
                    vertex.label.cex=1.1,
                    vertex.label.font=2, #bold
                    vertex.label.degree=lab.locs,
                    vertex.label.color="black",
                    #vertex.size=V(Net)$vertex.size*8,
                    layout=layout_in_circle(Net))
#evenness value
diversity(edges_highLUI3 [,"weight"])/log(length(V(highLUI.Net3))) #1.25... should be max 1??




# HUBS:  LOW LUI
edges_lowLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                     rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                            to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                   rep(c("ES1","ES2","ES3"),3)), 
                            weight = c(3.5,1.5,0,0,0,0,
                                       0,0,0,0,2,0,
                                       0,0,0,0,4,0.5,
                                       1,1,0,
                                       0,0,0.5,
                                       3,0,0.5))
#exclude 0s so it doesn't plot them!
edges_lowLUI3 <- edges_lowLUI3 [(edges_lowLUI3$weight!=0),]

#plot the network
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
Net <-lowLUI.Net3
V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
#lab.locs <- radian.rescale(x=1:length(V(lowLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

lowLUI.Net3<- plot(lowLUI.Net3, # only change this line
                   edge.width=E(Net)$weight*5,
                   vertex.color=V(Net)$vertex_colour, 
                   # vertex.label.dist=0.5,# set to 5 if using lab.locs
                   # vertex.label.cex=1.1,
                   # vertex.label.font=2, #bold
                   # vertex.label.degree=lab.locs,
                   vertex.label.color="transparent", # "black" # vertex.labels=NA,##no va
                   vertex.size=V(Net)$vertex.size*8,
                   layout=layout_in_circle(Net))




##Hubs: HIGH LUI
edges_highLUI3 <- data.frame(from = c(rep("BD1",6),rep("BD2",6),rep("BD3",6),
                                      rep("EF1",3),rep("EF2",3),rep("EF3",3)), 
                             to = c(rep(c("EF1","EF2","EF3","ES1","ES2","ES3"),3),
                                    rep(c("ES1","ES2","ES3"),3)), 
                             weight = c(0,0,2,0,0,0,
                                        0,1,0,0,1,0,
                                        1,1,0,0,0,3,
                                        1,0,0,
                                        0,0,0,
                                        1,1,0))
#exclude 0s so it doesn't plot them!
edges_highLUI3 <- edges_highLUI3 [(edges_highLUI3$weight!=0),]

#plot the network
highLUI.Net3 <- graph_from_data_frame(d=edges_highLUI3,vertices=vertex_info3, directed=F)
Net <-highLUI.Net3
V(Net)$vertex.size=strength(Net, vids = V(Net), mode = "all",loops = F, weights = abs(E(Net)$weight))
V(Net)$vertex_colour = c(rep("gold",3), rep ("lightblue",3), rep("indianred",3)) 
#lab.locs <- radian.rescale(x=1:length(V(highLUI.Net3)), direction=-1, start=0) #start = offset from 12 o'clock in radians; direction = 1 for clockwise; -1 for anti-clockwise.

highLUI.Net3<- plot(highLUI.Net3, # only change this line
                    edge.width=E(Net)$weight*5,
                    vertex.color=V(Net)$vertex_colour, 
                    # vertex.label.dist=0.5,# set to 4 if using lab.locs
                    # vertex.label.cex=1.1,
                    # vertex.label.font=2, #bold
                    # vertex.label.degree=lab.locs,
                    vertex.label.color="transparent", # "black" # vertex.labels=NA,##no va
                    vertex.size=V(Net)$vertex.size*8,
                    layout=layout_in_circle(Net))




# network metrics
lowLUI.Net3 <- graph_from_data_frame(d=edges_lowLUI3,vertices=vertex_info3, directed=F)
highLUI.Net3 <- graph_from_data_frame(d=edges_highLUI3,vertices=vertex_info3, directed=F)

lowLUI_density <-edge_density(lowLUI.Net3, loops=F) #0.28
highLUI_density <-edge_density(highLUI.Net3, loops=F) #0.25

lowLUI_modularity<-modularity(lowLUI.Net3,membership(cluster_walktrap(lowLUI.Net3)), weights=abs(E(lowLUI.Net3)$weight)) #0.48
highLUI_modularity<-modularity(highLUI.Net3,membership(cluster_walktrap(highLUI.Net3)), weights=abs(E(highLUI.Net3)$weight))#0.40

library(vegan)
H <- diversity(edges_lowLUI3 [,"weight"])
S <- length(V(lowLUI.Net3))
lowLUI_J <- H/log(S) #0.93

H <- diversity(edges_highLUI3 [,"weight"])
S <- length(V(highLUI.Net3))
highLUI_J <- H/log(S) #0.95



#### LUI VALUES IN A LINE
#LUI from BExIS
LUI2008_2015 <- read.table("...\\LUI_2008-2015.txt",header=TRUE, sep="\t", dec=",")
head(LUI2008_2015)

#ForMI from Kahl and Bauhus 2014
ForMI <- read.table("...\\16466_ForMI - Forest Management Intensity Index_1.3.2\\16466.txt",header=TRUE, sep="\t", dec=".")
head(ForMI)


library(ggplot2)
library(RColorBrewer)

figLUI_G<-data.frame(y=0,x=LUI2008_2015$LUI)
figLUI_F<-data.frame(y=0,x=ForMI$ForMI)

ggplot(figLUI_G, aes(x=x)) +
  ylim(0,0.5)+ xlim(0,3.5)+
  geom_point(aes(x,y,color=x), size = 2)+ #colour = CONTOUR, fill = INSIDE
  theme_classic()+
  scale_color_gradient(low="blue", high="red")

ggplot(figLUI_F, aes(x=x)) +
  ylim(0,0.5)+ xlim(0,3.5)+
  geom_point(aes(x,y,color=x), size = 2)+ #colour = CONTOUR, fill = INSIDE
  theme_classic()+
  scale_color_gradient(low="blue", high="red")
