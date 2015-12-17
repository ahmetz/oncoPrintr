#' oncoPrint
#'
#' @param df  this is the dataframe used to plot oncoprints. required colnames: Sample, Gene, VarClass
#' @param sort Boelean indicating whether genes should be sorted or now (default: True)
#' @param convert Boelean indicating whether varclasses should be converted to more standard names (default: True)
#' @param total_samples Total number of samples. If not given, total row numbers of input will be used
#' @param geneName IF given, this gene sits on top of the oncoprint, it won't be used in sorting
#' @param annotations If given, these will be used for pre-clustering
#' @param annotation_order a list of annotation items to be displayed
#' @param merge_scnas Bolean indicating whether to plot scnas and mutations together superimposed (default: False)
#' @param second_df if merge_scnas is True, this is the second set of events to be merged
#' @param alteration_score this list determines the relative importance of different genomic alterations. Amplification > Deletion > Mutations etc.
#' @param printSamples Bolean indicating whether sample names should be printed under the oncoprint
#'
#' @return oncoPrint image
#' @export
#'
#' @examples TODO

oncoPrint <- function(df, sort=TRUE, convert = TRUE, total_samples = NA, geneName = NA, annotation = NA, annotation_order = NA, merge_scnas = F, df2 = NA, colors = list(Amplification = "red", Deletion = "blue"), alteration_score = list(Amplification = 5, Deletion = 4, Nonsense = 2.8, Frameshift = 2.5, Splicing = 2.5, InFrame = 2, Promoter = 2, Mutation =1, Missense=1, Present = 1, NotTested = 0, None = 0, NotPresent = 0, del = 3, homodel = 2, LOH = 1.5, CNLOH = 1), printSamples = T, xpadding = 0.1, ypadding = 0.1) {
  # This is the plotting function
  library(reshape2)
  library(dplyr)

  colnames(df) <- c("Sample", "Gene", "VarClass")
  if(!is.na(df2)){
    colnames(df2) <- c("Sample", "Gene", "VarClass")
  }
  
  remove_duplicates <- function(df){
  cat("Filtering duplicates\n")
  choose_alteration_type <- function(x){
    for (i in 1:length(x)){
      #cat(length(x), "\n")
      #cat(x[i, ]$Gene, ":", x[i, ], ":", i)
      if(grepl("stopgain", x[i])){
        return(x[grep("stopgain", x[i])])
      }else if(grepl("insertion", x[i])){
        return(x[grep("insertion", x[i])])
      }else if(grepl("deletion", x[i])){
        return(x[grep("deletion", x[i]) ])
      }else if(grepl("splicing", x[i])){
        return(x[grep("splicing", x[i]) ])
      }else if(grepl("nonsynonymous", x[i])){
        return(x[grep("nonsynonymous", x[i]) ])
      }else if(grepl("upstream", x[i])){
        return(x[grep("upstream", x[i]) ])
      }else if(grepl("Amplification", x[i])){
        return(x[grep("Amplification", x[i]) ])
      }else if(grepl("IntragenicDeletion", x[i])){
        return(x[grep("IntragenicDeletion", x[i]) ])
      }else if(grepl("Deletion", x[i])){
        return(x[grep("Deletion", x[i]) ])
      }else if(grepl("Present", x[i])){
        return(x[grep("Present", x[i]) ])
      }else if(grepl("NotTested", x[i])){
        return(x[grep("NotTested", x[i]) ])
      }else if(grepl("del", x[i])){
        return(x[grep("del", x[i]) ])
      }else if(grepl("homo-del", x[i])){
        return(x[grep("homo-del", x[i]) ])
      }else if(grepl("CN-LOH", x[i])){
        return(x[grep("CN-LOH", x[i]) ])
      }else if(grepl("LOH", x[i])){
        return(x[grep("LOH", x[i]) ])
      }else if(grepl("None", x[i])){
        return(x[grep("None", x[i]) ])
      }
    }
  }
  df %>% group_by(Gene, Sample) %>% summarise(VarClass = choose_alteration_type(VarClass))
}
  
  
  # check if there are any samples with no alterations. If so, remove them to add later on. 
  # These will have VarClass = "None"
  not_altered_sample_num <- NA
  not_altered_sample_names <- NA
  if(!merge_scnas){
    if(nrow(df %>% filter(VarClass == "None")) > 0){
      not_altered_sample_num <- nrow(df %>% filter(VarClass == "None"))
      not_altered_sample_names <- df %>% filter(VarClass == "None") %>% select(Sample)
      df <- df %>% filter(VarClass != "None")
    }
  }
  
  cat("line86\n")
  #remove duplicates of gene events within the same sample.
  #TO-DO do not remove if a gene has both a copy number alteration and a mutation
  df <- remove_duplicates(df)
  if(merge_scnas){
    df2 <- remove_duplicates(df2)
  }
  
  if(convert){# change the variant types to more general names
    df <- convert_varclass(df)
    if(merge_scnas){
      df2 <- convert_varclass(df2)
      
    }
    
  }else{
    df[grep("splicing", df$VarClass), ]["VarClass"] <- "Mutation"
    df[grep("stop", df$VarClass), ]["VarClass"] <- "Mutation"
    df[grep("nonsynonymous", df$VarClass), ]["VarClass"] <- "Mutation"
    df[grep("insertion", df$VarClass), ]["VarClass"] <- "Mutation"
    df[grep("deletion", df$VarClass), ]["VarClass"] <- "Mutation"
    df[grep("upstream", df$VarClass), ]["VarClass"] <- "Mutation"
  }
  
  #remove duplicates of gene events within the same sample.
  #TO-DO do not remove if a gene has both a copy number alteration and a mutation
  
  #df <- df %>% group_by(Gene, Sample) %>% unique()
  
  # if there is an annotation data frame, then figure out how many samples there are with no mutations and add them to the alterations matrix
  if(merge_scnas && !is.na(annotation)){
    alts <- acast(df, Gene ~ Sample)
    alts2 <- acast(df2, Gene ~ Sample)
    alterations <- paste.matrix(alts, alts2)
    colnames(annotation) <- c("sample", "class")
    annotation.samples <- annotation$sample
    
    missing.sample <- annotation[which(!annotation.samples%in%colnames(alterations)), ]
    missing.matrix <- matrix(NA, nrow = nrow(alterations), ncol = nrow(missing.sample))
    colnames(missing.matrix) <- missing.sample$sample
    
    alterations <- cbind(alterations, missing.matrix)
    alterations.c <- matrix(as.numeric(!is.na(alterations)), ncol = ncol(alterations))
    colnames(alterations.c) <- colnames(alterations)
    row.names(alterations.c) <- row.names(alterations)
    
  }else if(!is.na(annotation)){
    alterations <- acast(df, Gene ~ Sample)
    colnames(annotation) <- c("sample", "class")
    annotation.samples <- annotation$sample
    
    missing.sample <- annotation[which(!annotation.samples%in%colnames(alterations)), ]
    missing.matrix <- matrix(NA, nrow = nrow(alterations), ncol = nrow(missing.sample))
    colnames(missing.matrix) <- missing.sample$sample
    
    alterations <- cbind(alterations, missing.matrix)
    alterations.c <- matrix(as.numeric(!is.na(alterations)), ncol = ncol(alterations))
    colnames(alterations.c) <- colnames(alterations)
    row.names(alterations.c) <- row.names(alterations)
  }else if(merge_scnas){
    alts <- acast(df, Gene ~ Sample)
    alts2 <- acast(df2, Gene ~ Sample)
    alterations <- paste.matrix(alts, alts2)
    alterations.c <- matrix(as.numeric(!is.na(alterations)), ncol = ncol(alterations))
    colnames(alterations.c) <- colnames(alterations)
    row.names(alterations.c) <- row.names(alterations)
  }else{
  cat("line153\n")
    alterations.c <- acast(df, Gene ~ Sample, fun.aggregate = length) # This is the 0 and 1 version of the matrix
    alterations <- acast(df, Gene ~ Sample)
  }

  #convert variant type matrix to numerical values
  for (i in 1:nrow(alterations)){
    for(j in 1:ncol(alterations)){
      altered <- alterations[i, j]
      
      if(!is.na(altered)){ # there is an alteration
        if(grepl("," ,altered)){ # alteration is a mix of two seperated by a comma
          alts <- unlist(str_split(altered, ",")) # split the alterations
          alt1 <- alts[1] #assign them individually
          alt2 <- alts[2]
          # first alteration
          alt <- alt1
          alterations.c[i, j] <- alteration_score[[alt]]
          alt <- alt2
          alterations.c[i, j] <- alterations.c[i, j] + alteration_score[[alt]]
          
        }else{
          alt <- altered
          alterations.c[i, j] <- alteration_score[[alt]]
        }
      }else{
        alterations.c[i, j] <- 0
      }
    }
  }
  
  # Order the samples
  if(is.na(geneName)){
    geneName <- row.names(alterations)
  }
  alterations.c <- memoSort(alterations.c, geneOrder = geneName, annotations = annotation, annotation_order = annotation_order)
  alterations <- alterations[row.names(alterations.c), colnames(alterations.c)]
  
  ngenes <- nrow(alterations);
  nsamples <- ncol(alterations);

  
  # if there any samples with no alterations, add them here: 
  if(!is.na(not_altered_sample_num) ){
    mat <- matrix(data = rep(NA, ngenes*not_altered_sample_num), ncol = not_altered_sample_num, nrow = ngenes)
    colnames(mat) <- not_altered_sample_names$Sample
    alterations <- cbind(alterations, mat)
  }
  
  # check the total_samples variable. If there is a value given, make sure the # of samples match that
  nsamples <- ncol(alterations)
  if(!is.na(total_samples)){
    if(total_samples != nsamples){
      diff <- total_samples - nsamples
      mat <- matrix(data = rep(NA, ngenes*diff), ncol = diff, nrow = ngenes)
      colnames(mat) <- paste(rep("MockSample", diff), "_", 1:diff, sep="")
      alterations <- cbind(alterations, mat)
      alterations.c <- cbind(alterations.c, mat)
    }
  }
  
  ngenes <- nrow(alterations)
  nsamples <- ncol(alterations)
  #cat("#ofGenes: ", ngenes, " # of Samples: ", nsamples)
  ### OncoPrint
  numOfOncos <- ngenes*nsamples
  oncoCords.base <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  oncoCords <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords.base) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  #crate a third matrix for copy number calls. this will be a second layer on top of grey background rects
  oncoCords.scna <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos );
  colnames(oncoCords.scna) <- c("xleft", "ybottom", "xright", "ytop", "altered");
  xpadding <- xpadding;
  ypadding <- ypadding;
  cnt <- 1;
  
  cat("nsamples: ", nsamples, " ngenes: ", ngenes, "\n")
  
  if (merge_scnas){ 

    for(i in 1:ngenes) {
      for(j in 1:nsamples) {
        #cat("i: ", i, " j: ", j, "\n")
        altered <- alterations[i, j]
        xleft <- j-1 + xpadding ;
        ybottom <- ((ngenes-i+1) -1) + ypadding;
        xright <- j - xpadding ;
        ytop <- (ngenes-i+1) -ypadding;
        oncoCords.base[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
        #browser()
        
        if(!is.na(altered)){ # there is an alteration
          if(grepl("," ,altered)){ # alteration is a mix of two seperated by a comma
            alts <- unlist(str_split(altered, ",")) # split the alterations
            alt1 <- alts[1] #assign them individually
            alt2 <- alts[2]
            # first alteration
            altered <- alt1  
            if(altered == "Mutation" || altered == "Missense" || altered == "Nonsense" ||altered == "Splicing" || altered == "Frameshift" || altered == "Promoter" || altered == "InFrame") {
              ytop2 <- ytop-0.25
              ybottom2 <- ybottom+0.25
              oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered);
            }else if( altered == "Amplification" || altered == "Deletion" || altered == "Present" || altered == "NotPresent" || altered == "NotTested"|| altered == "homodel" || altered == "del" || altered == "CNLOH" || altered == "LOH"){
              oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
            }
            # second alteration
            altered <- alt2
            if(altered == "Mutation" || altered == "Missense" || altered == "Nonsense" ||altered == "Splicing" || altered == "Frameshift" || altered == "Promoter" || altered == "InFrame") {
              ytop2 <- ytop-0.25
              ybottom2 <- ybottom+0.25
              oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered);
            }else if( altered == "Amplification" || altered == "Deletion" || altered == "Present" || altered == "NotPresent" || altered == "NotTested"|| altered == "homodel" || altered == "del" || altered == "CNLOH" || altered == "LOH"){
              oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
            }
            
            
          }else{ # alteration does not have a comma
            if(altered == "Mutation" || altered == "Missense" || altered == "Nonsense" ||altered == "Splicing" || altered == "Frameshift" || altered == "Promoter" || altered == "InFrame") {
              ytop2 <- ytop-0.25
              ybottom2 <- ybottom+0.25
              oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered);
            }else if( altered == "Amplification" || altered == "Deletion" || altered == "Present" || altered == "NotPresent" || altered == "NotTested" || altered == "homodel" || altered == "del" || altered == "CNLOH" || altered == "LOH"){
              oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
            }
          }
          
        }else{  # There is no alteration
          oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
        }
        cnt <- cnt+1;
      }
    }
    
  }else{
    for(i in 1:ngenes) {
      for(j in 1:nsamples) {
        altered <- alterations[i, j];
        xleft <- j-1 + xpadding ;
        ybottom <- ((ngenes-i+1) -1) + ypadding;
        xright <- j - xpadding ;
        ytop <- (ngenes-i+1) -ypadding;
        oncoCords.base[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
        #browser()
        if(!is.na(altered)){
          if(altered == "Mutation" || altered == "Missense" || altered == "Nonsense" ||altered == "Splicing" || altered == "Frameshift" || altered == "Promoter" || altered == "InFrame") {
            ytop2 <- ytop-0.25
            ybottom2 <- ybottom+0.25
            oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered);
          }else if( altered == "Amplification" || altered == "Deletion" || altered == "Present" || altered == "NotPresent" ||altered == "NotTested" || altered == "homodel" || altered == "del" || altered == "CNLOH" || altered == "LOH"){
            oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
          }else{
            oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
          }
        }else{
          oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered);
        }
        cnt <- cnt+1;
      }
    }
  }
  
  labels = rownames(alterations)
  gene_prcnt <- list()
  if(!is.na(total_samples)){
    labels=list()
    for(i in (1:nrow(alterations))){
      gene = row.names(alterations)[i]
      cnt <- sum(ifelse(test = is.na(alterations[i,]), yes = 0, no = ifelse(test =grepl(",", alterations[i,]), yes = 2, no = 1)))
      prcnt <- cnt/length(alterations.c[i,])*100
      gene_prcnt[gene] <- prcnt
      labels = c(labels, paste(gene, "  ",round(prcnt,1), "%", sep=""))
    }
  }
  
  cnt <- nsamples*ngenes
  colors <- rep(NA, cnt)
  colors.scna <- rep(NA, cnt)
  colors[ which(oncoCords[, "altered"] == "Mutation") ] <- "#26A818"
  colors[ which(oncoCords[, "altered"] == "Missense") ] <- "#26A818"
  colors[ which(oncoCords[, "altered"] == "Nonsense") ] <- "black"
  colors[ which(oncoCords[, "altered"] == "Splicing") ] <- "#A05E35"
  colors[ which(oncoCords[, "altered"] == "Frameshift") ] <- "#A05E35"
  colors[ which(oncoCords[, "altered"] == "Promoter") ] <- "#2986E2"
  colors[ which(oncoCords[, "altered"] == "InFrame") ] <- "#F26529"
  colors[ which(oncoCords[, "altered"] == "Present") ] <- "black"

  colors.scna[ which(oncoCords.scna[, "altered"] == "Present") ] <- "darkorchid2"
  colors.scna[ which(oncoCords.scna[, "altered"] == "NotPresent") ] <- "#DCD9D3"
  colors.scna[ which(oncoCords.scna[, "altered"] == "NotTested") ] <- "darkgrey"
  colors.scna[ which(oncoCords.scna[, "altered"] == "del") ] <- "red"
  colors.scna[ which(oncoCords.scna[, "altered"] == "LOH") ] <- "darkkhaki"
  colors.scna[ which(oncoCords.scna[, "altered"] == "homodel") ] <- "brown4"
  colors.scna[ which(oncoCords.scna[, "altered"] == "CNLOH") ] <- "deepskyblue"
  #colors.scna[ which(oncoCords.scna[, "altered"] == "Amplification") ] <- "#EA2E49"
  colors.scna[ which(oncoCords.scna[, "altered"] == "Amplification") ] <- "blue"
  #colors.scna[ which(oncoCords.scna[, "altered"] == "Deletion") ] <- "#174D9D"
  colors.scna[ which(oncoCords.scna[, "altered"] == "Deletion") ] <- "red"
  
  c48 <- c("#1d915c","#5395b4","#964a48","#2e3b42","#b14e72", "#402630","#f1592a","#81aa90","#f79a70","#b5ddc2","#8fcc8b","#9f1f63","#865444", "#a7a9ac","#d0e088","#7c885c","#d22628","#343822","#231f20","#f5ee31","#a99fce","#54525e","#b0accc","#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")
  
  #cat("\n", "samples*genes: ", cnt, "length of colors", length(colors))
  #change the 
  def.par <- par(no.readonly = TRUE)
  leftmargin = 1/nsamples*500 
  if(leftmargin > 20){leftmargin <- 20}
  bottommargin = 1/ngenes*500
  if(bottommargin > 5){bottommargin <- 1}
  if(!is.na(annotation)){
    
    split.screen(rbind(c(0.05, 0.95, 0.95, 0.99), c(0.05,0.95,0.15, 0.94), c(0.05, 0.95, 0.01, 0.15)))
    screen(1)
    par(mar=c(0,10,0, 0), mgp=c(3, 0.7, 0))
    plot(c(0, nsamples), c(0,1), type="n", main="", xlab="Samples", xaxt="n", ylab="", yaxt="n", frame.plot = F)
    counts <- data.frame(table(annotation$class))
    xleft <- 0
    axis.points <- list()
    subtype.labels <-list()
    for (i in 1:length(annotation_order)){
      subtype.labels <- c(subtype.labels, annotation_order[i])
      xright <- counts[which(counts$Var1 == annotation_order[i]), ]['Freq']
      rect(xleft, 0, xleft + xright$Freq, 1, col = c48[i])
      axis.points <- c(axis.points, (xleft + xleft+xright$Freq)/2)
      xleft <- xleft + xright$Freq
    }
    text(x=axis.points, y = 0.5, labels = subtype.labels)
    
    
    screen(2)
    par(mar=c(0,10,0,0), mgp=c(3, 0.7, 0))
    plot(c(0, nsamples), c(0, ngenes), type="n", main="", xlab="Samples", xaxt="n", ylab="", yaxt="n", frame.plot = F)
    rect(oncoCords.base[, "xleft"], oncoCords.base[, "ybottom"],oncoCords.base[, "xright"], oncoCords.base[, "ytop"], col="#DCD9D3", border=NA);
    rect(oncoCords.scna[, "xleft"], oncoCords.scna[, "ybottom"],oncoCords.scna[, "xright"], oncoCords.scna[, "ytop"], col=colors.scna, border=NA);
    rect(oncoCords[, "xleft"], oncoCords[, "ybottom"],oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border=NA);
    
    axis(2, at=(ngenes:1)-.5, labels=labels, las=2, lwd = 0);
    
        
    #add legend
    screen(3)
    par(mar=c(0,0,0,0))
    legend(x="topleft", c("Missense", "Nonsense", "Truncating", "In-Frame", "Promoter"), fill = c('#26A818', 'black',  '#A05E35', '#F26529', '#2986E2'), horiz=T, border = F, cex=0.9, bty = 'n')
    
    legend(x="bottomleft", c("Amplification", "Deletion", "LOH", "Present"), fill = c( '#EA2E49', '#174D9D', 'darkkhaki', 'darkorchid2'), horiz=T, border = F, cex=0.9, bty = 'n', x.intersp = 0.5)
    close.screen(all.screens = TRUE)
    
  }else{
    split.screen(rbind(c(0.05,0.95,0.15, 0.95), c(0.05, 0.95, 0.05, 0.15)))
    
    #add oncoprints
    screen(1)
    
    par(mar=c(1,10,0.25, 9), mgp=c(3, 0.7, 0))
    plot(c(0, nsamples), c(0, ngenes), type="n", main="", xlab="Samples", xaxt="n", ylab="", yaxt="n", frame.plot = F);
    rect(oncoCords.base[, "xleft"], oncoCords.base[, "ybottom"],oncoCords.base[, "xright"], oncoCords.base[, "ytop"], col="#DCD9D3", border=NA);
    rect(oncoCords.scna[, "xleft"], oncoCords.scna[, "ybottom"],oncoCords.scna[, "xright"], oncoCords.scna[, "ytop"], col=colors.scna, border=NA);
    rect(oncoCords[, "xleft"], oncoCords[, "ybottom"],oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border=NA);
    
    axis(2, at=(ngenes:1)-.5, labels=labels, las=2, lwd = 0);
    #printing samples or not
    if(printSamples){
      text((1:nsamples)-.5, par("usr")[3]+1,srt=45, adj = 1,  labels = colnames(alterations), xpd=T, cex = 0.6)
    }
    #add legend
    screen(2)
    par(mar=c(0,0,0,0))
    legend(x="topleft", c("Missense mutation", "Nonsense mutation", "Truncating mutation", "In-Frame mutation", "Promoter mutation"), fill = c('#26A818', 'black',  '#A05E35', '#F26529', '#2986E2'), horiz=T, border = F, cex=0.9, bty = 'n')
    legend(x="bottomleft", c( "Amplification", "Deletion", "Present", "LOH" ,"CNLOH"), fill = c('blue', 'red', 'darkorchid2', 'darkkhaki', 'deepskyblue'), horiz=T, border = F, cex=0.9, bty = 'n')
    close.screen(all.screens = TRUE)
  }
  par(def.par) 
  res <- list()
  res$sortedMatrix <- alterations
  res$sampleOrder <- colnames(alterations)
  res$geneOrder <- rownames(alterations)
  res$gene_prcnt <- gene_prcnt
  res$alterations.value <- alterations.c
  res
}

