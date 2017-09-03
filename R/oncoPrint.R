#' oncoPrint 
#' 
#' This function creates onco prints of genomic data based on mutations, copy number alterations, fusions, and other user defined categories. For each category 3 fields are required : Sample ID, Gene, and Variant Class. There are options to display binary events as well as continous variables. 
#'
#' @param data  A list of dataframes to be included in the oncoprint. For each data frame required columns are: Sample, Gene, Variant Type
#' @param sort Boelean indicating whether genes should be sorted or now (default: True)
#' @param convert Boelean indicating whether varclasses should be converted to more standard names (default: True)
#' @param total_samples Total number of samples. If not given, total row numbers of input will be used
#' @param geneName If given, genes are not automatically sorted. Instead this variable is used for sorting. Partial list of genes can be input
#' @param annotations If given, these will be used for creating sub-groups of samples. Clustering will happen within annotatoion groups and then merged
#' @param annotation_order This is the order of annotation categories to display. Required if annotations are given
#' @param alteration_score this list determines the relative importance of different genomic alterations. Amplification > Deletion > Mutations etc.
#' @param printSamples Bolean indicating whether sample names should be printed under the oncoprint
#' @param xpadding numerical value of the padding between two consecutive data columns
#' @param ypadding numerical value of the padding between two consecutive data rows

#' @return oncoPrint image and a list of values
#' @export
#'
#' @examples TODO
#' @import reshape2
#' @import dplyr
#' @import grid
#' @import stringr


oncoPrint <- function(data = NULL, sort=TRUE, convert = TRUE, total_samples = NULL, geneName = NULL, annotation = NULL, annotation_order = NULL, continuous_data = NULL, categorical_data = NULL, categorical_data_colors = NULL, annotation_colors = NULL, categorical_data_order = NULL, onco_colors = list(Mutation = "#26A818", Missense = "#26A818", Nonsense = "black", Splicing = "#ffaa00", Frameshift = "#A05E35" , Promoter = "#2986E2", InFrame = "#F26529", Present = "darkorchid2", NotPresent = "#DCD9D3", NotTested = "darkgrey", del = "red", LOH = "#D17878", homodel = "brown4", CNLOH =  "deepskyblue", Amplification = "#EA2E49", Deletion = "#174D9D", Yes = "#155B6B", No = "#12C8F9", Unknown = "azure1", Fusion =  "#B641F9", Pathogenic="white"), alteration_score = list(Amplification = 5, Fusion = 4.5, Deletion = 4, Pathogenic = 3,  Nonsense = 2.8, Frameshift = 2.6, Splicing = 2.5, InFrame = 2, Promoter = 2, Mutation =1, Missense=1, Present = 1, NotTested = 0, None = 0, NotPresent = 0, Yes = 0, No = 0, del = 3, homodel = 2, LOH = 1.5, CNLOH = 1), printSamples = F, xpadding = 0.1, ypadding = 0.1) {
  
  # The function here started with the gist from Arman Aksoy here https://gist.github.com/armish/564a65ab874a770e2c26 and developped into this
  
  
  # This is the plotting function
  require(reshape2)
  require(dplyr)
  require(stringr)
  
  if (is.data.frame(data)){
    merge_scnas = F
    df = data
    df2 = NULL
  }else{
    if (length(data) == 1) {
      merge_scnas = F
      df2 = NULL
      df = data[[1]]
    }else{
      merge_scnas = T
      df = data[[1]]
      df2 = data[2:length(data)]
    }
  }
  colnames(df) <- c("Sample", "Gene", "VarClass")
  df$Sample <- as.character(df$Sample)
  df$Gene <- as.character(df$Gene)
  df$VarClass <- as.character(df$VarClass)
  
  # check if there are any samples with no alterations. If so, remove them to add later on. 
  # These will have VarClass = "None"
  not_altered_sample_num <- NA
  not_altered_sample_names <- NA
  if(!merge_scnas){
    if(nrow(df %>% dplyr::filter(VarClass == "None")) > 0){
      not_altered_sample_num <- nrow(df %>% dplyr::filter(VarClass == "None"))
      not_altered_sample_names <- df %>% dplyr::filter(VarClass == "None") %>% select(Sample)
      df <- df %>% dplyr::filter(VarClass != "None")
    }
  }
  
  #remove duplicates of gene events within the same sample.
  events_in_data <- as.character()
  
  cat("Preparing input files\n")
  cat("Dim of df pre-process: ", dim(df), "\n")
  df <- remove_duplicates(df)
  cat("Dim of df after duplicate removal: ", dim(df), "\n")
  if (convert){
    df <- convert_varclass(df)
  }
  events_in_data <- c(events_in_data, unique(df[, 3]))
  
  
  cat("Finished preparing input files\n")
  cat("Dim of df after class conversion: ", dim(df), "\n")
  
  if(merge_scnas && !is.null(annotation)){
    alts <- acast(df, Gene ~ Sample)
    cat("There are ", length(df2), "additional data frames to process\n")
    
    for (dframe in df2){
      cat("Preparing additional data frames for input\n")
      colnames(dframe) <- c("Sample", "Gene", "VarClass")
      dframe$Sample <- as.character(dframe$Sample)
      dframe$Gene <- as.character(dframe$Gene)
      dframe$VarClass <- as.character(dframe$VarClass)
      if(convert){
        cat("Dimensions of df to convert : ", dim(dframe), "\n")
        dframe <- convert_varclass(dframe)
      }
      
      if(merge_scnas){
        cat("Dimensions of df to merge : ", dim(dframe), "\n")
        dframe <- remove_duplicates(dframe)
      }
      events_in_data <- c(events_in_data, unique(dframe[, 3]))
      
      alts2 <- acast(dframe, Gene ~ Sample)
      alts <- paste.matrix(alts, alts2)
      cat("Merged Matrix:\n")
      
    }
    
    cat("Finished additional data frames for input\n")
    alterations <- alts
    colnames(annotation) <- c("sample", "class")
    annotation.samples <- annotation$sample
    
    missing.sample <- annotation[which(!annotation.samples%in%colnames(alterations)), ]
    missing.matrix <- matrix(NA, nrow = nrow(alterations), ncol = nrow(missing.sample))
    colnames(missing.matrix) <- missing.sample$sample
    
    alterations <- cbind(alterations, missing.matrix)
    alterations.c <- matrix(as.numeric(!is.na(alterations)), ncol = ncol(alterations))
    colnames(alterations.c) <- colnames(alterations)
    row.names(alterations.c) <- row.names(alterations)
    
  }else if(!is.null(annotation)){
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
    for (dframe in df2){
      colnames(dframe) <- c("Sample", "Gene", "VarClass")
      dframe$Sample <- as.character(dframe$Sample)
      dframe$Gene <- as.character(dframe$Gene)
      dframe$VarClass <- as.character(dframe$VarClass)
      if(convert){
        dframe <- convert_varclass(dframe)
      }
      if(merge_scnas){
        dframe <- remove_duplicates(dframe)
      }
      events_in_data <- c(events_in_data, unique(dframe[, 3]))
      alts2 <- acast(dframe, Gene ~ Sample)
      alts <- paste.matrix(alts, alts2)
    }
    alterations <- alts
    alterations.c <- matrix(as.numeric(!is.na(alterations)), ncol = ncol(alterations))
    colnames(alterations.c) <- colnames(alterations)
    row.names(alterations.c) <- row.names(alterations)
  }else{
    alterations.c <- acast(df, Gene ~ Sample, fun.aggregate = length) # This is the 0 and 1 version of the matrix
    alterations <- acast(df, Gene ~ Sample)
  }
  #convert variant type matrix to numerical values
  for (i in 1:nrow(alterations)){
    for(j in 1:ncol(alterations)){
      altered <- alterations[i, j]
      
      if(!is.na(altered)){ # there is an alteration
        if(grepl("," ,altered)){ # alteration is a mix of two seperated by a comma
          alts <- unlist(stringr::str_split(altered, ",")) # split the alterations
          alterations.c[i, j] <- 0
          for (alt in alts){
            alterations.c[i, j] <- alterations.c[i, j] + alteration_score[[alt]]
          }
          
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
  if(is.null(geneName)){
    geneName <- row.names(alterations)
  }
  
  if(length(setdiff(geneName, row.names(alterations))) != 0){
    genes <- setdiff(geneName, row.names(alterations))
    empty_rows <- matrix(rep(NA, length(genes)*ncol(alterations)), nrow = length(genes))
    empty_rows.c <- matrix(rep(NA, length(genes)*ncol(alterations)), nrow = length(genes))
    row.names(empty_rows) <- genes
    message("missing genes: ", genes, "matrix dim: ", ncol(empty_rows), "-", nrow(empty_rows))
    alterations <- rbind(alterations, empty_rows)
    alterations.c <- rbind(alterations.c, empty_rows)
  }
  alterations.c <- sampleSort(alterations.c, geneOrder = geneName, annotations = annotation, annotation_order = annotation_order)
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
  
  if(!is.null(total_samples)){
    if(total_samples != nsamples){
      diff <- total_samples - nsamples
      mat <- matrix(data = rep(NA, ngenes*diff), ncol = diff, nrow = ngenes)
      colnames(mat) <- paste(rep("MockSample", diff), "_", 1:diff, sep="")
      alterations <- cbind(alterations, mat)
      alterations.c <- cbind(alterations.c, mat)
    }
  }
  if(!is.null(categorical_data)){ 
    if(length(unique(categorical_data[[1]])) != nsamples){
      diff <- length(unique(categorical_data[[1]])) - nsamples
      mat <- matrix(data = rep(NA, ngenes*diff), ncol = diff, nrow = ngenes)
      colnames(mat) <- setdiff(unique(categorical_data[[1]]), colnames(alterations))
      alterations <- cbind(alterations, mat)
      alterations.c <- cbind(alterations.c, mat)
    }
  }
  
  
  ### Set up the matrices that will hold the coordinates for the oncoprints for variety of alterations
  ngenes <- nrow(alterations)
  nsamples <- ncol(alterations)
  numOfOncos <- ngenes*nsamples
  oncoCords.base <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords.base) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  oncoCords <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  oncoCords.scna <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords.scna) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  oncoCords.fusion <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords.fusion) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  oncoCords.borders <- matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
  colnames(oncoCords.borders) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  xpadding <- xpadding
  ypadding <- ypadding
  cnt <- 1;
  
  
  message("nsamples: ", nsamples, " ngenes: ", ngenes, "\n")
  
  mutation_alterations <- c("Mutation" ,  "Missense" , "Nonsense" , "Splicing" ,  "Frameshift" ,  "Promoter" , "InFrame")
  scna_alterations <- c("Amplification" , "Deletion" , "homodel" , "del" , "CNLOH" , "LOH")
  misc_alterations <- c("Present" , "NotPresent" , "NotTested", "Yes" , "No", "Unknown")
  fusion_alterations <- c("Fusion")
  border_alterations <- c("Pathogenic")
  
  barplot_data <- matrix(rep(0, 3*ngenes), nrow = 3)
  colnames(barplot_data) <- row.names(alterations)
  row.names(barplot_data) <- c("Mutation", "SCNA", "Fusion")
  
  #create empty matrix for categorical variables
  emptyCat_data <- matrix(rep(0, 3*length(unique(categorical_data[, 2]))), nrow = 3)
  colnames(emptyCat_data) <- unique(categorical_data[, 2])
  barplot_data <- cbind(emptyCat_data, barplot_data)
  
  #adding data frame for sample total mutations/scna
  mutnum_data <- matrix(rep(0,3*length(unique(df$Sample))), nrow =3)
  colnames(mutnum_data) <- unique(df$Sample)
  row.names(mutnum_data) <- c("Mutation","SCNA", "Fusion")
  
  for ( i in 1:length(colnames(mutnum_data))){
    mutnum_data["Mutation",i] = length(which(df$Sample==colnames(mutnum_data)[i] & df$VarClass %in% mutation_alterations))
    mutnum_data["SCNA",i] = length(which(df$Sample==colnames(mutnum_data)[i] & df$VarClass %in% scna_alterations))
    mutnum_data["Fusion",i] = length(which(df$Sample==colnames(mutnum_data)[i] & df$VarClass %in% fusion_alterations))
  }
    cat("here\n")
  if (merge_scnas){ 
    
    for(i in 1:ngenes) {
      for(j in 1:nsamples) {
        gene <- row.names(alterations)[i]
        altered <- alterations[i, j]
        xleft <- j-1 + xpadding 
        ybottom <- ((ngenes-i+1) -1) + ypadding
        xright <- j - xpadding 
        ytop <- (ngenes-i+1) -ypadding
        oncoCords.base[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
        #browser()
        
        if(!is.na(altered)){ # there is an alteration
          if(grepl("," ,altered)){ # alteration is a mix of multiple seperated by a comma
            alts <- unlist(str_split(altered, ",")) # split the alterations
            for (altered in alts){
              if(altered %in% mutation_alterations) {
                ytop2 <- ytop-0.25
                ybottom2 <- ybottom+0.25
                oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
                barplot_data["Mutation", gene] <- barplot_data["Mutation", gene] + 1 
              }else if( altered %in% scna_alterations ){
                oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
                barplot_data["SCNA", gene] <- barplot_data["SCNA", gene] + 1
              }else if( altered %in% misc_alterations){
                oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
                
              }else if(altered %in% fusion_alterations){
                ytop2 <- ytop-0.1
                ybottom2 <- ybottom+0.1
                oncoCords.fusion[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
                barplot_data["Fusion", gene] <- barplot_data["Fusion", gene] + 1
              }else if(altered %in% border_alterations){
                oncoCords.borders[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
              }
            }
          }else{ # alteration does not have a comma
            if(altered %in% mutation_alterations) {
              ytop2 <- ytop-0.25
              ybottom2 <- ybottom+0.25
              oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
              barplot_data["Mutation", gene] <- barplot_data["Mutation", gene] + 1 
            }else if( altered %in% scna_alterations ){
              oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
              barplot_data["SCNA", gene] <- barplot_data["SCNA", gene] + 1
            }else if(  altered %in% misc_alterations){
              oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
            }else if(altered %in% fusion_alterations){
              ytop2 <- ytop-0.1
              ybottom2 <- ybottom+0.1
              oncoCords.fusion[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
              barplot_data["Fusion", gene] <- barplot_data["Fusion", gene] + 1
            }else if(altered %in% border_alterations){
              oncoCords.borders[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
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
        gene <- row.names(alterations)[i]
        altered <- alterations[i, j]
        xleft <- j-1 + xpadding 
        ybottom <- ((ngenes-i+1) -1) + ypadding
        xright <- j - xpadding 
        ytop <- (ngenes-i+1) -ypadding
        oncoCords.base[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
        #browser()
        if(!is.na(altered)){
          if(altered %in% mutation_alterations) {
            ytop2 <- ytop-0.25
            ybottom2 <- ybottom+0.25
            oncoCords[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
            barplot_data["Mutation", gene] <- barplot_data["Mutation", gene] + 1 
          }else if( altered %in% scna_alterations){
            oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
            barplot_data["SCNA", gene] <- barplot_data["SCNA", gene] + 1
          }else if(  altered %in% misc_alterations){
            oncoCords.scna[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
          }else if(altered %in% fusion_alterations){
            ytop2 <- ytop-0.1
            ybottom2 <- ybottom+0.1
            oncoCords.fusion[cnt, ] <- c(xleft, ybottom2, xright, ytop2, altered)
            barplot_data["Fusion", gene] <- barplot_data["Fusion", gene] + 1
          }else if(altered %in% border_alterations){
            oncoCords.borders[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
          }else{
            oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
          }
        }else{
          oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
        }
        cnt <- cnt+1
      }
    }
  }
  
  cnt <- 1
  oncoCords.catData <- matrix()

            
  if(!is.null(categorical_data)){
    if(length(categorical_data_colors) < length(unique(categorical_data[, 3]))){
      warning("You don't have sufficient amount of colors.\n")
    }
    message("Processing categorical data")
    ystart <- max(as.numeric(oncoCords.base[,4])) + ypadding
    ncategory <- length(unique(categorical_data[[2]]))
    oncoCords.catData <-  matrix( rep(0, nsamples * ncategory * 5), nrow= nsamples * ncategory)
    colnames(oncoCords.catData) <- c("xleft", "ybottom", "xright", "ytop", "altered")
    cat_data <- unname(unlist(unique(categorical_data[, 2])))
    if(!is.null(categorical_data_order)){
     cat_data <- categorical_data_order
    }
    for (j in 1:nsamples){
     for (i  in 1:length(cat_data)){
       sample <- colnames(alterations)[j]
       category <- cat_data[i]
       # idx <- which(categorical_data[[3]] == category & categorical_data[[1]] == sample)
       # altered <- ifelse(categorical_data[[3]][idx], categorical_data[[3]][idx], NA)
       altered <- unname(unlist(categorical_data[categorical_data[, 1] == sample & categorical_data[, 2] == category, ][, 3]))
       xleft <- j-1 + xpadding 
       ybottom <- ystart + i-1 + ypadding
       xright <- j - xpadding 
       ytop <- ystart + i - ypadding
       
       #message("cnt:", cnt, ", altered: ", altered, ", sample: ", sample, ", category: ", category) # ", idx: ", idx, ", sample.idx: ", idx, ", cat.idx: ", idx)  
       oncoCords.catData[cnt, ] <- c(xleft, ybottom, xright, ytop, altered)
       
       cnt <- cnt + 1
     }
    }
    
    oncoCords.base <- rbind(oncoCords.catData, oncoCords.base)
    message("Finished processing categorical data")
  }
  
  
  oncoCords.contData <- matrix()
  ystart <- max(as.numeric(oncoCords.base[,4])) + ypadding
  if(!is.null(continuous_data)){
    oncoCords.contData <-  matrix( rep(0, numOfOncos * 5), nrow=numOfOncos )
    colnames(oncoCords.contData) <- c("xleft", "ybottom", "xright", "ytop", "altered")
  }
  
  gene_prcnt <- list()
  if(!is.null(total_samples)){
    
    labels=list()
    for(i in (1:nrow(alterations))){
      gene = row.names(alterations)[i]
      cnt <- sum(ifelse(test = is.na(alterations[i,]), yes = 0, no = ifelse(test =grepl(",", alterations[i,]), yes = 2, no = 1)))
      prcnt <- cnt/length(alterations.c[i,])*100
      gene_prcnt[gene] <- prcnt
      labels = c(labels, paste(gene, "  ",round(prcnt,1), "%", sep=""))
    }
  }else{
    labels = rownames(alterations)
  }
  
  if (!is.null(categorical_data)){
    labels <- c(rev(unname(unlist(unique(categorical_data[, 2])))), labels)
  }
  barplot_data <- barplot_data[, match(rownames(alterations), colnames(barplot_data))]
  
  
  ## Set up color schema for different classes of alterations
  cnt <- nrow(oncoCords.base)
  colors <- rep(NA, cnt)
  colors.scna <- rep(NA, cnt)
  colors.fusion <- rep(NA, cnt)
  colors.cat <- rep(NA, cnt)
  colors.border <- rep(NA, cnt)
  for (alteration in mutation_alterations){
    colors[ which(oncoCords[, "altered"] == alteration) ] <- onco_colors[[alteration]]
  }
  
  for (alteration in scna_alterations){
    colors.scna[ which(oncoCords.scna[, "altered"] == alteration) ] <- onco_colors[[alteration]]
  }
  
  for (alteration in border_alterations){
    colors.border[ which(oncoCords.borders[, "altered"] == alteration) ] <- onco_colors[[alteration]]
  }
  c48 <- c("#1d915c","#5395b4","#964a48","#2e3b42","#b14e72", "#402630","#f1592a","#81aa90","#f79a70","#b5ddc2","#8fcc8b","#9f1f63","#865444", "#a7a9ac","#d0e088","#7c885c","#d22628","#343822","#231f20","#f5ee31","#a99fce","#54525e","#b0accc","#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")
  if(!is.null(categorical_data) & is.null(categorical_data_colors)){
    
    for (i in 1:length(unname(unlist(unique(categorical_data[, 3])))) ) {
      cat(i, "\n")
      categorical_data_colors[unname(unlist(unique(categorical_data[, 3])))[i]] <- c48[i] 
      message(categorical_data_colors)
    }
    for (alteration in unname(unlist(unique(categorical_data[, 3]))) ){
      colors.cat[ which(oncoCords.catData[, "altered"] == alteration) ] <- unname(unlist(categorical_data_colors[[alteration]] ) )
    }
  }
  
  if(!is.null(categorical_data) & !is.null(categorical_data_colors)){
    
    for (alteration in unname(unlist(unique(categorical_data[, 3]))) ){
      colors.cat[ which(oncoCords.catData[, "altered"] == alteration) ] <- unname(unlist(categorical_data_colors[[alteration]]))
    }
  }
  
  colors.fusion[ which(oncoCords.fusion[, "altered"] == "Fusion") ] <- onco_colors[["Fusion"]]
  
  
  ngenes <- nrow(alterations) + length(unique(categorical_data[[2]]))
  cat(ngenes)
  #change the 
  def.par <- par(no.readonly = TRUE)
  leftmargin = 1/nsamples*500 
  if(leftmargin > 20){leftmargin <- 20}
  bottommargin = 1/ngenes*500
  if(bottommargin > 5){bottommargin <- 1}

  #recommend output to pdf with 10x5" dimensions i.e - pdf("test.pdf", width = 10, height = 5, paper="special")
  split.screen(rbind(c(0.01, 0.85, 0.95, 0.99), 
                      c(0.01,0.85,0.15, 0.95),
                      c(0.01, 0.85, 0.01, 0.15),
                      c(0.85, 0.99, 0.15, 0.99),
                      c(0.86, 0.99, 0.01, 0.15)))
  
  # split.screen(rbind(c(0.01,0.8525,0.78, 0.99), #top sample bar
  #                    c(0.01,0.85,0.15, 0.8), #onco
  #                    c(0.1, 0.85, 0.01, 0.2), #legend
  #                    c(0.83, 0.99, 0.15, 0.8), # gene bar
  #                    c(0.84, 0.99, 0.01, 0.15))) #bar legend
  # 
  # screen(1)
  # par(mar=c(0,8.87,0.2,0))
  #barplot(mutnum_data, horiz = F, axisnames = F,col= c("#2986E2", "#F26529", "#619744"),cex.names = 0.5, cex.axis = 0.5, yaxs = "i", border=NA)
  
  if(!is.null(annotation)){
    screen(1)
    par(mar=c(0,10,0,0), mgp=c(3, 0.7, 0))
    plot(c(0, nsamples), c(0,1), type="n", main="", xlab="Samples", xaxt="n", ylab="", yaxt="n", frame.plot = F)
    counts <- data.frame(table(annotation$class))
    xleft <- 0
    axis.points <- list()
    subtype.labels <-list()
    for (i in 1:length(annotation_order)){
      subtype.labels <- c(subtype.labels, annotation_order[i])
      xright <- counts[which(counts$Var1 == annotation_order[i]), ]['Freq']
      if(!is.null(annotation_colors)){
        rect(xleft, 0, xleft + xright$Freq, 1, col = annotation_colors[i], border = "white")
       }else{
        rect(xleft, 0, xleft + xright$Freq, 1, col = c48[i], border = "white")
       }
      axis.points <- c(axis.points, (xleft + xleft+xright$Freq)/2)
      xleft <- xleft + xright$Freq
    }
    text(x=axis.points, y = 0.5, labels = subtype.labels)
  }
  
  
  screen(2)
  #bottom, left, top and right in lines of text
  par(mar=c(0,10,0,0), mgp=c(3, 0.7, 0))
  plot(c(0, nsamples), c(0, ngenes), type="n", main="", xlab="Samples", xaxt="n", ylab="", yaxt="n", frame.plot = F, xaxs = "i");
  rect(oncoCords.base[, "xleft"], oncoCords.base[, "ybottom"],oncoCords.base[, "xright"], oncoCords.base[, "ytop"], col="#DCD9D3", border=NA);
  rect(oncoCords.scna[, "xleft"], oncoCords.scna[, "ybottom"],oncoCords.scna[, "xright"], oncoCords.scna[, "ytop"], col=colors.scna, border=NA);
  rect(oncoCords.fusion[, "xleft"], oncoCords.fusion[, "ybottom"],oncoCords.fusion[, "xright"], oncoCords.fusion[, "ytop"], col=colors.fusion, border=NA);
  rect(oncoCords[, "xleft"], oncoCords[, "ybottom"],oncoCords[, "xright"], oncoCords[, "ytop"], col=colors, border=NA)
  rect(oncoCords.borders[, "xleft"], oncoCords.borders[, "ybottom"],oncoCords.borders[, "xright"], oncoCords.borders[, "ytop"], col=NA, border="black")
  if(!is.null(categorical_data)){
    rect(oncoCords.catData[, "xleft"], oncoCords.catData[, "ybottom"],oncoCords.catData[, "xright"], oncoCords.catData[, "ytop"], col=colors.cat, border=NA);
  }
  axis(2, at=(length(labels):1)-.5, labels=labels, las=2, lwd = 0, cex=0.8, cex.axis=0.7);
  #printing samples or not
  if(printSamples){
    text((1:nsamples)-.5, par("usr")[2]+.3,srt=45, adj = 1,  labels = colnames(alterations), xpd=T)
  }
  #add legend
  screen(3)
  par(mar=c(0,0,0,0))
  
  events_in_data <- unlist(events_in_data)
  if(length(mutation_alterations[mutation_alterations %in% events_in_data]) > 0){
    legend(x = 0, y = 1, names(onco_colors[names(onco_colors) %in% mutation_alterations[mutation_alterations %in% events_in_data]]), fill = unlist(onco_colors[names(onco_colors) %in% mutation_alterations[mutation_alterations %in% events_in_data]]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Mutations")
  }
  if (length(scna_alterations[scna_alterations %in% events_in_data]) > 0){
    legend(x = 0.15, y = 1, names(onco_colors[names(onco_colors) %in% scna_alterations[scna_alterations %in% events_in_data]]), fill = unlist(onco_colors[names(onco_colors) %in% scna_alterations[scna_alterations %in% events_in_data]]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "SCNA")
  }
  if(length(misc_alterations[misc_alterations %in% events_in_data])> 0){
    legend(x = 0.3, y = 1, names(onco_colors[names(onco_colors) %in% misc_alterations[misc_alterations %in% events_in_data]]), fill = unlist(onco_colors[names(onco_colors) %in% misc_alterations[misc_alterations %in% events_in_data]]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Misc Fetatures")
  }
  if(length(fusion_alterations[fusion_alterations %in% events_in_data]) >0 ){
    legend(x = 0.45, y=1, names(onco_colors[names(onco_colors) %in% fusion_alterations[fusion_alterations %in% events_in_data]]), fill = unlist(onco_colors[names(onco_colors) %in% fusion_alterations[fusion_alterations %in% events_in_data]]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Fusions")
  }
  if(!is.null(categorical_data)){
    if(length(categorical_data_colors) <= 7){
      legend(x = 0.6, y=1, names(categorical_data_colors), fill = unlist(categorical_data_colors), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Categorical Data")
    }else if(length(categorical_data_colors) > 5 & length(categorical_data_colors) <= 10){
      legend(x = 0.6, y=1, names(categorical_data_colors[1:5]), fill = unlist(categorical_data_colors[1:5]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Categorical Data")
      legend(x = 0.75, y=1, names(categorical_data_colors[6:10]), fill = unlist(categorical_data_colors[6:10]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "")
    }
    else if(length(categorical_data_colors) > 10 ){
      legend(x = 0.6, y=1, names(categorical_data_colors[1:5]), fill = unlist(categorical_data_colors[1:5]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "Categorical Data")
      legend(x = 0.75, y=1, names(categorical_data_colors[6:10]), fill = unlist(categorical_data_colors[6:10]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "")
      legend(x = 0.9, y=1, names(categorical_data_colors[11:length(categorical_data_colors)]), fill = unlist(categorical_data_colors[11:length(categorical_data_colors)]), horiz = F, border = F, cex = 0.7, bty = "n" , title = "")
    }
  }
  
  screen(4)
  par(mar=c(1.75,0.1,2.4,1))
  barplot(barplot_data[, rev(colnames(barplot_data))], horiz = T, axisnames = F, col= c("#21600A", "#602C0A", "#619744"), border = "white", xlab = paste("Total Samples = ", total_samples, sep=""), cex.names = 0.5, cex.axis = 0.5, yaxs ="i")
  
  screen(5)
  par(mar=c(0,0,0,0))
  legend(x = 0, y = 1,c("Mutations", "SCNA", "Fusion"), fill = c("#21600A", "#602C0A", "#619744"), bty="n", cex=0.75)
  close.screen(all.screens = TRUE)

  par(def.par) 
  res <- list()
  res$sortedMatrix <- alterations
  res$sampleOrder <- colnames(alterations)
  res$geneOrder <- rownames(alterations)
  res$gene_prcnt <- gene_prcnt
  res$alterations.value <- alterations.c
  res$barplot_data <- barplot_data
  res$sample_barplot_data <- mutnum_data
  res$oncoCords <- oncoCords.base
  res$catData <-  oncoCords.catData
  res$events <- events_in_data
  res
}

