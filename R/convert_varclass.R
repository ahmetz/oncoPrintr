#' convert_varclass
#' 
#' @description if conversion of variant classes is selected, this function will do that
#'
#' @param df dataframe of genomic alterations
#'
#' @return data frame with variant classes converted
#' @details variant classes are converted based on this schema:\cr
#' splicing -> Splicing\cr
#' stopgain or stoploss -> Nonsense\cr
#' nonsynonymous -> Missense\cr
#' frameshift_deletion or insertion -> Frameshift\cr
#' nonframeshift_deletion or insertion -> InFrame\cr
#' upstream - Promoter\cr
#' IntragenicDeletion -> Frameshift\cr
#' @export 
#'
#' @examples TODO
convert_varclass <- function(df, class_lookup= c('splicing'="Splicing", 'stop'="Nonsense", 'nonsynonymous'="Missense", 'nonframeshift'="InFrame", 'frameshift'="Frameshift", 'upstream'="Promoter", 'IntragenicDeletion'="Deletion", 'Translation_Start_Site'="Frameshift")){
  cat("Converting variant classes\n")
  if (sum(grepl("splicing", df$VarClass) > 0)){
    df[grep("splicing", df$VarClass), ]["VarClass"] <- class_lookup['splicing']
  }
  if(sum(grepl("stop", df$VarClass) > 0)){
    df[grep("stop", df$VarClass), ]["VarClass"] <- class_lookup['stop']
  }
  if(sum(grepl("nonsynonymous", df$VarClass) > 0)){
    df[grep("nonsynonymous", df$VarClass), ]["VarClass"] <- class_lookup['nonsynonymous']
  }
  if(sum(grepl("nonframeshift", df$VarClass) > 0)){
    df[grep("nonframeshift", df$VarClass), ]["VarClass"] <- class_lookup['nonframeshift']
  }
  if(sum(grepl("frameshift", df$VarClass) > 0)){
    df[grep("frameshift", df$VarClass), ]["VarClass"] <- class_lookup['frameshift']
  }
  if(sum(grepl("upstream", df$VarClass) > 0)){
    df[grep("upstream", df$VarClass), ]["VarClass"] <- class_lookup['upstream']
  }
  if(sum(grepl("IntragenicDeletion", df$VarClass) > 0)){
    df[grep("IntragenicDeletion", df$VarClass), ]["VarClass"] <- class_lookup['IntragenicDeletion']
  }
  if(sum(grepl("silent", df$VarClass) > 0)){
    df[grep("silent_SNV", df$VarClass), ]["VarClass"] <- class_lookup['Silent']
  }
  if(sum(grepl("Translation_Start_Site", df$VarClass) > 0)){
    df[grep("Translation_Start_Site", df$VarClass), ]["VarClass"] <- class_lookup['Translation_Start_Site']
  }
  
  cat("Finished converting variant classes\n")
  return(df)
}
