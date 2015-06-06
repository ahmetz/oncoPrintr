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
convert_varclass <- function(df){
  if (sum(grepl("splicing", df$VarClass) > 0)){
    df[grep("splicing", df$VarClass), ]["VarClass"] <- "Splicing"
  }
  if(sum(grepl("stop", df$VarClass) > 0)){
    df[grep("stop", df$VarClass), ]["VarClass"] <- "Nonsense"
  }
  if(sum(grepl("nonsynonymous", df$VarClass) > 0)){
    df[grep("nonsynonymous", df$VarClass), ]["VarClass"] <- "Missense"
  }
  if(sum(grepl("frameshift", df$VarClass) > 0)){
    df[grep("frameshift", df$VarClass), ]["VarClass"] <- "Frameshift"
  }
  if(sum(grepl("nonframeshift", df$VarClass) > 0)){
    df[grep("nonframeshift", df$VarClass), ]["VarClass"] <- "InFrame"
  }
  if(sum(grepl("upstream", df$VarClass) > 0)){
    df[grep("upstream", df$VarClass), ]["VarClass"] <- "Promoter"
  }
  if(sum(grepl("IntragenicDeletion", df$VarClass) > 0)){
    df[grep("IntragenicDeletion", df$VarClass), ]["VarClass"] <- "Frameshift"
  }
  
  return(df)
}
