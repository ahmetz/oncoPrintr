#' paste.matrix 
#'
#' @description takes two genomic evetn matrices as input and merges them. the matrices do not have the have the same dimensions. If a gene and a sample have two different alterations, they will be merged seperated with a comma
#'
#' @param M1 matrix of genomic alterations
#' @param M2 matrix of genomic alterations
#'
#' @return merged matrix of M1 and M2
#' @export
#'
#' @examples TODO
paste.matrix <- function(M1, M2){
  rownames1 <- row.names(M1)
  rownames2 <- row.names(M2)
  colnames1 <- colnames(M1)
  colnames2 <- colnames(M2)
  row_names_all <- unique(c(rownames1, rownames2))
  col_names_all <- unique(c(colnames1, colnames2))
  
  n = length(row_names_all)
  m = length(col_names_all)
  ret <- matrix(rep(NA, n *m), ncol = m)
  colnames(ret) <- col_names_all
  row.names(ret) <- row_names_all
  
  for(row in row_names_all){
    for(col in col_names_all){
      x = NA
      y = NA
      if(row %in% rownames1 && col %in% colnames1){
        x = M1[row, col]
      }
      if(row %in% rownames2 && col %in% colnames2){
        y = M2[row, col]
      }
      
      if(is.na(x) && is.na(y)){
        ret[row, col] <- NA
      }else if(is.na(x) && !is.na(y)){
        ret[row, col] <- y
      }else if(!is.na(x) && is.na(y)){
        ret[row, col]<- x
      }else if(!is.na(x) && !is.na(y)){
        ret[row, col] <- paste(x, y, sep=",")
      }
    }
  }
  ret <- as.matrix(ret)
  colnames(ret) <- col_names_all
  row.names(ret) <- row_names_all
  return(ret)
}