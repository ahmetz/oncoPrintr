#' sampleSort
#'
#' @param M matrix of alterations.
#' @param geneName If certain genes needs to be on the top, this value should provide them in a list
#' @param annotations this holds the per sample annotation values. the matrix will be split into n sub matrices and sorted then merged at the end.
#' @param annotation_order (Required) specify the order of the annotation. If not present, R will sort it alphabetically
#'
#' @return matrix with samples and genes sorted
#' @export
#'
#' @examples TODO2
sampleSort <- function(M, geneOrder = geneOrder, annotations = NULL, annotation_order = NULL) {

  M <- M[geneOrder, ]
  if(!is.null(annotations)){
    colnames(annotations) <- c("sample", "class")
    classes <- unique(annotations$class)
    M2 <- matrix(nrow=nrow(M))
    colnames(M2) <- "Drop"
    rownames(M2) <- rownames(M)
    M2 <- M2[ , drop=F]
    for(class in annotation_order){
      if(! class %in% annotations){
        message("An annotation supplied in the annotation order does not exist in the annotation matrix! ")
        stop()
      }
      samples <- annotations[which(annotations[, 2] == class), ][, 1]
      sub_mat <- M[, which(colnames(M)%in%samples$sample), drop=FALSE]
      if (ncol(sub_mat) == 1){
        M2 <- cbind(M2, sub_mat[, 1,drop=F])
      }else{
        sub_mat.t <- t(sub_mat)
        sub_mat.t <-  sub_mat.t[do.call(order, as.data.frame(sub_mat.t)), ]
        sub_mat <- t(sub_mat.t)
        sub_mat <- sub_mat[, ncol(sub_mat):1]
        M2 <- cbind(M2, sub_mat)
      }
    }
    M <- M2[, 2:ncol(M2)]

    return(M)
  }else{
  M.t <- t(M)
  M.t <-  M.t[do.call(order, as.data.frame(M.t)), ]
  M <- t(M.t)
  M <- M[, ncol(M):1]
  return(M);
  }
}
