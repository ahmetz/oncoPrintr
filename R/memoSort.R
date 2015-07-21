#' memoSort
#'
#' @param M matrix of alterations.
#' @param geneName If certain genes needs to be on the top, this value should provide them in a list
#' @param annotations this holds the per sample annotation values. the matrix will be split into n sub matrices and sorted then merged at the end.
#' @param annotation_order specify the order of the annotation. If not present, R will sort it alphabetically
#'
#' @return matrix with samples and genes sorted
#' @export
#'
#' @examples TODO
memoSort <- function(M, geneName = NA, annotations = NA, annotation_order = NA) {
  M <- M[geneName, ]
#   scoreCol <- function(x) {
#     score <- 0;
#     for(i in 1:length(x)) {
#       if(x[i]) {
#         score <- score + 5000^(length(x)-i)
#       }
#     }
#     #score <- score + sum(x*(length(x):1))
#     return(score);
#   }
#   
  # try giving the gene order manually. based on number of samples that has that event...
  #if(!is.na(geneName)) {
  #  a <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)
  #  geneOrder <- a[[2]][-which(names(a$x) %in% geneName) ]
  #  for (gene in rev(geneName)){
  #    geneOrder <- c(which(rownames(M) == gene), geneOrder)
  #  }
  #  #geneOrder <- c(which(rownames(M)%in%geneName), geneOrder)
  #}else{
  #  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix
  #}
  
  if(!is.na(annotations)){
    colnames(annotations) <- c("sample", "class")
    classes <- unique(annotations$class)
    M2 <- matrix(nrow=nrow(M))
    colnames(M2) <- "Drop"
    rownames(M2) <- rownames(M)
    M2 <- M2[geneOrder, , drop=F]
    for(class in annotation_order){
      samples <- annotations[which(annotations[, 2] == class), ][, 1]
      sub_mat <- M[, which(colnames(M)%in%samples$sample), drop=FALSE]
      
      if (ncol(sub_mat) == 1){
        M2 <- cbind(M2, sub_mat[, 1,drop=F])
      }else{
        
        #scores <- apply(sub_mat[geneOrder, ], 2, scoreCol)
        #sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix
        sub_mat.t <- t(sub_mat)
        sub_mat.t <-  sub_mat.t[do.call(order, as.data.frame(sub_mat.t)), ]
        sub_mat <- t(sub_mat.t)
        sub_mat <- sub_mat[, ncol(sub_mat):1]
        M2 <- cbind(M2, sub_mat)
      }
      #cat(class, ", ", sampleOrder, "\n")
    }
    M <- M2[, 2:ncol(M2)]
    #print(M)
    return(M)
  }else{
    #scores <- apply(M[geneOrder, ], 2, scoreCol);
    #sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix
  }
  M.t <- t(M)
  M.t <-  M.t[do.call(order, as.data.frame(M.t)), ]
  M <- t(M.t)
  M <- M[, ncol(M):1]
  return(M);
  
}
