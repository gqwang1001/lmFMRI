#' @title Find the approximate set of neighbors of each voxels. 
#' @description Find the approximate set of neighhorhood of each voxels, given the distance radius.The number of elements in each neighborhood is (2*radius)^3.
#'
#' @param mask Input Masked image which the masked voxels are labeled with zeros
#' @param radius Indices that dist(targetInd - queryInd) <= radius

#' @return The result from \code{\link{vcgKDtree}}
#' @export

findNearestNeighbors = function(maskImg, radius = 1){
  require(Rvcg)
  ind.mask <- which(maskImg>0, arr.ind = T)
  nn = Rvcg::vcgKDtree(target = ind.mask, query = ind.mask, k = radius^3*8, threads = 10)
  return(nn)
}