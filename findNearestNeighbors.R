findNearestNeighbors = function(maskImg, radius = 1){
  ## This function aims to find the approximate set of neighbors of each voxels. 
  ## The number of elements in each neighborhood is (2*radius)^3.
  
  ## maskImg: Input Masked image which masked voxels are labeled with zeros
  ## radius:  Indices that dist(targetInd - queryInd) <= radius
  require(Rvcg)
  
  ind.mask <- which(maskImg>0, arr.ind = T)
  length.mask <- nrow(ind.mask)
  
  nn = Rvcg::vcgKDtree(target = ind.mask, query = ind.mask, k = radius^3*8, threads = 10)
  return(nn)
}