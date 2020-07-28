localtstat = function(center, betamap, radius, nn, adjust.p.method = "BH"){
  ## This function computes the local t-statistics moderated by limma and psuedo t-statistics and p-values 
  ## Input:
  ## center.local: Integer: Index for target voxel
  ## betamap:      Activation map
  ## radius:       Radius for finding exact neighbor indices that dist(targetInd - queryInd) <= radius
  ## nn:           Results from findNearestNeighbor function
  ## adjust.p.method: Refer to the topTable function in limma package. 
  
  ## Output:
  ##tstat   :   local t-statistics modorated by limma
  ##pval    :   unadjusted p-value
  ##adj.pval:   adjusted  p-value
  ##psuedo.t :  psuedo t-statistics
  require(limma)

  idx = nn$index[center, nn$distance[center,]<= radius]
  fit = lmFit(betamap[idx,])
  eb.fit = eBayes(fit)
  tstat = eb.fit$t[1]
  pval = eb.fit$p.value[1]
  
  local.eb.fit = lapply(eb.fit, "[", 1)
  
  adj.pval = topTable(eb.fit, sort.by = 'none', adjust.method = adjust.p.method, number = length(idx))$adj.P.Val[1]
  psuedo.t = eb.fit$coefficients[1]/sqrt(mean((eb.fit$sigma)^2))/eb.fit$stdev.unscaled[1]
  
  return(c(tstat, pval, adj.pval, psuedo.t))
}

