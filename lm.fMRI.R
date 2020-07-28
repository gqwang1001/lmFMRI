lm.fMRI = function(betamap, mask, adj.p.method = "BH"){
  
  ## betamap: V-by-N matrix of activative maps. V: number of voxels, N: number of subjects
  ## mask: V-by-1 vector mask vector. Masked voxels have the value of zeros. 
  ## adjust.p.method: Refer to the topTable function in limma package. 
  
  require(limma)
  betamap.mask = betamap[mask>0, ]
  
  fit = lmFit(betamap)
  eb.fit = eBayes(fit)
  
  tvec.limma = eb.fit$t 
  pval.limma = eb.fit$p.value
  adj.pval.limma = topTable(eb.fit, 
                            sort.by = 'none', 
                            adjust.method = adj.p.method, 
                            number = nrow(betamap))$adj.P.Val
  
  # Restore the 3D image
  tmap.limma = mask
  tmap.limma[mask>0] = tvec.limma
  
  pvalmap.limma = mask
  pvalmap.limma[mask>0] = pval.limma
  
  ## adjust pvalue
  adj.pvalmap.limma = mask
  adj.pvalmap.limma[mask>0] = adj.pval.limma
  
  return(list(tmap.limma,
              pvalmap.limma,
              adj.pval.limma))
}