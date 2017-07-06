library(rhdf5)
library(qvalue)
library(dplyr)
# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/Gene_Mapping/"
# subFolderBase <- "OutGeneMapping.chr."
# chrSpecific = T
# range <- 1:22
# writeGlobalSig = T
# writeGeneSig = T
# threshold = 0.05
# multipleTestingGlobal = "ST"
#################

baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/On_PEER_Factors/"
subFolderBase <- "OutMappingPeerFactors."
chrSpecific = F
range <- 2:101
writeGlobalSig = T
writeGeneSig = F
threshold = 0.1
multipleTestingGlobal = "BF"


####
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
for(i in range){
  if(chrSpecific){
    tmp <- h5dump(file = paste(subFolderBase,i,"/qtl_results_",i,".h5",sep=""),)
  } else {
    tmp <- h5dump(file = paste(subFolderBase,i,"/qtl_results_all.h5",sep=""),)
  }
  for (j in names(tmp)) tmp[[j]][["feature"]] <- j
  observedFeatures = observedFeatures+length(tmp)
  df <- bind_rows(tmp)
  if(multipleTestingGlobal=="BF"){
    df <- df[df$corr_p_value<threshold,]
  }
  if(nrow(df)>0){
    results = rbind(results,df)  
  }
}

colnames(results)[which(colnames(results)=="corr_p_value")] <- "feature_corr_p_value"
if(multipleTestingGlobal=="ST"){
  results["global_corr_p_value"] <- qvalue(results$feature_corr_p_value)$qvalues  
} else if (multipleTestingGlobal=="BF"){
  results["global_corr_p_value"] <- results$feature_corr_p_value*observedFeatures
  results$global_corr_p_value[results$global_corr_p_value>1]<-1
}


results <- results[order(results$global_corr_p_value,decreasing = F),]

if(writeGlobalSig){
  write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = results[results$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeGeneSig){
  write.table(paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),x = results[results$gene_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}
