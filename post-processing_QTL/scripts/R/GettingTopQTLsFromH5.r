library(rhdf5)
library(qvalue)
library(dplyr)
# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithoutCorrection/Gene_Mapping/"
# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithCorrection/GeneLevel/"
# subFolderBase <- "OutGeneMapping.chr."

# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithoutCorrection/Transcript_Mapping/"
#baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithCorrection/TranscriptLevel/"
#subFolderBase <- "OutTranscriptMapping.chr."

# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithoutCorrection/Exon_Mapping/"
# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithCorrection/ExonLevel/"
# subFolderBase <- "OutExonMapping.chr."

# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithoutCorrection/Splicing_Mapping/"
# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/WithCorrection/Splicing/"
# subFolderBase <- "OutSplicingMapping.chr."
# subFolderBase <- "OutSplicingMapping.LowCovBlanked.chr."
# subFolderBase <- "OutSplicingMapping_qq.chr."

# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/GENES_IPS_QTL/GeneLevel_Effects_Corrected/"
# subFolderBase <- ""

baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/TRANS/"
subFolderBase <- "OutGeneMapping.Trans.ExamplarGenes.GWAS.chr."

chrSpecific = T
range <- 1:22
#range <- c(1:17,19:22)
writeGlobalSig = T
writeGlobalSigTop = T
writeFeatureSig = T
threshold = 0.05
peerCorrected = T
multipleTestingGlobal = "ST"
# multipleTestingGlobal = "BF"
#################

# baseFolder <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Effects_On_PEER_Factors/On_PEER_Factors/"
# subFolderBase <- "OutMappingPeerFactors."
# chrSpecific = F
# range <- 2:101
# writeGlobalSig = T
# writeGeneSig = F
# threshold = 0.1
# multipleTestingGlobal = "BF"


####
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
for(i in range){
  folderName = paste(subFolderBase,i,sep="")
  if(peerCorrected){
    folderName=paste(folderName,".PeerCorrected",sep="")
  }
  if(length(list.files(folderName))>2){
    if(chrSpecific){
      tmp <- h5dump(file = paste(folderName,"/qtl_results_",i,".h5",sep=""),)
    } else {
      tmp <- h5dump(file = paste(folderName,"/qtl_results_all.h5",sep=""),)
    }
  }
  if(length(tmp)>0){
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
}
colnames(results)[which(colnames(results)=="corr_p_value")] <- "feature_corr_p_value"
if(!length(is.na(results$feature_corr_p_value))==0){
 results <- results[-which(is.na(results$feature_corr_p_value)),]
}


if(multipleTestingGlobal=="ST"){
  results["global_corr_p_value"] <- qvalue(results$feature_corr_p_value)$qvalues
} else if (multipleTestingGlobal=="BF"){
  results["global_corr_p_value"] <- results$feature_corr_p_value*observedFeatures
  results$global_corr_p_value[results$global_corr_p_value>1]<-1
}


results <- results[order(results$global_corr_p_value,decreasing = F),]

if(writeGlobalSigTop){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig){
  write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = results[results$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeFeatureSig){
  write.table(paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),x = results[results$feature_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

