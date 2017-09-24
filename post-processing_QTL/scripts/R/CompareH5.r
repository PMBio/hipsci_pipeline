library(rhdf5)
library(qvalue)
library(dplyr)

#Set 1
baseFolder1 <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Output/HipSci/WithoutCorrection/Gene_Mapping/"
subFolderBase1 <- "OutGeneMapping.chr."

#Set 2
baseFolder2 <- "C:/OnlineFolders/BitSync/Current_Work/EBI_HipSci/QTL_Output/HipSci/WithoutCorrection/Gene_Mapping/EBI_Cluster/"
subFolderBase2 <- "OutGeneMapping.chr."


chrSpecific = T
range <- 21:22
#range <- c(1:17,19:22)
writeGlobalSig = T
writeGlobalSigTop = T
writeFeatureSig = T
threshold = 0.05
peerCorrected = F
# multipleTestingGlobal = "ST"
multipleTestingGlobal = "BF"
#################


##N features and SNPs is different!!!!!!
###WHUTTT?

observedFeatures1 <- 0
observedFeatures2 <- 0
results1 <- NULL
results2 <- NULL
for(i in range){
  folderName1 = paste(baseFolder1, subFolderBase1,i,sep="")
  folderName2 = paste(baseFolder2, subFolderBase2,i,sep="")
  if(peerCorrected){
    folderName1=paste(folderName1,".PeerCorrected",sep="")
    folderName2=paste(folderName2,".PeerCorrected",sep="")
  }
  
  if(length(list.files(folderName1))>2 && length(list.files(folderName2))>2){
    
    
    if(chrSpecific){
      tmp1 <- h5dump(file = paste(folderName1,"/qtl_results_",i,".h5",sep=""))
      tmp2 <- h5dump(file = paste(folderName2,"/qtl_results_",i,".h5",sep=""))
    } else {
      tmp1 <- h5dump(file = paste(folderName1,"/qtl_results_all.h5",sep=""))
      tmp2 <- h5dump(file = paste(folderName2,"/qtl_results_all.h5",sep=""))
    }
  }
  if(length(tmp1)>0 && length(tmp2)>0){
    for (j in names(tmp1)) tmp1[[j]][["feature"]] <- j
    for (j in names(tmp2)) tmp2[[j]][["feature"]] <- j
    observedFeatures1 = observedFeatures1+length(tmp1)
    observedFeatures2 = observedFeatures2+length(tmp2)
    df1 <- bind_rows(tmp1)
    df2 <- bind_rows(tmp2)
    if(nrow(df1)>0){
      results1 = rbind(results1,df1)  
    }
    if(nrow(df2)>0){
      results2 = rbind(results2,df2)  
    }
  }
}

colnames(results1)[which(colnames(results1)=="corr_p_value")] <- "feature_corr_p_value"
colnames(results2)[which(colnames(results2)=="corr_p_value")] <- "feature_corr_p_value"
if(length(which(is.na(results1$feature_corr_p_value)))!=0){
 results <- results[-which(is.na(results1$feature_corr_p_value)),]
}

if(length(which(is.na(results2$feature_corr_p_value)))!=0){
  results <- results2[-which(is.na(results2$feature_corr_p_value)),]
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

