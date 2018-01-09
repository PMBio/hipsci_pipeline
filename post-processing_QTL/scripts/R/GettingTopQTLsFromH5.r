library(rhdf5)
library(qvalue)
library(dplyr)

filterFile <-  ""
filterFeatures <- NULL
if(!(filterFile=="" || is.na(filterFile))){
 filterFeatures = read.delim(filterFile,as.is=T,header=F)[,1]
}

##Settings
baseFolder <- "./"

writeGlobalSig = T
writeGlobalSigTop = T
writeFeatureSig = T
writeNominalSig = F
topResultBased = T
threshold = 0.05
threshold2 = 1.0
multipleTestingGlobal = "ST"
# multipleTestingGlobal = "BF"
#################
##Read files.
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
snpAnnotation <- NULL
filesToRead <- list.files("./",pattern = ".h5",full.names = T)
if(file.exists("./qtl_results_all.h5")){
  tmp <- h5dump(file = "./qtl_results_all.h5")
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
  snpAnnotation <- read.delim("./snp_metadata_all.txt",as.is=T)
} else {
  for(i in filesToRead){
    
    baseName <- gsub(gsub(i,pattern = ".h5",replacement = ""),pattern = "./qtl_results",replacement = "")
    
    if(length(list.files("./",pattern = paste(baseName,"\\.",sep="")))==3 | length(list.files("./",pattern = paste(baseName,"\\.",sep="")))==4){
      tmp <- h5dump(file = i)
      snpAnnotationTmp <- read.delim(paste("./snp_metadata",baseName,".txt",sep=""),as.is=T)
      snpAnnotation <- rbind(snpAnnotation,snpAnnotationTmp)
    } else {
      print(paste("Skipping",i,"because not necessary data is available or to many files with the same name."))
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
}
snpAnnotation <- unique(snpAnnotation)
#colnames(results)[which(colnames(results)=="corr_p_value")] <- "feature_corr_p_value"
if(length(which(is.na(results$empirical_feature_p_value)))!=0){
  results <- results[-which(is.na(results$empirical_feature_p_value)),]
}

if(!is.null(filterFeatures)){
  results <- results[which(results$feature%in%filterFeatures),]
}

##Add assessed allele
results["Assesed Allele"] <- ""
snpAnnotation <- snpAnnotation[which(snpAnnotation$snp_id %in% results$snp_id),]
for(i in 1:nrow(snpAnnotation)){
  results$`Assesed Allele`[which(results$snp_id==snpAnnotation$snp_id[i])] <- snpAnnotation$assessed_allele[i]
}

results <- results[order(results$empirical_feature_p_value, results$p_value,decreasing = F),]

if(topResultBased){
  resultsFull <- results
  results <- results[which(!duplicated(results$feature)),]
}

##Multiple testing
if(multipleTestingGlobal=="ST"){
  results["global_corr_p_value"] <- qvalue(results$empirical_feature_p_value)$qvalues
} else if (multipleTestingGlobal=="BF"){
  results["global_corr_p_value"] <- results$empirical_feature_p_value*observedFeatures
  results$global_corr_p_value[results$global_corr_p_value>1]<-1
}

results <- results[order(results$global_corr_p_value,results$empirical_feature_p_value, results$p_value,decreasing = F),]

if(writeGlobalSigTop){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}


if(writeGlobalSigTop & topResultBased){
  write.table(paste(baseFolder,"top_results_global_level_",threshold,".txt",sep=""),x = results[intersect(which(results$global_corr_p_value<threshold),which(!duplicated(results$feature))),],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig & !topResultBased){
  write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = results[results$global_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeFeatureSig & !topResultBased){
  write.table(paste(baseFolder,"results_gene_level_",threshold,".txt",sep=""),x = results[results$feature_corr_p_value<threshold,],sep="\t",row.names=F,quote=F)
}

if(writeNominalSig){
  write.table(paste(baseFolder,"results_nominal_level_",threshold2,".txt",sep=""),x = resultsFull[resultsFull$p_value<threshold2,],sep="\t",row.names=F,quote=F)
}

if(writeGlobalSig & topResultBased){
  if(length(which(results$global_corr_p_value<threshold))>0){
    minFeatureCorrected <- max(results$empirical_feature_p_value[which(results$global_corr_p_value<threshold)])
    write.table(paste(baseFolder,"results_global_level_",threshold,".txt",sep=""),x = resultsFull[resultsFull$empirical_feature_p_value<minFeatureCorrected,],sep="\t",row.names=F,quote=F)    
  }
}
