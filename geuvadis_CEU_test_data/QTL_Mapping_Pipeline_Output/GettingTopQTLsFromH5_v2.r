library(rhdf5)
library(qvalue)
library(dplyr)

##Settings
baseFolder <- "C:/Users/mjbonder/Documents/GitHub/hipsci_pipeline/geuvadis_CEU_test_data/QTL_Mapping_Pipeline_Output/"

multipleTestingGlobal = "ST"

#################
##Read files.
setwd(baseFolder)
observedFeatures <- 0
results <- NULL
snpAnnotation <- NULL
featureAnnotation <- NULL
filesToRead <- list.files(".",pattern = "qtl_results",full.names = T)
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
      df <- df[order(df$empirical_feature_p_value, df$p_value,decreasing = F),]
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
      featureAnnotationTmp <- read.delim(paste("./feature_metadata",baseName,".txt",sep=""),as.is=T)
      featureAnnotation <- rbind(featureAnnotation,featureAnnotationTmp)
    } else {
      print(paste("Skipping",i,"because not necessary data is available or to many files with the same name."))
    }
    if(length(tmp)>0){
      for (j in names(tmp)) tmp[[j]][["feature"]] <- j
      observedFeatures = observedFeatures+length(tmp)
      df <- bind_rows(tmp)
      if(nrow(df)>0){
        df <- df[order(df$empirical_feature_p_value, df$p_value,decreasing = F),]
        results = rbind(results,df)  
      }
    }
  }
}

rm(snpAnnotationTmp, featureAnnotationTmp, df)

snpAnnotation <- unique(snpAnnotation)
featureAnnotation <- unique(featureAnnotation)

##Add information
t <- match(results$snp_id, snpAnnotation$snp_id)
resultsAnnotated <- cbind(results, snpAnnotation[t,])

t <- match(results$feature , featureAnnotation$feature_id )
resultsAnnotated <- cbind(resultsAnnotated, featureAnnotation[t,])

results <- resultsAnnotated
results <- results[,-c(5,6)]
rm(resultsAnnotated, snpAnnotation, featureAnnotation)
##

resultsFull <- results
results <- resultsFull[which(!duplicated(resultsFull$feature)),]

save(results, resultsFull, file = "../Trans_20180503.Rds")
#load("../Trans_20180503.Rds")


if(length(which(is.na(results$empirical_feature_p_value)))!=0){
  results <- results[-which(is.na(results$empirical_feature_p_value)),]
}

##Multiple testing

results["global_corr_p_value"] <- qvalue(results$empirical_feature_p_value)$qvalues
results <- results[order(results$global_corr_p_value,results$empirical_feature_p_value, results$p_value,decreasing = F),]

write.table(paste(baseFolder,"top_results_global_level.txt",sep=""),x = results,sep="\t",row.names=F,quote=F)
write.table(paste(baseFolder,"results_global_level.txt",sep=""),x = resultsFull,sep="\t",row.names=F,quote=F)

