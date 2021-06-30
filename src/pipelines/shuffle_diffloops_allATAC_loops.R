setwd("/oak/stanford/groups/khavari/users/mguo123/psych_project/data/GWAS_enrichment")
SNP.all <- read.delim("./resources/LinkedSNP-all-sorted.bed",header=F,stringsAsFactors = F)
DiseaseName <- names(table(SNP.all$V4))

for (tissue in c("Astro" , "H9_D0", "H9_D10", "H9_D28", "SL_D0",  "SLC_D0")){
print(tissue)
file.shuffle <- list.files('./result_diffloop_atac/',pattern=paste0("^",tissue,"(.*?)\\_overlap_all.bed$"),full.names=T)
                                                                                                                                   
                                                                                                                     

tab.shuffle <- matrix(NA,ncol=length(DiseaseName),nrow=length(file.shuffle))

for (i in 1:length(file.shuffle)){
  info <- file.info(file.shuffle[i])
  
  if (info$size != 0){
    
    tmp <- read.delim(file.shuffle[i],header=F,stringsAsFactors = F)
    disease_SNP <- unique(tmp[,4:5])
    tab.tmp <- table(disease_SNP$V4)
    tab.shuffle[i,] <- tab.tmp[match(DiseaseName,names(tab.tmp))]}
}
colnames(tab.shuffle) <- DiseaseName
tab.shuffle <- data.frame(File=file.shuffle,tab.shuffle,check.names = F,stringsAsFactors = F)
write.table(tab.shuffle,file=paste0("./result_diffloop_atac/",tissue,"_shuffle_100Times_actuallastline_diffloop_allATAC_loops.txt"),sep='\t',quote=F,row.names=F)

}

## psycmpra
SNP.psych <- read.delim("./resources/psych_mpra_lib.bed",header=F,stringsAsFactors = F)
SNPs_arr <- names(table(SNP.psych$V4))

for (tissue in c("Astro" , "H9_D0", "H9_D10", "H9_D28", "SL_D0",  "SLC_D0")){
print(tissue)
file.shuffle <- list.files('./result_diffloop_atac/',pattern=paste0("^",tissue,"(.*?)\\_overlap_psychmpra.bed$"),full.names=T)
                                                                                                                                                                                                                                                   
tab.shuffle <- matrix(NA,ncol=length(SNPs_arr),nrow=length(file.shuffle))

for (i in 1:length(file.shuffle)){
  info <- file.info(file.shuffle[i])
 print(i) 
  if (info$size != 0){
    
    tmp <- read.delim(file.shuffle[i],header=F,stringsAsFactors = F)
    disease_SNP <- unique(tmp[,3:4])
    tab.tmp <- table(disease_SNP$V4)
    tab.shuffle[i,] <- tab.tmp[match(SNPs_arr,names(tab.tmp))]}
}

colnames(tab.shuffle) <- SNPs_arr 
tab.shuffle <- data.frame(File=file.shuffle,tab.shuffle,check.names = F,stringsAsFactors = F)
write.table(tab.shuffle,file=paste0("./result_diffloop_atac/",tissue,"_shuffle_100Times_actuallastline_diffloop_allATAC_loops_psych.txt"),sep='\t',quote=F,row.names=F)

}#eend tissue





# ### Need to merge disease names
#   [3] "Alzheimer Dementia"                                                                                                                                                   
#   [4] "Alzheimer Dementia,Cardiovascular Disease,Coronary artery disease,Myocardial Infarction"                                                                              
#   [5] "Alzheimer Dementia,Cardiovascular Disease,Longevity"                                                                                                                  
#   [6] "Alzheimer Dementia,Coronary artery disease,Longevity,Type II diabetes"                                                                                                
#   [7] "Alzheimer Dementia,Crohn’s Disease,Inflammatory Bowel Disease"                                                                                                        
#   [8] "Alzheimer Dementia,Myocardial Infarction"                                                                                                                             
#   [9] "Alzheimer Dementia,Parental Lifespan,Longevity"                                                                                                                       
#  [10] "Alzheimer Dementia,Schizophrenia"                                                                                                                                     
#  [11] "Alzheimer Dementia,Stroke"                                                                                                                                            
#  [12] "Alzheimer Dementia,Systemic Lupus Erythematosus"                                                                                                                      
#  [13] "Alzheimer Dementia,Type II diabetes"                                                                                                                                  
#  [14] "Anxiety"                                                                                                                                                                                                                                                                                                 
#  [25] "Asthma,Longevity"                                                                                                                                                                                                                                                                                                      
#  [31] "Attention Deficit Hyperactivity Disorder"                                                                                                                             
#  [32] "Attention Deficit Hyperactivity Disorder,Bipolar"                                                                                                                     
#  [33] "Attention Deficit Hyperactivity Disorder,Bipolar,Hypertension,Schizophrenia"                                                                                          
#  [34] "Attention Deficit Hyperactivity Disorder,Bipolar,Sarcoidosis,Schizophrenia"                                                                                           
#  [35] "Attention Deficit Hyperactivity Disorder,Bipolar,Schizophrenia"                                                                                                       
#  [36] "Attention Deficit Hyperactivity Disorder,Bipolar,Schizophrenia,Ulcerative colitis"                                                                                    
#  [37] "Attention Deficit Hyperactivity Disorder,Hypertension"                                                                                                                
#  [38] "Attention Deficit Hyperactivity Disorder,Rheumatoid_Arthritis"                                                                                                                                                                                                                                                          
#  [41] "Bipolar"                                                                                                                                                              
#  [42] "Bipolar,Depression"                                                                                                                                                   
#  [43] "Bipolar,Depression,Schizophrenia"                                                                                                                                     
#  [44] "Bipolar,Ischemic Stroke,Stroke"                                                                                                                                       
#  [45] "Bipolar,Schizophrenia"                                                                                                                                                
#  [46] "Bipolar,Schizophrenia,Type I diabetes"                                                                                                                                                                                                                                                                                     
#  [49] "Borderline Personality Disorder"                                                                                                                                                                                                                              
#  [80] "Coronary artery disease,Parental Lifespan,Longevity"                                                                                                                                                                                                                                                 
# [117] "Depression"                                                                                                                                                           
# [118] "Depression,Schizophrenia"                                                                                                                                                                                                                                                                                                   
# [130] "Hypertension,Longevity"                                                                                                                                                                                                                                                                            
# [132] "Hypertension,Schizophrenia"                                                                                                                                                                                                                  
# [138] "Inflammatory Bowel Disease,Schizophrenia"                                                                                                                                                                                                                                                                                   
# [151] "Longevity"                                                                                                                                                                                                                                                                                                
# [160] "Multiple Sclerosis,Schizophrenia"                                                                                                                                                                                                                                                                                   
# [164] "Obsessive Compulsive Disorder"                                                                                                                                        
# [165] "OCD"                                                                                                                                                                  
# [170] "Panic Disorder"                                                                                                                                                       
# [171] "Parental Lifespan,Longevity"                                                                                                                                          
# [172] "Parental Lifespan,Longevity,Myocardial Infarction"                                                                                                                    
# [173] "Post-traumatic Stress Disorder"                                                                                                                                                                                                                                                                  
# [184] "Schizophrenia"                                                                                                                                                        
# [185] "Schizophrenia,Type II diabetes"    