library(data.table)
library(tidyverse)
library(ggplot2)
library(knitr)
library(patchwork)
library(stringr)
library(dplyr)
library(readxl)

# read in files
annotation<-"../VEP-filtered_200k_MC4R.txt"
curated<-"../prioritizedVariants_withGR38pos-UCSC_plink.txt"
exomeSeq<-"../Obesity/final_Obesity.raw"
diseaseCarriers<-"../carriers_by_obesity-variants.txt"
phenotypes<-"../phenotypes/20211213/ukbbExomeSeq_diabetesByDoc-collapsed_ageDiabetesDiag-collapsed_HbA1cCollapsedConverted_allDiabetes-collapsed.txt"
esm1b<-"../ESM1b/MC4R/MC4R.xlsx"
pcs<-"../phenotypes/basic_demographics.tab"
plottingRanges<-"../variant_prediction_methods_ranges.csv"
annotation<-data.frame(annotation,header=TRUE)
curated<-read.table(curated,header=TRUE,sep ="\t")
exomeSeq<-data.frame(exomeSeq)
diseaseCarriers<-read.table(diseaseCarriers,header = TRUE)
phenotypes<-read.table(phenotypes,header=TRUE)
esm1b<-read_excel(esm1b)
pcs<-data.frame(fread(pcs),header = TRUE)
plottingRanges<-read.csv(plottingRanges)

#file clean up
condition<-"Obesity"
gene<-"MC4R"
diseaseCarriers<-diseaseCarriers$f.eid
colnames(phenotypes)
phenKeep<-"maxBMI"
units<-"(kg/m^2)"
thresholdPhen<-30
phenotypes<-phenotypes[,c("f.eid",phenKeep)]
phenotypes<-phenotypes[!is.na(phenotypes[,c(phenKeep)]),]
pcs$maxAge<-apply(pcs[,c("age_visit_0","age_visit_1","age_visit_2","age_visit_3")], 1, max, na.rm=TRUE)
phenotypes<-merge(x=phenotypes,y=pcs,by.x="f.eid",by.y="f.eid",all.x=TRUE)
esm1b$aaChange<-gsub(".*_","",esm1b$variant)
colnames(exomeSeq)[7:ncol(exomeSeq)]<-gsub(x = colnames(exomeSeq)[7:ncol(exomeSeq)],pattern = "\\.",replacement = ":")
colnames(exomeSeq)[7:ncol(exomeSeq)]<-gsub(x = colnames(exomeSeq)[7:ncol(exomeSeq)],pattern = "\\_.*",replacement = "")
curated$vep<-gsub(x = curated$plink,pattern = "\\.",replacement = ":")
curated$vep<-gsub("\\_.*","",curated$vep)
curated<-curated%>%
  dplyr::filter(Gene==gene)
transcript<-"ENST00000299766"
annotation<-annotation[annotation$Feature==transcript,]
annotation$refAA<-unlist(lapply(str_split(annotation$Amino_acids, "/"), `[`, 1))
annotation$altAA<-unlist(lapply(str_split(annotation$Amino_acids, "/"), `[`, 2))
annotation$aaChange<-paste0(annotation$refAA,annotation$Protein_position,annotation$altAA)
annotation$function_change<-"Unknown"
annotation$function_change[!grepl(pattern = paste(c("benign","conflicting","not","uncertain"),collapse = "|"),annotation$CLIN_SIG)&!is.na(annotation$CLIN_SIG)]<-"ClinVar LP/P"
annotation$function_change[annotation$X.Uploaded_variation%in%(curated$vep[curated$Condition==condition])]<-condition
GoF<-c("T11S", "V103I","F201L","G231S","I251L","I289L","L304F","I317V","Y332C") # From Lotta et al. 2019
annotation$function_change[annotation$aaChange%in%GoF]<-"Protective"
annotation<-merge(x=annotation,y=esm1b,by="aaChange", all.x=TRUE)
annotation$gnomADe_AF[is.na(annotation$gnomADe_AF)]<-0

# filter for individuals with one missense variant, and no other coding variants in monogenic gene
check<-exomeSeq[,colnames(exomeSeq)%in%(annotation$X.Uploaded_variation[!annotation$Consequence%in%c("missense_variant,splice_region_variant","missense_variant","synonymous_variant","intron_variant","3_prime_UTR_variant","5_prime_UTR_variant")])]
ncol(check)
check<-cbind(exomeSeq[exomeSeq$FID%in%ids,1:6],check)
check<-check%>%
  filter(rowSums(check[7:(ncol(check))],na.rm = TRUE)==0)
nrow(check)
check<-exomeSeq[exomeSeq$FID%in%check$FID,]
ids<-check$FID
check<-check[,colnames(check)%in%(annotation$X.Uploaded_variation[annotation$Consequence%in%c("missense_variant,splice_region_variant","missense_variant")])]
check<-cbind(exomeSeq[exomeSeq$FID%in%ids,1:6],check)
numMissense<-data.frame(matrix(nrow = nrow(check),ncol = 2))
colnames(numMissense)<-c("FID","number_missense")
missenseCount<-function(x){
  sum(x[7:ncol(check)][!is.na(x[7:ncol(check)])]!=0)
}
test<-apply(check,1,missenseCount)
numMissense$FID<-check$FID
numMissense$number_missense<-test
oneMiss<-exomeSeq[exomeSeq$FID%in%numMissense$FID[numMissense$number_missense==1],]
oneMiss<-oneMiss[,colnames(oneMiss)%in%(annotation$X.Uploaded_variation[annotation$Consequence%in%c("missense_variant,splice_region_variant","missense_variant")])]
oneMiss<-cbind(exomeSeq[exomeSeq$FID%in%numMissense$FID[numMissense$number_missense==1],1:6],oneMiss)
whichOneMiss<-data.frame(oneMiss$FID)
colnames(whichOneMiss)<-"FID"
whichOneMiss$missensePlink<-NA
for(i in 1:nrow(oneMiss)){
  test<-oneMiss[i,]
  for (j in 7:ncol(test)) {
    if(is.na(test[1,j])){next}
    else if(test[1,j]!=0){
      whichOneMiss[i,1]<-test[1,"FID"]
      whichOneMiss[i,2]<-colnames(test)[j]
      break
    }
  }
}
whichOneMiss<-merge(x=whichOneMiss,y=annotation,by.x = "missensePlink",by.y="X.Uploaded_variation")
whichOneMiss<-merge(x=whichOneMiss,y=phenotypes,by.x="FID",by.y="f.eid")
by_maxMiss<-whichOneMiss%>%
  group_by(missensePlink)%>%
  summarise(mean_of_phen=mean(!!sym(phenKeep),na.rm=TRUE))
by_maxMiss<-merge(x=by_maxMiss,y=annotation,by.x="missensePlink",by.y="X.Uploaded_variation")
by_maxMiss<-by_maxMiss[!is.na(by_maxMiss$mean_of_phen),]

# plotting mean phenotype by variant score 
for(i in plottingRanges$method){
  if(i=="EVE"){
    cor<-cor.test(by_maxMiss[,c(plottingRanges[plottingRanges$method==i,"score_name"])],by_maxMiss$mean_of_phen,use="complete.obs")
    print(
      ggplot(by_maxMiss,aes(x=!!sym(plottingRanges[plottingRanges$method==i,"score_name"]),y=mean_of_phen))+
        geom_point(aes(color=factor(function_change),shape=EVE_CLASS),size=3,alpha = 0.6)+
        theme_classic(base_size = 12)+
        labs(y=paste(phenKeep,units),
             x=paste0(i," score"),
             subtitle =  paste0("1 missense & only intron/syn/UTR,n=",nrow(whichOneMiss),",R=", round(cor$estimate,digits = 4),",cor p-test=",formatC(cor$p.value,digits = 4,format="e")),
             title = paste0("mean(phen of variant_i carriers)~",i," score of variant_i"),
             color="Variant annotation",
             caption = paste0("at least 1 carrier of ",sum(!is.na(by_maxMiss[,plottingRanges[plottingRanges$method==i,"score_name"]]))," variants"))+
        geom_hline(yintercept=thresholdPhen, linetype="dashed", color = "red",size=1.2,alpha=.5)+
        geom_smooth(method=lm, se=FALSE,size=0.8,color="black",alpha=.3)+
        scale_color_manual(values=c("deepskyblue","navy","deeppink3","grey73"))+
        xlim(plottingRanges[plottingRanges$method==i,"range_min"],plottingRanges[plottingRanges$method==i,"range_max"])
    )
  }
  else if(i=="RAW CADD"){
    cor<-cor.test(by_maxMiss[,c(plottingRanges[plottingRanges$method==i,"score_name"])],by_maxMiss$mean_of_phen,use="complete.obs")
    print(
      ggplot(by_maxMiss,aes(x=!!sym(plottingRanges[plottingRanges$method==i,"score_name"]),y=mean_of_phen))+
        geom_point(aes(color=factor(function_change)),size=3,alpha = 0.6)+
        theme_classic(base_size = 12)+
        labs(y=paste(phenKeep,units),
             x=paste0(i," score"),
             subtitle =  paste0("1 missense & only intron/syn/UTR,n=",nrow(whichOneMiss),",R=", round(cor$estimate,digits = 4),",cor p-test=",formatC(cor$p.value,digits = 4,format="e")),
             title = paste0("mean(phen of variant_i carriers)~",i," score of variant_i"),
             color="Variant annotation",
             caption = paste0("at least 1 carrier of ",sum(!is.na(by_maxMiss[,plottingRanges[plottingRanges$method==i,"score_name"]]))," variants"))+
        geom_hline(yintercept=thresholdPhen, linetype="dashed", color = "red",size=1.2,alpha=.5)+
        geom_smooth(method=lm, se=FALSE,size=0.8,color="black",alpha=.3)+
        scale_color_manual(values=c("deepskyblue","navy","deeppink3","grey73"))
    )
  }
  else{
    cor<-cor.test(by_maxMiss[,c(plottingRanges[plottingRanges$method==i,"score_name"])],by_maxMiss$mean_of_phen,use="complete.obs")
    print(
      ggplot(by_maxMiss,aes(x=!!sym(plottingRanges[plottingRanges$method==i,"score_name"]),y=mean_of_phen))+
        geom_rect(aes(xmin = plottingRanges[plottingRanges$method==i,"likelyPath_min"], xmax = plottingRanges[plottingRanges$method==i,"likelyPath_max"], ymin = -Inf, ymax = Inf),
                  fill = "lightcyan", alpha = 0.1)+
        geom_point(aes(color=factor(function_change)),size=3,alpha = 0.6)+
        theme_classic(base_size = 12)+
        labs(y=paste(phenKeep,units),
             x=paste0(i," score"),
             subtitle =  paste0("1 missense & only intron/syn/UTR,n=",nrow(whichOneMiss),",R=", round(cor$estimate,digits = 4),",cor p-test=",formatC(cor$p.value,digits = 4,format="e")),
             title = paste0("mean(phen of variant_i carriers)~",i," score of variant_i"),
             color="Variant annotation",
             caption = paste0("at least 1 carrier of ",sum(!is.na(by_maxMiss[,plottingRanges[plottingRanges$method==i,"score_name"]]))," variants"))+
        geom_hline(yintercept=thresholdPhen, linetype="dashed", color = "red",size=1.2,alpha=.5)+
        geom_smooth(method=lm, se=FALSE,size=0.8,color="black",alpha=.3)+
        scale_color_manual(values=c("deepskyblue","navy","deeppink3","grey73"))
    )
    print(
      ggplot(by_maxMiss,aes(x=!!sym(plottingRanges[plottingRanges$method==i,"score_name"]),y=mean_of_phen))+
        geom_point(aes(color=factor(function_change)),size=3,alpha = 0.6)+
        theme_classic(base_size = 12)+
        labs(y=paste(phenKeep,units),
             x=paste0(i," score"),
             subtitle =  paste0("1 missense & only intron/syn/UTR,n=",nrow(whichOneMiss),",R=", round(cor$estimate,digits = 4),",cor p-test=",formatC(cor$p.value,digits = 4,format="e")),
             title = paste0("mean(phen of variant_i carriers)~",i," score of variant_i"),
             color="Variant annotation",
             caption = paste0("at least 1 carrier of ",sum(!is.na(by_maxMiss[,plottingRanges[plottingRanges$method==i,"score_name"]]))," variants"))+
        geom_hline(yintercept=thresholdPhen, linetype="dashed", color = "red",size=1.2,alpha=.5)+
        geom_smooth(method=lm, se=FALSE,size=0.8,color="black",alpha=.3)+
        scale_color_manual(values=c("deepskyblue","navy","deeppink3","grey73"))
    )
  }
}

# adjusting for covariates
regressOutCorr<-function(missenseInd,pheno){
  variables<-c(paste0("PC",c(1:10)),"factor(genetic_sex)","maxAge","PRS_center_scale")
  f<-as.formula(paste(
    pheno,
    paste(variables,collapse=" + "),
    sep=" ~ ")
  )
  regression_reg<-lm(f,data=missenseInd,na.action=na.exclude)
  missenseInd$residual<-resid(regression_reg)
  by_maxMissResid<-missenseInd%>%
    group_by(missensePlink)%>%
    summarise(mean_of_resid=mean(residual,na.rm=TRUE))
  nrow(by_maxMissResid)
  sum(is.na(by_maxMissResid$mean_of_resid))#cant use these variants
  by_maxMissResid<-merge(x=by_maxMissResid,y=annotation,by.x="missensePlink",by.y="X.Uploaded_variation")
  by_maxMissResid<-by_maxMissResid[!is.na(by_maxMissResid$mean_of_resid),]
  cor<-cor.test(by_maxMissResid$mean_of_resid,by_maxMissResid$score)
  return(c(wantPRS,cor$estimate,cor$p.value))
}
saveResults<-data.frame(matrix(ncol=7,nrow=3))
colnames(saveResults)<-c("individuals_included","MAF_max","total_individuals","total_missense","PRS_regressed_out","pearson_esm1b_meanResid","p_corr")
iterate<-c(0.01,0.001,0.0001,0.00001,0.000001)
for(i in iterate){
  test<-df[df$gnomADe_AF<i,]
  test<-whichOneMiss[whichOneMiss$gnomADe_AF<i,]
  saveResults[nrow(saveResults)+1,]<-c("everyone",i,nrow(test),length(unique(test$missensePlink)),
                                       regressOutCorr(missenseInd = test,pheno=phenKeep,wantPRS = FALSE))
}

# additional carrier sets based on different definitions
annotation$function_change<-"Unknown"
annotation$function_change[grepl(pattern = paste(c("pathogenic$","pathogenic\\,"),collapse = "|"),annotation$CLIN_SIG)&!is.na(annotation$CLIN_SIG)]<-"ClinVar LP/P"
check<-exomeSeq[,colnames(exomeSeq)%in%annotation$X.Uploaded_variation[annotation$function_change=="ClinVar LP/P"]]
check<-cbind(exomeSeq[,1:6],check)
check<-check%>%
  filter(rowSums(check[7:(ncol(check))],na.rm = TRUE)>=1)
forGenes<-data.frame(check$FID)
colnames(forGenes)<-"FEID"
forGenes$monogenicGene<-gene
write.table(forGenes,file="../loose_LP-P/carriers_by_obesity-clinVar-genes.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
annotation$function_change<-"Unknown"
annotation$function_change[!grepl(pattern = paste(c("benign","conflicting","not","uncertain","protective"),collapse = "|"),annotation$CLIN_SIG)&!is.na(annotation$CLIN_SIG)]<-"ClinVar LP/P"
check<-exomeSeq[,colnames(exomeSeq)%in%annotation$X.Uploaded_variation[annotation$function_change=="ClinVar LP/P"]]
check<-cbind(exomeSeq[,1:6],check)
check<-check%>%
  filter(rowSums(check[7:(ncol(check))],na.rm = TRUE)>=1)
nrow(check)
forGenes<-data.frame(check$FID)
colnames(forGenes)<-"FEID"
forGenes$monogenicGene<-gene
write.table(forGenes,file="../strict_LP-P/carriers_by_obesity-clinVar-genes.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
VUS_esm1b_LP<-whichOneMiss%>%
  filter(score<(-7.5)&function_change=="Unknown")
write.table(VUS_esm1b_LP, "missense_carriers_ESM1b-LP-P_MC4R.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)
