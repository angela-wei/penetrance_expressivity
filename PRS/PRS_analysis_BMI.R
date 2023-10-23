# read in files & clean
library(data.table)
library(ggplot2)
library(tidyverse)
library(patchwork)
prs<-"../bmi_prs_khera.sum"
phenotypes<-"../phenotypes/20211213/ukbbExomeSeq_diabetesByDoc-collapsed_ageDiabetesDiag-collapsed_HbA1cCollapsedConverted_allDiabetes-collapsed.txt"
pcs<-"../phenotypes/generalCovariatesEUR_unrel.csv"
allCurCarriers<-"../carriers_by_obesity-variants.txt"
clinVarLooseCarriers<-"../loose_LP-P/carriers_by_obesity-clinVar-genes.txt"
clinVarstrictCarriers<-"../strict_LP-P/carriers_by_obesity-clinVar-genes.txt"
esm1b<-"../missense_carriers_ESM1b-LP-P_MC4R.csv"
prs<-data.frame(fread(prs))
phenotypes<-data.frame(fread(phenotypes))
allCurCarriers<-data.frame(fread(allCurCarriers))
clinVarLooseCarriers<-read.delim(clinVarLooseCarriers,header=TRUE,sep=" ")
clinVarstrictCarriers<-read.delim(clinVarstrictCarriers,header=TRUE,sep=" ")
esm1b<-data.frame(fread(esm1b,header=TRUE,fill=TRUE))
pcs<-data.frame(fread(computerpath,pcs))
colnames(prs)<-c("FID","PRS")
prs<-prs[prs$FID%in%pcs$f.eid,] #filter for eur and unrel
phenKeep<-"maxBMI"
phenotypes<-phenotypes[!is.na(phenotypes[,phenKeep]),]
everything<-merge(x=prs,y=phenotypes,by.x="FID",by.y="f.eid")
everything<-everything[!is.na(everything[,phenKeep]),]
everything<-merge(x=everything,y=pcs,by.x="FID",by.y="f.eid")
prs_mean<-mean(everything$PRS,na.rm = TRUE)
prs_std<-sd(everything$PRS, na.rm = TRUE)
everything$PRS_center_scale<-(everything$PRS-prs_mean)/prs_std

# Figure 5: PRS within carriers
carrier_prs<-everything[everything$FID%in%allCurCarriers$f.eid,]
variables<-c(paste0("PC",c(1:10)),"factor(sexMale)","age","PRS_center_scale")
carrier_prs_regress<-as.formula(paste(phenKeep,paste(variables,collapse=" + "),sep=" ~ "))
carrier_prs_regress<-lm(carrier_prs_regress,data=carrier_prs,na.action=na.exclude)
summary(carrier_prs_regress)
cor.test(x = carrier_prs[,"PRS_center_scale"],y=carrier_prs[,phenKeep])
ggplot(carrier_prs,aes(x=PRS_center_scale,y=!!sym(phenKeep)))+
  geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  theme_classic()+
  labs(x="BMI PRS, centered & scaled",
       y="BMI, kg/m^2",
       subtitle=paste0(nrow(carrier_prs)," unrel, eur, obesity carriers in UKB200K"),
       title="Obesity carrier BMI correlates with carrier BMI PRS")

# Figure 4: carrier phenotype
carriers_mean_phen<-data.frame(matrix(nrow=1,ncol = 5))
colnames(carriers_mean_phen)<-c("group","total_individuals","meanPhen","L95CI","U95CI")
summaryInfo<-function(ids,groupName){
  phenToCheck<-phenotypes[phenotypes$f.eid%in%ids,]
  carriers_ci<-t.test(phenToCheck[,phenKeep])
  return(c(groupName,sum(!is.na(phenToCheck[,phenKeep])),mean(phenToCheck[,phenKeep]),carriers_ci$conf.int[1],carriers_ci$conf.int[2]))
}
carriers_mean_phen[1,]<-summaryInfo(ids=allCurCarriers$f.eid,groupName = "cur")
carriers_mean_phen[nrow(carriers_mean_phen)+1,]<-summaryInfo(ids=clinVarLooseCarriers$FEID,groupName = "CV_L")
carriers_mean_phen[nrow(carriers_mean_phen)+1,]<-summaryInfo(ids=clinVarstrictCarriers$FEID,groupName = "CV_S")
carriers_mean_phen[nrow(carriers_mean_phen)+1,]<-summaryInfo(ids=esm1b$FID,groupName = "unk-7.5")
carriers_mean_phen$meanPhen<-as.numeric(carriers_mean_phen$meanPhen)
carriers_mean_phen$L95CI<-as.numeric(carriers_mean_phen$L95CI)
carriers_mean_phen$U95CI<-as.numeric(carriers_mean_phen$U95CI)
diff_carriers_plot<-ggplot(carriers_mean_phen,aes(x=group,y=meanPhen))+
  geom_point()+
  theme_classic(base_size = 15)+
  ylab("")+
  geom_errorbar(aes(ymin=L95CI, ymax=U95CI), width=.1)+
  ylim(23,33)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = margin(0, 0, 0, 0, "pt"))
diff_carriers_plot

# Figure 4: non-carriers
noncarriers<-everything[!everything$FID%in%c(allCurCarriers$f.eid,clinVarLooseCarriers$FEID,clinVarstrictCarriers$FEID,esm1b$FID),]
noncarriers<-noncarriers[order(noncarriers$PRS_center_scale),]
noncarriers$PRS_center_scale_percentile<-rank(noncarriers$PRS_center_scale)/length(noncarriers$PRS_center_scale)
makePercentiles<-function(numPercentile,df){
  print(paste(nrow(df)/numPercentile," indiviudals per bin"))
  percentiles<-data.frame(c(0:(numPercentile-1)))
  colnames(percentiles)<-"percentile"
  percentiles$meanPhen<-NA
  percentiles$L95CI<-NA
  percentiles$U95CI<-NA
  for (i in percentiles$percentile) {
    check<-df%>%
      filter(PRS_center_scale_percentile>=i/numPercentile & PRS_center_scale_percentile<(i+1)/numPercentile)
    percentiles[i+1,"meanPhen"]<-mean(check[,phenKeep],na.rm=TRUE)
    ci<-t.test(check$maxBMI)
    percentiles[i+1,"L95CI"]<-ci$conf.int[1]
    percentiles[i+1,"U95CI"]<-ci$conf.int[2]
  }
  percentiles$comparedToCarriers<-percentiles$meanPhen>=carriers_mean_phen[carriers_mean_phen$group=="cur","meanPhen"]
  return(percentiles)
}
for (i in c(100,250,500,1000)) {
  percentiles<-makePercentiles(i,noncarriers)
  noncarriersPercentile<-ggplot(percentiles,aes(x=percentile,y=meanPhen))+
    geom_errorbar(aes(ymin=L95CI, ymax=U95CI,color=comparedToCarriers), width=.1,show.legend = FALSE)+
    geom_point(show.legend = FALSE,aes(color=comparedToCarriers))+
    theme_classic(base_size = 15)+
    ylab("Mean BMI, kg/m^2")+
    xlab(paste0("Percentile BMI PRS noncarriers (",i," bins)"))+
    ylim(23,33)+
    scale_color_manual(values=c("black","red")) +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"))
  p<-noncarriersPercentile+plot_spacer()+diff_carriers_plot+plot_layout(widths=c(5,0.005,3))
  p
  print(p)
}
noncarriersGreater<-percentiles[percentiles$meanPhen>=carriers_mean_phen[carriers_mean_phen$group=="cur","meanPhen"],]
