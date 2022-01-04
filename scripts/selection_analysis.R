#####**** libraries and functions ****#####
#####*** libraries ***#####
library("data.table")
library("gdata")
library("ape")
library("org.Hs.eg.db")
library("ROCR")
library("RColorBrewer") 
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3")) 
library("plyr")
#####*** function: change the format of deltaSVM data ***#####
dataMod<-function(deltaSVM) {
  ## change the first column into bed format
  deltaSVM<-as.data.frame(deltaSVM)
  names(deltaSVM)<-c("ID","deltaSVM","varNumb","pValue")
  splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
  splString<-data.frame(unlist(splString))
  ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
  ID.bed<-ID.bed[,c(1:3)]
  ID.bed[,2]<-as.numeric(ID.bed[,2])-1
  deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:4)])
  colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
  deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
  deltaSVM$end<-as.numeric(as.character(deltaSVM$end))
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  return(deltaSVM)
}


#####**** mouse analysis ****#####
#####*** take CEBPA binding sites as an example ***#####
#####** check the performance of trained model to distinguish binding sites and random sequences **#####
## five fold cross validation 
cv<-fread("data/mouse_SVM_model/CEBPA/bl6-CEBPA_gkmtrain.cvpred.txt")
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
## plot
plot(perf,lwd = 3,cex=1.4)
# calculate AUC
auc_result <-performance( pred, measure = "auc")
unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,"(AUC = 0.968)",bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####** deltaSVM summary **#####
deltaSVMcons<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_conserved_deltaSVM_highertailTest.txt")
deltaSVMgain<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_specific_gain_deltaSVM_highertailTest.txt")
deltaSVMloss<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_specific_loss_deltaSVM_lowertailTest.txt")
## change the format of deltaSVM files
deltaSVMcons<-dataMod(deltaSVMcons)
deltaSVMgain<-dataMod(deltaSVMgain)
deltaSVMloss<-dataMod(deltaSVMloss)
## plot deltaSVM and pvalue, use conserved binding sites as an example 
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
n=2
hist(deltaSVMcons$deltaSVM,breaks = 60,main="Conserved",xlab="deltaSVM",xlim=c(-30,30),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])
hist(deltaSVMcons$pValue,breaks = 60,main="Conserved",xlab="Pvalue",xlim=c(0,1),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[3])
## calculate proportion of positive selection
par(mfrow=c(1,1))
par(mar=c(8,7,4,2))
posCons<-deltaSVMcons[deltaSVMcons[,6]<0.01,]
posGain<-deltaSVMgain[deltaSVMgain[,6]<0.01,]
posLoss<-deltaSVMloss[deltaSVMloss[,6]<0.01,]
barplot(c(nrow(posCons)/nrow(deltaSVMcons),nrow(posGain)/nrow(deltaSVMgain),nrow(posLoss)/nrow(deltaSVMloss)),
        ylim=c(0,0.3),col=pal[2],ylab="Proportion of positive binding sites \n (p<0.01)",main="CEBPA binding sites",
        cex.lab=2,cex.main=2,cex.axis=2)
text(x=c(1,2,3),y=-0.02,cex=2,srt = 45,adj = 1,labels = c("Conserved","Gain","Loss"),xpd = TRUE)

#########** validate positive selection based on binding intensity **########
## use conserved binding sites as an example 
binding_intensity<-fread("data/mouse_ChIP-Seq/bl6_CEBPA.txt")
binding_intensity$chr<-paste0(rep("chr",nrow(binding_intensity)),binding_intensity$chr)
colnames(binding_intensity)[c(2,3)]<-c("start","end")
deltaSVM_intensity<-merge(deltaSVMcons,binding_intensity,by=c("chr","start","end"))
pos<-subset(deltaSVM_intensity,deltaSVM_intensity$pValue<0.01)
nonPos<-subset(deltaSVM_intensity,deltaSVM_intensity$pValue>=0.01)
## plot the intensity between positive binding sites and non-positive binding sites
par(mfrow=c(1,1))
par(mar=c(7,5,4,2))
boxplot(pos$intensity,nonPos$intensity,col=c(pal[1],pal[2]),
        ylim= c(-100,2000),ylab="Binding intensity in C57BL/6J",xaxt = "n",
        notch=T,pch=16,outcex=0.5,main="Conserved",cex.lab=1.2,cex.main=1.2
)
text(x=c(1,2),y=-100-2000/15,cex.lab=1.2,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
text(x=c(1,2), y=-100, labels=paste0("n=", c(nrow(pos),nrow(nonPos))))
wtest<-wilcox.test(pos$intensity,nonPos$intensity)
legend("topleft",legend=paste("p=",signif(wtest$p.value,3)),bty = 'n')


#####**** human analysis ****#####
#####*** take CEBPA binding sites as an example ***#####
#####** check the performance of trained model to distinguish binding sites and random sequences **#####
## five fold cross validation 
cv<-fread("data/human_SVM_model/CEBPA/hsap_CEBPA_gkmtrain.cvpred.txt")
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
## plot
plot(perf,lwd = 3,cex=1.4)
# calculate AUC
auc_result <-performance( pred, measure = "auc")
unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,"(AUC = 0.986)",bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####** deltaSVM summary **#####
deltaSVM<-fread("data/human_deltaSVM/hsap_CEBPA_deltaSVM_highertailTest.txt")
## change the format of deltaSVM files
deltaSVM<-dataMod(deltaSVM)
## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
n=2
hist(deltaSVM$deltaSVM,breaks = 60,main="CEBPA binding sites",xlab="deltaSVM",xlim=c(-30,30),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])
hist(deltaSVM$pValue,breaks = 60,main="CEBPA binding sites",xlab="Pvalue",xlim=c(0,1),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[3])

#########** validation for positive selection  **########
#####* compare the ratio of #substitutions and #polymorphisms betwween positive and non-positive sites  *#####
sub_poly_Numb<-read.table("data/human_substitutions_polymorphisms/CEBPA_sub_poly.txt",header = T)
sub_poly_Numb$subNumb<-as.numeric(as.character(sub_poly_Numb$subNumb))
sub_poly_Numb$polyNumb<-as.numeric(as.character(sub_poly_Numb$polyNumb))
sub_poly_Numb<-na.omit(sub_poly_Numb)
deltaSVM_sub_poly_Numb<-merge(deltaSVM,sub_poly_Numb,by=c("chr","start","end"))
## positive and non-positive sites
pos<-subset(deltaSVM_sub_poly_Numb,deltaSVM_sub_poly_Numb$pValue<0.01)
nonPos<-subset(deltaSVM_sub_poly_Numb,deltaSVM_sub_poly_Numb$pValue>=0.01)
subNumb<-c(sum(pos$subNumb),sum(nonPos$subNumb))
polyNumb<-c(sum(pos$polyNumb),sum(nonPos$polyNumb))
## plot
par(mfrow=c(1,1))
par(mar=c(9,5,4,2))
bp<-barplot(subNumb/polyNumb,ylim=c(0,1.5),main="Human TFBS (CEBPA)",cex.lab=1.5,cex.main=1.5,cex.axis =1.5,ylab="# Substitutions / # polymorphisms",col=c(pal[1],pal[2]))
text(x=bp,y=0-1.5/15,cex=1.5,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
text(x=bp, y=0+1.5/15, labels=paste0("n=", c(nrow(pos),nrow(nonPos))),cex=1.5)
# fisher exact test
fTest<-fisher.test(matrix(c(subNumb,polyNumb),nrow = 2,ncol = 2))
legend("topleft",legend=paste("p=",signif(fTest$p.value,3)),bty = 'n',cex=1.5)

#####* compare expression variation across populations betwween the putative target genes of positive and non-positive sites  *#####
## target genes of TFBS
TFBSgene<-fread("data/human_TFBS_target_genes/hsap_CEBPA_targetGene.txt")
TFBSgene<-TFBSgene[,c(1:3,7)]
colnames(TFBSgene)[c(1:4)]<-c("chr","start","end","geneID")
TFBSgene<-unique(TFBSgene)
## expression count
geneCount<-fread("data/human_gtex/getx_liver_gene_tpm.txt")
colnames(geneCount)[1]<-"geneID"
geneCount$geneID<-gsub("\\..*","",geneCount$geneID)
geneCount[,c(2:176)]<-log2(geneCount[,c(2:176)]+1)
geneCount$var<-apply(geneCount[,c(2:176)],1,function(x) var(x))
geneCount$mean<-apply(geneCount[,c(2:176)],1,function(x) mean(x))
deltaSVMGene<-merge(deltaSVM,TFBSgene,by=c("chr","start","end"))
deltaSVMGeneExp<-merge(deltaSVMGene,geneCount,by="geneID")
## positive and non-positive sites
pos<-subset(deltaSVMGeneExp,deltaSVMGeneExp$pValue<0.01)
nonPos<-subset(deltaSVMGeneExp,deltaSVMGeneExp$pValue>=0.01)
posFilter<-pos[!pos$geneID%in%nonPos$geneID,]
nonPosFilter<-nonPos[!nonPos$geneID%in%pos$geneID,]
pos<-posFilter
nonPos<-nonPosFilter
## plot
par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(pos$var,nonPos$var,xaxt = "n",
        ylim=c(-0.2,2),notch=T,pch=16,outcex=0.5,main="Human CEBPA TFBS",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
        ylab="Expression variance across populations",col=c(pal[1],pal[2]))
text(x=c(1,2),y=-0.36,cex=1.5,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
text(x=c(1,2), y=-0.2, cex=1.5,labels=paste0("n=", c(nrow(pos),nrow(nonPos))))

wtest<-wilcox.test(pos$var,nonPos$var)
legend("topleft",legend=paste("p=",signif(wtest$p.value,3)),bty = 'n',cex=1.5)

#####**** fly analysis ****#####
#####*** check the performance of trained model to distinguish binding sites and random sequences ***#####
## five fold cross validation 
cv<-fread("data/fly_SVM_model/dmel_CTCF_gkmtrain.cvpred.txt")
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
## plot
plot(perf,lwd = 3,cex=1.4)
# calculate AUC
auc_result <-performance( pred, measure = "auc")
unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,"(AUC = 0.911)",bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####*** deltaSVM summary ***#####
deltaSVMcons<-fread("data/fly_deltaSVM/dmel_CTCF_conserved_deltaSVM_highertailTest.txt")
deltaSVMgain<-fread("data/fly_deltaSVM/dmel_CTCF_gain_deltaSVM_highertailTest.txt")
deltaSVMloss<-fread("data/fly_deltaSVM/dmel_CTCF_loss_deltaSVM_lowertailTest.txt")
## change the format of deltaSVM files
deltaSVMcons<-dataMod(deltaSVMcons)
deltaSVMgain<-dataMod(deltaSVMgain)
deltaSVMloss<-dataMod(deltaSVMloss)
## plot deltaSVM and pvalue, use conserved binding sites as an example 
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
n=2
hist(deltaSVMcons$deltaSVM,breaks = 60,main="Conserved",xlab="deltaSVM",xlim=c(-20,20),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])
hist(deltaSVMcons$pValue,breaks = 60,main="Conserved",xlab="Pvalue",xlim=c(0,1),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[3])
## calculate proportion of positive selection
par(mfrow=c(1,1))
par(mar=c(8,7,4,2))
posCons<-deltaSVMcons[deltaSVMcons[,6]<0.01,]
posGain<-deltaSVMgain[deltaSVMgain[,6]<0.01,]
barplot(c(nrow(posCons)/nrow(deltaSVMcons),nrow(posGain)/nrow(deltaSVMgain)),
        ylim=c(0,0.5),col=pal[2],ylab="Proportion of positive binding sites",main="CTCF binding sites",
        cex.lab=2,cex.main=2,cex.axis=2)
text(x=c(1,2),y=-0.02,cex=2,srt = 45,adj = 1,labels = c("Conserved","Gain"),xpd = TRUE)

#####**** huamn CTCF adaptive evolution between tissues ****#####
#####*** check the performance of trained model to distinguish binding sites and random sequences ***#####
## five fold cross validation 
cv<-fread("data/human_CTCF_adaptation/human_SVM_model/all_merged_ctcf_gkmtrain.cvpred.txt")
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
## plot
plot(perf,lwd = 3,cex=1.4)
# calculate AUC
auc_result <-performance( pred, measure = "auc")
unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,"(AUC = 0.879)",bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####*** deltaSVM summary ***#####
deltaSVM<-fread("data/human_CTCF_adaptation/human_deltaSVM/ctcf_deltaSVM_highertailTest.txt")
deltaSVM<-dataMod(deltaSVM)
## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
n=1.5
hist(deltaSVM$deltaSVM,breaks = 200,main=" ",xlab="deltaSVM",xlim=c(-30,30),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])
hist(deltaSVM$pValue,breaks = 50,main=" ",xlab="Pvalue",xlim=c(0,1),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[3])

#####*** positive selection and pleitropy ***#####
organNumb<-fread("data/human_CTCF_adaptation/human_ctcf_binding/all_merged_annotated.bed")
## from the fourth column, if the value > 0, indicating this binding site is expressed in this tissue
colnames(organNumb)<-c("chr","start","end","adrenal gland","B cell","esophagus muscularis mucosa","retinal pigment epithelial cell",
                       "omental fat pad","gastrocnemius medialis","astrocyte of the cerebellum","astrocyte of the spinal cord",
                       "brain microvascular endothelial cell","choroid plexus epithelial cell","heart left ventricle","liver","upper lobe of left lung",
                       "neural cell (in vitro differentiated)","ovary", "pancreas","foreskin fibroblast","peyer patch","prostate gland","lower leg skin",
                       "spleen","stomach","testis","thoracic aorta","thyroid gland","tibial nerve","transverse colon","uterus","vagina")
organNumb$orgNumb<-apply(organNumb[,c(4:32)],1,function(x) length(x[x>0]))
deltaSVMorgNumb<-merge(deltaSVM,organNumb,by=c("chr","start","end"))
posTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue<0.01)
nonPosTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue>=0.01)
## plot
par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(posTFBS$orgNumb,nonPosTFBS$orgNumb,ylim=c(0,40),ylab="Number of tissues / cell types",main="Human CTCF binding sites",xaxt="n",
        col=c(pal[1],pal[2]),notch=T,pch=16,outcex=0.5,cex.lab=1.2,cex.main=1.2)
text(x=c(1,2),y=-2.5,cex.lab=1.2,srt = 45,adj = 1,labels = c("Positive sites", "Non-positive sites"),xpd = TRUE)
wTest<-wilcox.test(posTFBS$orgNumb,nonPosTFBS$orgNumb)
legend("topleft",legend=paste("p=",signif(wTest$p.value,3)),bty = 'n')

#####*** positive selection and different tissues ***#####
propPos<-c()
posNumb<-c()
allNumb<-c()
for (i in c(7:35)) {
  subData<-deltaSVMorgNumb[,c(1:6,i)]
  tissueBasedData<-subData[subData[,7]>0,]
  propPos[i-6]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])/nrow(tissueBasedData)
  posNumb[i-6]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])
  allNumb[i-6]<-nrow(tissueBasedData)
}
organProp<-data.frame("organ" = rep(NA,29), "prop" = rep(NA,29),"posNumb" = rep(NA,29),"allNumb" = rep(NA,29))
organProp$organ<-colnames(organNumb)[c(4:32)]
organProp$prop<-propPos
organProp$posNumb<-posNumb
organProp$allNumb<-allNumb
organProp<-organProp[order(organProp$prop),]
## fisher test for proportions
ftest<-data.frame("odds_ratio" = rep(NA,29), "pvalue" = rep(NA,29))
for (i in c(1:29)) {
  temp<-fisher.test(matrix(c(organProp[i,3],organProp[i,4] ,nrow(deltaSVMorgNumb[deltaSVMorgNumb$pValue<0.01,]),nrow(deltaSVMorgNumb)),nrow = 2,ncol = 2))
  ftest[i,1]<-temp$estimate
  ftest[i,2]<-temp$p.value
  
}
##plot
par(mfrow=c(1,1))
par(mar=c(15, 5, 2, 15) + 0.1)
plotColor<-c(pal[6],pal[2],"black",pal[10],pal[5],pal[4], pal[2],pal[2],pal[1],pal[4],pal[1],pal[7],pal[4],pal[7],
             pal[10],pal[8],pal[6],pal[9],pal[9],pal[9],pal[9],pal[10],pal[5],pal[1],pal[5],pal[1],pal[1],pal[1],pal[1])
plotLegend<-c("Nervous system","Male reproductive system","Immune system","Endocrine system","Integumentary system","Respiratory system",
              "Cardiovascular system","Digestive system","Female reproductive system","Skeletomuscular system")
plot(c(1:29),organProp$prop, ylab="Proportion of positive binding sites", main="Human CTCF binding sites",xlab="", 
     pch=16,cex=3,xaxt="n",cex.lab=1.5,cex.main=1.5, cex.axis=1.5,type="p",col=plotColor, ylim=c(0.02,0.04))
legend("topright",legend=plotLegend,col = unique(rev(plotColor)),cex=1.2,
       pch=rep(16,times=10),pt.cex=2,bty="n",inset=c(-0.38,0),xpd = TRUE)
abline(v=c(1:29),col="grey", lty=4)
axis(side = 1, at = c(1:29), labels=F)
text(c(1:29), 0.0185, srt = 45,cex=1.3, adj = 1,labels = organProp$organ, xpd = TRUE)

#####*** number of substitutions for non positive sites in different tissues ***#####
humanChimp_subsNumb<-fread("data/human_CTCF_adaptation/human_chimp_substitutions/humanChimp_subNumb.txt")
colnames(humanChimp_subsNumb)<-c("chr","start","end","subNumb")
deltaSVMorgNumb_subsNumb<-merge(deltaSVMorgNumb,humanChimp_subsNumb,by=c("chr","start","end"))
numbSub<-c()
for (i in c(7:35)) {
  subData<-deltaSVMorgNumb_subsNumb[,c(1:6,i,37)]
  tissueBasedData<-subData[subData[,7]>0,]
  nonPosData<-tissueBasedData[tissueBasedData$pValue>=0.01,]
  numbSub[i-6]<-mean(nonPosData$subNumb/(nonPosData$end-nonPosData$start))
}
organSub<-data.frame("organ" = rep(NA,29), "numb" = rep(NA,29))
organSub$organ<-colnames(organNumb)[c(4:32)]
organSub$numb<-numbSub
organSub<-organSub[order(organSub$numb),]

##plot
par(mfrow=c(1,1))
par(mar=c(15, 5, 2, 15) + 0.1)
plotColor<-c(pal[6],pal[2],"black",pal[4],pal[1],pal[2],pal[7],pal[1],pal[5],pal[10],
             pal[2],pal[1],pal[7],pal[5],pal[9],pal[9],pal[10],pal[1],pal[4],pal[9],
             pal[10],pal[8],pal[6],pal[4],pal[5],pal[9],pal[1], pal[1], pal[1])
plotLegend<-c("Nervous system","Endocrine system","Male reproductive system","Digestive system",
              "Integumentary system","Respiratory system","Immune system","Cardiovascular system",
              "Female reproductive system","Skeletomuscular system")
plot(c(1:29),organSub$numb, ylab="Average number of substitutions  \n Non-positive sites", main="Human CTCF binding sites",xlab="",
     pch=16,cex=3,xaxt="n",col=plotColor,cex.lab=1.5,cex.main=1.5, cex.axis=1.5, type="p", ylim=c(0.013,0.014))

legend("topright",legend=plotLegend,col = unique(rev(plotColor)),cex=1.2,
       pch=rep(16,times=10),pt.cex=2,bty="n",inset=c(-0.38,0),xpd = TRUE)
abline(v=c(1:29),col="grey", lty=4)
axis(side = 1, at = c(1:29), labels=F)
text(c(1:29), 0.0129, srt = 45,cex=1.3, adj = 1,labels = organSub$organ, xpd = TRUE)

#####**** mouse CTCF adaptive evolution between tissues ****#####
#####*** check the performance of trained model to distinguish binding sites and random sequences ***#####
## five fold cross validation 
cv<-fread("data/mouse_CTCF_adaptation/mouse_SVM_model/all_merged_ctcf__gkmtrain.cvpred.txt")
colnames(cv)<-c("position","prediction","real_state","cv_number")
pred <- prediction(cv$prediction, cv$real_state) 
perf <- performance( pred, "tpr", "fpr" )
## plot
plot(perf,lwd = 3,cex=1.4)
# calculate AUC
auc_result <-performance( pred, measure = "auc")
unlist(slot(auc_result, "y.values"))
legend(0.3,0.6,"(AUC = 0.919)",bty="n",lwd = 3,cex=1.4,col = c("black", "blue","white")) 

#####*** deltaSVM summary ***#####
deltaSVM<-fread("data/mouse_CTCF_adaptation/mouse_deltaSVM/ctcf_deltaSVM_highertailTest.txt")
deltaSVM<-dataMod(deltaSVM)

## plot deltaSVM and pvalue
par(mfrow=c(1,2))
par(mar=c(7,5,4,2))
n=1.5
hist(deltaSVM$deltaSVM,breaks = 50,main=" ",xlab="deltaSVM",xlim=c(-30,30),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[4])
hist(deltaSVM$pValue,breaks = 50,main=" ",xlab="Pvalue",xlim=c(0,1),
     cex.lab=n,cex.axis=n,cex.main=2,col=pal[3])

#####*** positive selection and pleitropy ***#####
organNumb<-fread("data/mouse_CTCF_adaptation/mouse_ctcf_binding/all_merged_annotated.bed")
## from the fourth column, if the value > 0, indicating this binding site is expressed in this tissue
colnames(organNumb)<-c("chr","start","end","bone marrow","cerebellum","cortical plate","heart","kidney","liver",
                       "lung","olfactory bulb","small intestine","testis","thymus")
organNumb$orgNumb<-apply(organNumb[,c(4:14)],1,function(x) length(x[x>0]))
deltaSVMorgNumb<-merge(deltaSVM,organNumb,by=c("chr","start","end"))
posTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue<0.01)
nonPosTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue>=0.01)
## plot
par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(posTFBS$orgNumb,nonPosTFBS$orgNumb,ylim=c(0,12),ylab="Number of tissues / cell types",main= "Mouse CTCF binding sites",xaxt="n",
        col=c(pal[1],pal[2]),notch=T,pch=16,outcex=0.5,cex.lab=1.2,cex.main=1.2)
text(x=c(1,2),y=-1,cex.lab=1.2,srt = 45,adj = 1,labels = c("Positive sites", "Non-positive sites"),xpd = TRUE)
wTest<-wilcox.test(posTFBS$orgNumb,nonPosTFBS$orgNumb)
legend("topleft",legend=paste("p=",signif(wTest$p.value,3)),bty = 'n')

#####*** positive selection and different tissues ***#####
propPos<-c()
posNumb<-c()
allNumb<-c()
for (i in c(7:17)) {
  subData<-deltaSVMorgNumb[,c(1:6,i)]
  tissueBasedData<-subData[subData[,7]>0,]
  propPos[i-6]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])/nrow(tissueBasedData)
  posNumb[i-6]<-nrow(tissueBasedData[tissueBasedData$pValue<0.01,])
  allNumb[i-6]<-nrow(tissueBasedData)
}
organProp<-data.frame("organ" = rep(NA,11), "prop" = rep(NA,11))
organProp$organ<-colnames(organNumb)[c(4:14)]
organProp$prop<-propPos
organProp$posNumb<-posNumb
organProp$allNumb<-allNumb
organProp<-organProp[order(organProp$prop),]
## plot
par(mar=c(15, 5, 2, 15) + 0.1)
plotColor<-c(pal[10],pal[1],pal[4],pal[1],pal[1],pal[5],pal[10],pal[7],pal[9],"blue",pal[8])
plotLegend<-c("Respiratory system","Excretory system","Endocrine system","Cardiovascular system","Immune system","Male reproductive system","Nervous system","Digestive system")
plot(c(1:11),organProp$prop, ylab="Proportion of positive binding sites", main="Mouse CTCF binding sites",xlab="", pch=16,cex=3,
     xaxt="n", type="p",col=plotColor,cex.lab=1.5,cex.main=1.5, cex.axis=1.5, ylim=c(0.02,0.04))
legend("topright",legend=plotLegend,col = unique(rev(plotColor)),cex=1.2,
       pch=rep(16,times=11),pt.cex=2,bty="n",inset=c(-0.38,0),xpd = TRUE)
abline(v=c(1:11),col="grey", lty=4)
axis(side = 1, at = c(1:11), labels=F)
text(c(1:11), 0.0185, srt = 45, cex=1.3,adj = 1,labels = organProp$organ, xpd = TRUE)

#####*** number of substitutions for non positive sites in different tissues ***#####
mouseSpret_subsNumb<-fread("data/mouse_CTCF_adaptation/mouse_spret_substitutions/mouseSpret_subNumb.txt")
colnames(mouseSpret_subsNumb)<-c("chr","start","end","subNumb")
deltaSVMorgNumb_subsNumb<-merge(deltaSVMorgNumb,mouseSpret_subsNumb,by=c("chr","start","end"))
numbSub<-c()
for (i in c(7:17)) {
  subData<-deltaSVMorgNumb_subsNumb[,c(1:6,i,19)]
  tissueBasedData<-subData[subData[,7]>0,]
  nonPosData<-tissueBasedData[tissueBasedData$pValue>=0.01,]
  numbSub[i-6]<-mean(nonPosData$subNumb/(nonPosData$end-nonPosData$start))
}
organSub<-data.frame("organ" = rep(NA,11), "numb" = rep(NA,11))
organSub$organ<-colnames(organNumb)[c(4:14)]
organSub$numb<-numbSub
organSub<-organSub[order(organSub$numb),]
## plot
plotColor<-c(pal[10],pal[1],pal[4],pal[5],pal[10],pal[7],pal[8],pal[1],"blue",pal[9],pal[1])
plotLegend<-c("Nervous system","Endocrine system","Excretory system","Respiratory system","Cardiovascular system",
              "Immune system","Male reproductive system","Digestive system")

par(mar=c(15, 5, 2, 15) + 0.1)
plot(c(1:11),organSub$numb, ylab="Average number of substitutions  \n Non-positive sites", main="Mouse CTCF binding sites",xlab="",
     pch=16,cex=3,xaxt="n",cex.lab=1.5,cex.main=1.5, cex.axis=1.5, type="p",ylim=c(0.01,0.014),col=plotColor)

legend("topright",legend=plotLegend,col = unique(rev(plotColor)),cex=1.2,
       pch=rep(16,times=10),pt.cex=2,bty="n",inset=c(-0.38,0),xpd = TRUE)
abline(v=c(1:11),col="grey", lty=4)
axis(side = 1, at = c(1:11), labels=F)
text(c(1:11), 0.0097, srt = 45, adj = 1,cex=1.3,labels = organSub$organ, xpd = TRUE)
