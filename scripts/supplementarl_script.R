#####** function: change the format of deltaSVM data **#####
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
  
  ## qvalue
  deltaSVM<-deltaSVM[order(deltaSVM$pValue),]
  deltaSVM$qValue<-qvalue(deltaSVM$pValue,pi0=1)$qvalues
  return(deltaSVM)
}


#####** Figure S1**#####
deltaSVMnullCons<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_conserved_deltaSVM_nullDis.txt")
deltaSVMnullCons<-deltaSVMnullCons[c(10001:20000),]
par(mfrow=c(1,1))
par(mar=c(5,5,4,2))
hist(deltaSVMnullCons$V2,breaks = 100,xlab="deltaSVM",main="Null distribution of one CEBPA conserved binding site \n(chr12: 73750975-73751322)")

#####** Figure S2**#####

nullDis<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_specific_loss_deltaSVM_nullDis.txt")
colnames(nullDis)<-c("ID","deltaSVM")

countLengt<-function(x) {return(length(x[x<0]))}
nullDisMean <- aggregate(deltaSVM ~ ID, FUN=mean, data=nullDis) 
splString<-strsplit(nullDisMean$ID,"_",fixed=TRUE)
splString<-data.frame(unlist(splString))
ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
ID.bed<-ID.bed[,c(1:3)]
ID.bed[,2]<-as.numeric(ID.bed[,2])-1
nullDisMean<-cbind(ID.bed, nullDisMean[,2])
colnames(nullDisMean)<-c("chr","start","end","mean")
nullDisMean<-data.frame(nullDisMean)


deltaSVMcons<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_conserved_deltaSVM_highertailTest.txt")
deltaSVMgain<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_specific_gain_deltaSVM_highertailTest.txt")
deltaSVMloss<-fread("data/mouse_deltaSVM/CEBPA/bl6-CEBPA_specific_loss_deltaSVM_lowertailTest.txt")

## take deltaSVMloss as an example
deltaSVMloss<-dataMod(deltaSVMloss)

## positive 
a<-subset(deltaSVMloss,deltaSVMloss$pValue<0.01)
a<-merge(a,nullDisMean,by=c("chr","start","end"))
## negative 
b<-subset(deltaSVMloss,deltaSVMloss$pValue>=0.01)
b<-merge(b,nullDisMean,by=c("chr","start","end"))

boxplot(as.numeric(as.character(a$mean)),as.numeric(as.character(b$mean)))
wilcox.test(as.numeric(as.character(a$mean)),as.numeric(as.character(b$mean)))
par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(as.numeric(as.character(a$mean)),as.numeric(as.character(b$mean)),xaxt = "n",
        ylim=c(-6,2),notch=T,pch=16,outcex=0.5,main="Lineage specific loss",cex.lab=1.5,cex.main=1.5,cex.axis=1.5,
        ylab="Mean deltaSVM of null distribution",col=c(pal[1],pal[2]))
text(x=c(1,2),y=-6.5,cex=1.5,srt = 45,adj = 1,labels = c("Positive sites","Non-positive sites"),xpd = TRUE)
text(x=c(1,2), y=-5.5, cex=1.5,labels=paste0("n=", c(nrow(a),nrow(b))))

wtest<-wilcox.test(as.numeric(as.character(a$mean)),as.numeric(as.character(b$mean)))
legend("topleft",legend=paste("p=",signif(wtest$p.value,3)),bty = 'n',cex=1.5)

#####** Figure S3, S4, S5, S6, S8 **#####
## please check the code of Figure 2 

#####** Figure S7 **#####
## please check the code of Figure 3

#####** Figure S9 **#####
## please check the code of Figure 2B

#####** Figure S10, S11 **#####
## please check the code of Figure 3C and 3D respectively

#####** Figure S12 **#####
deltaSVM<-fread("data/human_CTCF_adaptation/human_deltaSVM/ctcf_deltaSVM_highertailTest.txt")
deltaSVM<-dataMod(deltaSVM)
organNumb<-fread("data/human_CTCF_adaptation/human_ctcf_binding/all_merged_annotated.bed")
colnames(organNumb)<-c("chr","start","end","adrenal gland","B cell","esophagus muscularis mucosa","retinal pigment epithelial cell",
                       "omental fat pad","gastrocnemius medialis","astrocyte of the cerebellum","astrocyte of the spinal cord",
                       "brain microvascular endothelial cell","choroid plexus epithelial cell","heart left ventricle","liver","upper lobe of left lung",
                       "neural cell (in vitro differentiated)","ovary", "pancreas","foreskin fibroblast","peyer patch","prostate gland","lower leg skin",
                       "spleen","stomach","testis","thoracic aorta","thyroid gland","tibial nerve","transverse colon","uterus","vagina")
organNumb$orgNumb<-apply(organNumb[,c(4:32)],1,function(x) length(x[x>0]))
deltaSVMorgNumb<-merge(deltaSVM,organNumb,by=c("chr","start","end"))
posTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue<0.01)
nonPosTFBS<-subset(deltaSVMorgNumb,deltaSVMorgNumb$pValue>=0.01)

organName<-colnames(organNumb)[c(4:32)]
systemName<-c("Endocrine system","Immune system","Digestive system","Nervous system","Integumentary system","Skeletomuscular system",
              "Nervous system","Nervous system","Nervous system","Nervous system","Cardiovascular system","Endocrine system",
              "Respiratory system","Nervous system","Female reproductive system","Endocrine system","Male reproductive system","Immune system",
              "Male reproductive system","Integumentary system","Immune system","Digestive system","Male reproductive system","Cardiovascular system",
              "Endocrine system","Nervous system","Digestive system","Female reproductive system","Female reproductive system")

SystemColor<-cbind(c("Nervous system","Male reproductive system","Immune system","Endocrine system","Integumentary system","Respiratory system",
                     "Cardiovascular system","Digestive system","Female reproductive system","Skeletomuscular system"),
                   c(1,5,10,9,6,8,7,4,2,"black")
)
colnames(SystemColor)<-c("systemName","colorname")

organSystem<-cbind(organName,systemName)
organSystemColorName<-merge(organSystem,SystemColor,by="systemName")
organSystemColorName<-data.frame(organSystemColorName)

par(mfrow=c(1,1))
par(mar=c(8,5,4,2))
boxplot(posTFBS$orgNumb,nonPosTFBS$orgNumb,ylim=c(0,40),ylab="Number of tissues / cell types",main="Human CTCF binding sites",xaxt="n",
        col=c(pal[1],pal[2]),notch=T,pch=16,outcex=0.5,cex.lab=1.2,cex.main=1.2)
text(x=c(1,2),y=-2.5,cex.lab=1.2,srt = 45,adj = 1,labels = c("Positive sites", "Non-positive sites"),xpd = TRUE)
wTest<-wilcox.test(posTFBS$orgNumb,nonPosTFBS$orgNumb)
legend("topleft",legend=paste("p=",signif(wTest$p.value,3)),bty = 'n')

#####** Figure S13 and  S14 **#####
## please check the code of Figure 4

#####** Figure S15 **#####
deltaSVM<-fread("data/human_CTCF_adaptation/human_deltaSVM/ctcf_each_substitution_deltaSVM.txt")
## change the first column into bed format
deltaSVM<-as.data.frame(deltaSVM)
names(deltaSVM)<-c("ID","position","upRef","ref","nextRef","sub","deltasvm","contentA","contentT","contentC","contentG","contentCG")
## "ID" is the peak ID, "position" is which nucleotide in the peak with substitution, 
## "ref" is the ancestor sequence, "sub" is the substitution sequence, "upRef" is the upstream sequence of "ref",
## "nextRef" is the downstream sequence of  "ref", "deltasvm" is the predicted affinity change of the substitution,
## "contentA","contentT","contentC","contentG" is the number of A, T, C, G in a peak 
## "contentCG" is the number of CG in a peak

splString<-strsplit(deltaSVM$ID,"_",fixed=TRUE)
splString<-data.frame(unlist(splString))
ID.bed<-matrix(splString$unlist.splString., ncol=5, byrow=TRUE)
ID.bed<-ID.bed[,c(1:3)]
ID.bed[,2]<-as.numeric(ID.bed[,2])-1
deltaSVM<-cbind(ID.bed, deltaSVM[,c(2:12)])
colnames(deltaSVM)[c(1:3)]<-c("chr","start","end")
deltaSVM$start<-as.numeric(as.character(deltaSVM$start))
deltaSVM$end<-as.numeric(as.character(deltaSVM$end))

## calculate the number of A T C G
tempdata<-deltaSVM[,c(1:3,10:14)]
tempdata<-unique(tempdata)
## remove C and G of CpG from content C and content G
tempdata$contentC<-tempdata$contentC-tempdata$contentCG
tempdata$contentG<-tempdata$contentG-tempdata$contentCG

sumContent<-colSums(tempdata[,c(4:8)])

## number of 12 possible substitutions (not for CpG sequences)
deltaSVM_nonCpG<-subset(deltaSVM, !(deltaSVM$upRef=="C" & deltaSVM$ref=="G") & !(deltaSVM$ref=="C" & deltaSVM$nextRef=="G"))

# AC
deltaAC<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="A" & deltaSVM_nonCpG$sub=="C")
# AG
deltaAG<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="A" & deltaSVM_nonCpG$sub=="G")
# AT
deltaAT<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="A" & deltaSVM_nonCpG$sub=="T")

# TA
deltaTA<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="T" & deltaSVM_nonCpG$sub=="A")
# TC
deltaTC<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="T" & deltaSVM_nonCpG$sub=="C")
# TG
deltaTG<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="T" & deltaSVM_nonCpG$sub=="G")

# CA
deltaCA<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="C" & deltaSVM_nonCpG$sub=="A")
# CG
deltaCG<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="C" & deltaSVM_nonCpG$sub=="G")
# CT
deltaCT<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="C" & deltaSVM_nonCpG$sub=="T")

# GA
deltaGA<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="G" & deltaSVM_nonCpG$sub=="A")
# GC
deltaGC<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="G" & deltaSVM_nonCpG$sub=="C")
# GT
deltaGT<-subset(deltaSVM_nonCpG, deltaSVM_nonCpG$ref=="G" & deltaSVM_nonCpG$sub=="T")

# calculate mutation rate and affinity change

mutationType_nonCG<-c(nrow(deltaAC),nrow(deltaAG),nrow(deltaAT),
                      nrow(deltaTA),nrow(deltaTC),nrow(deltaTG),
                      nrow(deltaCA),nrow(deltaCG),nrow(deltaCT),
                      nrow(deltaGA),nrow(deltaGC),nrow(deltaGT))
median_affinity_nonCG<-lapply(affinity_change_nonCG, median)

## for CpG sequences
deltaSVM_CpG<-subset(deltaSVM, (deltaSVM$upRef=="C" & deltaSVM$ref=="G") | (deltaSVM$ref=="C" & deltaSVM$nextRef=="G"))

# CA
deltaCG_CA<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="C" & deltaSVM_CpG$sub=="A")
# CG
deltaCG_CG<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="C" & deltaSVM_CpG$sub=="G")
# CT
deltaCG_CT<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="C" & deltaSVM_CpG$sub=="T")

# GA
deltaCG_GA<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="G" & deltaSVM_CpG$sub=="A")
# GC
deltaCG_GC<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="G" & deltaSVM_CpG$sub=="C")
# GT
deltaCG_GT<-subset(deltaSVM_CpG, deltaSVM_CpG$ref=="G" & deltaSVM_CpG$sub=="T")

# calculate mutation rate and affinity change 
mutationType_CG<-c(nrow(deltaCG_CA),nrow(deltaCG_CG),nrow(deltaCG_CT),
                   nrow(deltaCG_GA),nrow(deltaCG_GC),nrow(deltaCG_GT))

meadin_affinity_CG<-lapply(affinity_change_CG, median)

## put everything together
library(ggplot2)

mutationRate<-data.frame(c(mutationRate_nonCG, mutationRate_CG))

median_affinity_nonCG<-lapply(affinity_change_nonCG, median)
median_affinity_CG<-lapply(affinity_change_CG, median)
affinity_change<-data.frame(c(as.numeric(median_affinity_nonCG),as.numeric(median_affinity_CG)))

mutationAffinity<-cbind(mutationRate,affinity_change)
colnames(mutationAffinity)<-c("subRate","deltaSVM")
mutationAffinity$type <-c(rep("A",12),rep("B",6))
mutationAffinity$name<-labName<-c("A->C","A->G","A->T","T->A","T->C","T->G",
                                  "C->A","C->G","C->T","G->A","G->C","G->T",
                                  "C->A","C->G","C->T","G->A","G->C","G->T")


myplot<-ggplot(mutationAffinity, aes(x=subRate, y=deltaSVM, group=type,label=name)) +
  geom_point(aes(shape=type, color=type), size=3)+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  scale_size_manual(values=c(3,4))+
  ylim(-1,1)+xlim(-0.002,0.025)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.key = element_rect(colour = NA, fill = NA),legend.position="top",
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_text()

print(myplot)

#####** FIGURE 16 **#####
## plot biding affinity change

affinity_change_nonCG<-list(deltaAC$deltasvm,deltaAG$deltasvm,deltaAT$deltasvm,
                            deltaTA$deltasvm,deltaTC$deltasvm,deltaTG$deltasvm,
                            deltaCA$deltasvm,deltaCG$deltasvm,deltaCT$deltasvm,
                            deltaGA$deltasvm,deltaGC$deltasvm,deltaGT$deltasvm)

par(mar=c(4, 5, 2, 1) + 0.1)
labName<-c("A->C","A->G","A->T","T->A","T->C","T->G","C->A","C->G","C->T","G->A","G->C","G->T")
bp<-boxplot(affinity_change_nonCG,ylim=c(-50,50),ylab="deltaSVM",main="Without CpG sequence",xaxt="n",
            notch=T,pch=16,outcex=0.3,cex.lab=1.2,cex.main=1.2)
text(x=c(1:12),y=-55,cex.lab=1.2,srt = 45,adj = 1,labels =labName ,xpd = TRUE)

affinity_change_CG<-list(deltaCG_CA$deltasvm,deltaCG_CG$deltasvm,deltaCG_CT$deltasvm,
                         deltaCG_GA$deltasvm,deltaCG_GC$deltasvm,deltaCG_GT$deltasvm)

labName<-c("C->A","C->G","C->T","G->A","G->C","G->T")
bp<-boxplot(affinity_change_CG,ylim=c(-50,50),ylab="deltaSVM",main="CpG sequence",xaxt="n",
            notch=T,pch=16,outcex=0.3,cex.lab=1.2,cex.main=1.2)
text(x=c(1:6),y=-55,cex.lab=1.2,srt = 45,adj = 1,labels =labName ,xpd = TRUE)

