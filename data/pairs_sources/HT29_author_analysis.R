
###This is the code for the adhoc figure generation, paths for the figure location as well as input must be modified
###All input needed for this code is coming from the 'PILOT' sub folder system or from the 'ENCORE/input', 
###both paths should we modified as needed when running locally.
###RColorBrewer together with ineq, pROC and flextablelibraries  are needed to run these lines

###Load library and function for cohen D
library("pROC")
library("RColorBrewer")
library("ineq")
library('flextable')


dcohen<-function(set1,set2){
  func=(mean(set1,na.rm=T)-mean(set2,na.rm=T))/sqrt((sd(set1,na.rm=T)^2+sd(set2,na.rm=T)^2)/2)
  return(func)
}

####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE1
####### ####### ####### ####### ####### ####### ####### ####### 
####### FIG1F
stat_pil=as.matrix(read.delim('/PILOT/COUNTS/stats/stats_all.csv',sep=','))

count=as.numeric(stat_pil[stat_pil[,'Sample']=='PILOT_Plasmid' & 
                            stat_pil[,'type']%in%c("sgrna1_exists","sgrna2_exists","sgrna_match","swap","OEM1_exact","OEM2_exact"),'value'])
names(count)=stat_pil[stat_pil[,'Sample']=='PILOT_Plasmid' & 
                        stat_pil[,'type']%in%c("sgrna1_exists","sgrna2_exists","sgrna_match","swap","OEM1_exact","OEM2_exact"),'type']
count=count[c("sgrna1_exists","sgrna2_exists","sgrna_match","OEM1_exact","OEM2_exact","swap")]

pdf('/pathFortheFigure/Fig1F_stats_v2.pdf',width=5,height=5)
x=barplot(count,ylab='% total reads',las=3,names=c("sgrna1","sgrna2","sgrna1+sgrna2","sgrna1_only","sgrna2_only","swap"),
          ylim=c(0,100),border=NA)

text(x,count,round(count,2))
dev.off()

####### FIG1G
## path_stat_plasmid -> path where the stats for the COLO1 plasmid are located, in the github rep thwy are under ANALYSIS/ENCORE/input
#load table and modify 
stat=as.matrix(read.delim('path_stat_plasmid/lib-COLO-1.stats',sep=' '))
temp=stat[,1]
stat=apply(stat[,2:4],2,as.numeric)
rownames(stat)=temp

vect=c(stat['MATCH','lib.COLO.1.2']*100,
       stat['OEM1','lib.COLO.1.2']*100,
       stat['OEM2','lib.COLO.1.2']*100,
       stat['SWAP','lib.COLO.1.2']*100)

pdf('/pathFortheFigure/Fig1G_COLO1_plasmid_reads.pdf',width=5,height=5)

x=barplot(c(vect,0,0),ylim=c(0,100),ylab='GI library reads (% of reads)',
          names=c('gRNA1+gRNA2',"sgrna1_only","sgrna2_only",'swap','',''),las=3,border=NA)
text(x[,1],c(vect,0,0),round(c(vect,0,0),2),pos=4,srt=90)

dev.off()

###### FIG1-E
#Counts for the pilot plasmid
tab_exact=as.matrix(read.delim('PILOT/MERGED/PILOT_EXACT_COUNTS.txt'))

counts=apply(tab_exact[,9:ncol(tab_exact)],2,as.numeric)
rownames(counts)=tab_exact[,'ID']
gini_plm=counts[,'PILOT_Plasmid']

#package to calculate GINI index
library("ineq")

ineq::Gini(log(counts[,'PILOT_Plasmid']+1))

###lines for the curve

plm=cbind(counts[counts[,'PILOT_Plasmid']>0,'PILOT_Plasmid'],0,0)
plm=plm[order(plm[,1]),]
colnames(plm)=c('ori','read','sgRNA')

plm[,'read']=cumsum(plm[,'ori'])/sum(plm[,'ori'])
plm[,'sgRNA']=1:nrow(plm)
plm[,'sgRNA']=plm[,'sgRNA']/max(plm[,'sgRNA'])

pdf('/pathFortheFigure/Fig1E_plasmid1_reads.pdf',width=5,height=5)

plot(plm[,'sgRNA'],plm[,'read'],type='l',col=rgb(0.2,0.7,0.4),lwd=2,
     xlab='Fraction of total sgRNAs',ylab='Fraction of total reads')

points(c(0,1),c(0,1),type='l')
dev.off()


####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE2
####### ####### ####### ####### ####### ####### ####### ####### 

###### FIG2-A

anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

lfc_plm=readRDS('PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

#####Average replicas

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

pal1=c(rgb(108/255,116/255,134/255),rgb(105/255,156/255,201/255),rgb(205/255,232/255,236/255))

##Selection of positive and negative 
pos=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic')
neg=c('non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')

##Cohen D and t test calculations for the plot 
d_cohen_text=c(dcohen( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,2],
                mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,2]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,2],
                mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,2]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,2],
                mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,2]))
d_cohen_text=round(d_cohen_text,2)

x=c(t.test( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,2])$p.value,
    t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,2],
           mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,2])$p.value,
    t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,2],
           mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,2])$p.value)


pdf('pathFortheFigure/Fig2A_boxplot_500PX_D3.pdf',width=4,height=5)
boxplot(mat[anot[,'Notes']!='DistanceCut',2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,2],
        col=pal1[c(1,2,3,2,3,2,3)],xlab='sgRNA Scaffold',ylab='log2FC',las=3)

text(c(2.5,4.5,6.5),rep(2.5,3),d_cohen_text)
text(c(2.5,4.5,6.5),rep(3,3),x)

dev.off()

###### FIG2-C

####Aggregate gen level

sel1=anot[,'Scaffold']=='Wildtype'
sel2=anot[,'Scaffold']=='Modified6'
sel3=anot[,'Scaffold']=='Modified7'

temp=data.frame(G1=anot[sel1,'Gene1'],G2=anot[sel1,'Gene2'],lfc=mat[sel1,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2, data=temp, mean)
temp=as.matrix(temp)
wt=temp
temp=data.frame(G1=anot[sel2,'Gene1'],G2=anot[sel2,'Gene2'],lfc=mat[sel2,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m6=temp

temp=data.frame(G1=anot[sel3,'Gene1'],G2=anot[sel3,'Gene2'],lfc=mat[sel3,'PILOT_500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m7=temp

rownames(wt)=paste(wt[,1],wt[,2],sep='_')
temp=wt[,'lfc']
names(temp)=paste(wt[,2],wt[,1],sep='_')
wt=wt[rownames(wt)%in%names(temp),]
wt=cbind(wt,temp[rownames(wt)])

rownames(m6)=paste(m6[,1],m6[,2],sep='_')
temp=m6[,'lfc']
names(temp)=paste(m6[,2],m6[,1],sep='_')
m6=m6[rownames(m6)%in%names(temp),]
m6=cbind(m6,temp[rownames(m6)])

rownames(m7)=paste(m7[,1],m7[,2],sep='_')
temp=m7[,'lfc']
names(temp)=paste(m7[,2],m7[,1],sep='_')
m7=m7[rownames(m7)%in%names(temp),]
m7=cbind(m7,temp[rownames(m7)])

wt_text=wt[,1:2]
m6_text=m6[,1:2]
m7_text=m7[,1:2]

wt=cbind(as.numeric(wt[,3]),as.numeric(wt[,4]))
colnames(wt)=c('sgRNA1','sgRNA2')
m6=cbind(as.numeric(m6[,3]),as.numeric(m6[,4]))
colnames(m6)=c('sgRNA1','sgRNA2')
m7=cbind(as.numeric(m7[,3]),as.numeric(m7[,4]))
colnames(m7)=c('sgRNA1','sgRNA2')

###Remove redundancy

unique_ID=t(combn(unique(c(wt_text[,1],wt_text[,2])),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_wt=!rownames(wt_text)%in%unique_ID

unique_ID=t(combn(unique(c(m6_text[,1],m6_text[,2])),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_m6=!rownames(m6_text)%in%unique_ID

unique_ID=t(combn(unique(c(m7_text[,1],m7_text[,2])),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_m7=!rownames(m7_text)%in%unique_ID


###WT
data <- data.frame(sgRNA1 = wt[,1], sgRNA2 = wt[,2])
# Fit linear model
reg<-lm(sgRNA1 ~ sgRNA2, data = data)
# Scatter plot

pdf('pathFortheFigure/Fig2C_WT_500PX_D3.pdf',width=15,height=5)
par(mfrow=c(1,3))

data <- data.frame(sgRNA1 = wt[sel_wt,1], sgRNA2 = wt[sel_wt,2])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='Wildtype (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[sel_m6,1], sgRNA2 = m6[sel_m6,2])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='M6 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[sel_m7,1], sgRNA2 = m7[sel_m7,2])
reg<-lm(sgRNA1 ~ sgRNA2, data = data)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='M7 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))
dev.off()

###### FIG2-D and supplementary 2E 
###### FIG2-D corresponds to the 2 row of plots, supplementary2E is be the full set

######Selection of guides and colors

selection=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic',
            'non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')
colores=c("#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#3182bd", "#393b79","#756bb1")

###vectors to store values


AUCS_wt=rep(0,length(selection))
names(AUCS_wt)=selection
AUCS_m6=rep(0,length(selection))
names(AUCS_m6)=selection
AUCS_m7=rep(0,length(selection))
names(AUCS_m7)=selection
AUCS_ze=rep(0,length(selection))
names(AUCS_ze)=selection

AUCS_max=rep(0,9)
names(AUCS_max)=c('wt_100X','m6_100X','m7_100X',
                  'wt_500X','m6_500X','m7_500X',
                  'wt_PCR500X','m6_PCR500X','m7_PCR500X')

n_wt=AUCS_wt
n_m6=AUCS_m6
n_m7=AUCS_m7

n_max=AUCS_max

####RECALL CURVES

pdf('pathFortheFigure/Fig2Dsup_recall_D3_ROC_v2.pdf',width=10,height=10)
par(mfrow=c(3,3))

################ 100X

### WT
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_100X_D14_WT')


for (i in 1:4){
  stype_df=sort(mat[scaf1,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}


for (i in 5:7){
  stype_df=sort(mat[scaf1,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['wt_100X']=x$auc
n_max['wt_100X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

#####M6
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_100X_D14_M6')

for (i in 1:4){
  stype_df=sort(mat[scaf2,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf2,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m6_100X']=x$auc
n_max['m6_100X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

####M7

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_100X_D14_M7')

for (i in 1:4){
  stype_df=sort(mat[scaf3,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf3,'PILOT_100X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m7_100X']=x$auc
n_max['m7_100X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 


tab1=cbind(round(AUCS_wt,2),n_wt,
           round(AUCS_m6,2),n_m6,
           round(AUCS_m7,2),n_m7)
colnames(tab1)=c("AUC-WT", "n-WT", "AUC-M6", "n-M6", "AUC-M7", "n-M7" )

##################
################ 500X

### WT
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_WT')


for (i in 1:4){
  stype_df=sort(mat[scaf1,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}


for (i in 5:7){
  stype_df=sort(mat[scaf1,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['wt_500X']=x$auc
n_max['wt_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

#####M6
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_M6')

for (i in 1:4){
  stype_df=sort(mat[scaf2,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf2,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m6_500X']=x$auc
n_max['m6_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

####M7

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_M7')

for (i in 1:4){
  stype_df=sort(mat[scaf3,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf3,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m7_500X']=x$auc
n_max['m7_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 


tab2=cbind(round(AUCS_wt,2),n_wt,
           round(AUCS_m6,2),n_m6,
           round(AUCS_m7,2),n_m7)
colnames(tab2)=c("AUC-WT", "n-WT", "AUC-M6", "n-M6", "AUC-M7", "n-M7" )

##################
################ PCR500X

### WT
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_PCR500X_D14_WT')


for (i in 1:4){
  stype_df=sort(mat[scaf1,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}


for (i in 5:7){
  stype_df=sort(mat[scaf1,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['wt_PCR500X']=x$auc
n_max['wt_PCR500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

#####M6
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_PCR500X_D14_M6')

for (i in 1:4){
  stype_df=sort(mat[scaf2,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf2,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m6_PCR500X']=x$auc
n_max['m6_PCR500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

####M7

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_PCR500X_D14_M7')

for (i in 1:4){
  stype_df=sort(mat[scaf3,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf3,'PILOT_PCR500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m7_PCR500X']=x$auc
n_max['m7_PCR500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 


tab3=cbind(round(AUCS_wt,2),n_wt,
           round(AUCS_m6,2),n_m6,
           round(AUCS_m7,2),n_m7)
colnames(tab3)=c("AUC-WT", "n-WT", "AUC-M6", "n-M6", "AUC-M7", "n-M7" )

dev.off()

###Printing tabels with numbers 


ft <- flextable(data.frame(rbind(c(AUCS_max['wt_100X'],n_max['wt_100X'],
                                   AUCS_max['m6_100X'],n_max['m6_100X'],
                                   AUCS_max['m7_100X'],n_max['m7_100X']),tab1)))
ft <- add_header_row(ft,
                     colwidths = c(2,2,2),
                     values = c("Wildtype", "M6",'M7'))
ft <-   valign(ft,valign = "center", part = "header")

ft <- labelizor(
  x = ft, 
  part = "header", 
  labels = c("wt_100X" = "AUC", "wt_100X.1" = "n",
             "m6_100X" = "AUC", "m6_100X.1" = "n",
             "m7_100X" = "AUC", "m7_100X.1" = "n"))
ft1=ft


ft <- flextable(data.frame(rbind(c(AUCS_max['wt_500X'],n_max['wt_500X'],
                                   AUCS_max['m6_500X'],n_max['m6_500X'],
                                   AUCS_max['m7_500X'],n_max['m7_500X']),tab2)))
ft <- add_header_row(ft,
                     colwidths = c(2,2,2),
                     values = c("Wildtype", "M6",'M7'))
ft <-   valign(ft,valign = "center", part = "header")

ft <- labelizor(
  x = ft, 
  part = "header", 
  labels = c("wt_500X" = "AUC", "wt_500X.1" = "n",
             "m6_500X" = "AUC", "m6_500X.1" = "n",
             "m7_500X" = "AUC", "m7_500X.1" = "n"))
ft2=ft

ft <- flextable(data.frame(rbind(c(AUCS_max['wt_PCR500X'],n_max['wt_PCR500X'],
                                   AUCS_max['m6_PCR500X'],n_max['m6_PCR500X'],
                                   AUCS_max['m7_PCR500X'],n_max['m7_PCR500X']),tab3)))
ft <- add_header_row(ft,
                     colwidths = c(2,2,2),
                     values = c("Wildtype", "M6",'M7'))
ft <-   valign(ft,valign = "center", part = "header")

ft <- labelizor(
  x = ft, 
  part = "header", 
  labels = c("wt_PCR500X" = "AUC", "wt_PCR500X.1" = "n",
             "m6_PCR500X" = "AUC", "m6_PCR500X.1" = "n",
             "m7_PCR500X" = "AUC", "m7_PCR500X.1" = "n"))
ft3=ft

pdf('/Fig2Dsup_recall_D3_tab100X_v2.pdf',width=6,height=6)
plot(ft1, fit = "fixed", just = "center",main='100X')
dev.off()
pdf('/Fig2Dsup_recall_D3_tab500X_v2.pdf',width=6,height=6)
plot(ft2, fit = "fixed", just = "center",main='500X')
dev.off()
pdf('/Fig2Dsup_recall_D3_tabPCR500X_v2.pdf',width=6,height=6)
plot(ft3, fit = "fixed", just = "center",main='PCR500X')
dev.off()

####### FIG2B

##########SCORE comparaisson

score_ag=readRDS('/input/SCORE_median_lfc.rds')

lfc=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')
anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

mat=cbind(apply(lfc[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

###Getting genes ready

anot=cbind(anot,'')
colnames(anot)[ncol(anot)]='GEN'

sel=(anot[,'sgRNA1_Class']=='essential' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA1_Class']=='non-essential' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA1_Class']=='unclassified' & anot[,'sgRNA2_Class']%in%c('intergenic','non-targeting'))

anot[sel,'GEN']=anot[sel,'Gene1']

sel=(anot[,'sgRNA2_Class']=='essential' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA2_Class']=='non-essential' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting')) |
  (anot[,'sgRNA2_Class']=='unclassified' & anot[,'sgRNA1_Class']%in%c('intergenic','non-targeting'))

anot[sel,'GEN']=anot[sel,'Gene2']

#######
sel=anot[,'Scaffold']=='Wildtype' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''

scaff_wt=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_wt=as.matrix(aggregate(. ~ gen , data=scaff_wt, mean))

sel=anot[,'Scaffold']=='Modified6' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''
scaff_m6=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_m6=as.matrix(aggregate(. ~ gen , data=scaff_m6, mean))

sel=anot[,'Scaffold']=='Modified7' & anot[,'Scaffold']!='DistanceCut' & anot[,'GEN']!=''
scaff_m7=data.frame(gen=anot[sel,'GEN'],
                    P_100X_D14=mat[sel,1],
                    P_500X_D14=mat[sel,2],
                    P_PCR500X_D14=mat[sel,3])
scaff_m7=as.matrix(aggregate(. ~ gen , data=scaff_m7, mean))

temp=scaff_wt[,1]
scaff_wt=apply(scaff_wt[,2:4],2,as.numeric)
rownames(scaff_wt)=temp

temp=scaff_m6[,1]
scaff_m6=apply(scaff_m6[,2:4],2,as.numeric)
rownames(scaff_m6)=temp

temp=scaff_m7[,1]
scaff_m7=apply(scaff_m7[,2:4],2,as.numeric)
rownames(scaff_m7)=temp

score_ag=score_ag[names(score_ag)%in%c(rownames(scaff_wt),rownames(scaff_m6),rownames(scaff_m7))]

scaff_wt=cbind(scaff_wt,score_ag[rownames(scaff_wt)])
scaff_m6=cbind(scaff_m6,score_ag[rownames(scaff_m6)])
scaff_m7=cbind(scaff_m7,score_ag[rownames(scaff_m7)])

scaff_wt=na.omit(scaff_wt)
scaff_m6=na.omit(scaff_m6)
scaff_m7=na.omit(scaff_m7)

##500X Wildtype scaffold
data <- data.frame(DUAL = scaff_wt[,'P_500X_D14'], SCO = scaff_wt[,4])
# Fit linear model
reg<-lm(SCO ~ DUAL, data = data)

# Scatter plot

pdf('pathFortheFigure/Fig2B_WT_500PX_D3.pdf',width=5,height=5)

plot(data,xlim=c(-6,1),ylim=c(-6,1),
     main='Gene level correlation single KO',xlab='DualKO (log2FC)',ylab='SingleKO-SCORE (log2FC)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
text(-6,-1,paste('cor:',round(cor(data[,1],data[,2]),2)),pos=4)
text(-6,-1.5,paste('reg:',round(coef(reg)[2],2)),pos=4)
text(-6,-2,paste('n:',nrow(data)),pos=4)

dev.off()

####### ####### ####### ####### ####### ####### ####### ####### 
############## Fig2E
#reload anotation and project score table

anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']
anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

pos=unique(c(anot[anot[,'sgRNA1_Class']=='essential','Gene1'],
             anot[anot[,'sgRNA2_Class']=='essential','Gene2']))
neg=unique(c(anot[anot[,'sgRNA1_Class']=='non-essential','Gene1'],
             anot[anot[,'sgRNA2_Class']=='non-essential','Gene2']))

score_ag=readRDS('/inpout/SCORE_median_lfc.rds')

pdf('pathFortheFigure/Fig2E.pdf',width=6,height=6)

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='HT29-project Score')
stype_df=sort(score_ag)
binary=rep(0,length(stype_df))
binary[names(stype_df)%in%pos]=1
x=roc(	binary,stype_df,direction=">")
 points(x$specificities,x$sensitivities,type='l',col=rgb(0.75,0.2,0.2))

dev.off()
############
####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for FIGURE3
####### ####### ####### ####### ####### ####### ####### ####### 
############## Fig3A
##numeric values
colo1=readRDS('/input/colo1_plm.rds')
##character values
IDs=readRDS('/input/IDs_colo1.rds')
colo1=colo1[,grep('136',colnames(colo1))]
colo1=apply(colo1[,1:3],1,mean)
##Selecting conditions
sel1=c('AnchorCombinations','LibraryCombinations','GIControlsCombinations')
sel2=c('AnchorSingletons','LibrarySingletons','GIControlsSingletons')
sel3=c('ESSENTIAL-INTERGENIC','ESSENTIAL-NONTARGET',
       'INTERGENIC-INTERGENIC','NONESSENTIAL-NONESSENTIAL','NONTARGET-NONTARGET')

pdf('pathFortheFigure/Fig3A_boxplot_ENCORE_mean_v2.pdf',width=4,height=5)
boxplot(colo1[IDs[,'MyNote']%in%sel1],colo1[IDs[,'MyNote']%in%sel2],colo1[IDs[,'MyNote']==sel3[1]],
        colo1[IDs[,'MyNote']==sel3[2]],colo1[IDs[,'MyNote']==sel3[3]],colo1[IDs[,'MyNote']==sel3[4]],
        colo1[IDs[,'MyNote']==sel3[5]],names=c('DoubleKO','SingleKO',sel3),las=3,ylab='ENCORE-HT19 log(FC)',col='white')
dev.off()

dcohen(colo1[IDs[,'Note']=='PositiveControls'],
       colo1[IDs[,'Note']=='NegativeControls'])

t.test(colo1[IDs[,'Note']=='PositiveControls'],
       colo1[IDs[,'Note']=='NegativeControls'])$p.value


############# Fig3B 

IDs_c1=readRDS('/input/IDs_colo1.rds')
c1=as.matrix(readRDS('/input/colo1_plm_GEN.rds'))

row_mean_na_rm <- function(x) {
  mean(x, na.rm = TRUE)
}

cell=unique(unlist(strsplit(colnames(c1),'_'))[c(T,F)])
temp=matrix(0,nrow(c1),length(cell))
rownames(temp)=rownames(c1)
colnames(temp)=cell

for (i in 1:ncol(temp)){
  temp[,i]=apply(c1[,grep(colnames(temp)[i],colnames(c1))],1,row_mean_na_rm)
}
c1=temp

####Genetic interactions 
map_mat=c1[grep('GIControls',rownames(c1)),]
map_ids=IDs_c1[IDs_c1[,'Note']=='GIControls',]
map_ids=cbind(map_ids,paste(map_ids[,'Note'],map_ids[,'Gene1'],map_ids[,'Gene2'],sep='_'))
colnames(map_ids)[ncol(map_ids)]='pair'
comb=map_mat[rownames(map_mat)%in%map_ids[map_ids[,'MyNote']=='GIControlsCombinations','pair'],]
sing=map_mat[rownames(map_mat)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons','pair'],]
ctrl=c('ADAD1','CYLC2','KLK12')
sum(map_ids[map_ids[,'MyNote']=='GIControlsSingletons','Gene1']%in%ctrl)
anch=sing[rownames(sing)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons' & map_ids[,'Gene1']%in%ctrl,'pair'],]
libr=sing[rownames(sing)%in%map_ids[map_ids[,'MyNote']=='GIControlsSingletons' & map_ids[,'Gene2']%in%ctrl,'pair'],]

combi_c1=rbind(c1[c(grep('LibraryCombinations',rownames(c1)),grep('AnchorCombinations',rownames(c1))),],comb)
a_c1=rbind(c1[grep('AnchorSingletons',rownames(c1)),],anch)
l_c1=rbind(c1[grep('LibrarySingletons',rownames(c1)),],libr)

IDs_combi_c1=cbind(unlist(strsplit(rownames(combi_c1),'_'))[c(F,T,F)],
                   unlist(strsplit(rownames(combi_c1),'_'))[c(F,F,T)])

colnames(IDs_combi_c1)=c('lib','anc')

temp=data.frame(a_c1,gen=unlist(strsplit(rownames(a_c1),'_'))[c(F,F,T)])
temp=aggregate(. ~ gen , data=temp, median, na.action = NULL, na.rm=T)
a_c1=apply(as.matrix(temp[,2:ncol(temp)]),2,as.numeric)
rownames(a_c1)=temp$gen

temp=data.frame(l_c1,gen=unlist(strsplit(rownames(l_c1),'_'))[c(F,T,F)])
temp=aggregate(. ~ gen , data=temp, median, na.action = NULL, na.rm=T)
l_c1=apply(as.matrix(temp[,2:ncol(temp)]),2,as.numeric)
rownames(l_c1)=temp$gen
##REady to compare with paralogs
rownames(IDs_combi_c1)=paste(IDs_combi_c1[,1],IDs_combi_c1[,2],sep='_')

m_c1=c(a_c1[,grep('136',colnames(a_c1))],l_c1[,grep('136',colnames(a_c1))])

combi_c1=combi_c1[,grep('136',colnames(combi_c1))]

combi_c1=cbind(combi_c1,
               m_c1[IDs_combi_c1[,'lib']],
               m_c1[IDs_combi_c1[,'anc']])

rownames(combi_c1)=paste(unlist(strsplit(rownames(combi_c1),'_'))[c(F,T,F)],
                         unlist(strsplit(rownames(combi_c1),'_'))[c(F,F,T)],sep='_')
colnames(combi_c1)=c('pair','lib','anc')
###filtering NAs
temp=rep(0,nrow(combi_c1))
for (i in 1:nrow(combi_c1)){
  temp[i]=sum(is.na(combi_c1[i,]))
}
combi_c1=combi_c1[temp==0,]
####Observed minus expected (type bliss)
combi_c1=cbind(combi_c1,combi_c1[,'pair']-(combi_c1[,'anc']+combi_c1[,'lib']))
colnames(combi_c1)[ncol(combi_c1)]='sco'
####Observed-expected as high effect (HSA)
combi_c1=cbind(combi_c1,0)
colnames(combi_c1)[ncol(combi_c1)]='HSA'

for (i in 1:nrow(combi_c1)){
  
  combi_c1[i,'HSA']=combi_c1[i,'pair']-min(combi_c1[i,c('anc','lib')])
  
}

combi_c1=combi_c1[order(combi_c1[,'sco']),]
combi_c1[combi_c1[,'pair']<=-1 & combi_c1[,'sco']<=-1 ,]

###uploading table
#######PARALOGS from Ryan et al https://doi.org/10.1016/j.cels.2021.08.006
para=as.matrix(read.delim('/input/paralog_scores_ryan.csv',sep=','))
rownames(para)=para[,'sorted_gene_pair']
para=para[,c('prediction_rank','prediction_percentile','prediction_score','sorted_gene_pair','A1','A2' )]

###Uploading biogrid table to select genetic interactions
biogrid=rbind(readRDS('/input/biogrid_human-4.4.224.mitab_part1.rds'),
              readRDS('/input/biogrid_human-4.4.224.mitab_part2.rds'))
#selecting interactors
sel_type1=c('psi-mi:MI:0915(physical association)', 'psi-mi:MI:0407(direct interaction)')
sel_type2=c('psi-mi:MI:2373(negative genetic interaction (sensu biogrid))')

inte_1=biogrid[biogrid[,'Interaction.Types']%in%sel_type1 & 
                 biogrid[,'Interaction.Detection.Method']%in%c('psi-mi:MI:0004(affinity chromatography technology)','psi-mi:MI:0096(pull down)'),c("Alt.IDs.Interactor.A","Alt.IDs.Interactor.B")]
inte_2=biogrid[biogrid[,'Interaction.Types']%in%sel_type2,c("Alt.IDs.Interactor.A","Alt.IDs.Interactor.B")]

colnames(inte_1)=c('A','B')
colnames(inte_2)=c('A','B')

for (i in 1:nrow(inte_1)){
  
  temp1=unlist(strsplit(inte_1[i,'A'],'\\|'))
  temp2=unlist(strsplit(inte_1[i,'B'],'\\|'))
  temp1=temp1[grep('entrez gene/locuslink',temp1)]
  temp2=temp2[grep('entrez gene/locuslink',temp2)]
  temp1=gsub('entrez gene/locuslink:','',temp1)
  temp2=gsub('entrez gene/locuslink:','',temp2)
  inte_1[i,'A']=temp1[1]
  inte_1[i,'B']=temp2[1]
  
}

for (i in 1:nrow(inte_2)){
  
  temp1=unlist(strsplit(inte_2[i,'A'],'\\|'))
  temp2=unlist(strsplit(inte_2[i,'B'],'\\|'))
  temp1=temp1[grep('entrez gene/locuslink',temp1)]
  temp2=temp2[grep('entrez gene/locuslink',temp2)]
  temp1=gsub('entrez gene/locuslink:','',temp1)
  temp2=gsub('entrez gene/locuslink:','',temp2)
  inte_2[i,'A']=temp1[1]
  inte_2[i,'B']=temp2[1]
  
}

inte_1=c(paste(inte_1[,'A'],inte_1[,'B'],sep='_'),
         paste(inte_1[,'B'],inte_1[,'A'],sep='_'))
inte_2=c(paste(inte_2[,'A'],inte_2[,'B'],sep='_'),
         paste(inte_2[,'B'],inte_2[,'A'],sep='_'))

inte_1=unique(inte_1)
inte_2=unique(inte_2)


GI_1=!rownames(combi_c1)%in%inte_2
GI_2=rownames(combi_c1)%in%inte_2

sort_bliss=cbind(c(1:nrow(combi_c1)),
                 combi_c1[,'sco'])

####Boxplot final 

### 1.-Paralogos
sel1=as.numeric(para[,'prediction_percentile'])<=10 & as.numeric(para[,'prediction_score'])<0.2
sel2=as.numeric(para[,'prediction_score'])>=0.2

para_1=c(paste(para[sel1,'A1'],para[sel1,'A2'],sep='_'),
         paste(para[sel1,'A2'],para[sel1,'A1'],sep='_'))
para_2=c(paste(para[sel2,'A1'],para[sel2,'A2'],sep='_'),
         paste(para[sel2,'A2'],para[sel2,'A1'],sep='_'))


x1=sort_bliss[!rownames(combi_c1)%in%c(para_1,para_2),2]
x2=sort_bliss[rownames(combi_c1)%in%para_1,2]
x3=sort_bliss[rownames(combi_c1)%in%para_2,2]

x4=combi_c1[GI_1,'sco']
x5=combi_c1[GI_2,'sco']
x6=combi_c1[rownames(combi_c1)%in%inte_2[inte_2%in%inte_1],'sco']

pdf('/pathFortheFigure/fig3B.pdf',width=4,height = 6)

box=boxplot(x1,c(x2,x3),x4,x5,x6,ylim=c(-2,1.5),col='white',
            ylab='bliss score (observed-expected)',main='Bliss score based on logFC (HT29)')

points(jitter(rep(2,length(x3)),amount=0.2),x3,pch=16,col=rgb(0.7,0.2,0.2))

text(1.5,1,paste('p.val:',round(t.test(x1,c(x2,x3),exact=T)$p.value,4)),cex=0.7)
text(3.5,1,paste('p.val:',round(t.test(x4,x5,exact=T)$p.value,4)),cex=0.7)
text(4.5,1,paste('p.val:',round(t.test(x4,x6,exact=T)$p.value,4)),cex=0.7)

text(1:5,rep(-1.8,5),paste('n:',c(length(x1),length(c(x3,x2)),length(x4),length(x5),length(x6))),cex=0.7)
text(2,-2,paste('n:',length(x3)),cex=0.7)

dev.off()

#############

###########Figure 3C 

colo1=readRDS('/Users/inigo/projects/ENCORE/160823/rds_cas9_test/colo1_plm_scale.rds')
IDs=readRDS('/Users/inigo/projects/ENCORE/160823/rds_cas9_test/IDs_colo1.rds')
colo1=colo1[,grep('136',colnames(colo1))]
colo1=apply(colo1[,1:3],1,mean)

sel_pair1=c('BCL2L1','EGFR','AKT2','EGFR','EGFR')
sel_pair2=c('MCL1','BRAF','AKT1','PIK3CA','DNMT1')

h1=colo1[IDs[IDs[,'Gene1']=='HDAC1' & IDs[,'Note']=='LibrarySingletons','ID']]
h2=colo1[IDs[IDs[,'Gene1']=='HDAC1' & IDs[,'Gene2']=='HDAC2','ID']]
h3=colo1[IDs[IDs[,'Gene2']=='HDAC2' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
a1=colo1[IDs[IDs[,'Gene1']=='ASF1A' & IDs[,'Note']=='LibrarySingletons','ID']]
a2=colo1[IDs[IDs[,'Gene1']=='ASF1A' & IDs[,'Gene2']=='ASF1B','ID']]
a3=colo1[IDs[IDs[,'Gene2']=='ASF1B' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
c1=colo1[IDs[IDs[,'Gene1']=='CNOT7' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
c2=colo1[IDs[IDs[,'Gene1']=='CNOT7' & IDs[,'Gene2']=='CNOT8','ID']]
c3=colo1[IDs[IDs[,'Gene2']=='CNOT8' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
m1=colo1[IDs[IDs[,'Gene1']=='MAPK1' & IDs[,'MyNote']=='GIControlsSingletons','ID']]
m2=colo1[IDs[IDs[,'Gene1']=='MAPK1' & IDs[,'Gene2']=='MAPK3','ID']]
m3=colo1[IDs[IDs[,'Gene2']=='MAPK3' & IDs[,'MyNote']=='GIControlsSingletons','ID']]

gi1_1=colo1[IDs[IDs[,'Gene1']==sel_pair1[1] & IDs[,'MyNote']=='GIControlsSingletons','ID']]
gi1_2=colo1[IDs[IDs[,'Gene1']==sel_pair1[1] & IDs[,'Gene2']==sel_pair2[1],'ID']]
gi1_3=colo1[IDs[IDs[,'Gene2']==sel_pair2[1] & IDs[,'Note']=='AnchorSingletons','ID']]

gi2_1=colo1[IDs[IDs[,'Gene1']==sel_pair1[2] & IDs[,'MyNote']=='GIControlsSingletons','ID']]
gi2_2=colo1[IDs[IDs[,'Gene1']==sel_pair1[2] & IDs[,'Gene2']==sel_pair2[2],'ID']]
gi2_3=colo1[IDs[IDs[,'Gene2']==sel_pair2[2] & IDs[,'Note']=='AnchorSingletons','ID']]

gi3_1=colo1[IDs[IDs[,'Gene2']==sel_pair1[3] & IDs[,'Note']=='AnchorSingletons','ID']]
gi3_2=colo1[IDs[IDs[,'Gene1']==sel_pair1[3] & IDs[,'Gene2']==sel_pair2[3],'ID']]
gi3_3=colo1[IDs[IDs[,'Gene2']==sel_pair2[3] & IDs[,'Note']=='AnchorSingletons','ID']]

gi4_1=colo1[IDs[IDs[,'Gene1']==sel_pair1[4] & IDs[,'MyNote']=='GIControlsSingletons','ID']]
gi4_2=colo1[IDs[IDs[,'Gene1']==sel_pair1[4] & IDs[,'Gene2']==sel_pair2[4],'ID']]
gi4_3=colo1[IDs[IDs[,'Gene2']==sel_pair2[4] & IDs[,'Note']=='AnchorSingletons','ID']]

gi5_1=colo1[IDs[IDs[,'Gene1']==sel_pair1[5] & IDs[,'MyNote']=='GIControlsSingletons','ID']]
gi5_2=colo1[IDs[IDs[,'Gene1']==sel_pair1[5] & IDs[,'Gene2']==sel_pair2[5],'ID']]
gi5_3=colo1[IDs[IDs[,'Gene2']==sel_pair2[5] & IDs[,'Note']=='AnchorSingletons','ID']]

pal=c(rgb(170/255,222/255,135/255),
      rgb(255/255,230/255,128/255),
      rgb(170/255,222/255,135/255))

nombres=c(sel_pair1[1],paste(sel_pair1[1],sel_pair2[1],sep='+'),sel_pair2[1],
          sel_pair1[2],paste(sel_pair1[2],sel_pair2[2],sep='+'),sel_pair2[2],
          sel_pair1[3],paste(sel_pair1[3],sel_pair2[3],sep='+'),sel_pair2[3],
          sel_pair1[4],paste(sel_pair1[4],sel_pair2[4],sep='+'),sel_pair2[4],
          sel_pair1[5],paste(sel_pair1[5],sel_pair2[5],sep='+'),sel_pair2[5])

pdf('pathFortheFigure/Fig3C_boxplot_examples.pdf',width=10,height=6)


par(mfrow=c(1,2))


boxplot(h1,h2,h3,a1,a2,a3,c1,c2,c3,m1,m2,m3,1,1,1,
        names=c('HDAC1','HDAC1+HDAC2','HDAC2',
                'ASF1A','ASF1A+ASF1B','ASF1B',
                'CNOT7','CNOT7+CNOT8','CNOT8',
                'MAPK1','MAPK1+MAPK3','MAPK3','','',''),col='white',ylim=c(-3,2),outline=F,las=3,
        ylab='logFC(HT29-sgRNA) known paralogs')

points(jitter(rep(1,length(h1)),amount=0.3),h1,col=pal[1],pch=16)
points(jitter(rep(2,length(h2)),amount=0.3),h2,col=pal[2],pch=16)
points(jitter(rep(3,length(h3)),amount=0.3),h3,col=pal[3],pch=16)
points(jitter(rep(4,length(a1)),amount=0.3),a1,col=pal[1],pch=16)
points(jitter(rep(5,length(a2)),amount=0.3),a2,col=pal[2],pch=16)
points(jitter(rep(6,length(a3)),amount=0.3),a3,col=pal[3],pch=16)
points(jitter(rep(7,length(c1)),amount=0.3),c1,col=pal[1],pch=16)
points(jitter(rep(8,length(c2)),amount=0.3),c2,col=pal[2],pch=16)
points(jitter(rep(9,length(c3)),amount=0.3),c3,col=pal[3],pch=16)
points(jitter(rep(10,length(m1)),amount=0.3),m1,col=pal[1],pch=16)
points(jitter(rep(11,length(m2)),amount=0.3),m2,col=pal[2],pch=16)
points(jitter(rep(12,length(m3)),amount=0.3),m3,col=pal[3],pch=16)

t_h=c(round(t.test(h1,h2)$p.value,4),round(t.test(h3,h2)$p.value,4))
t_a=c(round(t.test(a1,a2)$p.value,4),round(t.test(a3,a2)$p.value,4))
t_c=c(round(t.test(c1,c2)$p.value,4),round(t.test(c3,c2)$p.value,4))
t_m=c(round(t.test(m1,m2)$p.value,4),round(t.test(m3,m2)$p.value,4))

text(c(1,3),c(1,1),paste('pval:',t_h),pos=4,srt=90)
text(c(4,6),c(1,1),paste('pval:',t_a),pos=4,srt=90)
text(c(7,9),c(1,1),paste('pval:',t_c),pos=4,srt=90)
text(c(10,12),c(1,1),paste('pval:',t_m),pos=4,srt=90)


boxplot(gi1_1,gi1_2,gi1_3,gi2_1,gi2_2,gi2_3,gi3_1,gi3_2,gi3_3,gi4_1,gi4_2,gi4_3,gi5_1,gi5_2,gi5_3,
        names=nombres,col='white',ylim=c(-3,2),outline=F,las=3,
        ylab='logFC(HT29-sgRNA) negative GI (biogrid)')

points(jitter(rep(1,length(gi1_1)),amount=0.3),gi1_1,col=pal[1],pch=16)
points(jitter(rep(2,length(gi1_2)),amount=0.3),gi1_2,col=pal[2],pch=16)
points(jitter(rep(3,length(gi1_3)),amount=0.3),gi1_3,col=pal[3],pch=16)
points(jitter(rep(4,length(gi2_1)),amount=0.3),gi2_1,col=pal[1],pch=16)
points(jitter(rep(5,length(gi2_2)),amount=0.3),gi2_2,col=pal[2],pch=16)
points(jitter(rep(6,length(gi2_3)),amount=0.3),gi2_3,col=pal[3],pch=16)
points(jitter(rep(7,length(gi3_1)),amount=0.3),gi3_1,col=pal[1],pch=16)
points(jitter(rep(8,length(gi3_2)),amount=0.3),gi3_2,col=pal[2],pch=16)
points(jitter(rep(9,length(gi3_3)),amount=0.3),gi3_3,col=pal[3],pch=16)
points(jitter(rep(10,length(gi4_1)),amount=0.3),gi4_1,col=pal[1],pch=16)
points(jitter(rep(11,length(gi4_2)),amount=0.3),gi4_2,col=pal[2],pch=16)
points(jitter(rep(12,length(gi4_3)),amount=0.3),gi4_3,col=pal[3],pch=16)
points(jitter(rep(13,length(gi5_1)),amount=0.3),gi5_1,col=pal[1],pch=16)
points(jitter(rep(14,length(gi5_2)),amount=0.3),gi5_2,col=pal[2],pch=16)
points(jitter(rep(15,length(gi5_3)),amount=0.3),gi5_3,col=pal[3],pch=16)


t_1=c(round(t.test(gi1_1,gi1_2)$p.value,4),round(t.test(gi1_2,gi1_3)$p.value,4))
t_2=c(round(t.test(gi2_1,gi2_2)$p.value,4),round(t.test(gi2_2,gi2_3)$p.value,4))
t_3=c(round(t.test(gi3_1,gi3_2)$p.value,4),round(t.test(gi3_2,gi3_3)$p.value,4))
t_4=c(round(t.test(gi4_1,gi4_2)$p.value,4),round(t.test(gi4_2,gi4_3)$p.value,4))
t_5=c(round(t.test(gi5_1,gi5_2)$p.value,4),round(t.test(gi5_2,gi5_3)$p.value,4))

text(c(1,3),c(1,1),paste('pval:',t_1),pos=4,srt=90)
text(c(4,6),c(1,1),paste('pval:',t_2),pos=4,srt=90)
text(c(7,9),c(1,1),paste('pval:',t_3),pos=4,srt=90)
text(c(10,12),c(1,1),paste('pval:',t_4),pos=4,srt=90)
text(c(13,15),c(1,1),paste('pval:',t_5),pos=4,srt=90)

dev.off()



####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE1
####### ####### ####### ####### ####### ####### ####### ####### 
#################################################################
##### Supplementary fig1B

tab_exact=as.matrix(read.delim('PILOT/MERGED/PILOT_EXACT_COUNTS.txt'))

counts=apply(tab_exact[,9:ncol(tab_exact)],2,as.numeric)
rownames(counts)=tab_exact[,'ID']

plasmid_count=counts[,"PILOT_Plasmid"]
plasmid_count=plasmid_count[order(plasmid_count)]
seq_IDs=seq(from=1,to=length(plasmid_count),by=83)

guide_means <- rowMeans(counts)          # mean across all input/plasmid replicates ideally
guide_sds   <- apply(counts, 1, sd)
guide_vars  <- apply(counts, 1, var)
guide_cv    <- guide_sds / (guide_means + 1e-8) # avoid divide by zero

summary_df <- data.frame(
  guide = rownames(counts),
  mean_count = guide_means,
  sd = guide_sds,
  var = guide_vars,
  cv = guide_cv
)
summary_df=as.matrix(summary_df)
temp=rownames(summary_df)
summary_df=apply(summary_df[,2:5],2,as.numeric)
rownames(summary_df)=temp

###
pdf("pathFortheFigure/Fig1suppB.pdf", height = 7,width = 8)
plot(counts[,"PILOT_Plasmid"],as.numeric(summary_df[,"cv"]),
     main="raw counts",ylab="counts CV accross all conditions per guide",
     xlab="plasmid counts",col=rgb(0.7,0.2,0.2,0.5),pch=16)
abline(v=10,col=rgb(0.2,0.5,0.7),lty=2,lwd=2)
dev.off()


#################################################################
##### Supplementary fig1C

###Columns to keep plus naming
pos=c(0,10,11,12,13,14,16,17,19,20,22,23,25,26,29,30,31,32,33)+1
nam=c('Readname','VectorID1','V1edit1R1','V1direction1R1','V1edit2R1','V1direction2R1','Scaffold','Scaffold_edit','Linker','Linker_edit','tRNA','tRNA_edit','Improved','Improved_edit','VectorID2','V2edit1R1','V2direction1R1','V2edit2R2','V2direction2R2')
##REsults table
tab_final=matrix(0,4,3)

colnames(tab_final)=c('wt','m6','m7')
rownames(tab_final)=c('sgrna1_exists','sgrna2_exists',
                      'sgrna_match','swap')
###plasmid lib
x=as.matrix(read.table("PILOT/MAPPING/WALK_280/lib_S5_walk280.mapping.out", sep="\t", skip=3, header=FALSE))
x=x[,pos]
colnames(x)=nam

for (j in 1:ncol(x)){
  x[,j]=gsub(' ','',x[,j])
}

####Los NA nos estan jodiendo, vamso a intentar otra estratagema

for (j in 1:ncol(x)){
  x[is.na( x[,j]),j]='No_Hay'
}

##scaffolds
wt=x[x[,'Scaffold']=='Wildtype',]
m6=x[x[,'Scaffold']=='Modified6',]
m7=x[x[,'Scaffold']=='Modified7',]

#WT

a=wt

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'wt']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

#m6

a=m6

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'m6']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

#m7

a=m7

##Misma estrategia, a para limpia de NAs, selectiones y b para counts
total_number_of_reads=nrow(a)

###########

sel=(a[,'VectorID1']!="OEM2" & a[,'VectorID1']!="OEM2_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM2" & a[,'VectorID2']!="OEM2_WO" & a[,'VectorID2']!="NOMATCH")

b = sum(as.numeric(a[sel,'V1edit1R1'])==0 & 
          as.numeric(a[sel,'V2edit1R1'])==0 & 
          a[sel,'V1direction1R1']=="FALSE" & 
          a[sel,'V2direction1R1']=="FALSE")

sgRNA1_exists=(b/total_number_of_reads)*100

###########

sel=(a[,'VectorID1']!="OEM1" & a[,'VectorID1']!="OEM1_WO" & a[,'VectorID1']!="NOMATCH") | 
  (a[,'VectorID2']!="OEM1" & a[,'VectorID2']!="OEM1_WO" & a[,'VectorID2']!="NOMATCH")

b= sum(a[sel,'V1edit2R1']==0 & 
         a[sel,'V2edit2R2']==0 & 
         a[sel,'V1direction2R1']=="FALSE" & 
         a[sel,'V2direction2R2']=="R")

sgRNA2_exists=(b/total_number_of_reads)*100

#################

b=sum(  a[,'VectorID1']==a[,'VectorID2'] & 
          a[,'V1edit2R1']==0 & 
          a[,'V2edit1R1']==0 &
          a[,'V1direction1R1']=="FALSE" & 
          a[,'V1direction2R1']=="FALSE" &
          a[,'V2edit1R1']==0 &
          a[,'V2edit2R2']==0 &
          a[,'V2direction1R1']=="FALSE" & 
          a[,'V2direction2R2']=="R")

sgRNA_match = (b/total_number_of_reads)*100

#################

b =sum(a[,'VectorID1']==a[,'VectorID2'] & 
         a[,'VectorID1']=="SWAP" & 
         a[,'V1edit2R1']==0 & 
         a[,'V2edit1R1']==0 &
         a[,'V2edit1R1']==0 & 
         a[,'V2edit2R2']==0)
swap = (b/total_number_of_reads)*100

########
tab_final[,'m7']=c(sgRNA1_exists,sgRNA2_exists,sgRNA_match,swap)

pdf('pathFortheFigure/Fig1suppC_overall_stats.pdf',width=7,height=5)
x=barplot(t(tab_final),beside=T,ylim=c(0,100),legend=c('WT','M6','M7'),col=c(rgb(233/255,164/255,49/255),
                                                                             rgb(205/255,228/255,241/255),
                                                                             rgb(67/255,82/255,114/255)),border=NA,
          main='Overall stats (Pilot library plasmid)',ylab='percentage')

text(x[,1],tab_final[1,]+2,round(tab_final[1,],1))
text(x[,2],tab_final[2,]+2,round(tab_final[2,],1))
text(x[,3],tab_final[3,]+2,round(tab_final[3,],1))
text(x[,4],tab_final[4,]+2,round(tab_final[4,],1))
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### 
######### Supp figure 1D

dat=as.matrix(read.delim('PILOT/MAPPING/WALK_280/lib_S5_walk280.mapping.out',sep='\t',skip=3,header = F))

###Select swaps
dat[,'V14']=gsub(' ','',dat[,'V14'])
dat[,'V31']=gsub(' ','',dat[,'V31'])
dat[,'V33']=gsub(' ','',dat[,'V33'])

selec=dat[,'V11']=='SWAP' & dat[,'V30']=='SWAP' & dat[,'V14']=='0' & dat[,'V31']=='0' & dat[,'V33']=='0' 
sum(selec)/nrow(dat)
dat=dat[selec,]

###heatmap for los counts
heat_mat=matrix(0,length(unique(dat[,'V2'])),length(unique(dat[,'V6'])))

row_ID=unique(dat[,'V2'])
col_ID=unique(dat[,'V6'])

###transformar a IDs numericos para facilitar la tarea
temp_dat=cbind(dat[,c('V2','V6')],'','')

for(i in 1:length(col_ID)){
  
  temp_dat[temp_dat[,2]==col_ID[i],4]=i

}

for(i in 1:length(row_ID)){
  
  temp_dat[temp_dat[,1]==row_ID[i],3]=i
  
}

###counts for each pair 
freq_temp=table(paste(temp_dat[,3],temp_dat[,4],sep='_'))

####IDs numericos para la matriz final
rownames(heat_mat)=1:nrow(heat_mat)
colnames(heat_mat)=1:ncol(heat_mat)

for (i in 1:nrow(heat_mat)){
  
  temp=freq_temp[unlist(strsplit(names(freq_temp),'_'))[c(T,F)]==i]
    names(temp)=unlist(strsplit(names(temp),'_'))[c(F,T)]
    heat_mat[i,names(temp)]=temp
}

for (i in 1:nrow(heat_mat)){
  
  heat_mat[i,heat_mat[i,]==0]=NA
  
}

pdf('pathFortheFigure/Fig1suppD.pdf',height = 100,width = 100)
pheatmap(heat_mat,color=viridis(1000),
         cluster_cols=F,cluster_rows=F,na_col='white',labels_row=F,labels_col=F,border_color = NA,
         cellwidth = 1, cellheight = 1,fontsize=1)
dev.off()




####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE2
####### ####### ####### ####### ####### ####### ####### ####### 
######### Supp figure 2A


#### Reload Anotation 

anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

table(anot[,'sgRNA1_Class'])
table(anot[,'sgRNA2_Class'])

pos=unique(c(anot[anot[,'sgRNA1_Class']=='essential','Gene1'],
             anot[anot[,'sgRNA2_Class']=='essential','Gene2']))
neg=unique(c(anot[anot[,'sgRNA1_Class']=='non-essential','Gene1'],
             anot[anot[,'sgRNA2_Class']=='non-essential','Gene2']))

pdf('pathFortheFigure/Fig2ASupp.pdf',width=5,height=5.5)
boxplot(score_ag,
        score_ag[names(score_ag)%in%pos],
        score_ag[names(score_ag)%in%neg],0,0,0,0,
        ylim=c(-6,2),
        names=c('library','essential','non-essential','','','',''),main='HT29 projetc Score',ylab='logFC',las=3)
text(2.5,1,paste('cohens D:',round(dcohen(score_ag[names(score_ag)%in%neg],
                                          score_ag[names(score_ag)%in%pos]),2)))
text(c(1,2,3),c(2,2,2),paste('n:',c(length(score_ag),sum(names(score_ag)%in%neg),sum(names(score_ag)%in%pos))))
dev.off()

######### Supp figure 2B

anot=as.matrix(read.delim('/PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']
lfc_plm=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut'
scaf2=anot[,'Scaffold']=='Modified6'& anot[,'Notes']!='DistanceCut'
scaf3=anot[,'Scaffold']=='Modified7'& anot[,'Notes']!='DistanceCut'

##We cut scaf WT to be same pairs as the others 
guides_ID=paste(anot[,'sgRNA1_ID'],anot[,'sgRNA2_ID'],sep='')
scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut' & guides_ID%in%guides_ID[scaf2]

#####HAcemos la media para los 3  

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

pal1=c(rgb(108/255,116/255,134/255),rgb(105/255,156/255,201/255),rgb(205/255,232/255,236/255))
colores=pal1[c(1,1,1,2,2,2,3,3,3,2,2,2,3,3,3,2,2,2,3,3,3)]
pos=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic')
neg=c('non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')

pdf('/pathFortheFigure/Fig2BSupp_boxplot_D3.pdf',width=4,height=5)
boxplot(mat[anot[,'Notes']!='DistanceCut',1],
        mat[anot[,'Notes']!='DistanceCut',2],
        mat[anot[,'Notes']!='DistanceCut',3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%neg,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%neg,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%pos,3],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,1],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,2],
        mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%neg,3],
        col=alpha(colores,1),xlab='sgRNA Scaffold',ylab='log2FC',las=3)
dev.off()

##Cohen´sD and ttest

d_cohen_text=c(dcohen( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,1],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,1]),
               dcohen( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,2],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,2]),
               dcohen( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,3],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,3]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,1],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,1]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,2],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,2]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,3],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,3]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,1],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,1]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,2],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,2]),
               dcohen(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,3],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,3]))
d_cohen_text=round(d_cohen_text,2)

ttest_text=c(t.test( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,1],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,1])$p.value,
             t.test( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,2],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,2])$p.value,
             t.test( mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype' & anot[,'vector_class']%in%neg,3],
                       mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Wildtype'  & anot[,'vector_class']%in%pos,3])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,1],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,1])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,2],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,2])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6' & anot[,'vector_class']%in%neg,3],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified6'  & anot[,'vector_class']%in%pos,3])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,1],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,1])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,2],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,2])$p.value,
             t.test(mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7' & anot[,'vector_class']%in%neg,3],
                      mat[anot[,'Notes']!='DistanceCut' & anot[,'Scaffold']=='Modified7'  & anot[,'vector_class']%in%pos,3])$p.value)

############################################################################################################################
############################################################################################################################
######### Supp figure 2C


scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut'
scaf2=anot[,'Scaffold']=='Modified6'& anot[,'Notes']!='DistanceCut'
scaf3=anot[,'Scaffold']=='Modified7'& anot[,'Notes']!='DistanceCut'

breaksList = seq(0.5, 1, length.out = 100)

colnames(lfc_plm)=c("100X_R1","100X_R2","100X_R3",
                    "500X_R1","500X_R2","500X_R3",
                    "PCR500X_R1","PCR500X_R2","PCR500X_R3")

pdf('/pathFortheFigure/Fig2Csup_corr_D3.pdf',width=5,height=5)

pheatmap(cor(lfc_plm),col=colorRampPalette(brewer.pal(5,'Reds'))(100),
         fontsize = 7,cellwidth = 7, cellheight = 7,
         breaks = breaksList)

dev.off()


#####Positional effect Supp fig 2D

anot=as.matrix(read.delim('/PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

lfc_plm=readRDS('/PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

#####HAcemos la media para los 3  (aunque usaremos 500X)


mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

####Agregamos con todos juntos

####Vamos a agregar

sel1=anot[,'Scaffold']=='Wildtype'
sel2=anot[,'Scaffold']=='Modified6'
sel3=anot[,'Scaffold']=='Modified7'

temp=data.frame(G1=anot[sel1,'Gene1'],G2=anot[sel1,'Gene2'],lfc_500x=mat[sel1,'PILOT_500X_D14'],lfc_100x=mat[sel1,'PILOT_100X_D14'],lfc_PCR500x=mat[sel1,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2, data=temp, mean)
temp=as.matrix(temp)
wt=temp
temp=data.frame(G1=anot[sel2,'Gene1'],G2=anot[sel2,'Gene2'],lfc_500x=mat[sel2,'PILOT_500X_D14'],lfc_100x=mat[sel2,'PILOT_100X_D14'],lfc_PCR500x=mat[sel2,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m6=temp

temp=data.frame(G1=anot[sel3,'Gene1'],G2=anot[sel3,'Gene2'],lfc_500x=mat[sel3,'PILOT_500X_D14'],lfc_100x=mat[sel3,'PILOT_100X_D14'],lfc_PCR500x=mat[sel3,'PILOT_PCR500X_D14'])
temp <- aggregate(. ~ G1 + G2 , data=temp, mean)
temp=as.matrix(temp)
m7=temp


rownames(wt)=paste(wt[,1],wt[,2],sep='_')
temp=wt[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(wt[,2],wt[,1],sep='_')
wt=wt[rownames(wt)%in%rownames(temp),]
wt=cbind(wt,temp[rownames(wt),])
wt=wt[,3:8]
colnames(wt)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(wt)
wt=apply(wt,2,as.numeric)
rownames(wt)=temp

rownames(m6)=paste(m6[,1],m6[,2],sep='_')
temp=m6[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(m6[,2],m6[,1],sep='_')
m6=m6[rownames(m6)%in%rownames(temp),]
m6=cbind(m6,temp[rownames(m6),])
m6=m6[,3:8]
colnames(m6)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(m6)
m6=apply(m6,2,as.numeric)
rownames(m6)=temp

rownames(m7)=paste(m7[,1],m7[,2],sep='_')
temp=m7[,c("lfc_500x","lfc_100x","lfc_PCR500x")]
rownames(temp)=paste(m7[,2],m7[,1],sep='_')
m7=m7[rownames(m7)%in%rownames(temp),]
m7=cbind(m7,temp[rownames(m7),])
m7=m7[,3:8]
colnames(m7)=c('lfc_500x_1','lfc_100x_1','lfc_PCR500x_1',
               'lfc_500x_2','lfc_100x_2','lfc_PCR500x_2')
temp=rownames(m7)
m7=apply(m7,2,as.numeric)
rownames(m7)=temp


unique_ID=t(combn(unique(unlist(strsplit(rownames(wt),'_'))),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_wt=!rownames(wt)%in%unique_ID

unique_ID=t(combn(unique(unlist(strsplit(rownames(m6),'_'))),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_m6=!rownames(m6)%in%unique_ID

unique_ID=t(combn(unique(unlist(strsplit(rownames(m7),'_'))),2))
unique_ID=paste(unique_ID[,1],unique_ID[,2],sep='_')
sel_m7=!rownames(m7)%in%unique_ID


pdf('pathFortheFigure/FigSupp2D_all_D3.pdf',width=10,height=10)
par(mfrow=c(3,3))

###WT
data <- data.frame(sgRNA1 = wt[sel_wt,'lfc_500x_1'], sgRNA2 = wt[sel_wt,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA2 ~ sgRNA1, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = wt[sel_wt,'lfc_100x_1'], sgRNA2 = wt[sel_wt,'lfc_100x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = wt[sel_wt,'lfc_PCR500x_1'], sgRNA2 = wt[sel_wt,'lfc_PCR500x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='Wildtype (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

###M6
data <- data.frame(sgRNA1 = m6[sel_m6,'lfc_500x_1'], sgRNA2 = m6[sel_m6,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA2 ~ sgRNA1, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[sel_m6,'lfc_100x_1'], sgRNA2 = m6[sel_m6,'lfc_100x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m6[sel_m6,'lfc_PCR500x_1'], sgRNA2 = m6[sel_m6,'lfc_PCR500x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M6 (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))


###M7
data <- data.frame(sgRNA1 = m7[sel_m7,'lfc_500x_1'], sgRNA2 = m7[sel_m7,'lfc_500x_2'])
# Fit linear model
reg<-lm(sgRNA2 ~ sgRNA1, data = data)
# Scatter plot
plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (X500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[sel_m7,'lfc_100x_1'], sgRNA2 = m7[sel_m7,'lfc_100x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (X100-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)
legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

data <- data.frame(sgRNA1 = m7[sel_m7,'lfc_PCR500x_1'], sgRNA2 = m7[sel_m7,'lfc_PCR500x_2'])
reg<-lm(sgRNA2 ~ sgRNA1, data = data)

plot(data,xlim=c(-5,2),ylim=c(-5,2),
     main='M7 (PCR500-D14/D3)',
     col=rgb(0.2,0.5,0.8,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.5),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))

dev.off()

#######  Supp fig 2F plasmid versus D3 control for LFC

anot=as.matrix(read.delim(PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

anot[anot[,'Gene2']=='','Gene2']=anot[anot[,'Gene2']=='','sgRNA2_Class']
anot[anot[,'Gene1']=='','Gene1']=anot[anot[,'Gene1']=='','sgRNA1_Class']

lfc_d3=readRDS('PILOT/MERGED/PILOC_EXACT_LFC_D3.rds')

mat_d3=cbind(apply(lfc_d3[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_d3[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_d3[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat_d3)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

lfc_pl=readRDS('PILOT/MERGED/PILOC_EXACT_LFC_PLM.rds')

mat_pl=cbind(apply(lfc_pl[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
             apply(lfc_pl[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
             apply(lfc_pl[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat_pl)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

pdf('/pathFortheFigure/FigSupp2F.pdf',width=5,height=5)

plot(mat_d3[,2],mat_pl[,2],ylim=c(-8,4),xlim=c(-8,4),main='PILOT_500X',xlab='D3',ylab='plasmid',
     col=rgb(0.2,0.5,0.7,0.5),pch=16)
dev.off()

#########
######### Supplementary figure 2 G

anot=as.matrix(read.delim('PILOT/LIB/PilotLib_Annot.txt'))
rownames(anot)=anot[,'ID']

lfc_plm=readRDS('PILOT/MERGED/PILOC_EXACT_LFC_PLM.rds')

scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut'
scaf2=anot[,'Scaffold']=='Modified6'& anot[,'Notes']!='DistanceCut'
scaf3=anot[,'Scaffold']=='Modified7'& anot[,'Notes']!='DistanceCut'

##We cut scaf WT to be same pairs as the others 
guides_ID=paste(anot[,'sgRNA1_ID'],anot[,'sgRNA2_ID'],sep='')
scaf1=anot[,'Scaffold']=='Wildtype'& anot[,'Notes']!='DistanceCut' & guides_ID%in%guides_ID[scaf2]

#####HAcemos la media para los 3  

mat=cbind(apply(lfc_plm[,c('PILOT_100X_D14_R1','PILOT_100X_D14_R2','PILOT_100X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_500X_D14_R1','PILOT_500X_D14_R2','PILOT_500X_D14_R3')],1,mean),
          apply(lfc_plm[,c('PILOT_PCR500X_D14_R1','PILOT_PCR500X_D14_R2','PILOT_PCR500X_D14_R3')],1,mean))
colnames(mat)=c('PILOT_100X_D14','PILOT_500X_D14','PILOT_PCR500X_D14')

###### ROC 

selection=c('essential + non-targeting','non-targeting + essential','intergenic + essential','essential + intergenic',
            'non-essential + non-essential','non-targeting + non-targeting','intergenic + intergenic')
colores=c("#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#3182bd", "#393b79","#756bb1")

AUCS_wt=rep(0,length(selection))
names(AUCS_wt)=selection
AUCS_m6=rep(0,length(selection))
names(AUCS_m6)=selection
AUCS_m7=rep(0,length(selection))
names(AUCS_m7)=selection
AUCS_ze=rep(0,length(selection))
names(AUCS_ze)=selection

AUCS_max=rep(0,9)
names(AUCS_max)=c('wt_100X','m6_100X','m7_100X',
                  'wt_500X','m6_500X','m7_500X',
                  'wt_PCR500X','m6_PCR500X','m7_PCR500X')

n_wt=AUCS_wt
n_m6=AUCS_m6
n_m7=AUCS_m7

n_max=AUCS_max

######## Lets ROC!
#####
pdf('/Fig2Gsup_recall_PLM_ROC_v2.pdf',width=10,height=4)
par(mfrow=c(1,3))

##################
################ 500X

### WT
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_WT')


for (i in 1:4){
  stype_df=sort(mat[scaf1,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf1,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_wt[i]=x$auc
  n_wt[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['wt_500X']=x$auc
n_max['wt_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

#####M6
plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_M6')

for (i in 1:4){
  stype_df=sort(mat[scaf2,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf2,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m6[i]=x$auc
  n_m6[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m6_500X']=x$auc
n_max['m6_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 

####M7

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='PILOT_500X_D14_M7')

for (i in 1:4){
  stype_df=sort(mat[scaf3,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[5:7],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

for (i in 5:7){
  stype_df=sort(mat[scaf3,'PILOT_500X_D14'])
  stype_df=stype_df[names(stype_df)%in%anot[anot[,'vector_class']%in%c(selection[1:4],selection[i]),'ID']]
  binary=rep(0,length(stype_df))
  binary[names(stype_df)%in%anot[anot[,'vector_class']==selection[i],'ID']]=1
  
  x=roc(	binary,stype_df,direction=">")
  
  AUCS_m7[i]=x$auc
  n_m7[i]=sum(binary==1)
  
  points(x$specificities,x$sensitivities,type='l',col=colores[i])
  
}

binary=rep(0,length(stype_df))
binary[names(stype_df)%in%anot[anot[,'vector_class']%in%selection[1:4],'ID']]=1
x=roc(	binary,stype_df,direction=">")

AUCS_max['m7_500X']=x$auc
n_max['m7_500X']=sum(binary==1)
points(x$specificities,x$sensitivities,type='l',col=rgb(0.5,0.5,0.5),lty=2) 


tab2=cbind(round(AUCS_wt,2),n_wt,
           round(AUCS_m6,2),n_m6,
           round(AUCS_m7,2),n_m7)
colnames(tab2)=c("AUC-WT", "n-WT", "AUC-M6", "n-M6", "AUC-M7", "n-M7" )

dev.off()

####### ####### ####### ####### ####### ####### ####### ####### 
####### CODE for SUPP FIGURE3
####### ####### ####### ####### ####### ####### ####### ####### 

#######  Supp fig 3A A Gini index

colo1=as.matrix(read.delim('/input/COLO1_RUNMERGED_EXACT_ANNOTATED.txt'))

selec=c("lib.COLO.1",
        "SIDM00136_CPID2437","SIDM00136_CPID2440","SIDM00136_CPID2443",
        "SIDM00136_CPID1020","SIDM00136_CPID1023","SIDM00136_CPID1026")

gini_plm=counts[,'lib.COLO.1']

library("ineq")

##curve only plasmid
plm=cbind(counts[counts[,'lib.COLO.1']>0,'lib.COLO.1'],0,0)
plm=plm[order(plm[,1]),]
colnames(plm)=c('ori','read','sgRNA')

plm[,'read']=cumsum(plm[,'ori'])/sum(plm[,'ori'])
plm[,'sgRNA']=1:nrow(plm)
plm[,'sgRNA']=plm[,'sgRNA']/max(plm[,'sgRNA'])

pdf('/pathFortheFigure/FigSupp3A_plasmidCOLO_reads.pdf',width=5,height=5)

plot(plm[,'sgRNA'],plm[,'read'],type='l',col=rgb(0.2,0.7,0.4),lwd=2,
     xlab='Fraction of total sgRNAs',ylab='Fraction of total reads')
legend('topleft',paste('GINI:',round(ineq::Gini(counts[counts[,'lib.COLO.1']>0,'lib.COLO.1']),2)))
points(c(0,1),c(0,1),type='l')
dev.off()

##### Supplementary 3B
stat=as.matrix(read.delim('ANALYSIS/ENCORE/input/All.stats'))
plmd=as.matrix(read.delim('ANALYSIS/ENCORE/input/lib-COLO-1.stats',sep=' '))
colnames(plmd)=c('fet','plasmid','plasmid_total','plasmid_rel')

colnames(stat)=gsub('\\.1','_total',colnames(stat))
colnames(stat)=gsub('\\.2','_rel',colnames(stat))


temp=stat[,1]
stat=cbind(apply(stat[,grep("SIDM00136_CPID2437",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID2440",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID2443",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1020",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1023",colnames(stat))],2,as.numeric),
           apply(stat[,grep("SIDM00136_CPID1026",colnames(stat))],2,as.numeric))
rownames(stat)=temp

vect=cbind(as.numeric(plmd[,'plasmid_rel']),
           stat[,grep('_rel',colnames(stat))])


pdf('pathFortheFigure/FigSupp3B_plasmidCOLO_reads.pdf',width=10,height=5)

x=barplot(vect*100,beside=T,ylim=c(0,100),las=3,legend=c('gRNA1+gRNA2','SWAP','sgRNA1 only','sgRNA2 only','No match'),
          names=c('Plasmid','cas9_neg (r1)','cas9_neg (r2)','cas9_neg (r3)',
                  'HT29 (r1)','HT29 (r2)','HT29 (r3)'),ylab='HT29 reads (% total reads)')
text(x[1,],vect[1,]*100,round(vect[1,]*100,2),pos=4,srt=90)
text(x[2,],vect[2,]*100,round(vect[2,]*100,2),pos=4,srt=90)
text(x[3,],vect[3,]*100,round(vect[3,]*100,2),pos=4,srt=90)
text(x[4,],vect[4,]*100,round(vect[4,]*100,2),pos=4,srt=90)
text(x[5,],vect[5,]*100,round(vect[5,]*100,2),pos=4,srt=90)

dev.off()


############
#### Supplementary figure 3C

score_ag=readRDS('/inpout/SCORE_median_lfc.rds')

colo1=readRDS('/input/colo1_c91.rds')

IDs_colo1=readRDS('/input/IDs_colo1.rds')
sel=IDs_colo1[,'Note']%in%c('AnchorSingletons','LibrarySingletons','NegativeControls','PositiveControls')
IDs_colo1=IDs_colo1[sel,]
colo1=colo1[sel,]

#####Agregamos ahora nuestros datos 

colo1_ag=data.frame(gene=c(IDs_colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),'Gene2'],
                           IDs_colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),'Gene1']),
                    lfc1=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),1],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),1]),
                    lfc2=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),2],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),2]),
                    lfc3=c(colo1[IDs_colo1[,'Note']%in%c('AnchorSingletons'),3],
                           colo1[IDs_colo1[,'Note']%in%c('LibrarySingletons'),3]),row.names=NULL)
colo1_ag= aggregate(. ~ gene, data=colo1_ag, median)
colo1_ag=as.matrix(colo1_ag)
temp=colo1_ag[,1]
colo1_ag=apply(colo1_ag[,2:4],2,as.numeric)
rownames(colo1_ag)=temp

colo1_ag=colo1_ag[rownames(colo1_ag)%in%names(score_ag),]
score_ag=score_ag[names(score_ag)%in%c(rownames(colo1_ag))]

colo1_ag=apply(colo1_ag,1,median)

colo1_ag=cbind(colo1_ag,score_ag[names(colo1_ag)])

x=c(colo1_ag[,1])
y=c(colo1_ag[,2])

##500X Wildtype scaffold
data <- data.frame(X = x, Y = y)
# Fit linear model
reg<-lm(Y ~ X, data = data)

pdf("/pathFortheFigure/FigSupp3C_COLOvsSCORE.pdf",height=5, width=5)

plot(data,xlim=c(-4,2),ylim=c(-4,2),
     main='Gene level correlation single KO',xlab='logFC(ENCORE-HT29)',ylab='logFC(Score-HT29)',
     col=rgb(0.2,0.5,0.8),pch=16)
abline(h=0,v=0,col=rgb(0.5,0.5,0.5),lty=2)
abline(reg, col =rgb(0.7,0.5,0.2,0.7),lwd=2)

legend('topleft',c(paste('cor:',round(cor(data[,1],data[,2]),2)),
                   paste('reg:',round(coef(reg)[2],2)),
                   paste('n:',nrow(data))))
dev.off()

############
#### Supplementary figure 3D

colo1=readRDS('/ENCORE/input/colo1_plm.rds')
IDs=readRDS('/ENCORE/input/IDs_colo1.rds')
colo1=colo1[,grep('136',colnames(colo1))]
colo1=apply(colo1[,1:3],1,mean)

pos=names(colo1)[IDs[,'Note']=='PositiveControls']
neg=names(colo1)[IDs[,'Note']=='NegativeControls']

pdf('/pathFortheFigure/FigSupp3D',width=5,height=5.5)

plot(0,type='n',ylim=c(0,1),xlim=c(1,0),xlab=c("Sensitivity"),ylab=c('Specificity'),
     main='HT29-project ENCORE (pos. & neg. ctrl)')
stype_df=sort(colo1[c(pos,neg)])
binary=rep(0,length(stype_df))
binary[names(stype_df)%in%pos]=1

x=roc(	binary,stype_df,direction=">")

points(x$specificities,x$sensitivities,type='l',col=rgb(0.75,0.2,0.2))

text(0.6,0.6,'AUC: 0.9313')

dev.off()














