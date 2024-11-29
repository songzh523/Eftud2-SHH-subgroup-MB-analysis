args<-commandArgs(T)
datloc=args[1] ###location of the rmats output files 

filenam="SE.MATS.JC.txt"
totalfiles<-list.files(path=datloc,pattern = "^output*")
totalfilenams<-sub("output-","",totalfiles)

#################Fig1a,c,e,g#####################
for(ff in 1:length(totalfiles)){
  
datloc1=paste(datloc,totalfiles[ff],"/",sep="")
f<-read.table(paste(datloc,totalfiles[ff],"/",filenam,sep=""),header = T)
sigevent<-f[f$PValue<0.05,] 
wtpsi<-unlist(strsplit(as.character(sigevent$IncLevel1),",",fixed = T))
kopsi<-unlist(strsplit(as.character(sigevent$IncLevel2),",",fixed = T))
df<-as.data.frame(cbind(wtpsi,kopsi))
colnames(df)<-c("Control","KO")
df$Control<-as.numeric(as.character(df$Control))
df$KO<-as.numeric(as.character(df$KO))
df1<-reshape2::melt(df)

my_comparisons<-list(c("Control","KO"))
p<-ggboxplot(df1,"variable","value",fill="variable",xlab="",ylab = "PSI")+theme(legend.position="none")
print(p+stat_compare_means(comparisons = my_comparisons,label = "p.signif"))
}

#################Fig1b,d,f,h,Fib5a#####################
datloc1=paste(datloc,totalfiles[ff],"/",sep="")
filenams=list.files(datloc1,"*.MATS.JC.txt") 
type<-sub("\\..*","",filenams)
keps<-c("GeneID","geneSymbol","ID","PValue","FDR","IncLevel1","IncLevel2","IncLevelDifference")
for(i in 1:length(filenams)){
f0<-read.table(paste(datloc1,filenams[i],sep = ""),header = T)
if(i==1){
  f=f0[,keps]
  f$type<-type[i]
}else{
  tmp<-f0[,keps]
  tmp$type<-type[i]
  f=rbind(f,tmp)
  rm(tmp)
}
}
sigevent<-f[f$PValue<0.05,] 
sigevent$lab<-ifelse(sigevent$IncLevelDifference>0,"WT",ifelse(sigevent$IncLevelDifference<0,"KO","other"))

inp1<-round(prop.table(table(sigevent$type[sigevent$lab=="WT"])),4)*100
inp2<-round(prop.table(table(sigevent$type[sigevent$lab=="KO"])),4)*100

inp1<-as.data.frame(inp1)
inp1$Var2<-"WT"
inp2<-as.data.frame(inp2)
inp2$Var2<-"KO"
percents<-rbind(inp1,inp2)
percents$Var2<-factor(percents$Var2,levels = c("WT","KO"))

stackplot(percents,ord,"significantASevent",brewer.pal(5, "Set3"))


