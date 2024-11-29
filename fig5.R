args<-commandArgs(T)
datloc=args[1] ###location of the rmats output files 

filenam="SE.MATS.JC.txt"
totalfiles<-list.files(path=datloc,pattern = "^output*")
totalfilenams<-sub("output-","",totalfiles)
####function####
pieplt<-function(x,labelx,mainlab){
  piepercent<- round(100*x/sum(x), 1)
  label<-paste(labelx,"(",piepercent,"%)",sep="")
  print(pie(x,labels = label,main = mainlab, col=brewer.pal(5, "Set3")))
}
mousegene2KEGG<-function(v){
  DEG.entrez_id = mapIds(x = org.Mm.eg.db, keys = v, keytype = "SYMBOL", column = "ENTREZID")
  DEG.entrez_id = na.omit(DEG.entrez_id)
  enrich.kegg = enrichKEGG(gene = DEG.entrez_id,
                           organism = "mouse",
                           pAdjustMethod = "BH",
                           qvalueCutoff  = 0.05)
  enrichkegg1<-as.data.frame(enrich.kegg)
  if(dim(enrichkegg1)[1]>0){ 
    enrichkegg1$geneID<-as.character(enrichkegg1$geneID)
    for(g in 1:length(enrichkegg1$geneID)){ 
      gene0 <- enrichkegg1$geneID[g]
      gene1<-strsplit(gene0,split="/",fixed = T)
      gene2<-unlist(gene1)
      gensymbol0 = mapIds(x = org.Mm.eg.db, keys = gene2, keytype = "ENTREZID", column = "SYMBOL" )
      gensymbol =paste(as.vector(gensymbol0),collapse = "/")
      enrichkegg1$geneID[g]=as.character(gensymbol)
    }
  }
  return(enrichkegg1)
}


####fig5b####
for(ff in 1:length(totalfiles)){
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
sigevent<-f[f$PValue<0.05,] #f$PValue<0.05
sigevent$lab<-ifelse(sigevent$IncLevelDifference>0,"WT",ifelse(sigevent$IncLevelDifference<0,"KO","other"))

ord<-c("SE","A5SS","A3SS","MXE","RI")
inp<-table(sigevent$type[sigevent$lab=="WT"])
pieplt(inp[ord],ord,"WT")
inp<-table(sigevent$type[sigevent$lab=="KO"])
pieplt(inp[ord],ord,"KO")

####fig5c####
koevent<-sigevent[sigevent$IncLevelDifference<0,] 
koevents<-unique(koevent$geneSymbol)
koeventskegg<-mousegene2KEGG(as.character(koevents))
}




