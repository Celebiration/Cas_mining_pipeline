#!/usr/bin/Rscript
library(ggplot2)
library(Biostrings)
library(jsonlite)
Args<-commandArgs(trailingOnly=T)
dir<-Args[1]
file<-Args[2]
cut_off<-abs(as.integer(Args[3]))
lines1<-'【------------------------------'
lines2<-'------------------------------】'
print(paste0(lines1,"整合结果...",lines2))
data<-read.table(file,stringsAsFactors = F,header = F,quote = '"')
data<-data[,c(1:5,7,9:25)]
names(data)<-c("id","hmm_start","hmm_end","hmm_strand","contig_len",'hmm_query_name',"full_E-value","full_score","full_bias","best_E-value","best_score","best_bias","exp","reg","clu","ov","env","dom","rep","inc","target_description","contig_seq","protein_seq")
data[data$hmm_strand=="+","hmm_strand"]<-1
data[data$hmm_strand=="-","hmm_strand"]<-(-1)
plotdata1<-cbind(data,integraity=!grepl("partial",data$'target_description'))
crispr<-read_json(paste0(dir,"/hitted_CRISPR_out/result.json"),simplifyVector = T)
crispr2<-crispr$Sequences[unlist(lapply(crispr$Sequences$Crisprs, length))>=1,]
crispr3<-do.call("rbind",crispr2$Crisprs)
mytrim<-function(x,sep="_"){
  tmp<-unlist(strsplit(x,sep))
  return(paste0(tmp[1:length(tmp)-1],collapse = sep))
}
tmp<-unlist(lapply(crispr3$Name, mytrim))
crispr3<-cbind(tmp,crispr3)
if (is.null(crispr3)) {
  ssss<-c("tmp","Name","Start","End","DR_Consensus","Repeat_ID","DR_Length","Spacers","Potential_Orientation","CRISPRDirection","Evidence_Level","Conservation_DRs","Conservation_Spacers","Regions")
  crispr3<-data.frame(matrix(nrow=0,ncol=length(ssss)))
  names(crispr3)<-ssss
}
final<-merge(plotdata1,crispr3,by.x="id",by.y="tmp",all.x=T)
hmm_seq<-DNAStringSet(unlist(apply(final[,c('contig_seq','hmm_start','hmm_end','hmm_strand')],1,function(x){if (as.integer(x[4])==1){return(DNAString(substr(x[1],as.integer(x[2]),as.integer(x[3]))))}else{return(reverseComplement(DNAString(substr(x[1],as.integer(x[2]),as.integer(x[3])))))}})))
final<-cbind(final[,1:21],hmm_seq=hmm_seq,final[,22:ncol(final)])
print(paste0(lines1,"保存结果...",lines2))
save(final,file=paste0(dir,"/all_result.RData"))
final<-cbind(final[,1:4],aa_length=(final$hmm_end-final$hmm_start+1)/3,final[,5:ncol(final)])
write.csv(subset(final,select=-which(names(final) %in% c("contig_seq","Regions"))),paste0(dir,"/all_summary.csv"),row.names = F,quote = T)
print(paste0(lines1,"画图...",lines2))
maxx<-round(max(-log10(plotdata1$'full_E-value'[plotdata1$'full_E-value' != 0])))
maxx<-((maxx+9) %/% 10)*10
maxx<-max(maxx,50)
if (any(plotdata1$'full_E-value' == 0)){
  breaks = c(seq(0,maxx,10),maxx+40)
  labels = c(as.character(seq(0,maxx,10)),"E-value=0")
}else{
  breaks = c(seq(0,maxx,10))
  labels = as.character(seq(0,maxx,10))
}
graphics.off()
pdf(paste0(dir,"/cas12bij hmmsearch E-value distribution.pdf"),width = 10,height = 6)
ggplot(data=plotdata1,aes(x=-log10(`full_E-value`+1e-300),group=hmm_query_name,fill=hmm_query_name))+
  geom_histogram(binwidth = 1,position = "stack")+
  coord_cartesian(ylim = c(0,50))+
  theme_bw()+
  labs(x="-log10(E-value)")+
  ggtitle("cas12bij hmmsearch E-value distribution")+
  scale_x_continuous(breaks = breaks,labels = labels)+
  scale_fill_manual(values = c("#8AC486","#91ADDB","#EF5526","#F7C961"),name=NULL)+
  geom_vline(xintercept=cut_off)
dev.off()

