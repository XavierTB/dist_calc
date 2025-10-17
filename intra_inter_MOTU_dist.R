# This script computes distances intra- and inter- MOTUs after SWARM
# 
# Only MOTUs with more than one ESV are used
# For inter-MOTU distances it uses the representative sequence of each MOTU
# 
# The input file must be a .csv with all ESVs in rows. Columns can contain any information (typically values of abundance per samples), but MUST include a column "MOTU" with the code of the MOTU to which
# each ESV is assigned, and a column "seq" with the sequence of the ESV (all sequences should be aligned). 
# ESVs must be ordered by abundance, so the first ESV of each MOTU has the representative sequence


library(ape)


datas<-read.csv("input_file.csv",stringsAsFactors = F,check.names = F)

dim(datas)

length(levels(factor(datas$MOTU)))


#Select those with more than 1 ESV

pp<-tapply(datas$MOTU,factor(datas$MOTU,levels=unique(datas$MOTU)),length)

pp<-pp[pp>1]
noms<-names(pp)
  
datass<-datas[datas$MOTU%in%noms,]

dim(datass)

length(levels(factor(datass$MOTU)))

MOTU<-levels(factor(datass$MOTU,levels=unique(datass$MOTU)))

# compute intra-MOTU distances
intra<-NULL
dist<-0
for (i in 1:length(MOTU))
{
  seqs<-datass$seq[datass$MOTU==MOTU[i]]
  r<-strsplit(seqs,"", tolower)
  m<-as.DNAbin(r)
  d<-dist.dna(m,model="N",as.matrix=TRUE)
  dd<-mean(d[upper.tri(d)])
  dist<-dist+dd/length(MOTU)
  intra<-c(intra,dd)
} 
#mean intra_MOTU distance
dist

#in percent
dist*100/313

#quantiles
quantile(intra,probs=c(0.05,0.1,0.9,0.95))*100/313

length(intra)
mean(intra)


##For inter_MOTU distances we compare only the representative sequences
inter<-NULL
dist<-0
contador<-0
for (i in 1:(length(MOTU)-1))
{
  if (i/100-ceiling(i/100)==0) message("comparing MOTU ",i," of ",length(MOTU)," with the rest",Sys.time())
  for (j in (i+1):length(MOTU))
  {
  
  seqs<-rbind(datass$seq[datass$MOTU==MOTU[i]][1],datass$seq[datass$MOTU==MOTU[j]][1])
  r<-strsplit(seqs,"", tolower)
  m<-as.DNAbin(r)
  d<-dist.dna(m,model="N",as.matrix=TRUE)
  dd<-d[1,2]
  dist<-dist+dd
  contador<-contador+1
  inter<-c(inter,dd)
  } 
}
#mean inter-MOTU distance
distt<-dist/contador
distt

#in percent
distt*100/313

#quantiles
quantile(inter,probs=c(0.05,0.1,0.9,0.95))*100/313

length(inter)
mean(inter)

message("that's it")


