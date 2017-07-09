#intake for covar, of clustal aligned fasta or standard fasta.
coveIn=function(filepath, filename){
if (missing(filename)){
  file=dir(filepath)

}else{
file = dir(filepath, pattern=filename)
}
  
preseq = scan(paste0(filepath, file),what=character(),sep="\n")
presplit=sapply(preseq, tolower)
seqs=strsplit(presplit,split="")

seqList=list()
count=0
for (i in 1:length(seqs)){
  if (">" %in% seqs[[i]]){
  count=count+1

  }
}


seqList=list()
count=0

for (i in 1:length(seqs)){

  if (">" %in% seqs[[i]])
  {
    count=count+1
  }
  else{
    if (count==(length(seqList)+1))
    {
    seqList[[count]]=seqs[[i]] 
    }
    else{
    seqList[[count]]=c(seqList[[count]],seqs[[i]])
    }
  }

}
seqList[is.na(seqList)]='-'
seqList

}