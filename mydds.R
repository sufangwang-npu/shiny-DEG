library(DESeq2)
mydds <- function(inFile,condition){
 
  #header is too long,substitue gene.results with nothing
  ##colnames( inFile ) <- sub("(*).genes.results","\\1",colnames(inFile))
  #make all numbers into integer
  mycountData <- round(inFile,digits=0)
  data1Design <- data.frame(row.names = colnames( mycountData ),condition = condition)
  mycolData <- data1Design
  mydds <- DESeqDataSetFromMatrix(countData = mycountData, colData = mycolData,design = ~condition)
  #only keep expression level >100, minimize variance
  mydds <- mydds[ rowSums(counts(mydds)) > 10, ]
  mydds$condition<-relevel(mydds$condition,ref="ctrl")
  #mydds$condition <- factor(mydds$condition, levels=c("control","exp"))
  mydds <- DESeq(mydds)
  myres <- results(mydds)
  return(myres)
}