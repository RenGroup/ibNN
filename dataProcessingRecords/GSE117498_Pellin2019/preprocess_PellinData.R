#Preprocess of raw MPP data
pub.mpp<-read.table("/Users/XXX/neuralNetwork/Pellin2019/GSM3305360_MPP.raw_counts.tsv",
                    sep="\t",header=T)
pub.mpp<-pub.mpp[2:dim(pub.mpp)[1],] #trim the first row

library(biomaRt)
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes = pub.mpp$Barcode
genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
               values = genes ,mart = human, 
               attributesL = c("entrezgene_id"), 
               martL = human, uniqueRows=T)
count.symbol<-table(genes$HGNC.symbol) #count the number of appearances of each gene symbol
rm.symbol<-names(count.symbol)[count.symbol > 1]
genes.uniq<-genes[!genes$HGNC.symbol %in% rm.symbol,] #remove genes with multiple gene symbols
count.id<-table(genes.uniq$NCBI.gene..formerly.Entrezgene..ID)
rm.id<-names(count.id)[count.id > 1]
genes.uniq<-genes.uniq[!genes.uniq$NCBI.gene..formerly.Entrezgene..ID %in% rm.id,] #remove genes with multiple gene ids

pub.mpp.uniq<-merge(genes.uniq,pub.mpp,by.x = "HGNC.symbol",by.y = "Barcode")
pub.mpp.uniq<-na.omit(pub.mpp.uniq)

pub.mpp.uniq[,3:dim(pub.mpp.uniq)[2]]<-log2(pub.mpp.uniq[,3:dim(pub.mpp.uniq)[2]]+1)
write.table(pub.mpp.uniq,file = "/Users/chix/neuralNetwork/Pellin2019/exprMa_mpp_rmDup_pellin_github.txt",
            sep="\t",col.names = T,row.names = F,quote = F)

#Preprocess of raw HSC data
pub.hsc<-read.table("/Users/chix/neuralNetwork/Pellin2019/GSM3305359_HSC.raw_counts.tsv",
                    sep="\t",header=T)
pub.hsc<-pub.hsc[2:dim(pub.hsc)[1],] #trim the first row

genes = pub.hsc$Barcode
genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
               values = genes ,mart = human, 
               attributesL = c("entrezgene_id"), 
               martL = human, uniqueRows=T)
count.symbol<-table(genes$HGNC.symbol)
rm.symbol<-names(count.symbol)[count.symbol > 1]
genes.uniq<-genes[!genes$HGNC.symbol %in% rm.symbol,] #remove genes with multiple gene symbols
count.id<-table(genes.uniq$NCBI.gene..formerly.Entrezgene..ID)
rm.id<-names(count.id)[count.id > 1]
genes.uniq<-genes.uniq[!genes.uniq$NCBI.gene..formerly.Entrezgene..ID %in% rm.id,] #remove genes with multiple gene ids

pub.hsc.uniq<-merge(genes.uniq,pub.hsc,by.x = "HGNC.symbol",by.y = "Barcode")
pub.hsc.uniq<-na.omit(pub.hsc.uniq)

pub.hsc.uniq[,3:dim(pub.hsc.uniq)[2]]<-log2(pub.hsc.uniq[,3:dim(pub.hsc.uniq)[2]]+1)
write.table(pub.hsc.uniq,file = "/Users/chix/neuralNetwork/Pellin2019/exprMa_hsc_rmDup_pellin_github.txt",
            sep="\t",col.names = T,row.names = F,quote = F)