#install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedTCGAData")
BiocManager::install("RaggedExperiment")
install.packages("purrr")
install.packages("stringr")

library(curatedTCGAData)
library(RaggedExperiment)
library(purrr)
library(dplyr)
library(stringr)

curatedTCGAData(deseaseCode="PAAD", version="1.1.38")#これでassaysの種類を確認

PAADdata <- curatedTCGAData("PAAD",
                            assays=c("CNASNP","GISTIC_ThresholdedByGene","Methylation","Mutation","RNASeq2GeneNorm"),
                            dry.run = FALSE, version="1.1.38")

dim(colData(PAADdata))  #[1] 185 979

## col_data for graph annotation
dfcolData <- data.frame(matrix(nrow = nrow(colData(PAADdata))))
for(i in 1:ncol(colData(PAADdata))){dfcolData <- cbind(dfcolData, PAADdata@colData@listData[i])}
dfcolData <- dfcolData[,-1] #一番目の列を削除
#dim(dfcolData) [1] 185 979
rownames(dfcolData) <- PAADdata@colData@rownames
colnames(dfcolData) <- gsub("\\-|\\_", ".", colnames(dfcolData)) #列名に含まれる-と_を.に変換。gsubはそういう関数。
rownames(dfcolData) <- gsub("\\-|\\_", ".", rownames(dfcolData)) #行名に対しても同様の操作をする

toBeRemoved <- which(data.frame(map(dfcolData, ~sum(is.na(.))))> nrow(dfcolData)*0.3)# Rの~はlambda関数みたいなもん
dfcolData <- dfcolData[,-(toBeRemoved)]
#dim(dfcolData) [1] 185 369

View(dfcolData)

#特殊なケースを除くために、patient.histological.typeで絞る
unique(dfcolData$patient.histological.type)
dfcolData <- filter(dfcolData, patient.histological.type %in% c("pancreas-adenocarcinoma-other subtype","pancreas-adenocarcinoma ductal type"))
#dim(dfcolData) [1] 179 369 ちゃんと行数減ってる

#assayごとにデータを分類
assay_lst <- assays(PAADdata)

##Copy number alteration
CNA <- assay_lst[["PAAD_CNASNP-20160128"]] #assay_lst[[1]]と同じ。CNAほとんどのセルがNAなんだけど、大丈夫なん？
colnames(CNA) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:digit:]]+)$", "", colnames(CNA))
  #--->CNAのコラム名の"-01A-11D-A40V-01"とかの部分を""に置き換える
CNA <- log2(CNA +1)
  #--->CNAのすべてのセルに対して対数変換を行う
dim(CNA) #[1] 203871    368

length(unique(colnames(CNA))) #[1] 185　複数のコラムが重複している 
CNA <- CNA[,unique(colnames(CNA))]
dim(CNA) #[1] 203871    185

colnames(CNA) <- gsub("\\-", ".", colnames(CNA))
  #--->"TCGA-YY-A8LH"→"TCGA.YY.A8LH"

CNA <- CNA[, intersect(colnames(CNA), rownames(dfcolData))]
  #--->今の状態のCNAデータには、全部のhistological typeのデータが含まれているので
       #"pancreas-adenocarcinoma-other subtype","pancreas-adenocarcinoma ductal type"
       #のみに絞る

##

##

##

##
