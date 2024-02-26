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

curatedTCGAData(diseaseCode="PAAD", version="2.1.1")#これでassaysの種類を確認

PAADdata <- curatedTCGAData(diseaseCode="PAAD",
                            assays=c("CNASNP","CNVSNP","GISTIC_ThresholdedByGene","Methylation","Mutation","RNASeq2GeneNorm"),
                            dry.run = FALSE, version="2.1.1")
#View(colData(PAADdata)) どういう階層になってるのか見てみた。
dim(colData(PAADdata))  #[1] 185 979

table(PAADdata@colData@listData["histological_type"])

## make data frame
dfcolData <- data.frame(matrix(nrow = nrow(colData(PAADdata)))) #これに列方向にjoinしていく
for(i in 1:ncol(colData(PAADdata))){dfcolData <- cbind(dfcolData, PAADdata@colData@listData[i])}
dfcolData <- dfcolData[,-1] #一番目の列を削除
#dim(dfcolData) [1] 185 979
rownames(dfcolData) <- PAADdata@colData@rownames
colnames(dfcolData) <- gsub("\\-|\\_", ".", colnames(dfcolData)) #列名に含まれる-と_を.に変換。gsubはそういう関数。
rownames(dfcolData) <- gsub("\\-|\\_", ".", rownames(dfcolData)) #行名に対しても同様の操作をする

toBeRemoved <- which(data.frame(map(dfcolData, ~sum(is.na(.))))> nrow(dfcolData)*0.3)# Rの~はlambda関数みたいなもん
  #--->からのデータが全行の30%以上を占めるコラムを削除
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
sum(is.na(CNA))
  #--->空のセル多すぎるけど、、、
#colnames(CNA) <- gsub("(\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:alnum:]]+\\-[[:digit:]]+)$", "", colnames(CNA))
  #--->CNAのコラム名の"-01A-11D-A40V-01"とかの部分を""に置き換える
#CNA <- log2(CNA +1)
  #--->CNAのすべてのセルに対して対数変換を行う
#dim(CNA) #[1] 203871    368

#length(unique(colnames(CNA))) #[1] 185　複数のコラムが重複している 
#CNA <- CNA[,unique(colnames(CNA))]
#dim(CNA) #[1] 203871    185

#colnames(CNA) <- gsub("\\-", ".", colnames(CNA))
  #--->"TCGA-YY-A8LH"→"TCGA.YY.A8LH"

#CNA <- CNA[, intersect(colnames(CNA), rownames(dfcolData))]
  #--->今の状態のCNAデータには、全部のhistological typeのデータが含まれているので
       #"pancreas-adenocarcinoma-other subtype","pancreas-adenocarcinoma ductal type"
       #のみに絞る
#CNA <- t(CNA)


##Copy number variation
CNV <- assay_lst[["PAAD_CNVSNP-20160128"]]
sum(is.na(CNV))
  #--->こっちもほぼ空だ
#toBeRemoved_CNV <- which(data.frame(map(CNV, ~sum(is.na(.))))> nrow(dfcolData)*0.3)
#CNV_new <- CNV[, toBeRemoved_CNV]

##Gene expression
Exp <- assay_lst[["PAAD_RNASeq2GeneNorm-20160128"]]
sum(is.na(Exp))
  #--->発現量のデータはNA少ないね
colnames(Exp) <- gsub("\\-[[:digit:]][[:digit:]]", "", colnames(Exp))
  #--->Expのコラム名の"-01"とかの部分を""に置き換える。まあほとんど-01だけどたまに-11がある。
Exp <- log2(Exp + 1)
  #--->Expのすべてのセルに対して対数変換を行う。警告出るのはNAのデータがあるからだから、無視
dim(Exp) #[1] 18465   183
length(unique(colnames(Exp))) #[1] 177
Exp <- Exp[, unique(colnames(Exp))]
colnames(Exp) <- gsub("\\-", ".", colnames(Exp))
Exp <- Exp[, intersect(colnames(Exp), rownames(dfcolData))]
  #--->今の状態のExpデータには、全部のhistological typeのデータが含まれているので
       #"pancreas-adenocarcinoma-other subtype","pancreas-adenocarcinoma ductal type"
       #のみに絞る
Exp <- t(Exp)
#sum(is.na(Exp[, 'TP53']))

#発現変動が小さすぎる遺伝子のデータは取り除く
SDs <- apply(Exp,2,sd)　#2番目の引数で2を指定しているので、各列に対して関数を適用することになる。
orderedGenes <- names(SDs[order(SDs, decreasing = TRUE)[1:2000]])
Exp <- Exp[, orderedGenes]
  #--->なんかよく知られてる遺伝子めっちゃ消えちゃったんだけど、、、
       #でも元の論文のやつも有名な遺伝子ないね

#遺伝子発現量で降順に並べてみると？



##

##
