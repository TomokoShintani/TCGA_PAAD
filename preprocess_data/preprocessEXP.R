library(curatedTCGAData)
library(RaggedExperiment)
library(purrr)
library(dplyr)
library(stringr)

curatedTCGAData(diseaseCode="PAAD", version="2.1.1")

PAADdata <- curatedTCGAData(diseaseCode="PAAD",
                            assays=c("RNASeq2GeneNorm"),
                            dry.run = FALSE, version="2.1.1")

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

#View(dfcolData)

#特殊なケースを除くために、patient.histological.typeで絞る
unique(dfcolData$patient.histological.type)
dfcolData <- filter(dfcolData, patient.histological.type %in% c("pancreas-adenocarcinoma-other subtype","pancreas-adenocarcinoma ductal type"))
#dim(dfcolData) [1] 179 369 ちゃんと行数減ってる

#assayごとにデータを分類
assay_lst <- assays(PAADdata)


Exp <- assay_lst[["PAAD_RNASeq2GeneNorm-20160128"]]
sum(is.na(Exp))

colnames(Exp) <- gsub("\\-[[:digit:]][[:digit:]]", "", colnames(Exp))
Exp <- log2(Exp + 1)

dim(Exp) #[1] 18465   183
length(unique(colnames(Exp))) #[1] 177
Exp <- Exp[, unique(colnames(Exp))]
colnames(Exp) <- gsub("\\-", ".", colnames(Exp))
Exp <- Exp[, intersect(colnames(Exp), rownames(dfcolData))]
Exp <- t(Exp)

SDs <- apply(Exp,2,sd)　#2番目の引数で2を指定しているので、各列に対して関数を適用することになる。
orderedGenes <- names(SDs[order(SDs, decreasing = TRUE)[1:3000]])
Exp <- Exp[, orderedGenes]

#export
write.csv(Exp, "C:/Users/tomsh/OneDrive/デスクトップ/研究テーマ/TCGAPAADdata/preprocess_data/PAAD_Exp.csv", row.names = FALSE)