############# Load DEseq2 and data ##########


source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
biocLite("pcaExplorer")
library("pcaExplorer")
library("edgeR")
library("fdrtool")
countData <- read.table ("counts.txt", header=TRUE)
nrow(countData)
head(countData)

############# Pairwise DE ##########

spDesign = data.frame (
  row.names = colnames(countData),
  age = c("67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"67",	"32",	"32",	"32",	"32",	"32",	"32",	"32",	"32",	"32",	"32",	"32",	"32"),	
  tutor = c("tutor", "tutor",	"tutor",	"isolate",	"isolate",	"isolate",	"tutor",	"isolate",	"isolate",	"tutor",	"isolate",	"isolate",	"isolate",	"isolate",	"tutor",	"isolate",	"tutor",	"tutor",	"tutor",	"tutor",	"isolate",	"tutor",	"isolate",	"tutor",	"tutor",	"isolate",	"tutor",	"isolate",	"isolate",	"tutor",	"isolate",	"tutor",	"isolate",	"tutor",	"isolate",	"tutor"),
  song =  c("silence", "silence",	"song",	"song",	"silence",	"silence",	"song",	"song",	"silence",	"song",	"silence",	"song",	"song",	"song",	"silence",	"silence",	"silence",	"song",	"song"	, "silence",	"song",	"silence",	"silence",	"song",	"song",	"silence",	"song",	"song",	"song",	"song",	"silence",	"silence",	"silence",	"silence",	"song",	"silence"))
  spDesign

#P32 vs P67 (condition and song collapsed)
ddsAge<- DESeqDataSetFromMatrix (countData = countData,
                                     colData = spDesign,
                                     design = ~ age)

ddsAge <- DESeq(ddsAge)
resultsNames(ddsAge)

resAge <- results(ddsAge, name="age_67_vs_32")
pcaExplorer(dds = ddsAge)

head(resAge)
DESeq2::plotMA(resAge, alpha=0.5)
sum(resAge$padj < 0.1, na.rm=TRUE) ## 1371
sum(resAge$padj < 0.05, na.rm=TRUE) ## 895 genes
sum(resAge$pvalue < 0.01, na.rm=TRUE) ## 1392 genes

hist(resAge$pvalue, main="DESeq example data", xlab="p-values",cex.main=1.5)
#p value histogram looks good

write.csv(as.data.frame(resAge), file = "p32vp67.csv")

#Tutored vs Isolate (age and song collapsed)
ddsTutor<- DESeqDataSetFromMatrix (countData = countData,
                                 colData = spDesign,
                                 design = ~ tutor)

ddsTutor<- DESeq(ddsTutor)
resultsNames(ddsTutor)

resTutor <- results(ddsTutor, name="tutor_tutor_vs_isolate")
head(resTutor)
DESeq2::plotMA(resTutor, alpha=0.5)
sum(resTutor$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(resTutor$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(resTutor$pvalue < 0.01, na.rm=TRUE) ## 106 genes

hist(resTutor$pvalue, main="DESeq example data", xlab="p-values",cex.main=1.5)
write.csv(as.data.frame(resTutor), file = "TutorVsIsolate.csv")

#check pvalue histogram and correct with FDR tool

#removeNA
resTutor <- resTutor[ !is.na(resTutor$padj), ]
resTutor <- resTutor[ !is.na(resTutor$pvalue), ]
#removepadj for replacement
resTutor <- resTutor[, -which(names(resTutor) == "padj")]

FDR.resTutor<- fdrtool(resTutor$stat, statistic= "normal", plot = T)

FDR.resTutor$param[1, "sd"]

resTutor[,"padj"]  <- p.adjust(FDR.resTutor$pval, method = "BH")

hist(FDR.resTutor$pval, col = "royalblue4", 
     main = "Tutor vs Iso", xlab = "CORRECTED p-values")

table(resTutor[,"padj"] < 0.1)
DESeq2::plotMA(resTutor)
write.csv(as.data.frame(resTutor), file = "TutorVsIsolate.corrected.csv")


#Song vs Silence (age and tutor collapsed)
ddsSong<- DESeqDataSetFromMatrix (countData = countData,
                                   colData = spDesign,
                                   design = ~ song)

ddsSong<- DESeq(ddsSong)
resultsNames(ddsSong)
resSong <- results(ddsSong)
resSong <- results(ddsSong, name="song_song_vs_silence")
head(resSong)
DESeq2::plotMA(resSong, alpha=0.5)
sum(resSong$padj < 0.1, na.rm=TRUE) ## 1 genes
sum(resSong$padj < 0.05, na.rm=TRUE) ## 1 genes
sum(resSong$pvalue < 0.01, na.rm=TRUE) ## 95 genes

hist(resSong$pvalue, main="DESeq example data", xlab="p-values",cex.main=1.5)

write.csv(as.data.frame(resSong), file = "SongVsSilence.csv")

#check pvalue histogram and if skewed, use FDRtool

#removeNA
resSong <- resSong[ !is.na(resSong$padj), ]
resSong <- resSong[ !is.na(resSong$pvalue), ]
#removepadj for replacement
resSong <- resSong[, -which(names(resSong) == "padj")]

FDR.resSong<- fdrtool(resSong$stat, statistic= "normal", plot = T)

FDR.resSong$param[1, "sd"]

resSong[,"padj"]  <- p.adjust(FDR.resSong$pval, method = "BH")

hist(FDR.resSong$pval, col = "royalblue4", 
     main = "Song vs Silence", xlab = "CORRECTED p-values")

table(resSong[,"padj"] < 0.1)
DESeq2::plotMA(resSong)
write.csv(as.data.frame(resSong), file = "SongVsSilence.corrected.csv")


#### interaction of age and tutoring, controlling for song

ddsALLinteraction<- DESeqDataSetFromMatrix (countData = countData,
                                            colData = spDesign,
                                            design = ~ age + song + tutor + age:tutor)

ddsres <- DESeq(ddsALLinteraction)
design(ddsres)

resultsNames(ddsres)
resALLinteraction <- results(ddsres, name="age67.tutortutor")
write.csv(as.data.frame(resALLinteraction), file = "ALLage67.tutortutor.csv")

pcaExplorer(dds = ddsp67interaction)

###################
### p67 only #####
###################

countDatap67 <- read.table ("countsp67.txt", header=TRUE)
nrow(countDatap67)
head(countDatap67)

############# Pairwise DE ##########

spDesignp67 = data.frame (
  row.names = colnames(countDatap67),
  tutor = c("tutor", "tutor", "tutor", "isolate", "isolate", "isolate", "tutor", "isolate", "isolate","tutor","isolate", "isolate", "isolate", "isolate","tutor","isolate","tutor","tutor","tutor","tutor","isolate","tutor","isolate","tutor"),
  song =  c("silence","silence","song","song","silence","silence","song","song","silence","song","silence","song","song","song","silence","silence","silence","song","song","silence","song","silence","silence","song"))
spDesignp67

#p67 Tutored vs Isolate (age and song collapsed)
ddsp67tutor<- DESeqDataSetFromMatrix (countData = countDatap67,
                                   colData = spDesignp67,
                                   design = ~ tutor)

ddsp67tutor<- DESeq(ddsp67tutor)
resultsNames(ddsp67tutor)

pcaExplorer(dds = ddsp67tutorSong)


resp67tutor <- results(ddsp67tutor, name="tutor_tutor_vs_isolate")
head(resp67tutor)
DESeq2::plotMA(resp67tutor, alpha=0.5)
sum(resp67tutor$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(resp67tutor$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(resp67tutor$pvalue < 0.01, na.rm=TRUE) ## 108 genes

write.csv(as.data.frame(resp67tutor), file = "p67TutorVsIsolate.csv")

hist(resp67tutor$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
resp67tutor <- resp67tutor[ !is.na(resp67tutor$padj), ]
resp67tutor <- resp67tutor[ !is.na(resp67tutor$pvalue), ]
#removepadj for replacement
resp67tutor <- resp67tutor[, -which(names(resSong) == "padj")]

FDR.resp67tutor<- fdrtool(resp67tutor$stat, statistic= "normal", plot = T)

FDR.resp67tutor$param[1, "sd"]

resp67tutor[,"padj"]  <- p.adjust(FDR.resp67tutor$pval, method = "BH")

hist(FDR.resp67tutor$pval, col = "royalblue4", 
     main = "Tutor cs Isolate", xlab = "CORRECTED p-values")

table(resp67tutor[,"padj"] < 0.1)
DESeq2::plotMA(resp67tutor)
write.csv(as.data.frame(resp67tutor), file = "p67TutorVsIsolate.corrected.csv")

#p67 Tutored vs Isolate (age and song collapsed), controlling for song
ddsp67tutorSong<- DESeqDataSetFromMatrix (countData = countDatap67,
                                      colData = spDesignp67,
                                      design = ~ song + tutor)

ddsp67tutorSong<- DESeq(ddsp67tutorSong)
resultsNames(ddsp67tutorSong)


resp67tutorSong <- results(ddsp67tutorSong, name="tutor_tutor_vs_isolate")
head(resp67tutorSong)
DESeq2::plotMA(resp67tutorSong, alpha=0.5)
sum(resp67tutorSong$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(resp67tutorSong$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(resp67tutorSong$pvalue < 0.01, na.rm=TRUE) ## 108 genes

write.csv(as.data.frame(resp67tutorSong), file = "p67TutorVsIsolateSong.csv")

pcaExplorer(dds = ddsp67tutor)



###p67 Song vs Silent###

ddsp67song<- DESeqDataSetFromMatrix (countData = countDatap67,
                                      colData = spDesignp67,
                                      design = ~ song)

ddsp67song<- DESeq(ddsp67song)
resultsNames(ddsp67song)

resp67song <- results(ddsp67song, name="song_song_vs_silence")
head(resp67song)
DESeq2::plotMA(resp67song, alpha=0.5)
sum(resp67song$padj < 0.1, na.rm=TRUE) ## 3 genes
sum(resp67song$padj < 0.05, na.rm=TRUE) ## 3 genes
sum(resp67song$pvalue < 0.01, na.rm=TRUE) ## 105 genes

write.csv(as.data.frame(resp67song), file = "p67songVssilence.csv")

hist(resp67song$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
resp67song <- resp67song[ !is.na(resp67song$padj), ]
resp67song <- resp67song[ !is.na(resp67song$pvalue), ]
#removepadj for replacement
resp67song <- resp67song[, -which(names(resp67song) == "padj")]

FDR.resp67song<- fdrtool(resp67song$stat, statistic= "normal", plot = T)

FDR.resp67song$param[1, "sd"]

resp67song[,"padj"]  <- p.adjust(FDR.resp67song$pval, method = "BH")

hist(FDR.resp67song$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(resp67song[,"padj"] < 0.1)
DESeq2::plotMA(resp67song)
write.csv(as.data.frame(resp67song), file = "p67songVssilence.corrected.csv")


##testing interaction term, interaction between song response and tutoring

ddsp67interaction<- DESeqDataSetFromMatrix (countData = countDatap67,
                                     colData = spDesignp67,
                                     design = ~ song + tutor + song:tutor)

ddsres <- DESeq(ddsp67interaction)
design(ddsres)

resultsNames(ddsres)
resp67interaction <- results(ddsres, name="songsong.tutortutor")
write.csv(as.data.frame(resp67interaction), file = "p67interaction.csv")

pcaExplorer(dds = ddsp67interaction)

hist(resp67interaction$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
resp67interaction <- resp67interaction[ !is.na(resp67interaction$padj), ]
resp67interaction <- resp67interaction[ !is.na(resp67interaction$pvalue), ]
#removepadj for replacement
resp67interaction <- resp67interaction[, -which(names(resp67interaction) == "padj")]

FDR.resp67interaction <- fdrtool(resp67interaction$stat, statistic= "normal", plot = T)

FDR.resp67interaction$param[1, "sd"]

resp67interaction[,"padj"]  <- p.adjust(FDR.resp67interaction$pval, method = "BH")

hist(FDR.resp67interaction$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(resp67interaction[,"padj"] < 0.1)
DESeq2::plotMA(resp67interaction)
write.csv(as.data.frame(resp67interaction), file = "p67interaction.corrected.csv")


###p32 Tutor vs Isolate ###

p32countData <- read.table ("counts_p32.txt", header=TRUE)
head(p32countData)

############# Pairwise DE ##########

p32spDesign = data.frame (
  row.names = colnames(p32countData),
  tutor = c("tutor","isolate","tutor","isolate","isolate","tutor",  "isolate","tutor","isolate","tutor","isolate","tutor"),
  song =  c("song","silence","song",  "song",    "song",  "song",   "silence","silence","silence","silence","song","silence"))
p32spDesign

#tutor
p32ddstutor<- DESeqDataSetFromMatrix (countData = p32countData,
                                   colData = p32spDesign,
                                   design = ~ tutor)

#tutor, controlling for song
p32ddstutorSong<- DESeqDataSetFromMatrix (countData = p32countData,
                                      colData = p32spDesign,
                                      design = ~ song + tutor)


#keep <- rowSums(counts(ddstutor, normalized=TRUE) >= 5) >=6 

#ddstutor5<- ddstutor[keep,]
p32ddstutor <- DESeq(p32ddstutor)
resultsNames(p32ddstutor)
p32restutor <- results(p32ddstutor)
head(p32restutor)

pcaExplorer(dds = p32ddstutor)

sum(p32restutor$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(p32restutor$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(p32restutor$pvalue < 0.01, na.rm=TRUE) ## 56 genes

write.csv(as.data.frame(p32restutor), file = "p32tutorVSiso.csv")

hist(p32restutor$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32restutor <- p32restutor[ !is.na(p32restutor$padj), ]
p32restutor <- p32restutor[ !is.na(p32restutor$pvalue), ]
#removepadj for replacement
p32restutor <- p32restutor[, -which(names(p32restutor) == "padj")]

FDR.p32restutor<- fdrtool(p32restutor$stat, statistic= "normal", plot = T)

FDR.p32restutor$param[1, "sd"]

p32restutor[,"padj"]  <- p.adjust(FDR.p32restutor$pval, method = "BH")

hist(FDR.p32restutor$pval, col = "royalblue4", 
     main = "Tutor vs Iso", xlab = "CORRECTED p-values")

table(p32restutor[,"padj"] < 0.1)
DESeq2::plotMA(p32restutor)
write.csv(as.data.frame(p32restutor), file = "p32_tutorVSiso.corrected.csv")

#tutor controlling for song
p32ddstutorSong <- DESeq(p32ddstutorSong)
resultsNames(p32ddstutorSong)
p32restutorSong <- results(p32ddstutorSong)
head(p32restutorSong)

DESeq2::plotMA(p32restutorSong, alpha=0.5)
sum(p32restutorSong$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(p32restutorSong$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(p32restutorSong$pvalue < 0.01, na.rm=TRUE) ## 45 genes

write.csv(as.data.frame(p32restutorSong), file = "p32TutorVsIsolateSong.csv")

pcaExplorer(dds = ddsp67tutor)
###p32  Song vs Silent ###

#
p32ddssong<- DESeqDataSetFromMatrix (countData = p32countData,
                                      colData = p32spDesign,
                                      design = ~ song)

#keep <- rowSums(counts(ddssong, normalized=TRUE) >= 5) >=6 


p32ddssong <- DESeq(p32ddssong)
resultsNames(p32ddssong)
p32ressong <- results(p32ddssong)
head(p32ressong)

sum(p32ressong$padj < 0.1, na.rm=TRUE) ## 1 genes
sum(p32ressong$padj < 0.05, na.rm=TRUE) ## 1 genes
sum(p32ressong$pvalue < 0.01, na.rm=TRUE) ## 96 genes

write.csv(as.data.frame(p32ressong), file = "p32songVSsilence.csv")

hist(p32ressong$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32ressong <- p32ressong[ !is.na(p32ressong$padj), ]
p32ressong <- p32ressong[ !is.na(p32ressong$pvalue), ]
#removepadj for replacement
p32ressong <- p32ressong[, -which(names(p32ressong) == "padj")]

FDR.p32ressong<- fdrtool(p32ressong$stat, statistic= "normal", plot = T)

FDR.p32ressong$param[1, "sd"]

p32ressong[,"padj"]  <- p.adjust(FDR.p32ressong$pval, method = "BH")

hist(FDR.p32ressong$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p32ressong[,"padj"] < 0.1)
DESeq2::plotMA(p32ressong)
write.csv(as.data.frame(p32ressong), file = "p32_songVSiso.corrected.csv")


###p67tutored Song vs Silent ###

countDatap67tut <- read.table ("countsp67tutored.txt", header=TRUE)
nrow(countDatap67tut)
head(countDatap67tut)

############# Pairwise DE ##########

spDesignp67tut = data.frame (
  row.names = colnames(countDatap67tut),
  song =  c("song", "song", "song", "song", "song", "song", "silence", "silence", "silence", "silence", "silence", "silence"))
spDesignp67tut

p67tutddssong<- DESeqDataSetFromMatrix (countData = countDatap67tut,
                                      colData = spDesignp67tut,
                                      design = ~ song)


p67tutddssong <- DESeq(p67tutddssong)
resultsNames(p67tutddssong)
p67tutressong <- results(p67tutddssong)
head(p67tutressong)

sum(p67tutressong$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(p67tutressong$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(p67tutressong$pvalue < 0.01, na.rm=TRUE) ## 152 genes

write.csv(as.data.frame(p67tutressong), file = "p67TuTsongVSsilence.csv")

hist(p67tutressong$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool (this one is not terrible)

#removeNA
p67tutressong <- p67tutressong[ !is.na(p67tutressong$padj), ]
p67tutressong <- p67tutressong[ !is.na(p67tutressong$pvalue), ]
#removepadj for replacement
p67tutressong <- p67tutressong[, -which(names(p67tutressong) == "padj")]

FDR.p67tutressong<- fdrtool(p67tutressong$stat, statistic= "normal", plot = T)

FDR.p67tutressong$param[1, "sd"]

p67tutressong[,"padj"]  <- p.adjust(FDR.p67tutressong$pval, method = "BH")

hist(FDR.p67tutressong$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p67tutressong[,"padj"] < 0.1) #289
DESeq2::plotMA(p67tutressong)
write.csv(as.data.frame(p67tutressong), file = "p67TuT_songVSsilence.corrected.csv")

################################
###p67isolate Song vs Silent ###
################################

countDatap67iso <- read.table ("countsp67iso.txt", header=TRUE)
nrow(countDatap67iso)
head(countDatap67iso)

############# Pairwise DE ##########

spDesignp67iso = data.frame (
  row.names = colnames(countDatap67iso),
  song =  c("song", "song", "song", "song", "song", "song", "silence", "silence", "silence", "silence", "silence", "silence"))
spDesignp67iso

p67isoddssong<- DESeqDataSetFromMatrix (countData = countDatap67iso,
                                        colData = spDesignp67iso,
                                        design = ~ song)


p67isoddssong <- DESeq(p67isoddssong)
resultsNames(p67isoddssong)
p67isoressong <- results(p67isoddssong)
head(p67isoressong)

sum(p67isoressong$padj < 0.1, na.rm=TRUE) ## 0 genes
sum(p67isoressong$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(p67isoressong$pvalue < 0.01, na.rm=TRUE) ## 76 genes

write.csv(as.data.frame(p67isoressong), file = "p67IsosongVSsilence.csv")

hist(p67isoressong$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p67isoressong <- p67isoressong[ !is.na(p67isoressong$padj), ]
p67isoressong <- p67isoressong[ !is.na(p67isoressong$pvalue), ]
#removepadj for replacement
p67isoressong <- p67isoressong[, -which(names(p67isoressong) == "padj")]

FDR.p67isoressong<- fdrtool(p67isoressong$stat, statistic= "normal", plot = T)

FDR.p67isoressong$param[1, "sd"]

p67isoressong[,"padj"]  <- p.adjust(FDR.p67isoressong$pval, method = "BH")

hist(FDR.p67isoressong$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p67isoressong[,"padj"] < 0.1)
DESeq2::plotMA(p67isoressong) #230
write.csv(as.data.frame(p67isoressong), file = "p67iso_songVSsilence.corrected.csv")



####################################
###p67silence Tutored vs Isloate ###
####################################

countDatap67sil <- read.table ("countsp67silence.txt", header=TRUE)
nrow(countDatap67sil)
head(countDatap67sil)

spDesignp67sil = data.frame (
  row.names = colnames(countDatap67sil),
  condition =  c("tutored", "tutored", "tutored", "tutored", "tutored", "tutored", "isolate", "isolate", "isolate", "isolate", "isolate", "isolate"))
spDesignp67sil

p67silddstut<- DESeqDataSetFromMatrix (countData = countDatap67sil,
                                        colData = spDesignp67sil,
                                        design = ~ condition)


p67silddstut <- DESeq(p67silddstut)
resultsNames(p67silddstut)

p67silrestut <- results(p67silddstut)
head(p67silrestut)

sum(p67silrestut$padj < 0.1, na.rm=TRUE) ## 1 genes
sum(p67silrestut$padj < 0.05, na.rm=TRUE) ## 0 genes
sum(p67silrestut$pvalue < 0.01, na.rm=TRUE) ## 95 genes

write.csv(as.data.frame(p67silrestut), file = "p67sil_TutvsIso_2022.csv")

hist(p67silrestut$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#looks bad, hill shape ()

#removeNA
p67silrestut <- p67silrestut[ !is.na(p67silrestut$padj), ]
p67silrestut <- p67silrestut[ !is.na(p67silrestut$pvalue), ]
#removepadj for replacement
p67silrestut <- p67silrestut[, -which(names(p67silrestut) == "padj")]

FDR.p67silrestut<- fdrtool(p67silrestut$stat, statistic= "normal", plot = T)

FDR.p67silrestut$param[1, "sd"]

p67silrestut[,"padj"]  <- p.adjust(FDR.p67silrestut$pval, method = "BH")

hist(FDR.p67silrestut$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p67silrestut[,"padj"] < 0.1)
DESeq2::plotMA(p67silrestut) #230
write.csv(as.data.frame(p67silrestut), file = "p67sil_TutVsIso.corrected_2022.csv")

#######################################
####p67song - tutored vs isolate ######
#######################################

countDatap67song <- read.table ("countsp67song.txt", header=TRUE)
nrow(countDatap67song)
head(countDatap67song)

spDesignp67song = data.frame (
  row.names = colnames(countDatap67song),
  condition =  c("tutored", "tutored", "tutored", "tutored", "tutored", "tutored", "isolate", "isolate", "isolate", "isolate", "isolate", "isolate"))
spDesignp67song

p67songddstut<- DESeqDataSetFromMatrix (countData = countDatap67song,
                                       colData = spDesignp67song,
                                       design = ~ condition)

p67songddstut <- DESeq(p67songddstut)
resultsNames(p67songddstut)

p67songrestut <- results(p67songddstut)
head(p67songrestut)

sum(p67songrestut$padj < 0.1, na.rm=TRUE) ## 3 genes
sum(p67songrestut$padj < 0.05, na.rm=TRUE) ## 3 genes
sum(p67songrestut$pvalue < 0.01, na.rm=TRUE) ## 127 genes

write.csv(as.data.frame(p67songrestut), file = "p67song_TutvsIso_2022.csv")

hist(p67songrestut$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#looks bad (not terrible), hill shape ()

#removeNA
p67songrestut <- p67songrestut[ !is.na(p67songrestut$padj), ]
p67songrestut <- p67songrestut[ !is.na(p67songrestut$pvalue), ]
#removepadj for replacement
p67songrestut <- p67songrestut[, -which(names(p67songrestut) == "padj")]

FDR.p67songrestut<- fdrtool(p67songrestut$stat, statistic= "normal", plot = T)

FDR.p67songrestut$param[1, "sd"]

p67songrestut[,"padj"]  <- p.adjust(FDR.p67songrestut$pval, method = "BH")

hist(FDR.p67songrestut$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p67songrestut[,"padj"] < 0.1)
DESeq2::plotMA(p67songrestut) #2359
write.csv(as.data.frame(p67songrestut), file = "p67song_TutVsIso.corrected_2022.csv")


###p32 song and silence ###
>>
  
  ####################################
###p32silence Tutored vs Isloate ###
####################################

countDatap32sil <- read.table ("countsp32silence.txt", header=TRUE)
nrow(countDatap32sil)
head(countDatap32sil)

spDesignp32sil = data.frame (
  row.names = colnames(countDatap32sil),
  condition =  c("tutored", "tutored", "tutored",  "isolate", "isolate", "isolate"))
spDesignp32sil

p32silddstut<- DESeqDataSetFromMatrix (countData = countDatap32sil,
                                       colData = spDesignp32sil,
                                       design = ~ condition)


p32silddstut <- DESeq(p32silddstut)
resultsNames(p32silddstut)

p32silrestut <- results(p32silddstut)
head(p32silrestut)

sum(p32silrestut$padj < 0.1, na.rm=TRUE) ## 2 genes
sum(p32silrestut$padj < 0.05, na.rm=TRUE) ## 1 genes
sum(p32silrestut$pvalue < 0.01, na.rm=TRUE) ## 65 genes

write.csv(as.data.frame(p32silrestut), file = "p32sil_TutvsIso.csv")

hist(p32silrestut$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32silrestut <- p32silrestut[ !is.na(p32silrestut$padj), ]
p32silrestut <- p32silrestut[ !is.na(p32silrestut$pvalue), ]
#removepadj for replacement
p32silrestut <- p32silrestut[, -which(names(p32silrestut) == "padj")]

FDR.p32silrestut<- fdrtool(p32silrestut$stat, statistic= "normal", plot = T)

FDR.p32silrestut$param[1, "sd"]

p32silrestut[,"padj"]  <- p.adjust(FDR.p32silrestut$pval, method = "BH")

hist(FDR.p32silrestut$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p32silrestut[,"padj"] < 0.1)
DESeq2::plotMA(p32silrestut) #188
write.csv(as.data.frame(p32silrestut), file = "p32sil_TutVsIso.corrected.csv")

#######################################
####p32song - tutored vs isolate ######
#######################################

countDatap32song <- read.table ("countsp32song.txt", header=TRUE)
nrow(countDatap32song)
head(countDatap32song)

spDesignp32song = data.frame (
  row.names = colnames(countDatap32song),
  condition =  c("tutored", "tutored", "tutored", "isolate", "isolate", "isolate"))
spDesignp32song

p32songddstut<- DESeqDataSetFromMatrix (countData = countDatap32song,
                                        colData = spDesignp32song,
                                        design = ~ condition)

p32songddstut <- DESeq(p32songddstut)
resultsNames(p32songddstut)

p32songrestut <- results(p32songddstut)
head(p32songrestut)

sum(p32songrestut$padj < 0.1, na.rm=TRUE) ## 2 genes
sum(p32songrestut$padj < 0.05, na.rm=TRUE) ## 1 genes
sum(p32songrestut$pvalue < 0.01, na.rm=TRUE) ## 78 genes

write.csv(as.data.frame(p32songrestut), file = "p32song_TutvsIso.csv")

hist(p32songrestut$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32songrestut <- p32songrestut[ !is.na(p32songrestut$padj), ]
p32songrestut <- p32songrestut[ !is.na(p32songrestut$pvalue), ]
#removepadj for replacement
p32songrestut <- p32songrestut[, -which(names(p32songrestut) == "padj")]

FDR.p32songrestut<- fdrtool(p32songrestut$stat, statistic= "normal", plot = T)

FDR.p32songrestut$param[1, "sd"]

p32songrestut[,"padj"]  <- p.adjust(FDR.p32songrestut$pval, method = "BH")

hist(FDR.p32songrestut$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p32songrestut[,"padj"] < 0.1)
DESeq2::plotMA(p32songrestut) #296
write.csv(as.data.frame(p32silrestut), file = "p32song_TutVsIso.corrected.csv")


#####################################
######################################
###p32 isolate song vs silence
######################################
######################################
countDatap32iso <- read.table ("countsp32isolate.txt", header=TRUE)
nrow(countDatap32iso)
head(countDatap32iso)

############# Pairwise DE ##########

spDesignp32iso = data.frame (
  row.names = colnames(countDatap32iso),
  song =  c("song", "song", "song", "silence", "silence", "silence"))
spDesignp32iso

p32isoddssong<- DESeqDataSetFromMatrix (countData = countDatap32iso,
                                        colData = spDesignp32iso,
                                        design = ~ song)


p32isoddssong <- DESeq(p32isoddssong)
resultsNames(p32isoddssong)
p32isoressong <- results(p32isoddssong)
head(p32isoressong)

sum(p32isoressong$padj < 0.1, na.rm=TRUE) ## 3 genes
sum(p32isoressong$padj < 0.05, na.rm=TRUE) ## 3 genes
sum(p32isoressong$pvalue < 0.01, na.rm=TRUE) ## 121 genes

write.csv(as.data.frame(p32isoressong), file = "p32IsosongVSsilence.csv")

hist(p32isoressong$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32isoressong <- p32isoressong[ !is.na(p32isoressong$padj), ]
p32isoressong <- p32isoressong[ !is.na(p32isoressong$pvalue), ]
#removepadj for replacement
p32isoressong <- p32isoressong[, -which(names(p32isoressong) == "padj")]

FDR.p32isoressong<- fdrtool(p32isoressong$stat, statistic= "normal", plot = T)

FDR.p32isoressong$param[1, "sd"]

p32isoressong[,"padj"]  <- p.adjust(FDR.p32isoressong$pval, method = "BH")

hist(FDR.p32isoressong$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p32isoressong[,"padj"] < 0.1)
DESeq2::plotMA(p32isoressong) #445
write.csv(as.data.frame(p32isoressong), file = "p32iso_songVSsilence.corrected.csv")

################################
###p67tutored Song vs Silent ###
################################

countDatap32tut <- read.table ("countsp32tutored.txt", header=TRUE)
nrow(countDatap32tut)
head(countDatap32tut)

############# Pairwise DE ##########

spDesignp32tut = data.frame (
  row.names = colnames(countDatap32tut),
  song =  c("song", "song", "song",  "silence", "silence", "silence"))
spDesignp32tut

p32tutddssong<- DESeqDataSetFromMatrix (countData = countDatap32tut,
                                        colData = spDesignp32tut,
                                        design = ~ song)


p32tutddssong <- DESeq(p32tutddssong)
resultsNames(p32tutddssong)
p32tutressong <- results(p32tutddssong)
head(p32tutressong)

sum(p32tutressong$padj < 0.1, na.rm=TRUE) ## 4 genes
sum(p32tutressong$padj < 0.05, na.rm=TRUE) ## 3 genes
sum(p32tutressong$pvalue < 0.01, na.rm=TRUE) ## 74 genes

write.csv(as.data.frame(p32tutressong), file = "p32TutSongVsSilence.csv")

hist(p32tutressong$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool

#removeNA
p32tutressong <- p32tutressong[ !is.na(p32tutressong$padj), ]
p32tutressong <- p32tutressong[ !is.na(p32tutressong$pvalue), ]
#removepadj for replacement
p32tutressong <- p32tutressong[, -which(names(p32tutressong) == "padj")]

FDR.p32tutressong<- fdrtool(p32tutressong$stat, statistic= "normal", plot = T)

FDR.p32tutressong$param[1, "sd"]

p32tutressong[,"padj"]  <- p.adjust(FDR.p32tutressong$pval, method = "BH")

hist(FDR.p32tutressong$pval, col = "royalblue4", 
     main = "song vs silence", xlab = "CORRECTED p-values")

table(p32tutressong[,"padj"] < 0.1) #305
DESeq2::plotMA(p32tutressong)
write.csv(as.data.frame(p32tutressong), file = "p32TuT_songVSsilence.corrected.csv")

############################
# Tutored Only #
####################

###tutored p32 vs p67 ###

countData_tut <- read.table ("counts_tutored.txt", header=TRUE)
nrow(countData_tut )
head(countData_tut )

############# Pairwise DE ##########

spDesign_tut = data.frame (
  row.names = colnames(countData_tut),
  age =  c("p32", "p32","p32","p32","p32","p32", "p67","p67","p67","p67","p67","p67","p67","p67","p67","p67","p67","p67" ))
spDesign_tut

TutddsAge<- DESeqDataSetFromMatrix (countData = countData_tut,
                                        colData = spDesign_tut,
                                        design = ~ age)


TutddsAge <- DESeq(TutddsAge)
resultsNames(TutddsAge)
tutresAge <- results(TutddsAge)
head(tutresAge)

sum(tutresAge$padj < 0.1, na.rm=TRUE) ## 412 genes
sum(tutresAge$padj < 0.05, na.rm=TRUE) ## 275 genes
sum(tutresAge$pvalue < 0.01, na.rm=TRUE) ## 708 genes

write.csv(as.data.frame(tutresAge), file = "Tutored_p32vsp67.csv")

hist(tutresAge$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)

#check pvalue histogram and if skewed, use FDRtool (this one is actually not bad)

#############################
# isolate only ##############
#############################
###tutored p32 vs p67 ###

countData_iso <- read.table ("counts_isolate.txt", header=TRUE)
nrow(countData_iso )
head(countData_iso )

############# Pairwise DE ##########

spDesign_iso = data.frame (
  row.names = colnames(countData_iso),
  age =  c("p32", "p32","p32","p32","p32","p32", "p67","p67","p67","p67","p67","p67","p67","p67","p67","p67","p67","p67" ))
spDesign_iso

IsoddsAge<- DESeqDataSetFromMatrix (countData = countData_iso,
                                    colData = spDesign_iso,
                                    design = ~ age)


IsoddsAge <- DESeq(IsoddsAge)
resultsNames(IsoddsAge)
IsoresAge <- results(IsoddsAge)
head(IsoresAge)

sum(IsoresAge$padj < 0.1, na.rm=TRUE) ## 267 genes
sum(IsoresAge$padj < 0.05, na.rm=TRUE) ## 138 genes
sum(IsoresAge$pvalue < 0.01, na.rm=TRUE) ## 570 genes

write.csv(as.data.frame(IsoresAge), file = "Isolate_p32vsp67.csv")

hist(tutresAge$pvalue, main="DESeq", xlab="p-values",cex.main=1.5)










>>>>>>>>>>>

  #tutor, controlling for song treatment
ddstutor<- DESeqDataSetFromMatrix (countData = countData,
                                       colData = spDesign,
                                       design = ~ song + tutor)
                                       #design = ~ song + tutor + song:tutor)

keep <- rowSums(counts(ddstutor)) >= 5
ddstutor5<- ddstutor[keep,]
ddstutor5 <- DESeq(ddstutor5)
resultsNames(ddstutor5)

#transformation and data vis
vsd <- vst(ddstutor5, blind = FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(ddstutor5)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( ddstutor5$dex, ddstutor10$cell, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

plotPCA(vsd, ntop= 100, intgroup = c("age"))


pcaExplorer(dds = ddstutor5, dst = vsd)

#ddstutor <- ddstutor[which(mcols(ddstutor)$betaConv),]
#restutor <- results(ddstutor, name="songsong.tutortutor")
restutor <- results(ddstutor, name="")
DESeq2::plotMA(restutor)
head(restutor)
sum(restutor$padj < 0.1, na.rm=TRUE) ## 2 genes
sum(restutor$padj < 0.05, na.rm=TRUE) ## 2 genes
sum(restutor$pvalue < 0.01, na.rm=TRUE) ## 128 genes

write.csv(as.data.frame(restutor), file = "p67_tutorVSisoINTERACTION.csv")

vst <- varianceStabilizingTransformation (ddstutor, blind=TRUE)
pcaExplorer(dds = ddstutor)

#song, controlling for tutor
ddssong<- DESeqDataSetFromMatrix (countData = countData,
                                   colData = spDesign,
                                   design = ~ tutor + song)


ddssong <- DESeq(ddssong)
resultsNames(ddssong)

#ddssong <- ddssong[which(mcols(ddssong)$betaConv),]
ressong <- results(ddssong, name="song_song_vs_silence")
head(ressong)
DESeq2::plotMA(ressong, alpha=0.5)
sum(ressong$padj < 0.1, na.rm=TRUE) ## 2 genes
sum(ressong$padj < 0.05, na.rm=TRUE) ## 2 genes
sum(ressong$pvalue < 0.01, na.rm=TRUE) ## 128 genes

write.csv(as.data.frame(ressong), file = "p67_songVSsilence.csv")

#histogram of p values
library(ggplot2)
stats <- read.csv ("p67_songVSsilence_trim.csv", header=TRUE,row.names=1)
head(stats)
ggplot(stats, aes(x=pvalue)) + geom_histogram()