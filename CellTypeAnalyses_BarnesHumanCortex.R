#Barnes data download and analysis (GSE21935 or GDS4522)
#Megan Hagenauer, 1/2017-7/2017
#Affymetrix Human Genome U133 Plus 2.0 Array

#installed GEOQuery using R

library(GEOquery)

#Scrap code - didn't work:
gse <- getGEO("GSE21935", GSEMatrix = TRUE)

# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21935/matrix/
#   Found 1 file(s)
# GSE21935_series_matrix.txt.gz
# trying URL 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE21nnn/GSE21935/matrix/GSE21935_series_matrix.txt.gz'
# ftp data connection made, file length 13667168 bytes
# ==================================================
#   downloaded 13.0 MB
# 
# Error in download.file(myurl, destfile, mode = mode, quiet = TRUE, method = getOption("download.file.method.GEOquery")) : 
#   cannot open URL 'http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL570&form=text&view=full'

gse<-getGEO(filename=system.file("extdata/GSE21935.soft.gz",package="GEOquery"))
# Error in fileOpen(fname) : File  does not appear to exist.
gse<-getGEO(filename=system.file("extdata/GSE21935.txt.gz",package="GEOquery"))
# Error in fileOpen(fname) : File  does not appear to exist.
gds<-getGEO(filename=system.file("extdata/GDS4522.soft.gz",package="GEOquery"))
# Error in fileOpen(fname) : File  does not appear to exist.
gse <- getGEO(filename=system.file("extdata/GSE21935_family.soft.gz",package="GEOquery"))
# Error in fileOpen(fname) : File  does not appear to exist.
gse<-getGEO("GSE21935_family", GSEMatrix = TRUE)
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE21935_FAMILY/GSE21935_FAMILY/matrix/
#   Error in function (type, msg, asError = TRUE)  : 
#   Server denied you to change to the given directory


#Hey - I think this worked!  So I guess I need to find the dataset accession number for each of my series.
gds<-getGEO("GDS4522", GSEMatrix = TRUE)
# File stored at: 
#   /var/folders/z3/lr898whn6wd1y76yv_8v3txc0000gn/T//RtmpwuL4rQ/GDS4522.soft.gz

head(Meta(gds))
Table(gds)[1:5,]
#That looks like RMA output.

Columns(gds)
#Well... that gets us the diagnosis, gender, age, and brain region, but we would still need to enter pH and PMI by hand.

Columns(gds)[,1:5]
str(Columns(gds))

# 'data.frame':	42 obs. of  5 variables:
#   $ sample       : Factor w/ 42 levels "GSM545725","GSM545726",..: 38 39 30 26 41 20 42 23 22 34 ...
# $ disease.state: Factor w/ 2 levels "control","schizophrenia": 2 2 2 2 2 2 2 2 2 2 ...
# $ gender       : Factor w/ 2 levels "female","male": 1 1 1 1 1 1 1 1 1 1 ...
# $ age          : Factor w/ 30 levels "25 years","28 years",..: 2 7 11 12 14 17 24 25 26 30 ...
# $ description  : chr  "Value for GSM545762: S028_Scz_F_28; src: Brain BA22 post-mortem schizophrenic" "Value for GSM545763: S029_Scz_F_53; src: Brain BA22 post-mortem schizophrenic" "Value for GSM545754: S013_Scz_F_65; src: Brain BA22 post-mortem schizophrenic" "Value for GSM545750: S007_Scz_F_67; src: Brain BA22 post-mortem schizophrenic" ...


#Hey - this worked too!
gse <- getGEO("GSE21935", GSEMatrix = FALSE)
# File stored at: 
#   /var/folders/z3/lr898whn6wd1y76yv_8v3txc0000gn/T//RtmpwuL4rQ/GSE21935.soft.gz
# Parsing....
# Found 43 entities...
# GPL570 (1 of 43 entities)
# GSM545725 (2 of 43 entities)
# GSM545726 (3 of 43 entities)
# GSM545727 (4 of 43 entities)
# GSM545728 (5 of 43 entities)
# GSM545729 (6 of 43 entities)
# GSM545730 (7 of 43 entities)
# GSM545731 (8 of 43 entities)
# GSM545732 (9 of 43 entities)
# GSM545733 (10 of 43 entities)
# GSM545734 (11 of 43 entities)
# GSM545735 (12 of 43 entities)
# GSM545736 (13 of 43 entities)
# GSM545737 (14 of 43 entities)
# GSM545738 (15 of 43 entities)
# GSM545739 (16 of 43 entities)
# GSM545740 (17 of 43 entities)
# GSM545741 (18 of 43 entities)
# GSM545742 (19 of 43 entities)
# GSM545743 (20 of 43 entities)
# GSM545744 (21 of 43 entities)
# GSM545745 (22 of 43 entities)
# GSM545746 (23 of 43 entities)
# GSM545747 (24 of 43 entities)
# GSM545748 (25 of 43 entities)
# GSM545749 (26 of 43 entities)
# GSM545750 (27 of 43 entities)
# GSM545751 (28 of 43 entities)
# GSM545752 (29 of 43 entities)
# GSM545753 (30 of 43 entities)
# GSM545754 (31 of 43 entities)
# GSM545755 (32 of 43 entities)
# GSM545756 (33 of 43 entities)
# GSM545757 (34 of 43 entities)
# GSM545758 (35 of 43 entities)
# GSM545759 (36 of 43 entities)
# GSM545760 (37 of 43 entities)
# GSM545761 (38 of 43 entities)
# GSM545762 (39 of 43 entities)
# GSM545763 (40 of 43 entities)
# GSM545764 (41 of 43 entities)
# GSM545765 (42 of 43 entities)
# GSM545766 (43 of 43 entities)
#There were 50 or more warnings (use warnings() to see the first 50)
head(Meta(gse))
str(Meta(gse))

# List of 30
# $ contact_address        : chr "Gunnels Wood Road"
# $ contact_city           : chr "Stevenage"
# $ contact_country        : chr "United Kingdom"
# $ contact_department     : chr "Computational Biology"
# $ contact_email          : chr "julie.x.huxley-jones@gsk.com"
# $ contact_institute      : chr "GlaxoSmithKline Pharmaceuticals"
# $ contact_name           : chr "Julie,,Huxley-Jones"
# $ contact_phone          : chr "01438 768416"
# $ contact_state          : chr "Hertfordshire"
# $ contact_zip/postal_code: chr "SG1 2NY"
# $ contributor            : chr [1:21] "Michael,R,Barnes" "Julie,,Huxley-Jones" "Peter,,Maycox" "Mark,,Lennon" ...
# $ email                  : chr "geo@ncbi.nlm.nih.gov"
# $ geo_accession          : chr "GSE21935"
# $ institute              : chr "NCBI NLM NIH"
# $ last_update_date       : chr "Dec 19 2016"
# $ name                   : chr "Gene Expression Omnibus (GEO)"
# $ overall_design         : chr "Post-mortem derived BA22 tissue from schizophrenic and control patients were compared. Age, gender, post-mortem delay and pH of"| __truncated__
# $ platform_id            : chr "GPL570"
# $ platform_taxid         : chr "9606"
# $ pubmed_id              : chr "21538462"
# $ relation               : chr "BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA126875"
# $ sample_id              : chr [1:42] "GSM545725" "GSM545726" "GSM545727" "GSM545728" ...
# $ sample_taxid           : chr "9606"
# $ status                 : chr "Public on Dec 31 2011"
# $ submission_date        : chr "May 20 2010"
# $ summary                : chr [1:2] "Transcriptional analysis of the superior temporal cortex (BA22) in schizophrenia: Pathway insight into disease pathology and dr"| __truncated__ "Schizophrenia is a highly debilitating psychiatric disorder which is known to have heritable genetic and environmental componen"| __truncated__
# $ supplementary_file     : chr "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/GSE21935/GSE21935_RAW.tar"
# $ title                  : chr "Comparison of post-mortem tissue from Brodman Brain BA22 region between schizophrenic and control patients"
# $ type                   : chr "Expression profiling by array"
# $ web_link               : chr "http://www.ncbi.nlm.nih.gov/geo"

GSMList(gse)[[1]]
# An object of class "GSM"
# channel_count 
# [1] "1"
# characteristics_ch1 
# [1] "gender: Female"                                           
# [2] "age: 87"                                                  
# [3] "post-mortem delay: 14.5h"                                 
# [4] "ph: 6.5"                                                  
# [5] "disease state: control"                                   
# [6] "tissue: superior temporal cortex (Brodmann Area 22, BA22)"
# contact_address 
# [1] "Gunnels Wood Road"
# contact_city 
# [1] "Stevenage"
# contact_country 
# [1] "United Kingdom"
# contact_department 
# [1] "Computational Biology"
# contact_email 
# [1] "julie.x.huxley-jones@gsk.com"
# contact_institute 
# [1] "GlaxoSmithKline Pharmaceuticals"
# contact_name 
# [1] "Julie,,Huxley-Jones"
# contact_phone 
# [1] "01438 768416"
# contact_state 
# [1] "Hertfordshire"
# contact_zip/postal_code 
# [1] "SG1 2NY"
# data_processing 
# [1] "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 150."
# data_row_count 
# [1] "54675"
# description 
# [1] "n/a"
# extract_protocol_ch1 
# [1] "Total RNA was extracted from frozen BA22 using a Polytron type homogenizer (YellowLine DI 25 Basic) and TriZol reagent (Invitrogen, Paisley, UK) in a ratio of 1 ml of TriZol to 20 mg of tissue. RNA was further purified using RNeasy mini-columns (Qiagen, Valencia, CA, USA) including on-column DNAse-1 step and elution in water. Although pH was analysed in brain lysates using a pH meter, this was not considered to be rigorous enough to exclude or include samples and instead the RNA integrity number (RIN) was used to assess the quality of the RNA as the primary inclusion criterion. The quantity of extracted RNA was determined by spectrophotometry and quality was assessed using an Agilent 2100 Bioanalyzer (South Plainfield, NJ, USA) to determine the RIN. Based on the RIN, samples were classified into three quality groups—pass (RIN >7.0); borderline (RIN 6.0–7.0); fail (RIN<6.0). Following classification, there were 41 pass (RIN range of samples 7.0–9.0: average=7.7), 16 borderline (RIN range of samples 6.0–6.9; average=6.4) and 5 fail samples. Samples in the fail category were excluded from the study, and the remaining samples were randomized into four batches, containing an equal number of schizophrenic/control and male/female samples, for target generation and hybridization."
# geo_accession 
# [1] "GSM545725"
# hyb_protocol 
# [1] "For each batch, 10ug total RNA was processed to Biotin-labeled cRNA and hybridized to HG-U133_Plus_2.0 GeneChips in accordance with the Affymetrix protocol (Affymetrix, Santa Clara, CA, USA,)."
# label_ch1 
# [1] "biotin"
# label_protocol_ch1 
# [1] "For each batch, 10ug total RNA was processed to Biotin-labeled cRNA and hybridized to HG-U133_Plus_2.0 GeneChips in accordance with the Affymetrix protocol (Affymetrix, Santa Clara, CA, USA,)."
# last_update_date 
# [1] "Dec 31 2011"
# molecule_ch1 
# [1] "total RNA"
# organism_ch1 
# [1] "Homo sapiens"
# platform_id 
# [1] "GPL570"
# scan_protocol 
# [1] "Arrays were scanned on a GeneChip Scanner 3000 and fluorescence intensity for each feature of the array was obtained by using GeneChip Operating Software (Affymetrix)."
# series_id 
# [1] "GSE21935"
# source_name_ch1 
# [1] "Brain BA22 post-mortem control"
# status 
# [1] "Public on Dec 31 2011"
# submission_date 
# [1] "May 20 2010"
# supplementary_file 
# [1] "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM545nnn/GSM545725/GSM545725_jrs_U133p2_DB2_54_C002_R97548.CEL.gz"
# taxid_ch1 
# [1] "9606"
# title 
# [1] "C002_Control_F_87"
# type 
# [1] "RNA"
# An object of class "GEODataTable"
# ****** Column Descriptions ******
#   Column             Description
# 1 ID_REF                        
# 2  VALUE MAS5.0 signal intensity
# ****** Data Table ******
#   ID_REF    VALUE
# 1 AFFX-BioB-5_at 312.4181
# 2 AFFX-BioB-M_at 411.2899
# 3 AFFX-BioB-3_at 251.2278
# 4 AFFX-BioC-5_at 640.8456
# 5 AFFX-BioC-3_at 630.6627
# 54670 more rows ...

str(GSMList(gse))
GSMList(gse)[1]$characteristics_ch1
str(GSMList(gse)[1])
GSMList(gse)[1]@header
GSMList(gse)$GSM545725
str(GSMList(gse)$GSM545725)
GSMList(gse)$GSM545725
#I can't get the characteristics out of the file. Bah...
#I'm being told the individual sample file info isn't subsettable???
str(gse)
names(GSMList(gse)$GSM545725)

Columns(GSMList(gse)$GSM545725)
#That really doesn't get us anywhere interesting

Meta(GSMList(gse)$GSM545725)

Meta(GSMList(gse)$GSM545725)$characteristics_ch1
#Win!

# [1] "gender: Female"                                           
# [2] "age: 87"                                                  
# [3] "post-mortem delay: 14.5h"                                 
# [4] "ph: 6.5"                                                  
# [5] "disease state: control"                                   
# [6] "tissue: superior temporal cortex (Brodmann Area 22, BA22)"

Meta(GSMList(gse)$GSM545725)$characteristics_ch1[1]
# [1] "gender: Female"

#I'm going to have to run that through a for loop and sub() that, and then double check data structure.

sub("gender: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[1])
as.numeric(sub("age: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[2]))
as.numeric(sub("h", "", sub("post-mortem delay: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[3])))
as.numeric(sub("ph: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[4]))
sub("disease state: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[5])
sub("tissue: ","", Meta(GSMList(gse)$GSM545725)$characteristics_ch1[6])

SampleID<-as.matrix(names(GSMList(gse)))
Gender<-matrix("a", nrow=42, ncol=1)
Age<-matrix(0, nrow=42, ncol=1)
BrainpH<-matrix(0, nrow=42, ncol=1)
PMI<-matrix(0, nrow=42, ncol=1)
Diagnosis<-matrix("a", nrow=42, ncol=1)
Tissue<-matrix("a", nrow=42, ncol=1)

length(GSMList(gse))
#[1] 42

for(i in c(1:42)){
  Gender[i]<-sub("gender: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[1])
  Age[i]<-as.numeric(sub("age: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[2]))
  PMI[i]<-as.numeric(sub("h", "", sub("post-mortem delay: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[3])))
  BrainpH[i]<-as.numeric(sub("ph: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[4]))
  Diagnosis[i]<-sub("disease state: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[5])
  Tissue[i]<-sub("tissue: ","", Meta(GSMList(gse)[[i]])$characteristics_ch1[6])
}

Barnes_SampleCharacteristics<-data.frame(SampleID, Gender, Age, BrainpH, PMI, Diagnosis, Tissue, stringsAsFactors=F)

getwd()
#[1] "/Users/mhh/Documents/Microarray Gen/Barnes_GSE21935"

write.csv(Barnes_SampleCharacteristics, "Barnes_SampleCharacteristics.csv")



#Now let's re-run RMA on their data:

library(org.Hs.eg.db)
library(plyr)
library(affy)

#This is where I obtained the updated custom .cdf for defining the probesets:
http://nmg-r.bioinformatics.nl/NuGO_R.html

# hgu133a2hsentrezg.db_19.0.2.tar.gz
# hgu133a2hsentrezgcdf_19.0.0.tar.gz
# hgu133a2hsentrezgprobe_19.0.0.tar.gz


install.packages(pkgs = c("hgu133a2hsentrezg.db", "hgu133a2hsentrezgcdf", "hgu133a2hsentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

#Changed working directory to where the .cel files are located

data2<-ReadAffy(cdfname ="hgu133a2hsentrezg")
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDF.txt")
RMAExpression_customCDF<-read.delim("data_customCDF.txt", sep="\t")
str(RMAExpression_customCDF)
#'data.frame':	12136 obs. of  43 variables:
write.csv(RMAExpression_customCDF, "RMAExpression_customCDF.csv")

##################################
library(AffyRNADegradation)

tongs <- GetTongs(data2, chip.idx = 4)
PlotTongs(tongs)
tongs <- GetTongs(data2, chip.idx = 5)
PlotTongs(tongs)

rna.deg<- RNADegradation(data2, location.type = "index")
RNADegradPerSample<-d(rna.deg)
##############################

ScanDate<-protocolData(data2)$ScanDate
#Yep, there are definitely different scan dates here (2), although they are pretty close.
library(reshape2)
ScanDate_Split<-colsplit(ScanDate, pattern=" ", c("ScanDate", "ScanTime"))
table(ScanDate_Split$ScanDate)


#Alright, now I need the additional annotation:

head(RMAExpression_customCDF)
RMAExpression_EntrezID<-sub("_at", "", RMAExpression_customCDF[,1])
head(RMAExpression_EntrezID)
RMAExpression_customCDFAnnotation<-data.frame(RMAExpression_customCDF[,1], RMAExpression_EntrezID, stringsAsFactors = F )
colnames(RMAExpression_customCDFAnnotation)<-c("ProbesetID", "EntrezGeneID")

x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))
table(lengths(xx))
# 1 
# 59887 

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

RMAExpression_customCDFAnnotation2<-join(RMAExpression_customCDFAnnotation, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")


sum(is.na(RMAExpression_customCDFAnnotation2[,3])==F)
#[1] 12064
dim(RMAExpression_customCDFAnnotation2)
#[1] 12136     3
#So almost all of the results have gene symbols.

write.csv(RMAExpression_customCDFAnnotation2, "RMAExpression_customCDFAnnotation2.csv")

#Cleaning up the workspace:
rm(GeneSymbol)
rm(EntrezGeneID)
rm(gds)
rm(gse)
rm(data2)
rm(x)
rm(xx)

#Note: I double checked, and the sample characteristics are in the same order as the new RMA file.
Barnes_SampleCharacteristics<-data.frame(Barnes_SampleCharacteristics, ScanDate_Split, RNADegradPerSample)
write.csv(Barnes_SampleCharacteristics, "Barnes_SampleCharacteristics.csv")

RMAExpression_customCDFMatrix<-as.matrix(RMAExpression_customCDF[,-1])
is.numeric(RMAExpression_customCDFMatrix)

colnames(RMAExpression_customCDFMatrix)<-SampleID

png("Boxplot_RMAExpression.png", width=2000, height=400)
boxplot(RMAExpression_customCDFMatrix)
dev.off()
#Wow - those are some pretty homogenous looking samples...RMA normalization worked, I guess.

png("CorMatrix_RMAExpression.png")
heatmap(cor(RMAExpression_customCDFMatrix))
dev.off()
#No super blatant outliers

write.csv(cor(RMAExpression_customCDFMatrix), "CorMatrix_RMAExpression.csv")
#It looks like there are a couple of samples that could be considered outliers, but they're not terrible. Let's double check:

#So that I can easily re-use code:
SignalSortedNoNA3<-RMAExpression_customCDFMatrix

png("10 Boxplot Sample Sample Correlations.png")
boxplot(data.frame(cor(SignalSortedNoNA3)), cex=0.25, las=3, par(cex.axis=0.75))
Median10thQuantile<-median(apply((cor(SignalSortedNoNA3)), 1, quantile, 0.1))
abline(a=Median10thQuantile, b=0, col=2)
dev.off()
#2 of the samples are borderline outliers. I'm going to leave them in.

#############################################

Gender<-relevel(as.factor(Barnes_SampleCharacteristics$Gender), ref="Male")
Age<-Barnes_SampleCharacteristics$Age
PMI<-Barnes_SampleCharacteristics$PMI
BrainpH<-Barnes_SampleCharacteristics$BrainpH
Diagnosis<-as.factor(Barnes_SampleCharacteristics$Diagnosis)
ScanDate<-Barnes_SampleCharacteristics$ScanDate
ScanDate<-as.factor(ScanDate)
RNADegradPerSample<-RNADegradPerSample

#Characterizing subjects and looking for relationships between subject variables:
#Changed Working directory
setwd("~/Documents/Microarray Gen/Barnes_GSE21935/SubjVarVsSubjVar")

SubjectFactorVariables<-cbind(Diagnosis, Gender, ScanDate)
colnames(SubjectFactorVariables)<-c("Diagnosis", "Gender", "ScanDate")

SubjectContinuousVariables<-cbind(BrainpH, PMI, Age, RNADegradPerSample)
colnames(SubjectContinuousVariables)<-c("BrainpH", "PMI", "Age", "RNADegrad")

for (i in 1:length(SubjectContinuousVariables[1,])){
  png(paste(paste("Histogram of", colnames(SubjectContinuousVariables)[i], sep="  "), "png", sep="."))	
  hist(SubjectContinuousVariables[, i], col=i+1)
  dev.off()		
}
#Wide age distribution, but many of the subjects are elderly (~70-80 years old)
#wide pH distribution (5.8-7), but typically low (~6.2)
#PMI tends to be less than 20 hrs, one outlier at 30 hrs, so analyses for PMI are going to be unreliable.


#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    RegressionLine<-lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}


#Using boxplots to visually examine the relationships between the continuous subject variables and categorical subject variables:
for (i in 1:length(SubjectContinuousVariables[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("14", paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectContinuousVariables)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}

#PH varies dramatically by diagnosis!  Ouch!

#Creating a text file of contingency tables to visually examine the relationships between categorical subject variables:
CrossTabsIV<-file("14 Cross Tabs Between Subject Factors.txt")
out<-c(
  capture.output(
    
    summary(Diagnosis),
    summary(Gender),
    
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i],SubjectFactorVariables[,j])
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(ContingencyTable)
        print(paste("p-value=", chisq.test(ContingencyTable)$p.value))	
      }		
    }
  )
)
cat(out, file="14 Cross Tabs Between Subject Factors.txt", sep="\n", append=TRUE)
close(CrossTabsIV)
rm(out)


library(car)

StatisticalRelationshipsIV<-file("14 Statistical Relationships between Subject Variables.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. 
    vif(lm(SignalSortedNoNA3[1,]~BrainpH+PMI+Diagnosis+Gender+Age))
    
  ),
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(chisq.test(ContingencyTable))		
      }		
    }
  )
  
)
cat(out, file="14 Statistical Relationships between Subject Variables.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIV)
rm(out)


#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIV<-file("14 Flagged Relationships Between Subject Variables.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectContinuousVariables[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectContinuousVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectContinuousVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectContinuousVariables[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  ),
  
  #Using chi-square to examine the statistical relationships between the categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectFactorVariables[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        ContingencyTable<-table(SubjectFactorVariables[,i], SubjectFactorVariables[,j])
        if(chisq.test(ContingencyTable)$p.value<0.05){
          print(paste(colnames(SubjectFactorVariables)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", chisq.test(ContingencyTable)$p.value, sep="  "))	
        }else{}		
      }		
    }
  )
)
cat(out, file="14 Flagged Relationships Between Subject Variables.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIV)
rm(out)


################################
# #Run principal components analysis (PCA) to determine which major gradients of sample-sample correlations exist in the data (i.e. who is similar to whom):
setwd("~/Documents/Microarray Gen/Barnes_GSE21935/PCA")


pcaNormFilterednoOutliers<-prcomp(t(SignalSortedNoNA3))
tmp<-pcaNormFilterednoOutliers$x[,1:4]
write.table(tmp, "PCA_1_4.txt", sep="\t")


PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, RMAExpression_customCDFAnnotation2plus2)
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]

PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]

# #Output a scree plot for the PCA (no outliers):
png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()



# #Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()

# #Output a scatterplot illustrating the relationship between Principal components 3 & 4 (PC3 & PC4):
png("10 PC3 vs PC4.png")
plot(PC3noOutliers~PC4noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off()


SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

PCAoutput<-cbind(SubjectFactorVariables, SubjectContinuousVariables, SubjectPCA)
write.csv(PCAoutput, "PCAoutput.csv")

#Using a scatterplot with best fit line to visually examine the relationships between the continuous subject variables and SubjectPCA:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectContinuousVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), "png", sep="."))	
    plot(SubjectPCA[,i]~SubjectContinuousVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "), xlab=colnames(SubjectContinuousVariables)[j], ylab=colnames(SubjectPCA)[i])
    RegressionLine<-lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])
    abline(RegressionLine, col=2)
    mtext(paste("p-value=", round(summary.lm(RegressionLine)$coefficients[8], digits=4)))
    dev.off()		
  }		
}



#Using boxplots to visually examine the relationships between the PCA and categorical subject variables:
for (i in 1:length(SubjectPCA[1,])){
  for(j in 1:length(SubjectFactorVariables[1,])){
    png(paste("15", paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), "png", sep="."))	
    boxplot(SubjectPCA[,i]~SubjectFactorVariables[,j], main=paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "), xlab=colnames(SubjectFactorVariables)[j], ylab=colnames(SubjectPCA)[i])
    mtext(paste("p-value=", round(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], digits=4)))
    dev.off()		
  }		
}




#Outputting a text file containing the statistical relationships between all of the subject variables and PCA:

StatisticalRelationshipsIVvsPCA<-file("15 Statistical Relationships between Subject Variables and PCA.txt")
out<-c(
  
  capture.output(
    #Calculating the variance inflation factor (vif) to determine which subject variables are highly related to other subject variables in the data set. Most important, of course, is whether any of the subject variables strongly correlate with Diagnosis. Note that "Location on Chip" has been removed as a variable because it is partially redundant with gender. 
    summary.lm(lm(PC1noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(PC2noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(PC3noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSample))
  ),
  
  capture.output(
    summary.lm(lm(PC4noOutliers~BrainpH+PMI+Diagnosis+Gender+Age+ScanDate+RNADegradPerSample))
  ),
  
  
  #Using linear regression to examine the statistical relationships between PCA and the continuous subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], sep="  "))
        print(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j])))
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between PCA and categorical subject variables:
  
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], sep="  "))
        print(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j])))		
      }		
    }
  )
  
)
cat(out, file="15 Statistical Relationships between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(StatisticalRelationshipsIVvsPCA)
rm(out)



#Flagging variables that are collinear with other subject variables:
FlaggedRelationshipsBetweenIVandPCA<-file("15 Flagged Relationships Between Subject Variables and PCA.txt")
out<-c(
  
  #Using linear regression to examine the statistical relationships between the continuous subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectContinuousVariables[1,])){
        if(summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectContinuousVariables)[j], "p-value=", summary.lm(lm(SubjectPCA[,i]~SubjectContinuousVariables[,j]))$coefficient[8], sep="  "))}else{}
      }		
    }
  ),
  
  #Using anova to examine the statistical relationships between the continuous subject variables and categorical subject variables:
  capture.output(
    for (i in 1:length(SubjectPCA[1,])){
      for(j in 1:length(SubjectFactorVariables[1,])){
        if(summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1]<0.05){
          print(paste(colnames(SubjectPCA)[i], "vs", colnames(SubjectFactorVariables)[j], "p-value=", summary(aov(SubjectPCA[,i]~SubjectFactorVariables[,j]))[[1]][["Pr(>F)"]][1], sep="  "))	
        }else{}		
      }		
    }
  )
)
cat(out, file="15 Flagged Relationships Between Subject Variables and PCA.txt", sep="\n", append=TRUE)
close(FlaggedRelationshipsBetweenIVandPCA)
rm(out)

#ScanDate really doesn't seem to matter in this dataset.
#In contrast, RNADegradation really does...

#####################

#Running Cell Type analyses:
setwd("~/Documents/Microarray Gen/Barnes_GSE21935/BrainInABlender")

Install.packages("devtools")

library("devtools")

install_github("hagenaue/BrainInABlender")

library("BrainInABlender")
temp<-data.frame(as.character(RMAExpression_customCDFAnnotation2plus2[,3]), SignalSortedNoNA3, stringsAsFactors=F)
write.csv(temp, "UserInput.csv")
str(temp)
table(table(RMAExpression_customCDFAnnotation2[,3]))
Sir_UnMixALot(userInput=temp, dataColumns=c(2:43), geneColumn=1, species="human")
#hit an error message. :(
#I'm going to try reinstalling and trying again.
#Nope, same error message. We've got a bug. :(

#To troubleshoot, I'm going to see if the problem continues to exist if I get rid of the NA gene symbols

UserInput<-read.table("UserInput.csv", header=T, stringsAsFactors = F, sep=",")
Sir_UnMixALot(userInput=UserInput, dataColumns=c(2:43), geneColumn=1, species="human")
#Nope, crashed with that too.

length(unique(UserInput[,1]))
[1] 12125

UserInput[duplicated(UserInput[,1]),1]
#Just a couple of NAs and two March genes

#Running it again just to doublecheck the output:
# [1] "I like big brains and I cannot lie..."
# [1] "Good: The userInput is a data frame"
# [1] "Good: The dataColumns are integer values"
# [1] "Warning: Your geneColumn is not an integer value"
# [1] "Good: You have specified only one geneColumn"
# [1] "Good: Your species selection is one of our included species."
# [1] "The dimensions of the matrix for the user's gene expression data:"
# [1] 12136    42
# [1] "The number of rows in the dataset that show no variability (sd=0):"
# [1] 0
# [1] "The percentage of rows in the dataset that show no variability (sd=0):"
# [1] 0
# [1] "These rows will be removed."
# [1] "Output added to working directory: ZscoreInput.csv"
# [1] "...ZscoreInput.csv is a data frame that includes a a z-scored version of your input (centered and scaled by row), with all rows removed that had zero variability."
# [1] "The number of unique gene symbols (rows) included in the user's input after filtering out genes with expression that completely lacked variability (sd=0):"
# [1] 12125
# [1] "The gene symbols are not unique, the z-scored data will be averaged by gene symbol and then re-z-scored"
# Error in temp[, i] <- tapply(ZscoreInput[, i], GeneNamesForJoinedInput_NoSD0[,  : 
#                                                                                number of items to replace is not a multiple of replacement length

#Alright, let's run it line by line and see where it is getting stuck:

userInput<-UserInput
dataColumns<-c(2:43)
geneColumn<-1
species<-"human"

print("I like big brains and I cannot lie...")

if(is.data.frame(userInput)==F){print("Error: Your userInput is not a data.frame")}else{print("Good: The userInput is a data frame")}

if(is.integer(dataColumns)==F){print("Warning: Your dataColumns are not integer values")}else{print("Good: The dataColumns are integer values")}

if(is.integer(geneColumn)==F){print("Warning: Your geneColumn is not an integer value")}else{print("Good: The geneColumn an integer value")}
#MH: perhaps I should just check to make sure it is numeric?

if(length(geneColumn)>1){print("Error: You have specified more than one geneColumn")}else{print("Good: You have specified only one geneColumn")}

if(species%in%c("Mouse", "mouse", "Human", "human")){print("Good: Your species selection is one of our included species.")}else{print("Error: You have specified a different species than the two included in our database (mouse and human).")}


#MH: I generalized this to specify the gene column number:
GeneNamesForJoinedInput <- userInput[,geneColumn]
GeneNamesForJoinedInput <- as.matrix(GeneNamesForJoinedInput)
#These look ok

#MH: I generalized this to specify the data column numbers:
TempJoinedInput_AsNum <- userInput[,dataColumns]

#MH: I added a little feedback to the user here:
print("The dimensions of the matrix for the user's gene expression data:")
print(dim(TempJoinedInput_AsNum))



#MH: This code determines how many of the rows (probes/genes) have zero variability and therefore have to be removed from the   calculations - I've made it so that the user gets some feedback about what is happening:
JoinedInput_StDev<-apply(TempJoinedInput_AsNum, 1, function(y) sd(y, na.rm=T)) 

print("The number of rows in the dataset that show no variability (sd=0):")
#MH: this originally said "return" but I'm not sure why. I changed it to print
print(sum(JoinedInput_StDev==0))

print("The percentage of rows in the dataset that show no variability (sd=0):")
print((sum(JoinedInput_StDev==0)/length(JoinedInput_StDev))*100)
print("These rows will be removed.")

JoinedInput_AsNumMatrix_Log2_NoSD0<-TempJoinedInput_AsNum[JoinedInput_StDev>0,]
temp<-GeneNamesForJoinedInput
GeneNamesForJoinedInput_NoSD0<-temp[JoinedInput_StDev>0]

GeneNamesForJoinedInput_NoSD0 <-as.matrix(GeneNamesForJoinedInput_NoSD0)

ZscoreInput<-t(scale(t(JoinedInput_AsNumMatrix_Log2_NoSD0), center=T, scale=T))#Zscores the data 
write.csv(ZscoreInput, "ZscoreInput.csv")
print("Output added to working directory: ZscoreInput.csv")
print("...ZscoreInput.csv is a data frame that includes a a z-scored version of your input (centered and scaled by row), with all rows removed that had zero variability.")

sum(is.na(ZscoreInput))
# ZERO WOOOOOOOOOO

#All good so far - this is where the bug is likely to be located

#MH: I added a little feedback to the user here:
print("The number of unique gene symbols (rows) included in the user's input after filtering out genes with expression that completely lacked variability (sd=0):")
print(length(unique(GeneNamesForJoinedInput_NoSD0)))

if(length(unique(GeneNamesForJoinedInput_NoSD0))==length((GeneNamesForJoinedInput_NoSD0[,1]))){print("All gene symbols are now unique.")}else{
  print("The gene symbols are not unique, the z-scored data will be averaged by gene symbol and then re-z-scored")
  
  temp<-matrix(0, length(unique(GeneNamesForJoinedInput_NoSD0)), ncol(ZscoreInput))
  
  #This is where the bug is - the number of rows in the output is 12124 instead of 12125. Huh.
  for(i in c(1:ncol(ZscoreInput))){
    temp[,i]<-tapply(ZscoreInput[,i], GeneNamesForJoinedInput_NoSD0[,1], mean)
  }
  
  length(unique(GeneNamesForJoinedInput_NoSD0[,1]))
  length(names(table(GeneNamesForJoinedInput_NoSD0[,1])))
  #Weird - why would these values differ?
  #Ah - NA values
  length(names(table(GeneNamesForJoinedInput_NoSD0[,1], useNA="always")))
  #[1] 12125
  
#Fixed it - tapply ignores NAs
#Hurray! It's working now!
  
 
  
#But the output looks like crap - none of the publication-specific indices really correlate with each other. Bug?
  #Huh - some of the samples just have generally really high levels or really low levels of everything.
  #...despite quantile normalization. Weird, right?  Perhaps these were lower RIN samples?
  #...and correcting for the overall low signal values for the cell type specific genes doesn't seem to help.
  
#Also - some of Doyle's indices have almost no representation after overlap removal (Corticostriatal and CCK). Very unstable. We really need to add some code to remove these.
  
  #In the meantime, let's just ignore the Projection Neuron and Interneuron results and take a peek at the others.
  
  
  UserInput<-data.frame(as.character(RMAExpression_customCDFAnnotation2plus2[,3]), SignalSortedNoNA3, stringsAsFactors=F)
  
SirUnMixALotOutput<-Sir_UnMixALot(userInput=UserInput, dataColumns=c(2:43), geneColumn=1, species="human")
  str(SirUnMixALotOutput)
  
  PublicationSpecific_CellTypeIndex<-SirUnMixALotOutput$PublicationSpecific_CellTypeIndex
  AveragePrimary_CellTypeIndex<-SirUnMixALotOutput$AveragePrimary_CellTypeIndex
  
  str(PublicationSpecific_CellTypeIndex)

  temp<-cbind(as.matrix(SubjectContinuousVariables), t(PublicationSpecific_CellTypeIndex), t(AveragePrimary_CellTypeIndex)) 
  str(temp)
  SubjectContinuousVariables<-temp
  
  #Then I re-ran the subjvar and PCA code.
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/CellTypeVsPCA")
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/CellTypeVsSubjVar")
  
 ################################################ 
  
  #These cell type indices really don't look like coherent concepts to me. Is it a problem with BrainInABlender or a problem with the dataset?
  
  #Could the really low expression genes be causing a problem? Although that should really only just add noise.
  #Ah -except the expression levels as a whole are *super* low for this dataset. Huh. Maybe it is just a bad set of arrays. (look at he boxplot - the 25 quantile is around 5 and the 75% quantile is around 6.5)
  
  StDevPerGene<-apply(SignalSortedNoNA3, 1, sd)
  MeanPerGene<-apply(SignalSortedNoNA3, 1, mean)
  
  png("StDevPerGene_vs_MeanSignalPerGene.png")
  plot(StDevPerGene~MeanPerGene)
  dev.off()
  
  #Interesting- not much evidence for heteroskedasticity at all. That seems a little suspect, since the lower level genes shouldn't be expressed at all...
  
  #I'm going to try reading in some of our affy data and comparing:
  
  dlpfc_5708_rma<-read.delim("dlpfc_5708_rma.txt", header=T, sep="\t")
  is.numeric( dlpfc_5708_rma)
  boxplot(dlpfc_5708_rma)
  boxplot(dlpfc_5708_rma[,c(2:30)])
  
  rm(dlpfc_5708_rma)
  
  #Yes, our data has a median near 6 too but quite a bit more spread (75% quantile is up near 8 instead of 6.5, and the whiskers reach up to 10.5)
  
  StDevPerGene<-apply(dlpfc_5708_rma[,-1], 1, sd)
  MeanPerGene<-apply(dlpfc_5708_rma[,-1], 1, mean)
  
  png("StDevPerGene_vs_MeanSignalPerGene.png")
  plot(StDevPerGene~MeanPerGene)
  dev.off()
  
  #Our data has clear heteroskedasticity. Probesets with signal below 4.5 or 5 are particularly not variable.
  
  #I wonder what the signal looked like for the original Barnes chips.  Perhaps some of the chips had very little signal, and that is why the quantile normalized data has such low expression?
  #In that case, the data may be super variable.
  #I wonder how I can view the data without quantile normalizing it. Ah ha - the rma command can be changed.
  
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/GSE21935_RAW")
  
  
  data2<-ReadAffy(cdfname ="hgu133a2hsentrezg")
  eset2 <- rma(data2, normalize=F)
  write.exprs(eset2,file="data_customCDF_NoQuantileNorm.txt")
  RMAExpression_customCDF_NoQuantileNorm<-read.delim("data_customCDF_NoQuantileNorm.txt", sep="\t")
  str(RMAExpression_customCDF_NoQuantileNorm)
  'data.frame':	12136 obs. of  43 variables:
    write.csv(RMAExpression_customCDF_NoQuantileNorm, "RMAExpression_customCDF_NoQuantileNorm.csv")
  
  
  png("Boxplot_RMAExpression_NoQuantileNorm.png", width=2000, height=400)
  boxplot(as.matrix(RMAExpression_customCDF_NoQuantileNorm[,-1]))
  dev.off()
  #It's not super variable. There are a couple of samples where 50% of the data falls below 5.5, but not strikingly so. And those aren't the same samples that looked weird in the cell type analyses.
  
  
  #Nevertheless, let's try re-running the analysis with only genes with expression>5 (signal of 32)
  
  
  Percentile75<-apply(SignalSortedNoNA3, 1, function(x) quantile(x, probs=c(0.75)))
  hist(Percentile75)
  sum(Percentile75>5)
  #[1] 10459
  sum(Percentile75>5.5)
  #[1] 8863
  SignalSortedNoNA3_Filtered<-SignalSortedNoNA3[Percentile75>5.5,]
  GeneSymbols_Filtered<-as.character(RMAExpression_customCDFAnnotation2[,3])[Percentile75>5.5]
  
  temp<-data.frame(GeneSymbols_Filtered, SignalSortedNoNA3_Filtered, stringsAsFactors=F)
 str(temp)
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/BrainInABlender_Filtered")
  
  SirUnMixALotOutput<-Sir_UnMixALot(userInput=temp, dataColumns=c(2:43), geneColumn=1, species="human")
  #Looks basically the same. 
  
  
  #I'm going to double check that the sample are who they are supposed to be:
  RMAExpression_customCDFAnnotation2[RMAExpression_customCDFAnnotation2[,3]=="XIST",]
  
  SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="XIST",]
  #Interesting - the NA values are popping up with this. Maybe that is part of the problem.
  #Nope -shouldn't be after averaging by gene symbol.
  boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="XIST",][11,]~Gender)
  #Wow-definitely not full separation.
  hist(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="XIST",][11,], breaks=40)
  #And not bimodal :(
  
  boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="RPS4Y1",][10,]~Gender)
 #evenly distributed between males and females.
  hist(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="RPS4Y1",][10,], breaks=40)
  #And not bimodal...very low expression values (3.7)
  
  boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="DDX3Y",][11,]~Gender)
  #evenly distributed, signal is down around 5.2
  
  boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="USP9Y",][11,]~Gender)
  #evenly distributed, but good signal values - 7.2
  hist(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2[,3]=="USP9Y",][11,], breaks=40)
  #Not particularly bimodal
  
  
  tail(SignalSortedNoNA3)
 sum(is.na(RMAExpression_customCDFAnnotation2[,3]))
 [1] 72
str(RMAExpression_customCDFMatrix)

#I wonder what their RMA data looks like - perhaps this is an annotation problem?


gds<-getGEO("GDS4522", GSEMatrix = TRUE)
# File stored at: 
#   /var/folders/z3/lr898whn6wd1y76yv_8v3txc0000gn/T//RtmpwuL4rQ/GDS4522.soft.gz

head(Meta(gds))
Table(gds)[1:5,]
#That looks like RMA output.

str(Table(gds))
# 'data.frame':	54675 obs. of  44 variables:
#   $ ID_REF    : Factor w/ 54676 levels "1007_s_at","1053_at",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ IDENTIFIER: Factor w/ 30807 levels "ADAM32","AFG3L1P",..: 43 58 36 48 34 44 67 53 14 23 ...
# $ GSM545762 : num  510.8 30.7 40.6 183.5 32 ...
# $ GSM545763 : num  477.4 41.2 24.7 161.5 33 ...
# $ GSM545754 : num  411.8 33.6 23.9 177 21.4 ...

write.csv(Table(gds), "OriginalRMAfromGEO.csv")

OriginalMas5<-Table(gds)
colnames(OriginalMas5)


OriginalMas5_log2<-log2(Table(gds)[,-c(1:2)])

boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[1,]~Gender)
#well that's weird... there is variability, but it is showing higher expression in Males...
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[2,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[3,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[4,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[5,]~Gender)

boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="RPS4Y1",])[1,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="DDX3Y",])[1,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="DDX3Y",])[2,]~Gender)
boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="DDX3Y",])[3,]~Gender)
#These all look reversed from what they are supposed to be. 


cbind(colnames(OriginalMas5_log2), SampleID, Barnes_SampleCharacteristics)
#Ah. Well, that would explain that. Maybe. 

Temp<-as.data.frame(colnames(OriginalMas5_log2))
colnames(Temp)<-"SampleID"
SampleCharacteristics_OriginalMas5<-join(Temp, Barnes_SampleCharacteristics)

boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[1,]~SampleCharacteristics_OriginalMas5$Gender)
#Working! Whoot!

hist(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="XIST",])[1,])
#Bimodal, as it should be.

boxplot(as.matrix(OriginalMas5_log2[OriginalMas5[,2]=="RPS4Y1",])[1,]~SampleCharacteristics_OriginalMas5$Gender)
#Yep - all good now. 


#So... why is the custom .cdf RMA data all messed up looking?
#Let's double check sample order first, even though that doesn't explain the lack of bimodality.
#Double checked. 

#next, let's double-check the alignment between gene symbol and entrez gene.
OriginalMas5[OriginalMas5[,2]=="XIST",1]
#[1] 214218_s_at 221728_x_at 224588_at   224589_at   224590_at   227671_at   235446_at
#[8] 243712_at 

OriginalMas5[OriginalMas5[,2]=="RPS4Y1",1]
#201909_at

#cdf and the chip.

# hgu133plus2hsentrezg.db_19.0.2.tar.gz
# hgu133plus2hsentrezgcdf_19.0.0.tar.gz
# hgu133plus2hsentrezgprobe_19.0.0.tar.gz

install.packages(pkgs = c("hgu133plus2hsentrezg.db", "hgu133plus2hsentrezgcdf", "hgu133plus2hsentrezgprobe"), repos = "http://nmg-r.bioinformatics.nl/bioc")

#Changed working directory to where the .cel files are located

data2<-ReadAffy(cdfname ="hgu133plus2hsentrezg")
eset2 <- rma(data2)
write.exprs(eset2,file="data_customCDFplus2.txt")
RMAExpression_customCDFplus2<-read.delim("data_customCDFplus2.txt", sep="\t")
str(RMAExpression_customCDF)
#'data.frame':	12136 obs. of  43 variables:
  write.csv(RMAExpression_customCDFplus2, "RMAExpression_customCDFplus2.csv")

  
  head(RMAExpression_customCDFplus2)
  RMAExpression_EntrezIDplus2<-sub("_at", "", RMAExpression_customCDFplus2[,1])
  head(RMAExpression_EntrezIDplus2)
  RMAExpression_customCDFAnnotationplus2<-data.frame(RMAExpression_customCDFplus2[,1], RMAExpression_EntrezIDplus2, stringsAsFactors = F )
  colnames(RMAExpression_customCDFAnnotationplus2)<-c("ProbesetID", "EntrezGeneID")
  
  x <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(x)
  xx <- as.list(x[mapped_genes])
  
  GeneSymbol<-unlist(xx, use.names=FALSE)
  EntrezGeneID<-rep(names(xx), lengths(xx))
  table(lengths(xx))
  # 1 
  # 59887 
  
  EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)
  
  RMAExpression_customCDFAnnotation2plus2<-join(RMAExpression_customCDFAnnotationplus2, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")
  
  
  sum(is.na(RMAExpression_customCDFAnnotation2[,3])==F)
  #[1] 12064
  dim(RMAExpression_customCDFAnnotation2)
  #[1] 12136     3
  #So almost all of the results have gene symbols.
  
  write.csv(RMAExpression_customCDFAnnotation2plus2, "RMAExpression_customCDFAnnotation2plus2.csv")
  
  RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",]

  RMAExpression_customCDFAnnotation2plus2[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]
  
  str(RMAExpression_customCDFplus2)
  SignalSortedNoNA3<-as.matrix(RMAExpression_customCDFplus2[,-1])
  
  png("Boxplot_RMAExpression_customCDFplus2.png", width=2000, height=400)
  boxplot(SignalSortedNoNA3)
  dev.off()
  #Looks quantile normalized and the signal values now have an appropriate range. Hurray!
  
  png("XIST_vs_Gender_customCDFplus2.png")
  boxplot(SignalSortedNoNA3[RMAExpression_customCDFAnnotation2plus2[,3]=="XIST",][173,]~Gender, col=2)
  dev.off()
  #And we now safely pass gender check. *phew*
  
  #Alright, time to re-do all of the previous analysis.
  
  
  #Hurray!  Our findings replicated!  The power is much weaker though.
  
  temp<-cbind(as.matrix(SubjectContinuousVariables), t(PublicationSpecific_CellTypeIndex), t(AveragePrimary_CellTypeIndex))
  SubjectContinuousVariables<-temp
  
  Gender<-relevel(as.factor(Gender), ref="Male")
  

  
  
  ###############################
  
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/LM4_basic")
  
  #I need to output results from a basic linear model for comparison:
  
  BrainpHcentered<-BrainpH-median(BrainpH)
  Agecentered<-Age-median(Age)
  PMIcentered<-PMI-median(PMI)
  
  GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 7)
  GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 7)
  colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "RNADeg")
  colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "RNADeg")
  row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
  row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
  head(GeneByCellTypeSubjVar2_Pvalues)
  
  
  for(i in c(1:length(SignalSortedNoNA3[,1]))){
    
    temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
    
    GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
    GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
    
  }
  
  GeneByCellTypeSubjVar2_Pvalues2<-cbind(GeneByCellTypeSubjVar2_Pvalues, RMAExpression_customCDFAnnotation2plus2)
  GeneByCellTypeSubjVar2_Betas2<-cbind(GeneByCellTypeSubjVar2_Betas, RMAExpression_customCDFAnnotation2plus2)
  
  
  write.csv(GeneByCellTypeSubjVar2_Pvalues2, "GeneByCellTypeSubjVar2_Pvalues.csv")
  write.csv(GeneByCellTypeSubjVar2_Betas2, "GeneByCellTypeSubjVar2_Betas.csv")
  
  
  GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 7)
  colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "RNADeg")
  
  for(i in c(1:length(SignalSortedNoNA3[,1]))){
    
    temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
    
    GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
    
  }
  
  GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)
  
  write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")
  
  
  
  for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
    png(paste(paste("17 Histogram of Raw Pvalues for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), "png", sep="."))	
    hist(GeneByCellTypeSubjVar2_Pvalues[,i], breaks=100, col=i, main=paste("Raw P-values for", colnames(GeneByCellTypeSubjVar2_Pvalues)[i], sep="  "), xlab="Unadjusted p-value", ylab="Count")
    abline(a=(length(GeneByCellTypeSubjVar2_Pvalues[,1])/100), b=0)
    dev.off()		
  }	
  
  
  GeneByCellTypeSubjVar2_PvaluesAdj<-matrix(0, length(GeneByCellTypeSubjVar2_Pvalues[,1]), length(GeneByCellTypeSubjVar2_Pvalues[1,]))
  colnames(GeneByCellTypeSubjVar2_PvaluesAdj)<-colnames(GeneByCellTypeSubjVar2_Pvalues)
  row.names(GeneByCellTypeSubjVar2_PvaluesAdj)<-row.names(GeneByCellTypeSubjVar2_Pvalues)
  
  
  library(multtest)
  
  for (i in c(1:length(GeneByCellTypeSubjVar2_Pvalues[1,]))){
    
    #Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):
    TempPvalAdj<-mt.rawp2adjp(GeneByCellTypeSubjVar2_Pvalues[,i], proc=c("BH"))
    GeneByCellTypeSubjVar2_PvaluesAdj[,i]<-TempPvalAdj$adjp[order(TempPvalAdj$index),2]
    
  }
  
  GeneByCellTypeSubjVar2_PvaluesAdj2<-cbind(GeneByCellTypeSubjVar2_PvaluesAdj, RMAExpression_customCDFAnnotation2plus2)
  write.csv(GeneByCellTypeSubjVar2_PvaluesAdj2, "GeneByCellTypeSubjVar2_PvaluesAdj.csv")
  
  GeneByCellTypeSubjVar2DF<-as.data.frame(cbind(GeneByCellTypeSubjVar2_Betas, GeneByCellTypeSubjVar2_Pvalues, GeneByCellTypeSubjVar2_PvaluesAdj))
  
  temp<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2DF)
  write.csv(temp, "GeneByCellTypeSubjVar2DF.csv" )
  
  ########################
  
  
  row.names(AveragePrimary_CellTypeIndex)
  Astrocyte_All<-AveragePrimary_CellTypeIndex[1,]
  Oligodendrocyte_All<-AveragePrimary_CellTypeIndex[8,]
  Microglia_All<-AveragePrimary_CellTypeIndex[3,]
  Endothelial_All<-AveragePrimary_CellTypeIndex[2,]
  RBC_All<-AveragePrimary_CellTypeIndex[10,]
  Neuron_All<-AveragePrimary_CellTypeIndex[5,]
  Neuron_Interneuron<-AveragePrimary_CellTypeIndex[6,]
  Neuron_Projection<-AveragePrimary_CellTypeIndex[7,]
  Oligodendrocyte_Immature<-AveragePrimary_CellTypeIndex[9,]
  Mural_All<-AveragePrimary_CellTypeIndex[4,]
  
  summary.lm(lm(PC1noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
  
  summary.lm(lm(PC2noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
  
  summary.lm(lm(PC3noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
  
  summary.lm(lm(PC4noOutliers~Astrocyte_All+Oligodendrocyte_All+Endothelial_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+Neuron_All+RBC_All+Oligodendrocyte_Immature+Mural_All+BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+RNADegradPerSample))
  
  
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/LM4_PrevelentCelltypes")
  
  GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 12)
  GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 12)
  colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")
  colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")
  row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
  row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
  head(GeneByCellTypeSubjVar2_Pvalues)

  GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 12)
  colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Oligodendrocyte", "Microglia", "Neuron_Projection", "Neuron_Interneuron", "RNADeg")  
  
  for(i in c(1:length(SignalSortedNoNA3[,1]))){
    
    temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Oligodendrocyte_All+Microglia_All+Neuron_Projection+Neuron_Interneuron+RNADegradPerSample))
    
    GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
    GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
    GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
    
  }
  

  
  GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)
  
  write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")
  
  #########################################
  
  setwd("~/Documents/Microarray Gen/Barnes_GSE21935/LM4_Everything")
  #I don't know if the sample size for this dataset can actually handle this.
  
  GeneByCellTypeSubjVar2_Pvalues<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
  GeneByCellTypeSubjVar2_Betas<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
  colnames(GeneByCellTypeSubjVar2_Pvalues)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg" )
  colnames(GeneByCellTypeSubjVar2_Betas)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg")
  row.names(GeneByCellTypeSubjVar2_Pvalues)<-row.names(SignalSortedNoNA3)
  row.names(GeneByCellTypeSubjVar2_Betas)<-row.names(SignalSortedNoNA3)
  head(GeneByCellTypeSubjVar2_Pvalues)
  
  GeneByCellTypeSubjVar2_Tstat<-matrix(0, length(SignalSortedNoNA3[,1]), 16)
  colnames(GeneByCellTypeSubjVar2_Tstat)<-c("Intercept", "BrainPH", "PMI", "Age", "Gender", "SCHIZ", "Astrocyte", "Endothelial","Microglia","Mural","Neuron_All",  "Neuron_Projection", "Neuron_Interneuron",   "Oligodendrocyte", "RBC", "RNADeg")
  
  for(i in c(1:length(SignalSortedNoNA3[,1]))){
    
    temp<-summary.lm(lm(SignalSortedNoNA3[i,]~BrainpHcentered + PMIcentered+ Agecentered+Gender+Diagnosis+Astrocyte_All+Endothelial_All+Microglia_All+Mural_All+Neuron_All+Neuron_Projection+Neuron_Interneuron+Oligodendrocyte_All+RBC_All+RNADegradPerSample))
    
    GeneByCellTypeSubjVar2_Betas[i,]<-temp$coefficients[,1]
    GeneByCellTypeSubjVar2_Tstat[i,]<-temp$coefficients[,3]
    GeneByCellTypeSubjVar2_Pvalues[i,]<-temp$coefficients[,4]
    
  }

  
  GeneByCellTypeSubjVar2_Tstat2<-cbind(RMAExpression_customCDFAnnotation2plus2, GeneByCellTypeSubjVar2_Tstat)
  
  write.csv(GeneByCellTypeSubjVar2_Tstat2, "GeneByCellTypeSubjVar2_Tstat.csv")
 
  ########################
