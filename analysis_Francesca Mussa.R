###############################################################################
############################EXPRESSION ANALYSIS STEPS###########################
############

#GOALS: TO SEE THE DIFFERENCE GENE EXPRESSION BETWEEN ERBB2 CANCER AND OTHER 
#SUBTYPES OF BREAST CANCER 

#####STEP 1: read the RNAseq file, Patient file and abberrations file

data_clinical_patient <- read.delim("data_clinical_patient.txt",
                                    comment.char="#")
View(data_clinical_patient)
#dataset contain two different ID per each patient, subtype of breast cancer, 
#demographic information....
#There are 1084 patients

data_mrna_seq_v2_rsem <- read.delim("data_mrna_seq_v2_rsem.txt")
View(data_mrna_seq_v2_rsem)
#Columns are patients while rows are genes ID. 
#Data represent the number of read for each patient and for each gene


data_cna <- read.delim("data_cna.txt")
View(data_cna)
#Columns are patients while rows are genes ID. 
#Data represent the number of abberrations for each gene
#in each blood sample. if the number is >0, the gene is amplified.  

################################################################################
#### STEP 2: DEMONGRAPHIC CHARACTERISTICS OF PATIENTS

require(skimr)
my_skim<-skim_with(base=sfl(),numeric=sfl(hist=NULL))
my_skim(data_clinical_patient$AGE)
table(data_clinical_patient$SEX)


################################################################################
#####STEP 3: CREATE METADATA 
#divide between ERBB2 amplificated and ERBB2 not amplificated among cancer subtypes


row_ERBB2<- as.integer(which(data_cna$Hugo_Symbol == 'ERBB2'))

data_only_mrnaseq<-as.matrix(data_mrna_seq_v2_rsem[,-c(1,2)])
data_aberration<-as.matrix(data_cna[,-c(1,2)])#dataset without 2 first columns; only with samples

#Create a list of patient ID 
pat_rna = colnames(data_only_mrnaseq) 
pat_cna = colnames(data_aberration)


#overlap two ID
index_overlap  = is.element(pat_rna, pat_cna) 


rna_subset = data_only_mrnaseq[, index_overlap] # get subset of rna patients that have cna data.

rna_pat_ids_subsets = colnames(rna_subset) # get the patients ids that have both rnaseq and cna.

# intialize memory for metadata.

metadata = matrix(0, dim(rna_subset)[2],1) #matrix with 1 column and rows= number of samples; only 0 values

dimnames(metadata) <- list(NULL, "Amplified")

j=1
for (i in 1:length(rna_pat_ids_subsets)) {
  index = which(pat_cna == rna_pat_ids_subsets[i])
  metadata[j,] = 1*(data_cna[row_ERBB2,index] > 0)
  j=j+1
}
  

metadata[is.na(metadata)] =0  # missing value are 0


hist(metadata, xlabel="Aplified", ylabel="Frequency", col="lightblue", 
   main="Frequency of Amplified" )


################################################################################
##### STEP 4: DIFFERENTIAL ANALYSIS

##STEP4.1: Create Deseqdataset
library(DESeq2)
library(dplyr)

rna_subset[is.na(rna_subset)] = 0  # Impute with zeros the NA
rna_subset[rna_subset<0] = 0

#Create an unique id for gene
id_gene<-paste(data_mrna_seq_v2_rsem$Hugo_Symbol, 
               data_mrna_seq_v2_rsem$Entrez_Gene_Id, sep = "_") 
rna_subset<-cbind(id_gene,rna_subset)
rna_subset<-data.frame(rna_subset)

if (any(duplicated(rna_subset$id_gene))) {
  warning("Duplicate row names. Adding a suffix to make them unique.")
  rna_subset$id_gene <- make.unique(rna_subset$id_gene)
}

#Insert the unique id_gene in rna data
rownames(rna_subset)<-rna_subset$id_gene
rna_subset[, -1] <- lapply(rna_subset[, -1], as.numeric)

dds <- DESeqDataSetFromMatrix(countData = round(rna_subset[,-1]), #n°reads without id_gene
                              colData = metadata, 
                              design = ~Amplified)


###############################################################################
##STEP 4.2: Filter

##Filter data where we have only more than 10 read count

smallestGroupSize <- 3 #reads need to be in at least  3 samples
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

################################################################################
##STEP 4.3: Normalize

dds <- DESeq(dds)
res<-results(dds)
summary(res)

################################################################################
##STEP 4.4: Do principal component analysis

# Transform the data to visualize; extraction of normalized data
rld <- vst(dds, blind=FALSE)



# Do Principal Components Analysis on normalized data
pc = prcomp(t(assay(rld)))  #exchange columns with rows; we want rows=patients

summary_pca <- summary(pc)
print(pc$rotation[, 1:2])


#Plot for first 2 PCA
plotPCA(rld, intgroup="Amplified")

#VISUALZIATION
plotMA(res, ylim=c(-2,2), main="Fold-change versus normalized mean counts ")
legend("topright", legend = c("Significant", "Not Significant"), col = c("blue", "grey"), pch = 19)


################################################################################
#####STEP 6: TOP 10 DIFFERENTAILLY EXPRESSED GENES


res <- res[order(-abs(res$pvalue)), ]
signif = which(res$padj<0.05) #selection genes significant expressed
deg = res[signif,]
nrow(deg) 
top10_genes <- head(deg, 10)
top10_genes


################################################################################
#####STEP 7: PATHWAY ENRIGHMENT

# Separate between gene upregulate and downregulate

dup = deg[deg[,2]>0.,]

dim(dup) #upregulated genes

ddown = deg[deg[,2]<0.,]

dim(ddown) 

#Take entrez_id because KEGG doesn't read a ID with also letters.
entrez_id_up=data_mrna_seq_v2_rsem[keep[signif[deg[,2]>0.]],2] #take gene upregualted; they have log2FoldChange>0
entrez_id_down=data_mrna_seq_v2_rsem[keep[signif[deg[,2]<0.]],2] #take gene downregulate; they have log2FoldChange<0

# Pathway Enrichment

BiocManager::install("clusterProfiler")

library(clusterProfiler)

up_paths = enrichKEGG(gene = entrez_id_up, organism = 'hsa', pvalueCutoff = 0.05)
head(up_paths)


down_paths = enrichKEGG(gene = entrez_id_down, organism = 'hsa', pvalueCutoff = 0.05)
head(down_paths)


#Visualization pathways
BiocManager::install("pathview")
library(pathview)

#I chose to analyze a upregualted pathway

pathview(gene.data = entrez_id_up, pathway.id = up_paths$ID, species = "hsa", 
         kegg.native = TRUE, keytype = "kegg", multi.gene = TRUE)



###############################################################################
#####STEP 8: ClUSTER ANALYSIS

hc <- hclust(dist(deg), method="complete") #hierarchical clustering on the distance matrix
hc1 <- hclust(dist(deg[1:10,]), method="complete") #CA only on top 10 genes

library(tibble)
cut_tree <- cutree(hc, k=6)
table(cut_tree)

cut_tre_1<-enframe(x=cut_tree, name="gene", value="cluster")
head(cut_tree)
cut_tree1 <- cutree(hc1, k=2)
table(cut_tree1)

# Plot the dendrogram with color-coded clusters for the first 10 differentially expressed genes
plot(hc1, hang = -1, col.dendrogram = cut_tree, main="Cluster Dendogram for 10 genes")
abline(h=3000, col="red", lwd=2)


################################################################################
####STEP 9: SURVIVAL ANALYSIS

library(survival)
library(survminer)
library(tidyverse)
library(ggplot2)

#Create a new datset with all elements for calculating the survival
data_clinical_patient$PATIENT_ID1 = gsub("\\-", ".",data_clinical_patient$PATIENT_ID)

pat_cna=substr(pat_cna, 1, 12)
index_overlap1  = is.element(data_clinical_patient$PATIENT_ID1, pat_cna) 


dataset = data_clinical_patient[index_overlap1,] # get subset of  patients in both dataset

for (i in 1:length(dataset$PATIENT_ID1)) {
  index11 <- which(pat_cna == data_clinical_patient$PATIENT_ID1[i])
  
  # Verify the index11
  if (length(index11) > 0) {
    dataset$amplified[i] <- 1 * (data_cna[row_ERBB2, index11] > 0)
  } else {
    # If index11 is empty, to assign NA
    dataset$amplified[i] <- NA  
  }
}


#status infornation 
dataset$deceased<-ifelse(dataset$OS_STATUS=="0:LIVING", FALSE, TRUE)

#time: for all patients we take the last follow up day
dataset$survival<-ifelse(dataset$OS_STATUS=="0:LIVING", dataset$OS_MONTHS, #i supposed that OS was overall survival 
                                      data_clinical_patient$OS_MONTHS)

#fitting survival curve
fit<-survfit(Surv(survival, deceased)~amplified, data=dataset)
fit 

ggsurvplot(fit, data=dataset, pval=T, risk.table = T)

#Cox model

model<-coxph(Surv(survival, deceased)~amplified, data=dataset)
model

#How many women and male have amplified ERBB2

femmine_amplified<-dataset$SEX=="Female" & dataset$amplified==1 #323 on 1058 female
male_amplified<-dataset$SEX=="Male" & dataset$amplified==1 #4 on 12
table(male_amplified)

# Which subtype is most common?
if (any(dataset$SEX == "Female" & dataset$amplified == 1)) {
  print(table(dataset$SUBTYPE[dataset$SEX == "Female" & dataset$amplified == 1]))
} #most common is luminal A


if (any(dataset$SEX == "Male" & dataset$amplified == 1)) {
  print(dataset$SUBTYPE[dataset$SEX == "Male" & dataset$amplified == 1])
} 

#Ethnicity
if (any(dataset$SEX == "Female" & dataset$amplified == 1)) {
  print(table(dataset$ETHNICITY[dataset$SEX == "Female" & dataset$amplified == 1]))
}


if (any(dataset$SEX == "Male" & dataset$amplified == 1)) {
  print(table(dataset$ETHNICITY[dataset$SEX == "Male" & dataset$amplified == 1]))
}

