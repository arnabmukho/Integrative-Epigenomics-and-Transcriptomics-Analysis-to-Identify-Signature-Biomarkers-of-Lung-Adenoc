getwd()

setwd("/data/AM/meth")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(EDASeq)
library(ConsensusClusterPlus)
library(dplyr)
library(survminer)
library(edgeR)
library(sesameData)
library(sesame)

#Parameters to search and access the datasets in the TCGA database
query.exp <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

#Download the files in chunks
GDCdownload(
  query = query.exp,
  files.per.chunk = 100
)

#Extract the downloaded files and process in .rda format
luad.exp <- GDCprepare(
  query = query.exp, 
  save = TRUE, 
  save.filename = "luadExp.rda"
)

# list of samples are Primary Tumor
samples.primary.tumour <- luad.exp$barcode[luad.exp$shortLetterCode == "TP"]

# list of samples are solid tissue normal
samples.solid.tissue.normal <- luad.exp$barcode[luad.exp$shortLetterCode == "NT"]

# Preprocessing the data

dataPrep <- TCGAanalyze_Preprocessing(
  object = luad.exp, 
  cor.cut = 0.6
)      

#Data normalization
dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

# filtering transcripts by selecting a threshold
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)   

# Differential expression analysis
dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,samples.solid.tissue.normal],
  mat2 = dataFilt[,samples.primary.tumour],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT",
  pipeline = "edgeR"
)

#Visualizng the DEGs
TCGAVisualize_volcano(
  x = dataDEGs$logFC,
  y = dataDEGs$FDR,
  filename = "volcanoexp.png",
  x.cut = 3,
  y.cut = 10^-5,
  names = rownames(dataDEGs),
  color = c("black","red","darkgreen"),
  names.size = 2,
  xlab = " Gene expression fold change (Log2)",
  legend = "State",
  title = "Volcano plot (Primary Tumor vs Solid Tissue Normal)",
  width = 10
)

#export the list of DEGs
write.csv(dataDEGs, "/data/AM/DEG.csv", row.names = TRUE)

#Enrichment analysis of the DEGs
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes Normal Vs Tumor",
  RegulonList = dataDEGs$gene_name
)  

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  GOBPTab = ansEA$ResBP,
  GOCCTab = ansEA$ResCC,
  GOMFTab = ansEA$ResMF,
  PathTab = ansEA$ResPat,
  nRGTab = dataDEGs$gene_name,
  nBar = 10,
  filename = "Enrichment_analysis.png"
)

datFilt <- dataNorm %>% 
  TCGAanalyze_Filtering(method = "varFilter") %>%
  TCGAanalyze_Filtering(method = "filter1") %>%  
  TCGAanalyze_Filtering(method = "filter2",foldChange = 1)

data_Hc2 <- TCGAanalyze_Clustering(
  tabDF = datFilt,
  method = "consensus",
  methodHC = "ward.D2"
) 
# Add  cluster information to Summarized Experiment
colData(luad.exp)$groupsHC <- paste0("EC",data_Hc2[[6]]$consensusClass)
TCGAanalyze_survival(
  data = colData(luad.exp),
  clusterCol = "groupsHC",
  main = "TCGA kaplan meier survival plot from consensus cluster",
  legend = "mRNA Group",
  height = 10,
  risk.table = T,
  conf.int = F,
  color = c("black","red","blue","green3", "violet", "cyan"),
  filename = "survival_luad_expression_subtypes.png"
)


TCGAvisualize_Heatmap(
  data = t(datFilt),
  col.metadata =  colData(luad.exp)[,
                                    c("barcode",
                                      "groupsHC")
  ],
  col.colors =  list(
    groupsHC = c(
      "EC1"="black",
      "EC2"="red",
      "EC3"="blue",
      "EC4"="green3",
      "EC5"="violet",
      "EC6"="cyan")
  ),
  sortCol = "groupsHC",
  type = "expression", # sets default color
  scale = "row", # use z-scores for better visualization. Center gene expression level around 0.
  title = "Heatmap from concensus cluster", 
  filename = "case2_Heatmap.png",
  extremes = seq(-2,2,1),
  color.levels = colorRampPalette(c("green", "black", "red"))(n = 5),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  width = 1280,
  height = 720
)

BiocManager::install("ComplexHeatmap")

#Methylation

query.met <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(
  query = query.met, 
  files.per.chunk = 20,
)

luad.met <- GDCprepare(
  query = query.met,
  save = FALSE,
)

luad.met <- luad.met[rowSums(is.na(assay(luad.met))) == 0,]

luad.met <- TCGAanalyze_DMC(luad.met, groupCol = "sample_type",
                            group1 = "Primary Tumor",
                            group2="Solid Tissue Normal",
                            p.cut = 10^-5,
                            diffmean.cut = 0.25,
                            legend = "State",
                            plot.filename = "metvolcano.png")

#starburst plot
starburst <- TCGAvisualize_starburst(
  met = luad.met, 
  exp = dataDEGs,
  genome = "hg19",
  group1 = "Primary Tumor",
  group2 = "Solid Tissue Normal",
  filename = "starburst.png",
  met.platform = "Illumina Human Methylation 450",
  met.p.cut = 10^-5,
  exp.p.cut = 10^-5,
  diffmean.cut = 0.25,
  logFC.cut = 3,
  names = FALSE, 
  height = 10,
  width = 15,
  dpi = 300)

BiocManager::install("ELMER")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(ELMER)
library(parallel)

query.met <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "DNA Methylation", 
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(
  query = query.met, 
  files.per.chunk = 20,
)

luad.met <- GDCprepare(
  query = query.met,
  save = TRUE, 
  save.filename = "luadDNAmet.rda",
  summarizedExperiment = TRUE
)

query.exp <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(
  query = query.exp,
  files.per.chunk = 100
)

luad.exp <- GDCprepare(
  query = query.exp, 
  save = TRUE, 
  save.filename = "luadExp.rda"
)

distal.probes <- get.feature.probe(genome = "hg38", met.platform = "450K")

library(MultiAssayExperiment)

mae <- createMAE(
  exp = luad.exp, 
  met = luad.met,
  save = FALSE,
  linearize.exp = TRUE,
  filter.probes = distal.probes,
  save.filename = "mae_luad.rda",
  met.platform = "450K",
  genome = "hg38",
  TCGA = TRUE
)

# Removing FFPE samples
mae <- mae[,!mae$is_ffpe]

# identify probes that are hypo and hypermethylated in tumor samples compared to the normal samples.
group.col <- "sample_type"
group1 <-  "Primary Tumor"
group2 <- "Solid Tissue Normal"
direction <- c("hypo", "hyper")
dir.out <- file.path("luad",direction)
dir.create(dir.out, recursive = TRUE)

#Differential methylation

sig.diff <- get.diff.meth(
  data = mae, 
  group.col = group.col,
  group1 =  group1,
  group2 = group2,
  minSubgroupFrac = 0.2,
  sig.dif = 0.3,
  diff.dir = direction,
  dir.out = dir.out, 
  pvalue = 0.01
)

#significant probe-gene pairs 

# determining nearby 20 genes for Significant probes
# 20 upstream and 20 dowstream genes
nearGenes <- GetNearGenes(
  data = mae, 
  probes = sig.diff$probe, 
  numFlankingGenes = 20 
)

#gene-pair
pair <- get.pair(
  data = mae,
  group.col = group.col,
  group1 =  group1,
  group2 = group2,
  nearGenes = nearGenes,
  minSubgroupFrac = 0.4,
  permu.dir = file.path(dir.out,"permu"),
  permu.size = 100000,
  raw.pvalue  = 0.05,   
  Pe = 0.001,
  filter.probes = TRUE,
  filter.percentage = 0.05,
  filter.portion = 0.3,
  dir.out = dir.out,
  cores = 1,
  label = direction
)

# Identify enriched motif for significantly hypomethylated probes which have putative target genes.

enriched.motif <- get.enriched.motif(
  data = mae,
  probes = pair$Probe, 
  dir.out = dir.out, 
  label = direction,
  min.incidence = 10,
  lower.OR = 1.1
)

#get list of Transcription factors

TF <- get.TFs(
  data = mae, 
  group.col = group.col,
  group1 =  group1,
  group2 = group2,
  minSubgroupFrac = 0.4,
  enriched.motif = enriched.motif,
  dir.out = dir.out, 
  cores = 1, 
  label = direction
)

#relationship between nearby 20 gene expression vs DNA methylation probe using a linear regression curve

scatter.plot(
  data = mae,
  byProbe = list(probe = sig.diff$probe[1], numFlankingGenes = 20), 
  category = "definition", 
  dir.out = "plots",
  lm = TRUE,
  save = TRUE
) 

#The anti-correlated pairs of gene and probes with sample_type and smoking status groups using a heat map.
heatmapPairs(
  data = mae, 
  group.col = "sample_type",
  group1 = "Primary Tumor", 
  annotation.col = c("paper_Smoking.Status"),
  group2 = "Solid Tissue Normal",
  pairs = pair,
  filename =  "heatmap_1.pdf"
)
