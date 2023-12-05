# This is a script to run KEGG pathway enrichment on module genes
# Created and modified by Sandra Hellberg, 2023-02-15
# Email: sandra.hellberg@liu.se

# Description: Run pathway analysis using clusterProfiler

################################################################################

# 1. Setup
# 2. Load and format input files
# 3. KEGG pathway analysis
# 4. Visualization


# 1. Setup ####
R_packages <- c("magrittr", "stringr", "DOSE", "clusterProfiler","org.Hs.eg.db", "ggplot2", "viridis")

for (i in 1:length(R_packages)) {
  library(R_packages[i], character.only = T)
}

# Set directory structure
setwd("P:/LiU/2023_LNB/Data/")
input_dir <- "P:/LiU/2023_LNB/RDS_files/"
data_dir <- "P:/LiU/2023_LNB/Data/"

# 2. Load and format input files ####
consensus_module <- readRDS(paste0(input_dir, "LNB_module_consensus.RDS"))
consensus_module <- data.frame(Gene=consensus_module)

background_illumina <- read.delim(paste0(data_dir, "genes_ProbeFeatures.txt"))
background_illumina <- bitr(background_illumina$Gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) # convert symbol to ENTREZID
background_illumina <- data.frame(ENTREZID= background_illumina[,2])


background_ppi <- read.delim("PPI_network.txt")
background_ppi <- data.frame(ENTREZID=background_ppi[,1])
background_ppi <- unique(background_ppi) # all genes included in the PPI

# create a combined background with probefeatures and the PPI
universe <- merge(background_illumina, background_ppi, by="ENTREZID", all=TRUE )
table(duplicated(universe$ENTREZID)) #check if duplicated entrez ID


# 3. KEGG pathway analysis
# input data
pathway_input <- bitr(consensus_module$Gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db" ) # convert symbol to ENTREZID
table(duplicated(pathway_input$ENTREZID)) #check if duplicated entrez ID

kegg1 <- enrichKEGG(gene= pathway_input$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg",
                    universe = universe$gene, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")

kegg1 <- setReadable(kegg1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

heatplot(kegg1, showCategory=10) # plot was exported using export in R and edited using Adobe Illustrator
dotplot(kegg1, showCategory=10, font.size=15) + viridis::scale_fill_viridis()

summary_enrich <- kegg1@result
write.table(summary_enrich, "KEGG_consensus_module.txt", sep= ",", quote=F, row.names = F) # results provided in Supplemental Table S4
