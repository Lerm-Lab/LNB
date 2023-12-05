# This is a script to generate the network and network attributes for visualization in Cytoscape
# Created and modified by Sandra Hellberg, 2023-03-08
# Email: sandra.hellberg@liu.se
# Description: Create a network and network attributes for visualization in Cytoscape

################################################################################

# 1. Setup
# 2. Load and format input files
# 3. Create network for Cytoscape
# 4. Network attributes


# 1. Setup ####
R_packages <- c("biomaRt")

for (i in 1:length(R_packages)) {
  library(R_packages[i], character.only = T)
}

set.seed(123)

# Set working directory
setwd("P:/LiU/2023_LNB/")
input_dir <- "P:/LiU/2023_LNB/Data/"


# 2. Load and format input files ####

# Load mart 
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Load seed genes
seed_genes <- readRDS("RDS_files/LNB_module_consensus_seed.RDS")
seed_genes <- data.frame(Gene=seed_genes)

# Load modules
consensus_genes <- readRDS("RDS_files/LNB__module_consensus.RDS")
consensus_genes <- data.frame(Gene=consensus_genes)

# Load and format PPI network
ppi_network <- read.delim(paste0(input_dir, "PPI_network.txt"))

ppi_high <- ppi_network[ppi_network$combined_score>900,] # threshold 900 for combined_score

# 3. Create network for Cytoscape #### 
symbol_entrez <- getBM(attributes <- c("entrezgene_id", "external_gene_name"),
                       filters= "entrezgene_id",
                       values = unique(c(ppi_high$entrez1, ppi_high$entrez2)),
                       mart = mart)

ord1 <- match(ppi_high$entrez1, symbol_entrez$entrezgene_id)
ord2 <- match(ppi_high$entrez2, symbol_entrez$entrezgene_id)

idx <- (!is.na(ord1)) & (!is.na(ord2))

symbol1 <- symbol_entrez$external_gene_name[ord1[idx]]
symbol2 <- symbol_entrez$external_gene_name[ord2[idx]]
ppi_high <- data.frame(symbol1, symbol2)

ppi_module <- ppi_high[(ppi_high$symbol1 %in% consensus_genes$Gene) & (ppi_high$symbol2 %in% consensus_genes$Gene),]
write.table(ppi_module, "ppi_consensus_module_900.txt", row.names = F, quote = F) # network of the module to be used as input in Cytoscape

# 4. Network attributes ####

# attribute node size
attribute_node_size <- data.frame(table((ppi_module$symbol1)))
write.table(attribute_node_size, "ppi_consensus_module_900_nodesize.txt", row.names = F, quote = F)

# attribute seed genes
seed_in_module <- data.frame(Gene=intersect(consensus_genes$Gene, seed_genes$Gene))
seed_in_module$seed <- "YES"

idx <- which(is.element(consensus_genes$Gene,seed_genes$Gene))
not_in_module <- data.frame(Gene=consensus_genes[-idx,])
not_in_module$seed <- "NO"

attribute_seed <- rbind(seed_in_module, not_in_module)
write.table(attribute_seed, "ppi_consensus_module_900_seed.txt", row.names = F, quote = F)


# attributes hyper and hypo-methylation

logFC <- read.delim(paste0(input_dir, "LNB_MMD0.2_logFC.txt"))
logFC <- logFC[,c(2,3)]

colnames(logFC)[1] <- "GeneID"

logFC <- logFC[logFC$GeneID %in% seed_in_module$Gene,] 

genes_unique <- unique(logFC$GeneID)

# Divide into hyper and hypomethylated genes
hyper_dmgs_all <- logFC[logFC$logFC > 0, ]
hypo_dmgs_all <- logFC[logFC$logFC < 0, ]

# All probes for a gene should indicate the same direction of methylation to be considered hyper/hypo
hyper_dmgs <- hyper_dmgs_all[!hyper_dmgs_all$GeneID %in% hypo_dmgs_all$GeneID, ]
hypo_dmgs <- hypo_dmgs_all[!hypo_dmgs_all$GeneID %in% hyper_dmgs_all$GeneID, ]

DMG_list <- data.frame(Gene = unique(logFC$GeneID))

#classify genes
DMG_list[, "Direction"] <- "Mixed"
for (i in 1:nrow(DMG_list)) {
  if (DMG_list[i, "Gene"] %in% hyper_dmgs$GeneID) {
    DMG_list[i, "Direction"] <- "Hyper"
  } else if (DMG_list[i, "Gene"] %in% hypo_dmgs$GeneID) {
    DMG_list[i, "Direction"] <- "Hypo"
  }
}

write.table(DMG_list, "ppi_consensus_module_900_direction.txt", row.names = F, quote = F)

