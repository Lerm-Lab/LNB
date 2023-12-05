# This is a script to generate the network modules using MODifieR
# Created and modified by Sandra Hellberg, 2023-01-12
# Email: sandra.hellberg@liu.se

# Description: Run module inference using the R package MODifieR (de Weerd et al., Bioinformatics, 2020; https://gitlab.com/Gustafsson-lab/MODifieR)

################################################################################

# 1. Setup
# 2. Load and format input files
# 3. Module inference
# 4. Consensus module



# 1. Setup ####
pack_R <- c("dplyr", "MODifieR", "VennDiagram", "org.Hs.eg.db", "clusterProfiler")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(405)


# Set directory structure 

setwd("P:/LiU/2023_LNB/Data/")
main_dir <- "P:/LiU/2023_LNB/"
input_dir <- "P:/LiU/2023_LNB/RDS_files/"

library(MODifieR)

# 2. Load and format input data ####
input_genes <- read.delim("LNB_DMGs.txt") # input list of differentially methylated genes derived from mapping the DMCs to genes based on the annotation provided by Illumina
class(input_genes)

genes_for_entrez <- input_genes[,1]

# translate symbol into entrez id using bitr
genes_entrez <- bitr(genes_for_entrez, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)  # convert symbol to entrez

genes_entrez$pval <- 0 # create a column of p values and set them all to 0

input_modifier <- genes_entrez[,c(2:3)]

ppi_network <- read.delim("PPI_network.txt") # Import the full PPI network


# 3. Module inference ####
genes_modifier <- create_custom_microarray_input_object(diff_genes = input_modifier)


# MCODE ####
module_mcode <- mcode(MODifieR_input = genes_modifier,
                      ppi_network = ppi_network, 
                      hierarchy = 1, 
                      vwp = 0.5,
                      haircut = FALSE,
                      fluff = FALSE,
                      fdt = 0.8,
                      loops = TRUE)

mcode_output <- module_mcode$module_genes

# check how many of the seeds that are kept in the module
seed_genes <- intersect(mcode_output, input_modifier$ENTREZID)
mcode_ignored <- setdiff(input_modifier$ENTREZID, seed_genes) # set diff to identify the ignored seed genes

mcode_output <- bitr(mcode_output, fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)

saveRDS(mcode_output,paste0(input_dir,"LNB_mcode.RDS"))


# Diamond #####
module_diamond <- diamond(MODifieR_input = genes_modifier, 
                          ppi_network = ppi_network,
                          deg_cutoff = 0.05,
                          n_output_genes = 100,
                          seed_weight = 10,
                          include_seed = TRUE)



diamond_output <- module_diamond$module_genes


# translate entrez into symbol id using bitr
diamond_output <- bitr(diamond_output, fromType = "ENTREZID",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
saveRDS(diamond_output,paste0(input_dir,"LNB_diamond.RDS"))


# 4. Consensus module ####

consensus_module <- intersect(diamond_output$SYMBOL,mcode_output$SYMBOL) # combine Diamond and MCODE to form a consensus module

common_seeds <- intersect(consensus_module,genes_entrez$SYMBOL) # check number of common seed genes
saveRDS(common_seeds, paste0(input_dir, "LNB_module_consensus_seed.RDS"))

saveRDS(consensus_module, paste0(input_dir, "LNB__module_consensus.RDS"))
write.table(consensus_module, paste0(data_dir, "LNB_module_consensus.txt"), row.names = F, quote = F)
