
# Useful resources -----
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# https://rpubs.com/pranali018/enrichmentAnalysis

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(clusterProfiler, enrichplot, dplyr, org.Hs.eg.db, ggplot2)


# Univariate results -----
# KEGG only significant ----
proteomics = read.csv("./2025_12_17_proteins_associated_GDM.csv")

universe = proteomics 

proteomics = proteomics %>%
  dplyr::filter(adj.P.Val < 0.05)

keytypes(org.Hs.eg.db)


converted = bitr(proteomics$Protein, 
                 fromType = "ALIAS", 
                 toType = "ENTREZID", 
                 OrgDb = 'org.Hs.eg.db')

universe = bitr(universe$Protein,
                fromType = "ALIAS",
                toType = "ENTREZID",
                OrgDb = 'org.Hs.eg.db')

ego_kegg = enrichKEGG(gene  = converted$ENTREZID,
                 keyType = 'ncbi-geneid',
                 organism = 'hsa',
                 universe = universe$ENTREZID,
                 minGSSize = 1, 
                 maxGSSize = 501,
                 pvalueCutoff = 0.5,
                 pAdjustMethod = "BH"
                 )

head(ego_kegg, n = 20)

x = setReadable(ego_kegg, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)

results = x2@result %>%
  dplyr::filter(Count > 1)

write.csv(results, "./ora_kegg_univariate_results.csv", row.names = FALSE)

dotplot(x2)
upsetplot(x2)
heatplot(x2)
cnetplot(x2)
emapplot(x2)
treeplot(x2)


# GO CC only significant ----
ego_go = enrichGO(gene  = converted$ENTREZID,
                      OrgDb = 'org.Hs.eg.db',
                      keyType = 'ENTREZID',
                      ont = "CC", # ALL
                      minGSSize = 2, 
                      maxGSSize = 800,
                      readable = TRUE,
                      pvalueCutoff = 0.05
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(x2)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2)
cnetplot(x2)
emapplot(x2)
# treeplot(x2)

# GO BP only significant ----

ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "BP", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.07
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(ego_go)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2)
cnetplot(x2)
emapplot(x2)
treeplot(x2)

# GO MF only significant ----

ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "MF", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.05
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(ego_go)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2, label_format = 70)
cnetplot(x2)
emapplot(x2)
treeplot(x2)


# Multivariate results -----
# KEGG only significant ----

proteomics = read.csv("./2025_04_22_proteins_associated_GDM.csv")

universe = proteomics 

proteomics = readr::read_rds("./coef_proteomic_signature.RDS")

proteomics = proteomics %>%
  dplyr::rename(Protein = Coefficient, Coefficient = s1)

keytypes(org.Hs.eg.db)


converted = bitr(proteomics$Protein, 
                 fromType = "ALIAS", 
                 toType = "ENTREZID", 
                 OrgDb = 'org.Hs.eg.db')

universe = bitr(universe$Protein,
                fromType = "ALIAS",
                toType = "ENTREZID",
                OrgDb = 'org.Hs.eg.db')

ego_kegg = enrichKEGG(gene  = converted$ENTREZID,
                      keyType = 'ncbi-geneid',
                      organism = 'hsa',
                      # universe = universe$ENTREZID,
                      minGSSize = 2, 
                      maxGSSize = 800,
                      pvalueCutoff = 0.42,
                      pAdjustMethod = "fdr"
)

head(ego_kegg, n = 20)

x = setReadable(ego_kegg, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)

results = x2@result %>%
  dplyr::filter(Count > 6)

write.csv(results, "./ora_kegg_multivariate_results.csv", row.names = FALSE)

# results = readr::read_csv("./ora_kegg_multivariate_results.csv")

dotplot(x2, showCategory = results$Description)
upsetplot(x2)
heatplot(x2, showCategory = results$Description)
cnetplot(x2, showCategory = results$Description)
emapplot(x2, showCategory = results$Description)
treeplot(x2, showCategory = results$Description)


# GO CC only significant ----
ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "ALL", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.05
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(x2)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2)
cnetplot(x2)
emapplot(x2)
# treeplot(x2)

# GO BP only significant ----

ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "BP", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.07
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(ego_go)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2)
cnetplot(x2)
emapplot(x2)
treeplot(x2)

# GO MF only significant ----

ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "MF", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.05
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(ego_go)

# dotplot(x2, showCategory = results$Description, 
#         split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") 

# upsetplot(x2)
heatplot(x2, label_format = 70)
cnetplot(x2)
emapplot(x2)
treeplot(x2)


# GO ALL only significant ----

ego_go = enrichGO(gene  = converted$ENTREZID,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'ENTREZID',
                  ont = "ALL", 
                  minGSSize = 2, 
                  maxGSSize = 800,
                  readable = TRUE,
                  pvalueCutoff = 0.05
)

x = setReadable(ego_go, "org.Hs.eg.db", "ENTREZID")
x2 = pairwise_termsim(x)
results = x2@result %>%
  dplyr::filter(Count > 1)

dotplot(ego_go)

dotplot(x2, showCategory = 5,
        split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

# upsetplot(x2)
heatplot(x2, label_format = 70)
cnetplot(x2)
emapplot(x2)
treeplot(x2)

