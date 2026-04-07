

library(fgsea)
library(Biobase)
library(dplyr)
library(tibble)

# limma ----
metabolites = read.csv("/mnt/project/results/1_limma_analysis/2025_12_15_metabolites_associated_GDM.csv")

ranking = metabolites %>%
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::select(Metabolite, logFC) %>%
  dplyr::arrange(desc(logFC)) %>%
  tibble::deframe()



met_classes = list(
  Cholesterol = c("total_cholesterol", "total_cholesterol_minus_hdl_c", "remnant_cholesterol_non_hdl_non_ldl_cholesterol", 
                  "vldl_cholesterol",  "cholesterol_in_idl", "clinical_ldl_cholesterol", "ldl_cholesterol", "hdl_cholesterol",
                  "cholesterol_in_chylomicrons_and_extremely_large_vldl", "cholesterol_in_very_large_vldl", "cholesterol_in_large_vldl", "cholesterol_in_medium_vldl", "cholesterol_in_small_vldl", "cholesterol_in_very_small_vldl",
                  "cholesterol_in_large_ldl", "cholesterol_in_medium_ldl", "cholesterol_in_small_ldl",
                  "cholesterol_in_very_large_hdl",  "cholesterol_in_large_hdl", "cholesterol_in_medium_hdl", "cholesterol_in_small_hdl",
                  "cholesterol_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage", 
                  "cholesterol_to_total_lipids_in_very_large_vldl_percentage", "cholesterol_to_total_lipids_in_large_vldl_percentage", "cholesterol_to_total_lipids_in_medium_vldl_percentage", "cholesterol_to_total_lipids_in_small_vldl_percentage", "cholesterol_to_total_lipids_in_very_small_vldl_percentage", 
                  "cholesterol_to_total_lipids_in_idl_percentage", 
                  "cholesterol_to_total_lipids_in_large_ldl_percentage", "cholesterol_to_total_lipids_in_medium_ldl_percentage", "cholesterol_to_total_lipids_in_small_ldl_percentage",  
                  "cholesterol_to_total_lipids_in_very_large_hdl_percentage", "cholesterol_to_total_lipids_in_large_hdl_percentage", "cholesterol_to_total_lipids_in_medium_hdl_percentage", "cholesterol_to_total_lipids_in_small_hdl_percentage"),
  Triglycerides = c("total_triglycerides", 
                    "triglycerides_in_vldl", "triglycerides_in_idl", 
                    "triglycerides_in_ldl", "triglycerides_in_hdl", 
                    "triglycerides_in_chylomicrons_and_extremely_large_vldl", "triglycerides_in_very_large_vldl", "triglycerides_in_large_vldl", "triglycerides_in_medium_vldl", "triglycerides_in_small_vldl", "triglycerides_in_very_small_vldl", 
                    "triglycerides_in_large_ldl", "triglycerides_in_medium_ldl", "triglycerides_in_small_ldl", 
                    "triglycerides_in_very_large_hdl", "triglycerides_in_large_hdl",  "triglycerides_in_medium_hdl",  "triglycerides_in_small_hdl", 
                    "triglycerides_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage", 
                    "triglycerides_to_total_lipids_in_very_large_vldl_percentage", "triglycerides_to_total_lipids_in_large_vldl_percentage", "triglycerides_to_total_lipids_in_medium_vldl_percentage", "triglycerides_to_total_lipids_in_small_vldl_percentage", "triglycerides_to_total_lipids_in_very_small_vldl_percentage", 
                    "triglycerides_to_total_lipids_in_idl_percentage", 
                    "triglycerides_to_total_lipids_in_large_ldl_percentage", "triglycerides_to_total_lipids_in_medium_ldl_percentage", "triglycerides_to_total_lipids_in_small_ldl_percentage", 
                    "triglycerides_to_total_lipids_in_very_large_hdl_percentage", "triglycerides_to_total_lipids_in_large_hdl_percentage", "triglycerides_to_total_lipids_in_medium_hdl_percentage", "triglycerides_to_total_lipids_in_small_hdl_percentage"),
  Phospholipids = c("total_phospholipids_in_lipoprotein_particles", 
                    "phospholipids_in_vldl",  "phospholipids_in_idl", "phospholipids_in_ldl", "phospholipids_in_hdl", 
                    "phospholipids_in_chylomicrons_and_extremely_large_vldl", "phospholipids_in_very_large_vldl",  "phospholipids_in_large_vldl", "phospholipids_in_medium_vldl", "phospholipids_in_small_vldl", "phospholipids_in_very_small_vldl",
                    "phospholipids_in_large_ldl", "phospholipids_in_medium_ldl", "phospholipids_in_small_ldl",
                    "phospholipids_in_very_large_hdl", "phospholipids_in_large_hdl", "phospholipids_in_medium_hdl", "phospholipids_in_small_hdl", 
                    "phospholipids_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage", 
                    "phospholipids_to_total_lipids_in_very_large_vldl_percentage", "phospholipids_to_total_lipids_in_large_vldl_percentage", "phospholipids_to_total_lipids_in_medium_vldl_percentage", "phospholipids_to_total_lipids_in_small_vldl_percentage", "phospholipids_to_total_lipids_in_very_small_vldl_percentage", 
                    "phospholipids_to_total_lipids_in_idl_percentage", 
                    "phospholipids_to_total_lipids_in_large_ldl_percentage", "phospholipids_to_total_lipids_in_medium_ldl_percentage", "phospholipids_to_total_lipids_in_small_ldl_percentage", 
                    "phospholipids_to_total_lipids_in_very_large_hdl_percentage", "phospholipids_to_total_lipids_in_large_hdl_percentage", "phospholipids_to_total_lipids_in_medium_hdl_percentage", "phospholipids_to_total_lipids_in_small_hdl_percentage"),
  Cholesteryl_esters = c("total_esterified_cholesterol", 
                         "cholesteryl_esters_in_vldl", "cholesteryl_esters_in_idl", 
                         "cholesteryl_esters_in_ldl", "cholesteryl_esters_in_hdl", 
                         "cholesteryl_esters_in_chylomicrons_and_extremely_large_vldl",  "cholesteryl_esters_in_very_large_vldl", "cholesteryl_esters_in_large_vldl", "cholesteryl_esters_in_medium_vldl", "cholesteryl_esters_in_small_vldl", "cholesteryl_esters_in_very_small_vldl",
                         "cholesteryl_esters_in_large_ldl", "cholesteryl_esters_in_medium_ldl", "cholesteryl_esters_in_small_ldl",
                         "cholesteryl_esters_in_very_large_hdl", "cholesteryl_esters_in_large_hdl", "cholesteryl_esters_in_medium_hdl",   "cholesteryl_esters_in_small_hdl", 
                         "cholesteryl_esters_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage", 
                         "cholesteryl_esters_to_total_lipids_in_very_large_vldl_percentage", "cholesteryl_esters_to_total_lipids_in_large_vldl_percentage", "cholesteryl_esters_to_total_lipids_in_medium_vldl_percentage", "cholesteryl_esters_to_total_lipids_in_small_vldl_percentage", "cholesteryl_esters_to_total_lipids_in_very_small_vldl_percentage", 
                         "cholesteryl_esters_to_total_lipids_in_idl_percentage",
                         "cholesteryl_esters_to_total_lipids_in_large_ldl_percentage", "cholesteryl_esters_to_total_lipids_in_medium_ldl_percentage", "cholesteryl_esters_to_total_lipids_in_small_ldl_percentage",
                         "cholesteryl_esters_to_total_lipids_in_very_large_hdl_percentage", "cholesteryl_esters_to_total_lipids_in_large_hdl_percentage", "cholesteryl_esters_to_total_lipids_in_medium_hdl_percentage", "cholesteryl_esters_to_total_lipids_in_small_hdl_percentage"),
  Free_cholesterol = c("total_free_cholesterol",  
                       "free_cholesterol_in_vldl", "free_cholesterol_in_idl", 
                       "free_cholesterol_in_ldl", "free_cholesterol_in_hdl", 
                       "free_cholesterol_in_chylomicrons_and_extremely_large_vldl", "free_cholesterol_in_very_large_vldl", "free_cholesterol_in_large_vldl", "free_cholesterol_in_medium_vldl", "free_cholesterol_in_small_vldl", "free_cholesterol_in_very_small_vldl", 
                       "free_cholesterol_in_large_ldl", "free_cholesterol_in_medium_ldl", "free_cholesterol_in_small_ldl", 
                       "free_cholesterol_in_very_large_hdl", "free_cholesterol_in_large_hdl",  "free_cholesterol_in_medium_hdl",  "free_cholesterol_in_small_hdl", 
                       "free_cholesterol_to_total_lipids_in_chylomicrons_and_extremely_large_vldl_percentage", "free_cholesterol_to_total_lipids_in_very_large_vldl_percentage", "free_cholesterol_to_total_lipids_in_large_vldl_percentage",
                       "free_cholesterol_to_total_lipids_in_medium_vldl_percentage", "free_cholesterol_to_total_lipids_in_small_vldl_percentage", "free_cholesterol_to_total_lipids_in_very_small_vldl_percentage", 
                       "free_cholesterol_to_total_lipids_in_idl_percentage", 
                       "free_cholesterol_to_total_lipids_in_large_ldl_percentage", "free_cholesterol_to_total_lipids_in_medium_ldl_percentage",  "free_cholesterol_to_total_lipids_in_small_ldl_percentage", 
                       "free_cholesterol_to_total_lipids_in_very_large_hdl_percentage", "free_cholesterol_to_total_lipids_in_large_hdl_percentage", "free_cholesterol_to_total_lipids_in_medium_hdl_percentage", "free_cholesterol_to_total_lipids_in_small_hdl_percentage" 
  ),
  Total_lipids = c("total_lipids_in_lipoprotein_particles",  
                   "total_lipids_in_vldl", "total_lipids_in_idl", "total_lipids_in_ldl", "total_lipids_in_hdl",
                   "total_lipids_in_chylomicrons_and_extremely_large_vldl", "total_lipids_in_very_large_vldl",  "total_lipids_in_large_vldl", "total_lipids_in_medium_vldl", "total_lipids_in_small_vldl", "total_lipids_in_very_small_vldl", 
                   "total_lipids_in_large_ldl", "total_lipids_in_medium_ldl", "total_lipids_in_small_ldl",
                   "total_lipids_in_very_large_hdl", "total_lipids_in_large_hdl",  "total_lipids_in_medium_hdl", "total_lipids_in_small_hdl"),
  LP_concentration = c("total_concentration_of_lipoprotein_particles",  "concentration_of_chylomicrons_and_extremely_large_vldl_particles", "concentration_of_vldl_particles", "concentration_of_idl_particles",  "concentration_of_ldl_particles", "concentration_of_hdl_particles",
                       "concentration_of_very_large_vldl_particles", "concentration_of_large_vldl_particles",  "concentration_of_medium_vldl_particles", "concentration_of_small_vldl_particles", "concentration_of_very_small_vldl_particles", 
                       "concentration_of_large_ldl_particles", "concentration_of_medium_ldl_particles", "concentration_of_small_ldl_particles",
                       "concentration_of_very_large_hdl_particles", "concentration_of_large_hdl_particles", "concentration_of_medium_hdl_particles",  "concentration_of_small_hdl_particles"),
  Other_lipids = c("total_cholines", "phosphatidylcholines", "phosphoglycerides", "sphingomyelins", "triglycerides_to_phosphoglycerides_ratio"),
  ApoLP_size = c("apolipoprotein_b", "apolipoprotein_a1", "apolipoprotein_b_to_apolipoprotein_a1_ratio", "average_diameter_for_vldl_particles", "average_diameter_for_ldl_particles", "average_diameter_for_hdl_particles"),
  Fatty_acids = c("total_fatty_acids", "degree_of_unsaturation", "omega_3_fatty_acids", "omega_6_fatty_acids", "polyunsaturated_fatty_acids", "monounsaturated_fatty_acids", "saturated_fatty_acids", "linoleic_acid", "docosahexaenoic_acid", "omega_3_fatty_acids_to_total_fatty_acids_percentage", "omega_6_fatty_acids_to_omega_3_fatty_acids_ratio", "polyunsaturated_fatty_acids_to_total_fatty_acids_percentage", "monounsaturated_fatty_acids_to_total_fatty_acids_percentage", "saturated_fatty_acids_to_total_fatty_acids_percentage", "linoleic_acid_to_total_fatty_acids_percentage", "docosahexaenoic_acid_to_total_fatty_acids_percentage", "polyunsaturated_fatty_acids_to_monounsaturated_fatty_acids_ratio", "omega_6_fatty_acids_to_omega_3_fatty_acids_ratio"),
  Amino_acids = c("spectrometer_corrected_alanine", "alanine", "glutamine", "glycine", "histidine", "isoleucine", "leucine", "phenylalanine", "tyrosine", "valine", "total_concentration_of_branched_chain_amino_acids_leucine_isoleucine_valine"),
  Glycolysis = c("glucose", "lactate", "glucose_lactate", "pyruvate", "citrate"),
  Ketone_bodies = c( "acetate", "acetoacetate", "acetone", "x3_hydroxybutyrate"),
  Fluid_balance = c("creatinine", "albumin"),
  Inflammation = c("glycoprotein_acetyls")
)

fgsea_results = fgseaMultilevel(
  pathways = met_classes,
  stats = ranking,
  minSize = 1,
  maxSize = 50
)

fgsea_results = fgsea_results[order(fgsea_results$padj), ]



fgsea_results$leadingEdge = sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))


# Reemplazar "_" con espacio en los nombres de las clases
fgsea_results <- fgsea_results %>%
  mutate(pathway = gsub("_", " ", pathway))  # Sustituir "_" por espacios

# Crear una columna para marcar significancia con un asterisco
fgsea_results <- fgsea_results %>%
  mutate(Significance = ifelse(padj < 0.05, "*", ""))

# Crear el gráfico
library(ggplot2)
# msea_plot <-ggplot(fgsea_results, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity") +
#   geom_text(
#     aes(
#       label = Significance,
#       y = ifelse(NES > 0, NES + 0.2, NES - 0.2)  # Posición ajustada según NES
#     ),
#     size = 4, color = "black"
#   ) +
#   coord_flip() +
#   labs(
#     title = "Metabolite Set Enrichment Analysis",
#     x = "",
#     y = "Normalized Enrichment Score (NES)"
#   ) +
#   scale_fill_gradient2(
#     low = "#A7C7E7", mid = "white", high = "#F7A1A1",
#     midpoint = 0, space = "Lab", limits = c(-2.5, 2.5), oob = scales::squish
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centrar el título
#     axis.title.y = element_text(size = 12),
#     axis.title.x = element_text(size = 12),
#     axis.text.y = element_text(size = 10),
#     axis.text.x = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 9)
#   )

fgsea_results %>%  
  dplyr::mutate(direction = if_else(NES > 0, "Up-Regulated", "Down-Regulated")) %>% 
  ggplot(aes(y = reorder(pathway, NES), x = NES, color = pval, size = size)) + 
  facet_grid(. ~ direction, scales = "free_x", space = "free_x") + 
  scale_color_gradientn(name = "P-value", colors = c("red", "purple", "blue")) +
  scale_size(name = "Count", range = c(5, 12)) + 
  geom_point() +
  labs(y = "", x = "Enrichment score") + 
  theme_bw() + 
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"), 
    strip.text = element_text(size = 20),     
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(size = 15), 
    axis.title.x = element_text(size = 16)
  )



# elasticnet ----

metabolites = read.csv("/mnt/project/results/1_limma_analysis/2025_04_22_metabolites_associated_GDM.csv")

elasticnet_under = readRDS("/mnt/project/results/2_elasticnet/final_elasticnet_metabolomics.RDS")

coef_elasticnet = as.matrix(coef(elasticnet_under$finalModel, elasticnet_under$bestTune$lambda))
coef_elasticnet = coef_elasticnet %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "Coefficient") %>% 
  dplyr::filter(s1 != 0)

metabolites = coef_elasticnet %>% 
  dplyr::filter(Coefficient %in% metabolites$Metabolite)

ranking = metabolites %>%
  dplyr::arrange(desc(s1)) %>%
  tibble::deframe()

# met_classes = list(
#   Cholesterol = c("total_cholesterol", "total_cholesterol_minus_hdl_c", "remnant_cholesterol_non_hdl_non_ldl_cholesterol", "vldl_cholesterol", "clinical_ldl_cholesterol", "ldl_cholesterol", "hdl_cholesterol"),
#   Triglycerides = c("total_triglycerides", "triglycerides_in_vldl", "triglycerides_in_ldl", "triglycerides_in_hdl", "triglycerides_to_phosphoglycerides_ratio"),
#   Phospholipids = c("total_phospholipids_in_lipoprotein_particles", "phospholipids_in_vldl", "phospholipids_in_ldl", "phospholipids_in_hdl"),
#   Cholesteryl_esters = c("cholesteryl_esters_in_vldl", "cholesteryl_esters_in_ldl", "cholesteryl_esters_in_hdl"),
#   Free_Cholesterol = c("total_free_cholesterol", "free_cholesterol_in_vldl", "free_cholesterol_in_ldl", "free_cholesterol_in_hdl"),
#   Total_Lipids = c("total_lipids_in_lipoprotein_particles", "total_lipids_in_vldl", "total_lipids_in_ldl", "total_lipids_in_hdl"),
#   # Lipoprotein_particle_concentration = c("total_concentration_of_lipoprotein_particles", "concentration_of_vldl_particles", "concentration_of_ldl_particles", "concentration_of_hdl_particles"),
#   # Lipoprotein_particle_size = c("average_diameter_for_vldl_particles", "average_diameter_for_ldl_particles", "average_diameter_for_hdl_particles"),
#   # Other_lipids = c("phosphoglycerides", "total_cholines", "phosphatidylcholines", "sphingomyelins"),
#   # Apolipoproteins = c("apolipoprotein_b", "apolipoprotein_a1", "apolipoprotein_b_to_apolipoprotein_a1_ratio"),
#   Fatty_acids = c("total_fatty_acids", "degree_of_unsaturation", "omega_3_fatty_acids", "omega_6_fatty_acids", "polyunsaturated_fatty_acids", "monounsaturated_fatty_acids", "saturated_fatty_acids", "linoleic_acid", "cocosahexaenoic_acid", "omega_3_fatty_acids_to_total_fatty_acids_percentage", "omega_6_fatty_acids_to_omega_3_fatty_acids_ratio", "polyunsaturated_fatty_acids_to_total_fatty_acids_percentage", "monounsaturated_fatty_acids_to_total_fatty_acids_percentage", "saturated_fatty_acids_to_total_fatty_acids_percentage", "linoleic_acid_to_total_fatty_acids_percentage", "docosahexaenoic_acid_to_total_fatty_acids_percentage", "polyunsaturated_fatty_acids_to_monounsaturated_fatty_acids_ratio", "omega_6_fatty_acids_to_omega_3_fatty_acids_ratio"),
#   Branched_chain_amino_acids = c("total_concentration_of_branched_chain_amino_acids_leucine_isoleucine_valine", "isoleucine", "leucine", "valine"),
#   Other_amino_acids = c("alanine", "glutamine", "glycine", "histidine", "phenylalanine", "tyrosine"),
#   Glycolysis_related_metabolites = c("glucose", "lactate", "pyruvate", "citrate"),
#   # Ketone_bodies = c("x3_hydroxybutyrate", "acetate", "acetoacetate", "acetone"),
#   Fluid_balance = c("creatinine", "albumin"),
#   # Inflammation = c("glycoprotein_acetyls"),
#   Chylomicrons_and_very_large_VLDL = c("concentration_of_chylomicrons_and_extremely_large_vldl_particles", "total_lipids_in_chylomicrons_and_extremely_large_vldl", "phospholipids_in_chylomicrons_and_extremely_large_vldl", "cholesterol_in_chylomicrons_and_extremely_large_vldl", "cholesteryl_esters_in_chylomicrons_and_extremely_large_vldl", "free_cholesterol_in_chylomicrons_and_extremely_large_vldl", "triglycerides_in_chylomicrons_and_extremely_large_vldl"),
#   Very_large_VLDL = c("concentration_of_very_large_vldl_particles", "total_lipids_in_very_large_vldl", "phospholipids_in_very_large_vldl", "Cholesterol in Very Large VLDL", "cholesterol_in_very_large_vldl", "free_cholesterol_in_very_large_vldl", "triglycerides_in_very_large_vldl"),
#   Large_VLDL = c("concentration_of_large_vldl_particles", "total_lipids_in_large_vldl", "phospholipids_in_large_vldl", "cholesterol_in_large_vldl", "cholesteryl_esters_in_large_vldl", "free_cholesterol_in_large_vldl", "triglycerides_in_large_vldl"),
#   Medium_VLDL = c("concentration_of_medium_vldl_particles", "total_lipids_in_medium_vldl", "phospholipids_in_medium_vldl", "cholesterol_in_medium_vldl", "cholesteryl_esters_in_medium_vldl", "free_cholesterol_in_medium_vldl", "triglycerides_in_medium_vldl"),
#   Small_VLDL = c("concentration_of_small_vldl_particles", "total_lipids_in_small_vldl", "phospholipids_in_small_vldl", "cholesterol_in_small_vldl", "cholesteryl_esters_in_small_vldl", "free_cholesterol_in_small_vldl", "triglycerides_in_small_vldl"),
#   Very_small_VLDL = c("concentration_of_very_small_vldl_particles", "total_lipids_in_very_small_vldl", "phospholipids_in_very_small_vldl", "cholesterol_in_very_small_vldl", "cholesteryl_esters_in_very_small_vldl", "free_cholesterol_in_very_small_vldl", "triglycerides_in_very_small_vldl"),
#   IDL = c("concentration_of_idl_particles", "total_lipids_in_idl", "phospholipids_in_idl", "cholesterol_in_idl", "cholesteryl_esters_in_idl", "free_cholesterol_in_idl", "triglycerides_in_idl"),
#   Large_LDL = c("concentration_of_large_ldl_particles", "total_lipids_in_large_ldl", "phospholipids_in_large_ldl", "cholesterol_in_large_ldl", "cholesteryl_esters_in_large_ldl", "free_cholesterol_in_large_ldl", "triglycerides_in_large_ldl"),
#   Medium_LDL = c("concentration_of_medium_ldl_particles", "total_lipids_in_medium_ldl", "phospholipids_in_medium_ldl", "cholesterol_in_medium_ldl", "cholesteryl_esters_in_medium_ldl", "free_cholesterol_in_medium_ldl", "triglycerides_in_medium_ldl"),
#   Small_LDL = c("concentration_of_small_ldl_particles", "total_lipids_in_small_ldl", "phospholipids_in_small_ldl", "cholesterol_in_small_ldl", "cholesteryl_esters_in_small_ldl", "free_cholesterol_in_small_ldl", "triglycerides_in_small_ldl"),
#   Very_large_HDL = c("concentration_of_very_large_hdl_particles", "total_lipids_in_very_large_hdl", "phospholipids_in_very_large_hdl", "cholesterol_in_very_large_hdl", "cholesteryl_esters_in_very_large_hdl", "free_cholesterol_in_very_large_hdl", "triglycerides_in_very_large_hdl"),
#   Large_HDL = c("concentration_of_large_ldl_particles", "total_lipids_in_large_hdl", "phospholipids_in_large_hdl", "cholesterol_in_large_hdl", "cholesteryl_esters_in_large_hdl", "free_cholesterol_in_large_hdl", "triglycerides_in_large_hdl"),
#   Medium_HDL = c("concentration_of_medium_hdl_particles", "total_lipids_in_medium_hdl", "phospholipids_in_medium_hdl", "cholesterol_in_medium_hdl", "cholesteryl_esters_in_medium_hdl", "free_cholesterol_in_medium_hdl", "triglycerides_in_medium_hdl"),
#   Small_HDL = c("concentration_of_small_hdl_particles", "total_lipids_in_small_hdl", "phospholipids_in_small_hdl", "cholesterol_in_small_hdl", "cholesteryl_esters_in_small_hdl", "free_cholesterol_in_small_hdl", "triglycerides_in_small_hdl")
# )


fgsea_results = fgseaMultilevel(
  pathways = met_classes,
  stats = ranking,
  minSize = 1,
  maxSize = 23
)

fgsea_results = fgsea_results[order(fgsea_results$padj), ]



fgsea_results$leadingEdge = sapply(fgsea_results$leadingEdge, function(x) paste(x, collapse = ";"))


# Reemplazar "_" con espacio en los nombres de las clases
fgsea_results <- fgsea_results %>%
  mutate(pathway = gsub("_", " ", pathway))  # Sustituir "_" por espacios

# Crear una columna para marcar significancia con un asterisco
fgsea_results <- fgsea_results %>%
  mutate(Significance = ifelse(padj < 0.05, "*", ""))

# Crear el gráfico
library(ggplot2)
# msea_plot <-ggplot(fgsea_results, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
#   geom_bar(stat = "identity") +
#   geom_text(
#     aes(
#       label = Significance,
#       y = ifelse(NES > 0, NES + 0.2, NES - 0.2)  # Posición ajustada según NES
#     ),
#     size = 4, color = "black"
#   ) +
#   coord_flip() +
#   labs(
#     title = "Metabolite Set Enrichment Analysis",
#     x = "",
#     y = "Normalized Enrichment Score (NES)"
#   ) +
#   scale_fill_gradient2(
#     low = "#A7C7E7", mid = "white", high = "#F7A1A1", 
#     midpoint = 0, space = "Lab", limits = c(-2.5, 2.5), oob = scales::squish
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),  # Centrar el título
#     axis.title.y = element_text(size = 12),
#     axis.title.x = element_text(size = 12),
#     axis.text.y = element_text(size = 10),
#     axis.text.x = element_text(size = 10),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 9)
#   )
# 
# msea_plot

fgsea_results %>%  
  dplyr::mutate(direction = if_else(NES > 0, "Up-Regulated", "Down-Regulated")) %>% 
  ggplot(aes(y = reorder(pathway, NES), x = NES, color = pval, size = size)) + 
  facet_grid(. ~ direction, scales = "free_x", space = "free_x") + 
  scale_color_gradientn(name = "P-value", colors = c("red", "purple", "blue")) +
  scale_size(name = "Count", range = c(4, 10)) + 
  geom_point() +
  labs(y = "", x = "Enrichment score") + 
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "cm"), 
    strip.text = element_text(size = 20),     
    axis.text.y = element_text(size = 15), 
    axis.text.x = element_text(size = 15), 
    axis.title.x = element_text(size = 16)
  )
