


# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(dplyr, readr, ggplot2, ggrepel)

# Impot the results ------

# metab = read_csv("/mnt/project/results/1_limma_analysis/2025_04_22_metabolites_associated_GDM.csv")
prote = read_csv("/mnt/project/results/1_limma_analysis/2025_12_17_proteins_associated_GDM.csv")

# metab$logFC %>% summary()
# metab$adj.P.Val %>% summary()
# 
# metab$logFC %>% hist()
# metab$adj.P.Val %>% hist()

cols = c("Increased Levels" = "#FB9A99",
         "Decreased Levels" = "#A6CEE3",
         "ns" = "grey")
sizes = c("Increased Levels" = 3,
          "Decreased Levels" = 3,
          "ns" = 1.5)
alphas = c("Increased Levels" = 1,
           "Decreased Levels" = 1,
           "ns" = 0.15)


# metab %>%  
#   dplyr::mutate(direction = ifelse(logFC > 0 & adj.P.Val < 0.05 & abs(logFC) > 0.3, "Increased Levels", 
#                             ifelse(logFC < 0  & adj.P.Val < 0.05 & abs(logFC) > 0.3, "Decreased Levels", 
#                                           "ns")), 
#                 predictor_label = ifelse(direction != "ns", Metabolite, NA_character_)
#   ) %>% 
#   ggplot(aes(x = logFC, y = -log10(adj.P.Val), 
#              fill = direction, 
#              size = direction, 
#              alpha = direction, 
#              label = predictor_label)) + 
#   geom_point(shape = 21, colour = "black") + 
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.7, color = "darkgrey") + 
#   geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", size = 0.7, color = "darkgrey") +
#   scale_fill_manual(values = cols) + # point colour
#   scale_size_manual(values = sizes) + # point size
#   scale_alpha_manual(values = alphas) + # transparency 
#   geom_label_repel(max.overlaps = Inf, size = 3.5, box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50') +
#   theme_minimal(base_size = 14) +  # Tamaño base de 14 para mejor legibilidad
#   theme(
#     plot.caption = element_text(size = 10),
#     panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey90"),  # Líneas de cuadrícula principales suaves
#     panel.grid.minor = element_blank(),  # Elimina líneas de cuadrícula menores
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Borde negro para el panel
#     axis.title = element_text(face = "bold"),  # Ejes en negrita para destacar títulos
#     axis.text = element_text(color = "black"),  # Texto del eje en negro para mejor contraste
#     legend.position = "none",  # Elimina la leyenda
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Título centrado, negrita y tamaño grande
#     plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12)  # Subtítulo centrado y en negrita para equilibrio visual
#   ) +
#   labs(title = "Metabolites Associated with GDM", 
#        # subtitle = "Fórmula: HOMA-IR z-score ~ Metabolito + Hospital + Punto Temporal + (1 | Id) (N = 278)",
#        x = "log(Fold Change)",
#        y = "-log(FDR)", 
#        fill = "Direction")  #+ 
#   # scale_y_continuous(limits = c(0, 2.5),
#   #                    breaks = seq(0, 2.5, by = 0.5)
#   # ) + 
#   # scale_x_continuous(
#   #   limits = c(-1.5, 1.5)
#   # )

prote %>%  
  dplyr::mutate(direction = ifelse(logFC > 0 & adj.P.Val < 0.05 & abs(logFC) > 0, "Increased Levels", 
                                   ifelse(logFC < 0 & adj.P.Val < 0.05 & abs(logFC) > 0.0, "Decreased Levels", 
                                          "ns")), 
                predictor_label = ifelse(direction != "ns", Protein, NA_character_)
  ) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), 
             fill = direction, 
             size = direction, 
             alpha = direction, 
             label = predictor_label)) + 
  geom_point(shape = 21, colour = "black") + 
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.7, color = "darkgrey") + 
  # geom_vline(xintercept = c(-0.40, 0.40), linetype = "dashed", size = 0.7, color = "darkgrey") + 
  scale_fill_manual(values = cols) + # point colour
  scale_size_manual(values = sizes) + # point size
  scale_alpha_manual(values = alphas) + # transparency 
  geom_label_repel(max.overlaps = Inf, size = 5, box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50') +
  theme_minimal(base_size = 14) +  # Tamaño base de 14 para mejor legibilidad
  theme(
    axis.text.x = element_text(size = 13), 
    axis.text.y = element_text(size = 16),
    # plot.caption = element_text(size = 10),
    panel.grid.major = element_line(size = 0.5, linetype = 'solid', color = "grey90"),  # Líneas de cuadrícula principales suaves
    panel.grid.minor = element_blank(),  # Elimina líneas de cuadrícula menores
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Borde negro para el panel
    axis.title = element_text(face = "bold"),  # Ejes en negrita para destacar títulos
    # axis.text = element_text(size = 15, color = "black"),  # Texto del eje en negro para mejor contraste
    legend.position = "none",  # Elimina la leyenda
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  # Título centrado, negrita y tamaño grande
    plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 12)  # Subtítulo centrado y en negrita para equilibrio visual
  ) +
  labs(
       # title = "Volcano plt", 
       # subtitle = "Fórmula: HOMA-IR z-score ~ Metabolito + Hospital + Punto Temporal + (1 | Id) (N = 278)",
       x = "log(Fold Change)",
       y = "-log(FDR)", 
       fill = "Direction")  #+ 
# scale_y_continuous(limits = c(0, 2.5),
#                    breaks = seq(0, 2.5, by = 0.5)
# ) + 
# scale_x_continuous(
#   limits = c(-1.5, 1.5)
# )

