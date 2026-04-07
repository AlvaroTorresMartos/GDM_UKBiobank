
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  ggplot2, ggtext, patchwork)

results = readr::read_rds("/mnt/project/2026_01_12_cox_models.RDS") %>% 
  dplyr::mutate(omics = c(rep("metabolomics", 6), rep("proteomics", 6)), 
                HR = case_when(term == "quantile_met1" ~ "Reference", 
                               TRUE ~ sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high)), 
                term = case_when(term == "quantile_met1" ~ "Quintile 1", 
                                 term == "quantile_met2" ~ "Quintile 2", 
                                 term == "quantile_met3" ~ "Quintile 3", 
                                 term == "quantile_met4" ~ "Quintile 4", 
                                 term == "quantile_met5" ~ "Quintile 5", 
                                 term == "z_scaled" ~ "Continuous score", 
                                 )
                ) 

# metabolomics ------


forest_prote = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous( breaks = c(0, 1, 2, 3, 4, 5)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
    axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Metabolomics signature**") 

table1



table2 = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_t2dm, perc_t2dm)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::filter(omics == "metabolomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value < 0.05 ~ "< 0.05", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)





# proteomics -----


forest_prote = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(x = estimate, y = term)) +
  geom_point(size = 4, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0, color = "black") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  # scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5)) +
  theme_minimal(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  
    axis.text.y = element_text(size = 0, hjust = 0, color = "black"), 
        axis.text.x = element_text(color = "black"))  + 
  labs(x = "", y = "")

forest_prote





table1 = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = paste0("**", term, "**"), 
                term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = term)) +
  ggtext::geom_richtext(hjust = 0, fill = NA, label.color = NA, text.color = "black") +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown()
  )  +
  labs(title = "**Proteomics signature**") 

table1



table2 = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = count_total)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**N**") 

table2

table3 = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                count_percentage = sprintf("%.0f (%.2f%%)", count_t2dm, perc_t2dm)) %>% 
  ggplot(aes(y = term, x = 1, label = count_percentage)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 5),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**Events**") 

table3

table4 = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term)) %>% 
  ggplot(aes(y = term, x = 1, label = HR)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +  # Ancho para alinear
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**HR (95% CI)**") 

table4

table5 = results %>% 
  dplyr::filter(omics == "proteomics") %>% 
  dplyr::mutate(term = forcats::fct_rev(term), 
                p.value = as.numeric(p.value), 
                p.value = case_when(is.na(p.value)  ~ "",
                                    p.value < 0.001 ~ "< 0.001", 
                                    p.value < 0.05 ~ "< 0.05", 
                                    p.value > 0.001 ~ sprintf("%.2f", p.value)
                )) %>% 
  ggplot(aes(y = term, x = 1, label = p.value)) +
  geom_text(hjust = 0) +
  xlim(1, 1.5) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(5, 5, 5, 20),  # Ajusta para alineación
    axis.text.y = element_blank(), 
    plot.title = ggtext::element_markdown())  +
  labs(title = "**P-value**") 

table5





# Combinar horizontalmente
final_plot <- table1 + table2 + table3 + forest_prote + table4 + table5 + plot_layout(widths = c(2, 0.5, 1, 4, 1.3, 1))
print(final_plot)



