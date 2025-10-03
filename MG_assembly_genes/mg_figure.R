# plotting of figure with gene counts identified in MGs assemblies
# Requires to run code for MG assembly and gene annotation

# Load required libraries
library(tidyverse)
library(tidylog)
library(readr)
library(readxl)
library(openxlsx2)
library(ggbeeswarm)

# Prepare individual data 
individual_data <- select(all_trace_genes_renamed_diamond, X1, X2, X11, X12) |> # all_trace_genes_renamed_diamond is an output of external script
  tidyr::separate_wider_delim(X1, delim = "|", names = c("qseqid", "treatment")) |> 
  tidyr::separate_wider_delim(X2, delim = "|", names = c("accession_helper", "gene_id")) |>
  left_join(tax_merged, by = "accession_helper") |> 
  right_join(distinct_full, "qseqid") |> 
  tidyr::separate_wider_delim(X2, delim = ";", names = c("domain", "phylum", "class","order", "family", "genus", "species")) |>
  filter(function_diamond == "Carbon monoxide oxidation" | 
           function_diamond == "Methane oxidation" | 
           function_diamond == "hydrogenase") |>
  group_by(gene, class, treatment) |>  
  summarise(n = n(), .groups = "drop") |> 
  left_join(readxl::read_excel("gene_numbers.xlsx"), by = "treatment") |>
  mutate(perc = (n/gene_count)*1e6) |> # genes per mio of genes
  filter(!is.na(class)) |> # remove no hits from diamond tax
  left_join(readxl::read_excel("distinct_full_taxonomy_helper.xlsx"),
            by = c("gene", "class")) |> # join helper table to plot only most abundant
  mutate(class_helper = if_else(is.na(class_helper), true = "Others", false = class_helper)) |>
  mutate(class_helper = str_remove(class_helper, "c__")) |> 
  mutate(temp = fct_recode(temp, `Long−term warmed` = "warm")) |>
  mutate(temp = fct_recode(temp, Unwarmed = "ambient")) %>% 
  add_row(individual_data_2016) %>% # add individual data for 2016 MGs calculated in an external script, important row, adds two more samples to a plot
  mutate(gene = factor(gene, levels = c("NiFe", "CoxL", "PmoA")))

# Create summary for means and error bars
summary_data <- individual_data |>
  group_by(gene, class_helper, temp) |> 
  summarise(
    perc_mean = mean(perc), 
    se = sd(perc)/sqrt(n()), 
    n_points = n(),
    .groups = "drop"
  ) |>
  group_by(class_helper) |> 
  mutate(sort_help = mean(perc_mean)) |>
  ungroup() |>
  mutate(wd = if_else(gene == "CoxL", 1, 
                      if_else(gene == "NiFe", 1, 0.3))) %>% 
  mutate(gene = factor(gene, levels = c("NiFe", "CoxL", "PmoA")))

# Add sorting to individual data
individual_data_sorted <- individual_data |>
  left_join(summary_data |> select(class_helper, sort_help) |> distinct(),
            by = "class_helper") %>% 
  mutate(point_shape = ifelse(treatment == "A6" | treatment == "W6", "triangle", "circle"))

# Plot
plot_individual_points <- ggplot() +
  # Individual points layer with beeswarm to avoid overlapping
  geom_quasirandom(
    data = individual_data_sorted,
    aes(x = fct_reorder(class_helper, sort_help, .desc = T), 
        y = perc, 
        colour = temp,
        shape = point_shape),
    dodge.width = 0.8,
    size = 1.5,
    alpha = 1
  ) +
  # Mean points layer
  geom_point(
    data = summary_data,
    aes(x = fct_reorder(class_helper, sort_help, .desc = T),
        y = perc_mean,
        group = temp),
    position = position_dodge(width = 0.8),
    size = 4,
    shape = 95,
    colour = "black"
  ) +
  # Error bars layer
  geom_errorbar(
    data = summary_data,
    aes(x = fct_reorder(class_helper, sort_help, .desc = T), 
        ymin = perc_mean - se, 
        ymax = perc_mean + se,
        group = temp),
    width = 0.5,
    position = position_dodge(width = 0.8),
    linewidth = 1,
    colour = "black"
  ) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c("circle" = 16, "triangle" = 17)) +
  guides(shape = "none") +  # Hide shape legend
  guides(colour = "none") +  # Hide temp legend
  facet_wrap(~ gene, scales = "free") +
  labs(y = "Genes per 1M assembled genes") +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    legend.title = element_blank(),
    strip.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0, face = "italic", size = 9, color = "black"),
    axis.title.x = element_blank()
  )
plot_individual_points
# ggsave this

# export supplementary table
export_summary_data <- summary_data %>% 
  select(-sort_help, -wd) %>% 
  rename(
    `Gene`      = gene,
    `Taxon (class)` = class_helper,
    `Soil regime` = temp,
    `Genes per 1M assembled genes`    = perc_mean,
    `Std. Error`       = se,
    `N`                = n_points
  ) %>% 
  arrange(`Gene`, `Taxon (class)`, `Soil regime`, desc(`Genes per 1M assembled genes`)) %>% 
  mutate(
    `Genes per 1M assembled genes` = sprintf("%.2f%%", `Genes per 1M assembled genes`),
    `Std. Error`    = sprintf("%.3f", `Std. Error`),
    `N`             = formatC(`N`, format = "d", big.mark = ",")
  ) %>% 
  glimpse
colnames(export_summary_data) <- c(
  "Gene",
  "Taxon_Class",
  "Soil_Regime",
  "Genes_per_1M_assembled",
  "Std_Error",
  "N_Observations"
)
# openxlsx2::write_xlsx(export_summary_data, file = "Supplementary_Table.xlsx", overwrite = TRUE)

