library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(stringr)
library(RColorBrewer)

depmap = read.csv("InputData/Benchmarks/processed/ground_truth_depmap_hit_processed.csv")
koferle = read.csv("InputData/Benchmarks/processed/ground_truth_koferle_processed.csv") %>% 
  distinct() 
koferle <- koferle %>%
  group_by(sorted_gene_pair) %>%
  mutate(ground_truth = ifelse(any(ground_truth == 1), 1, ground_truth)) %>%
  distinct() %>%
  ungroup()

base_size = 8
text_size = 2#.25

pairs = rbind(read.csv("zdLFC Scripts\\zdLFC Output\\Parrish_PC9.csv") %>% select(X) %>% mutate(study = "Parrish"),
              read.csv("zdLFC Scripts\\zdLFC Output\\DeDe_zdLFC.csv")%>% select(X)%>% mutate(study = "Dede"),
              read.csv("zdLFC Scripts\\zdLFC Output\\Thompson_zdLFC.csv")%>% select(X)%>% mutate(study = "Thompson"),
              read.csv("zdLFC Scripts\\zdLFC Output\\ChymeraHAP1.csv")%>% select(X)%>% mutate(study = "CHyMErA"),
              read.csv("zdLFC Scripts\\zdLFC Output\\ITO.csv")%>% select(X)%>% mutate(study = "Ito"))


desired_order<- c("Not common with any benchmark",
                  "Common with Köferle benchmark only",
                  "Common with De Kegel bechmark only",
                  "Common with both benchmarks")


labels <- c("depmap" = "(b) De Kegel Benchmark", "koferle" = "(c) Köferle Benchmark")

Fig2bc = pairs %>% group_by(study) %>%
  left_join(depmap, by = c( "X" = "sorted_gene_pair")) %>%
  select(X, ground_truth, study) %>% 
  left_join(koferle, by = c("X" =  "sorted_gene_pair")) %>%
  rename(sorted_gene_pair = X,
         depmap = ground_truth.x,
         koferle = ground_truth.y) %>%
  pivot_longer(cols = c(depmap, koferle), names_to = "source", values_to = "label") %>%
  count(study, source, label) %>% 
  drop_na() %>%
  group_by(source, study) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label_text = paste0(round(prop * 100,1), "%\n(", n, ")")   ) %>%


  ggplot(
    aes(x = study, y = n, 
        fill = factor(label, levels = c(0, 1), labels = c("not-SL", "SL"))
    )) +
  geom_col(position = "fill") +  # Fill makes it proportional
  geom_text(
    aes(
      label = label_text
    ),
    position = position_fill(vjust = 0.5),
    hjust = 0.5,
    size = text_size, 
    color = "#404040"
  )+
  #facet_wrap(~source) +
  facet_grid(~source, labeller = labeller(source = labels)) +
  
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = "Proportion",
    x = "Study",
    fill = "Label (ground truth)"
  ) +
  theme_minimal(base_size = base_size) +
  theme(
    strip.text = element_text(size = base_size),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = base_size-2),        # shrink legend text
    legend.title = element_blank(),#element_text(size = base_size-2),       # optional, for the legend title
    legend.key.size = unit(0.3, "cm"),#    
    plot.subtitle = element_text(
    hjust = 0.5,        # Horizontal justification: 0 = left, 0.5 = center, 1 = right
    margin = margin(b = 10)  # Add space below subtitle
  ))+
  scale_fill_manual(
    values = c("not-SL" = "#E0E0E0", "SL" = "#4393C3"),
    breaks = c("SL", "not-SL")  # ✅ this controls the legend order
  )+
  
  ylab("Percentage of SL pairs") + xlab("") +
  coord_flip() 

common = pairs %>% group_by(study) %>%
  left_join(depmap, by = c( "X" = "sorted_gene_pair")) %>%
  select(X, ground_truth, study) %>% 
  left_join(koferle, by = c("X" =  "sorted_gene_pair")) %>%
  rename(sorted_gene_pair = X,
         depmap = ground_truth.x,
         koferle = ground_truth.y) %>%
  distinct() %>%
  mutate(
    category = case_when(
      !is.na(depmap) & is.na(koferle) ~ desired_order[3],
      !is.na(koferle) & is.na(depmap) ~ desired_order[2],
      !is.na(koferle) & !is.na(depmap) ~ desired_order[4],
      TRUE ~ desired_order[1])
  ) %>%
  distinct() %>%
  mutate(category  =  factor(category, levels = desired_order)) %>%
  mutate(total = n()) %>%
  count(study, category, total, .drop = FALSE) %>%  # keep total in output    
  mutate(
    percent = n / sum(n) * 100,
    label = paste0(round(percent,1),"%", "\n(",n , ")"),
    
  )

# Create a single stac


Fig2a = common %>% ungroup()%>%
  
  ggplot(aes(x = study, y = n, fill = category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    aes(
      label = label,
      color =  ifelse(category == "Not common with any benchmark", "#404040", "white")),
    position = position_fill(vjust = 0.5),
    hjust = 0.5,
    size = text_size) +
  
  labs( x = NULL,y = NULL,fill = NULL,subtitle = "(a) Overlap with benchmarks"
  ) +
  scale_color_identity()+
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)) +  # Y-axis as percentage
  theme_minimal(base_size = base_size) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = base_size-2),        # shrink legend text
    legend.title = element_blank(),#element_text(size = base_size-2),       # optional, for the legend title
    legend.key.size = unit(0.3, "cm")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = base_size-2),        # shrink legend text
    legend.title = element_text(size = base_size-2),       # optional, for the legend title
    legend.key.size = unit(0.3, "cm"),  # or smaller if needed
    plot.subtitle = element_text(
      size = base_size,
     hjust = 0.5,        # Horizontal justification: 0 = left, 0.5 = center, 1 = right
      #margin = margin(b = 10)  # Add space below subtitle
    )
    
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Common with both benchmarks" = "#584f67",
      "Common with De Kegel bechmark only" = "#26677f",
      "Common with Köferle benchmark only" = "#89374f",
      "Not common with any benchmark" = "#E0E0E0"
    )
  ) + 
  guides(fill = guide_legend(reverse = TRUE))  # Optional, in case needed for consistency


combined = Fig2a/Fig2bc
combined
ggsave("Analysis/Figures/plots/Figure-2.png",plot =  combined, dpi = 1000, width = 20, height = 14, units = "cm",bg = "white")
