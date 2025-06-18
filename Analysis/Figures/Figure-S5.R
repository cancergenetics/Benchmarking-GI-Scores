library(readxl)
library(tidyverse)
library(tidylog)
library(RColorBrewer) # for heatmap colors
library(kableExtra) # for formatting kables
library(corrr)
library(tibble)
library(tidyr)
library(dplyr)
library(orthrus)
library(scales)

base_size = 8
text_size = 2.25

custom_labeller <- function(labels) {
  # Generate prefixes (a), (b), etc.
  prefixes <- paste0("(", letters[seq_along(labels)] ,")\n ")
  # Concatenate prefixes with original facet titles
  labels <- paste0(prefixes, labels)
  return(labels)
}

results_filtered = read.csv("Analysis/Output/Filtered/Compiled/Combined_Results.csv")
results_unfiltered = read.csv("Analysis/Output/Compiled/Combined_Results.csv")
# compare AUROC of Gemini on CHYMERA before and after pre-process
filtered = results_filtered %>% 
  filter( Score == "Gemini(Sens)" , Cell.line == "All" ) 

unfiltered = results_unfiltered %>% 
  filter(Score == "Gemini(Sens)", Cell.line == "All")

cor_df <- filtered %>%
  inner_join(unfiltered, by = c("Metric" ,  "Score", "Cell.line", "Validation.Set", "Study.name"),
             suffix = c(".filtered", ".unfiltered")) %>%
  select_if(~ n_distinct(.) > 1) %>%
  filter(Metric == "AUROC") %>%
  
  mutate(difference =  (value.filtered - value.unfiltered)) %>%
  mutate(Validation.Set = case_when(
    Validation.Set == "DepMap Hits" ~ "De Kegel benchmark",
    Validation.Set == "Köferle List" ~ "Köferle benchmark",
    TRUE ~ Validation.Set
  )) %>%
  mutate(direction = ifelse(difference < 0, "Negative", "Positive")) %>%
  mutate(dropped = (Common.samples.unfiltered - Common.samples.filtered) / Common.samples.unfiltered * 100) %>%
  group_by(Metric, Validation.Set) %>%
  summarize(
    PearsonR = cor(difference, dropped, use = "complete.obs", method = "pearson"),
    SpearmanR = cor(difference, dropped, use = "complete.obs", method = "spearman")
  ) %>%  
  mutate(label1 = paste0("Pearson R = ", round(PearsonR, 2)), label2 = paste0("Spearman R = ", round(SpearmanR, 2))) 

diff = filtered %>% inner_join(unfiltered, by = c("Metric" ,  "Score", "Cell.line", "Validation.Set", "Study.name"),
                               suffix = c(".filtered", ".unfiltered")) %>%
  select_if(~ n_distinct(.) > 1) %>%
  mutate(difference =  (value.filtered - value.unfiltered) ) %>%
  mutate(Validation.Set = ifelse(Validation.Set == "DepMap Hits", "De Kegel benchmark", Validation.Set)) %>%
  mutate(Validation.Set = ifelse(Validation.Set == "Köferle List", "Köferle benchmark", Validation.Set)) %>%
  mutate(direction = ifelse(difference < 0, "Negative", "Positive")) %>%
  select(-X.filtered, -X.unfiltered) %>%
  group_by(Study.name, Metric, Validation.Set) %>% 
  mutate(dropped = (Common.samples.unfiltered -Common.samples.filtered)/Common.samples.unfiltered * 100 )%>%
  ungroup() %>%
  filter(Metric == "AUROC") %>%
  
  select(Metric, Validation.Set, Study.name, dropped, difference, direction) %>%
  ggplot() +
  geom_point(aes(x = difference, y = dropped, color = direction), size = 3)+
  geom_text(
    aes(
      x = difference,
      y = dropped,
      label = Study.name ),
    vjust = -0.75, color = "grey40", size = text_size) +
  facet_grid(cols = vars(Validation.Set), scales = "free_x",
             labeller = labeller(Validation.Set = custom_labeller))+
  
  theme_minimal(base_size = base_size) +
  scale_y_continuous(labels = label_percent(scale = 1), limits = c(0, 30)) + 
  scale_color_manual(
    values = c("Negative" = "#D6604D", "Positive" = "#4393C3"), # Blue for positive, red for negative
    name = "")+
  labs(
    x = "Change in Performance (AUROC)",
    y = "Percentage of dropped gene pairs") + 
  theme(
    panel.border = element_rect(color = "darkgrey", fill = NA, size = 0.1), # Adds a border around each facet panel
    legend.position = "bottom",    
    panel.spacing = unit(0.50, "cm")
  )  + 
  coord_cartesian(clip = "off") +  
  geom_text(
    data = cor_df,
    aes(x = -Inf, y = Inf, label = label1),
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE,
    color = "black", size = text_size)
ggsave("Analysis/Figures/plots/Figure-S5.png",plot =  diff, dpi = 600, width = 18,height = 10, units = "cm",bg = "white")
