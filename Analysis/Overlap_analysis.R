library(tidyr)
library(dplyr)

depmap = read.csv("InputData/Benchmarks/processed/ground_truth_depmap_hit_processed.csv")
koferle = read.csv("InputData/Benchmarks/processed/ground_truth_koferle_processed.csv") %>% 
  distinct() 
koferle <- koferle %>%
  group_by(sorted_gene_pair) %>%
  mutate(ground_truth = ifelse(any(ground_truth == 1), 1, ground_truth)) %>%
  distinct() %>%
  ungroup()

common = depmap %>% 
  rename(DepMap = ground_truth) %>%
  inner_join(koferle %>% rename(Koferle = ground_truth), by = c("sorted_gene_pair" = "sorted_gene_pair"))  %>% 
  distinct(sorted_gene_pair, .keep_all = TRUE) 

common$DepMap <- factor(common$DepMap, levels = c(0, 1), labels = c("Depmap(Not-SL)", "Depmap(SL)"))
common$Koferle <- factor(common$Koferle, levels = c(0, 1), labels = c("Koferle(Not-SL)", "Koferle(SL"))

contingency_table = table(common$DepMap, common$Koferle)
fisher.test(contingency_table)
