library(ggplot2)
library(ggpubr)
library(tidyverse)
library(patchwork)
library(maditr)
library(forcats)

setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Bonito Lab/CONSTAX_fixes/CONSTAXv2_ms_code/scripts/classification/")

sil_reg <- read.csv("../../data/classification/silva/part_cv_metrics_silva_reg.csv")
sil_reg_c <- read.csv("../../data/classification/silva/part_cv_metrics_silva_reg_conservative.csv")
sil_mothur <- read.csv("../../data/classification/silva/part_cv_metrics_silva_mothur.csv")

uni_reg_utax <- read.csv("../../data/classification/unite/part_cv_metrics_unite_reg_utax.csv")
uni_reg_utax_c <- read.csv("../../data/classification/unite/part_cv_metrics_unite_reg_utax_conservative.csv")
uni_reg_blast <- read.csv("../../data/classification/unite/part_cv_metrics_unite_reg_blast.csv")
uni_reg_blast_c <- read.csv("../../data/classification/unite/part_cv_metrics_unite_reg_blast_conservative.csv")

uni_mothur <- read.csv("../../data/classification/unite/part_cv_metrics_unite_mothur.csv")
uni_qiime <- read.csv("../../data/classification/unite/part_cv_metrics_unite_qiime.csv")

levels(sil_reg$classifier)[levels(sil_reg$classifier)=="Consensus"] <- "CB"
levels(sil_reg_c$classifier)[levels(sil_reg_c$classifier)=="Consensus"] <- "CBC"

levels(uni_reg_utax$classifier)[levels(uni_reg_utax$classifier)=="Consensus"] <- "CU"
levels(uni_reg_blast$classifier)[levels(uni_reg_blast$classifier)=="Consensus"] <- "CB"
levels(uni_reg_utax_c$classifier)[levels(uni_reg_utax_c$classifier)=="Consensus"] <- "CUC"
levels(uni_reg_blast_c$classifier)[levels(uni_reg_blast_c$classifier)=="Consensus"] <- "CBC"


colnames(uni_mothur)
colnames(uni_reg_blast)

# Extract unique classifiers not in the uni_reg_blast dataframe
uni_reg_utax %>%
  filter(classifier %in% c("UTAX", "CU")) -> uru_fil
uni_reg_utax_c %>%
  filter(classifier == "CUC") -> uruc_fil
uni_reg_blast_c %>%
  filter(classifier == "CBC") %>%
  rbind(., uru_fil, uruc_fil, uni_reg_blast, uni_mothur, uni_qiime) -> comb_unite_df

head(comb_unite_df)
comb_unite_df %>% # Extract metrics
  pivot_longer(cols = sensitivity:EPQ,
             names_to = "Metric",
             values_to = "value") -> uni_reg_long
sil_reg_c %>% # Only difference here is the Consensus BLast classifier
  filter(classifier == "CBC") %>%
  rbind(., sil_reg, sil_mothur) -> comb_sil_df

comb_sil_df %>% # Extract metrics
  pivot_longer(cols = sensitivity:EPQ,
               names_to = "Metric",
               values_to = "value") -> sil_reg_long

metric.labs <- c("Errors per Query", "Misclassification", "Overclassification", "Sensitivity")
names(metric.labs) <- c("EPQ", "MC", "OC", "sensitivity")
lev.labs <- c("Family", "Genus")
names(lev.labs) <- c("fam", "gen")

sil_reg_long %>%
  ggplot(aes(x=classifier, y=value, color=region)) +
  geom_point(position = position_dodge(width=0.6), alpha = 0.5) +
  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("BLAST","RDP","SINTAX", "CB", "CBC", "mothur-wang", "mothur-knn")) +
  scale_color_viridis_d(end = 0.8) +
  labs(x=NULL, y="Metric Value", color = "Region", title = "Bacteria") -> p_sil
p_sil
ggsave("../../figures/region_classification_part_cv_silva.png", p_sil, width = 10, height = 8, units = "in", dpi = 400)

uni_reg_long %>%
  ggplot(aes(x=classifier, y=value, color=region)) +
  geom_point(position = position_dodge(width=0.6), alpha = 0.5) +
  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("BLAST","RDP","SINTAX", "UTAX", "CB", "CBC", "CU", "CUC", "mothur-wang", "mothur-knn", "qiime2-Naive-Bayes")) +
  scale_color_viridis_d(end = 0.8) +
  labs(x="Classifier", y="Metric Value", color = "Region", title = "Fungi") -> p_uni
p_uni
ggsave("../../figures/region_classification_part_cv_unite.png", p_uni, width = 10, height = 8, units = "in", dpi = 400)

p_sil / p_uni -> comb_plot

### Speed tests
# Training speed, 1 thread

n_vary_train <- read.csv("../../data/speed_tests/speed_training_t1_nvary.csv")

n_vary_train <- n_vary_train %>%
  filter(step == "train") %>% # extract the times associated with the training step
  mutate(seq_per_time=seq_count/time) # calculate per sequence time

sum_n_vary_t <- n_vary_train %>%
  group_by(algorithm, seq_count) %>% # Calculate mean and standard deviation of train time by algorithm and reference sequence count
  summarise(mean_seq_per_time = mean(seq_per_time),
            sd_seq_per_time = sd(seq_per_time))

ggplot(n_vary_train, aes(x = seq_count, y = seq_per_time, color = algorithm)) +
  geom_point()

g_t1 <- ggplot(sum_n_vary_t, aes(x = seq_count, y = mean_seq_per_time, color = algorithm)) +
  geom_line(lwd=1.5) +
  geom_point(data = n_vary_train,
             aes(x = seq_count, y = seq_per_time, color = algorithm),
             alpha=0.4) +
  geom_errorbar(aes(ymin=mean_seq_per_time-sd_seq_per_time,
                    ymax=mean_seq_per_time+sd_seq_per_time),
                width=200, linetype=1) +
  labs(x = "Sequence Count",
       y = "Sequences trained per second",
       color = "Algorithm") +
  scale_color_discrete(breaks=c("blast", "utax"), labels=c("BLAST", "UTAX")) +
  theme(legend.position = "none")
g_t1

# Classification speed, 1000, 2000, 4000 query set sizes

classify_dat <- read.csv("../../data/speed_tests/speed_classify_comb.csv")
classify_dat <- classify_dat %>%
  mutate(time_per_seq=seq_count/time)
head(classify_dat)

classify_dat %>%
  filter(step == "classify") -> cl_df
cl_df %>%
  group_by(threads, algorithm, seq_count) %>%
  summarise(mean_time=mean(time_per_seq), sd_time = sd(time_per_seq)) -> sum_cl_df

data.frame(sum_cl_df)

class_lp <- ggplot(sum_cl_df, aes(x = threads, y = mean_time, color = algorithm, linetype = as.factor(seq_count))) +
  geom_line(lwd=1.5) +
  geom_point() +
  labs(x = "Threads",
       y = "Sequences classified per second",
       color = "Algorithm",
       linetype = "Query sequence\ncount") +
  geom_errorbar(aes(ymin=mean_time-sd_time, ymax=mean_time+sd_time), width=2, linetype=1) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  scale_color_discrete(breaks=c("blast", "utax"), labels=c("BLAST", "UTAX"))
class_lp

g <- g_t1 + class_lp #+ plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
g
ggsave("../../figures/train_class_time_lp.png", g, width = 6, height = 3, units="in", dpi = 400)

(g_t1 + class_lp) / p_sil / p_uni + plot_annotation(tag_levels = 'A') -> comb_plot
comb_plot


ggsave("../../figures/region_classification_part_v_sil.png", p_sil + labs(x="Classifier"), width = 7, height = 4, units = "in", dpi = 400)
ggsave("../../figures/region_classification_part_v_uni.png", p_uni, width = 7, height = 4, units = "in", dpi = 400)
ggsave("../../figures/speed_and_region_classification_part_cv.png", comb_plot, width = 8.5, height = 12, units = "in", dpi = 400)

comb_reg_df <- rbind(uni_reg_long, sil_reg_long)

comb_reg_df %>%
  group_by(database, partition_level, region, classifier, Metric) %>%
  summarise(mean=mean(value), sd = sd(value)) %>%
  mutate(tbl_entry = paste(round(mean*100, 1), round(sd*100, 1), sep = "�")) ->
  res_tbl

res_tbl$classifier <- fct_relevel(res_tbl$classifier, c("BLAST", "RDP", "SINTAX", "UTAX", "CB", "CBC", "CU", "CUC"))
  
res_tbl %>% dcast(partition_level + classifier ~ database + region + Metric,
                  value.var = "tbl_entry") -> out_tbl
out_tbl$classifier <- fct_relevel(out_tbl$classifier, c("BLAST", "RDP", "SINTAX", "UTAX", "CB", "CBC", "CU", "CUC"))
out_tbl %>%
  arrange(partition_level, classifier) -> out_tbl
out_tbl
write.csv(out_tbl, "../../tables/region_partition_summary.csv")

# Parameter optimization
### UNITE
unite_params <- read.csv("../../data/classification/unite/part_cv_metrics_unite_params.csv")

unite_params %>%
  filter(classifier=="Consensus") %>%
  mutate(max_hits = as.character(max_hits)) %>%
  pivot_longer(cols = sensitivity:EPQ,
               names_to = "Metric",
               values_to = "value") -> uni_param_long

metric.labs <- c("Errors per Query", "Misclassification", "Overclassification", "Sensitivity")
names(metric.labs) <- c("EPQ", "MC", "OC", "sensitivity")
lev.labs <- c("Family", "Genus")
names(lev.labs) <- c("fam", "gen")

uni_param_long %>%
  filter(conf==0.8) %>%
ggplot(aes(x=max_hits, y=value)) +
  scale_x_discrete(limits=c("1", "3", "5", "10", "20")) +
  geom_point(position = position_dodge(width=0.6), alpha = 0.5) +

  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +

  scale_color_viridis_d(end = 0.8) +
  labs(x="Max Hits", y="Metric Value") -> p_uni_param_mh
p_uni_param_mh
ggsave("../../figures/unite_mh.png", p_uni_param_mh + labs(title="Fungi"), width = 10, height = 4, units = "in", dpi = 400)

uni_param_long %>%
  filter(max_hits%in%c(5, 20)) %>%
  ggplot(aes(x=conf, y=value, color=max_hits)) +
  geom_point(alpha = 0.5) +
  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +

  labs(x="Confidence Threshold", y="Metric Value",color = "Max Hits") -> p_uni_param_cf
p_uni_param_cf
ggsave("../../figures/unite_cf.png", p_uni_param_cf + labs(title="Fungi"), width = 10, height = 4, units = "in", dpi = 400)

g <- (p_uni_param_mh + labs(tag = "A")) / (p_uni_param_cf + labs(tag = "B"))

ggsave("../../figures/unite_param_comb.png", g, width = 10, height = 8, units = "in", dpi = 400)

### SILVA parameter optimization
silva_params <- read.csv("../../data/classification/silva/part_cv_metrics_silva_params.csv")

silva_params %>%
  filter(classifier=="Consensus") %>%
  mutate(max_hits = as.character(max_hits)) %>%
  pivot_longer(cols = sensitivity:EPQ,
               names_to = "Metric",
               values_to = "value") -> sil_param_long

sil_param_long %>%
  filter(conf==0.8) %>%
  ggplot(aes(x=max_hits, y=value)) +
  scale_x_discrete(limits=c("1", "3", "5", "10", "20")) +
  geom_point(position = position_dodge(width=0.6), alpha = 0.5) +

  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="Max Hits", y="Metric Value") -> p_sil_param_mh
p_sil_param_mh
ggsave("../../figures/silva_mh.png", p_sil_param_mh + labs(title = "Bacteria"), width = 10, height = 4, units = "in", dpi = 400)

sil_param_long %>%
  filter(max_hits %in% c(5, 20)) %>%
  ggplot(aes(x=conf, y=value, color=max_hits)) +
  geom_point(alpha = 0.5) +
  facet_grid(partition_level~Metric,
             labeller = labeller(Metric=metric.labs, partition_level=lev.labs)) +
  theme(axis.text.x.bottom = element_text(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0.5)) +

  labs(x="Confidence Threshold", y="Metric Value", color="Max Hits") -> p_sil_param_cf
p_sil_param_cf
ggsave("../../figures/silva_cf.png", p_sil_param_cf + labs(title = "Bacteria"), width = 10, height = 4, units = "in", dpi = 400)


g <- (p_sil_param_mh + labs(tag = "A")) / (p_sil_param_cf + labs(tag = "B"))

ggsave("../../figures/silva_param_comb.png", g, width = 10, height = 8, units = "in", dpi = 400)

# Comparison of classification counts

comb_tax_sb <- read.delim("../../data/classification/silva/combined_taxonomy_sil_blast.txt")

comb_tax_ub <- read.delim("../../data/classification/unite/combined_taxonomy_uni_blast.txt")

comb_tax_sb %>%
  sapply(function(x) sum(x != "")) %>%
  as.data.frame() %>%
  `colnames<-`(c("Value")) %>%
  mutate(Level = rownames(.)) %>%
  .[2:(dim(.)[1]-4),] %>%
  replace_na(list(Value=0)) %>%
  mutate(Rank = unlist(strsplit(Level, "_"))[c(F, T, F)],
         Classifier = unlist(strsplit(Level, "_"))[c(F, F, T)],
         Value = Value/dim(comb_tax_sb)[1]) %>%
  mutate(Classifier = factor(Classifier, levels = c("BLAST", "RDP", "SINTAX", "Consensus"))) %>%
  ggplot(aes(x = Rank, y = Value, fill=Classifier)) +
  geom_bar(position="dodge", stat="identity") +
  labs(y = "Proportion of OTUs classified",
       title = "16S - SILVA") -> silva_class_count_barplot
silva_class_count_barplot

comb_tax_ub %>%
  sapply(function(x) sum(x != "")) %>%
  as.data.frame() %>%
  `colnames<-`(c("Value")) %>%
  mutate(Level = rownames(.)) %>%
  .[2:dim(.)[1],] %>%
  replace_na(list(Value=0)) %>%
  mutate(Rank = unlist(strsplit(Level, "_"))[c(T, F)],
         Classifier = unlist(strsplit(Level, "_"))[c(F, T)],
         Value = Value/dim(comb_tax_ub)[1]) %>%
  mutate(Classifier = factor(Classifier, levels = c("BLAST", "RDP", "SINTAX", "Consensus"))) %>%
  ggplot(aes(x = Rank, y = Value, fill=Classifier)) +
  geom_bar(position="dodge", stat="identity") +
  scale_x_discrete(limits=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) +
  theme(axis.text.x.bottom = element_text(angle=-45, vjust=0.5, hjust=0)) + # Adjust the column labels
  labs(y = "Proportion of OTUs classified",
       title = "ITS - UNITE") -> unite_class_count_barplot
unite_class_count_barplot

g <- silva_class_count_barplot + unite_class_count_barplot + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
g
ggsave("../../figures/class_count_bps.png", g, units= "in", width = 8, height = 4)