library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

setwd("../../data/speed_tests/")
setwd("C:/Users/julia/OneDrive - Michigan State University/Documents/MSU/Bonito Lab/CONSTAX_fixes/CONSTAXv2_ms_code/data/speed_tests/")

classify_dat <- read.csv("speed_classify_comb.csv")
classify_dat <- classify_dat %>%
  mutate(time_per_seq=seq_count/time)
head(classify_dat)

classify_dat %>%
  filter(step == "classify") -> cl_df
cl_df %>%
  group_by(threads, algorithm, seq_count) %>%
  summarise(mean_time=mean(time_per_seq), sd_time = sd(time_per_seq)) -> sum_cl_df

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

ggplot(cl_df, mapping=aes(x = threads, y = 1000/time, color = algorithm)) +
  geom_point() +
  stat_smooth(method = 'nls',
              formula = 'y~a*x^b',
              method.args = list(start=c(a = 1., b=0.5)),
              se=FALSE)

n_vary_train <- read.csv("speed_training_t1_nvary.csv")

n_vary_train <- n_vary_train %>%
  filter(step == "train") %>%
  mutate(seq_per_time=seq_count/time)

sum_n_vary_t <- n_vary_train %>%
  group_by(algorithm, seq_count) %>%
  summarise(mean_seq_per_time = mean(seq_per_time),
            sd_seq_per_time = sd(seq_per_time))
sum_n_vary_t

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
       y = "Sequences trainied per second",
       color = "Algorithm") +
  scale_color_discrete(breaks=c("blast", "utax"), labels=c("BLAST", "UTAX")) #+

g_t1

g <- g_t1 + class_lp + plot_layout(guides = "collect")
g

ggsave("../../figures/speed_train_line_class_line.png", g, width = 12, height = 8, units = "in", dpi = 400)
