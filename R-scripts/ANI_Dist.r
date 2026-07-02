##set working directory
setwd("/Users/martinsenekowitsch/Documents/FH/4_Semester/Thesis/fastANI")

library(tidyr)
library(dplyr)
library(ggplot2)

##load data
# naming scheme: data_COMPLETENESS_CONTAMINATION_N50(x1000)
data_all <- read.delim("fastANI_all_vs_infantis_reference_output.txt", sep = "\t", header = FALSE)
data_90_2_50 <- read.delim("90_2_50_filtered_ANI_results.txt", sep = "\t", header = FALSE)
data_90_2_100 <- read.delim("90_2_100_filtered_ANI_results.txt", sep = "\t", header = FALSE)
data_90_5_50 <- read.delim("90_2_50_filtered_ANI_results.txt", sep = "\t", header = FALSE)
data_90_5_100 <- read.delim("90_2_50_filtered_ANI_results.txt", sep = "\t", header = FALSE)
data_90_10_100 <- read.delim("90_10_100_filtered_ANI_results.txt", sep = "\t", header = FALSE)
data_all_vs_all <- read.delim("fastANI_all_vs_all_output.txt", sep = "\t", header = FALSE)

data_all_filter <- data_all %>% filter(V3 > 99.6)
data_90_2_50_filter <- data_90_2_50 %>% filter(V3 > 99.6)
data_90_2_100_filter <- data_90_2_100 %>% filter(V3 > 99.6)
data_90_5_50_filter <- data_90_5_50 %>% filter(V3 > 99.6)
data_90_5_100_filter <- data_90_5_100 %>% filter(V3 > 99.6)
data_90_10_100_filter <- data_90_10_100 %>% filter(V3 > 99.6)
data_all_vs_all_filter <- data_all_vs_all %>% filter(V3 < 99.96)


## plot



# Histogram with density curve
ggplot(data_all_filter, aes(x = data_all_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# data_90_2_50_filter
ggplot(data_90_2_50_filter, aes(x = data_90_2_50_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# data_90_2_100_filter
ggplot(data_90_2_100_filter, aes(x = data_90_2_100_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# data_90_5_50_filter
ggplot(data_90_5_50_filter, aes(x = data_90_5_50_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# data_90_5_100_filter
ggplot(data_90_5_100_filter, aes(x = data_90_5_100_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# data_90_10_100_filter
ggplot(data_90_10_100_filter, aes(x = data_90_10_100_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

# combined
combined <- bind_rows(
  data_all_filter       %>% mutate(source = "all"),
  data_90_2_50_filter   %>% mutate(source = "90_2_50"),
  data_90_2_100_filter  %>% mutate(source = "90_2_100"),
  data_90_5_50_filter   %>% mutate(source = "90_5_50"),
  data_90_5_100_filter  %>% mutate(source = "90_5_100"),
  data_90_10_100_filter  %>% mutate(source = "90_10_100")
)
col3 <- names(combined)[3]

ggplot(combined, aes(x = .data[[col3]])) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 fill = "skyblue",
                 color = "white") +
  geom_density(alpha = 0.2, fill = "blue") +
  facet_wrap(~ source, ncol = 3) +
  labs(
    title = "Distribution of 3rd Column",
    x = col3,
    y = "Density"
  ) +
  theme_minimal()

ggplot(combined, aes(x = .data[[col3]], color = source, fill = source)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  labs(
    title = "Overlapping density distributions (3rd column)",
    x = col3,
    y = "Density"
  ) +
  theme_minimal()

# all vs all
ggplot(data_all_vs_all, aes(x = data_all_vs_all[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

ggplot(data_all_vs_all_filter, aes(x = data_all_vs_all_filter[[3]])) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white") +
  geom_density(alpha = .2, fill = "blue") +
  labs(title = "Distribution of 3rd Column", 
       x = names(data)[3], 
       y = "Density") +
  theme_minimal()

