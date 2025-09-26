library(tidyverse)
library(ggpubr)

# Experiment 1
df_1 <- read.csv("/Users/ascarpa/Transib1_Hybrid_Dysgenesis/Transib1_Hybrid_Dysgenesis_Experiment_1.csv", stringsAsFactors = FALSE) %>%
  mutate(cross = paste0(sub("x.*", "", ID), "_F x ", sub(".*x", "", sub("_.*", "", ID)), "_M"))

df_1 <- df_1 %>%
  mutate(normal_ovaries = (Normal*2) + Intermediate,
         dysgenic_ovaries = 100 - normal_ovaries)


cross_1 <- c("Tom008_F x Ten015_M", "Ten015_F x Tom008_M")
cross_2 <- c("Cor018_F x Aka018_M", "Aka018_F x Cor018_M")

temps <- c(29, 25, 21)
crosses <- list(cross_1, cross_2)

# Loop over combinations of cross and temperature
for (cross in crosses) {
  for (temp in temps) {
    # Filter the data
    df_filtered <- df_1 %>%
      filter(Temperature == temp, cross %in% cross)
    
    # Create the contingency table
    table_mat <- df_filtered %>%
      group_by(cross) %>%
      summarise(normal = sum(normal_ovaries),
                dysgenic = sum(dysgenic_ovaries)) %>%
      column_to_rownames("cross") %>%
      as.matrix()
    
    # Perform chi-squared test
    chisq_result <- chisq.test(table_mat)
    
    # Print the result
    print(paste("Chi-squared test for Temperature:", temp, "and Cross:", cross))
    print(chisq_result)
  }
}

ggplot(df_1, aes(x = cross, y = Prop_Dysgenic, fill = factor(Temperature))) +
  geom_boxplot() +
  labs(title = "Proportion of Dysgenic Ovaries for Different Crosses and Temperatures",
       x = "Cross",
       y = "Proportion of Dysgenic Ovaries",
       fill = "Temperature") +
  scale_fill_manual(values = c("dodgerblue", "darkorange1", "indianred1")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Experiment 1 old

df_29 <- read.csv("/Users/ascarpa/Transib1_Hybrid_Dysgenesis/Transib1_Hybrid_Dysgenesis_29.csv", stringsAsFactors = FALSE) %>%
  mutate(cross = paste0(sub("x.*", "", ID), "_F x ", sub(".*x", "", sub("_.*", "", ID)), "_M"))
df_25 <- read.csv("/Users/ascarpa/Transib1_Hybrid_Dysgenesis/Transib1_Hybrid_Dysgenesis_25.csv", stringsAsFactors = FALSE) %>%
  mutate(cross = paste0(sub("x.*", "", ID), "_F x ", sub(".*x", "", sub("_.*", "", ID)), "_M"))

df_29_filtered <- df_29 %>%
  filter(cross %in% c("Tom008_F x Ten015_M", "Ten015_F x Tom008_M", 
                      "Cor018_F x Aka018_M", "Aka018_F x Cor018_M"))
df_25_filtered <- df_25 %>%
  filter(cross %in% c("Tom008_F x Ten015_M", "Ten015_F x Tom008_M", 
                      "Cor018_F x Aka018_M", "Aka018_F x Cor018_M"))

ggplot(df_29_filtered, aes(x = cross, y = Prop_Dysgenic)) +
  geom_boxplot(fill = "indianred1") +
  ggtitle("Transib1 induced dysgenesis in D. melanogaster, crosses at 29°C") +
  labs(x = "Cross", y = "Proportion Dysgenic") +
  stat_compare_means(comparisons = list(c("Tom008_F x Ten015_M", "Ten015_F x Tom008_M"),
                                        c("Cor018_F x Aka018_M", "Aka018_F x Cor018_M")),
                     method = "t.test", label = "p.signif") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  theme_minimal()

ggplot(df_25_filtered, aes(x = cross, y = Prop_Dysgenic)) +
  geom_boxplot(fill = "darkorange1") +
  ggtitle("Transib1 induced dysgenesis in D. melanogaster, crosses at 25°C") +
  labs(x = "Cross", y = "Proportion Dysgenic") +
  stat_compare_means(comparisons = list(c("Tom008_F x Ten015_M", "Ten015_F x Tom008_M"),
                                        c("Cor018_F x Aka018_M", "Aka018_F x Cor018_M")),
                     method = "t.test", label = "p.signif") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  theme_minimal()

# Experiment 2

df_e2 <- read.csv("/Users/ascarpa/Transib1_Hybrid_Dysgenesis/Transib1_Hybrid_Dysgenesis-Experiment_2.csv", stringsAsFactors = FALSE) %>%
  mutate(cross = paste0(sub("x.*", "", ID), "_F x ", sub(".*x", "", sub("_.*", "", ID)), "_M"))

order <- c("Har_F x Dmel68_M", "Dmel68_F x Har_M", 
                   "Har_F x Tom008_M", "Tom008_F x Har_M", 
                   "Har_F x Ten015_M", "Ten015_F x Har_M")

df_e2 <- df_e2 %>%
  mutate(cross = factor(cross, levels = order))


ggplot(df_e2, aes(x = cross, y = Prop_Dysgenic, fill = factor(Pair))) +
  geom_boxplot(show.legend = FALSE) +
  ggtitle("Transib1 induced dysgenesis in D. melanogaster, crosses at 29°C") +
  labs(x = "Cross", y = "Proportion Dysgenic") +
  stat_compare_means(comparisons = list(c("Har_F x Dmel68_M", "Dmel68_F x Har_M"),
                                        c("Har_F x Tom008_M", "Tom008_F x Har_M"),
                                        c("Har_F x Ten015_M", "Ten015_F x Har_M")),
                     method = "t.test", label = "p.signif") +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  theme_minimal()
