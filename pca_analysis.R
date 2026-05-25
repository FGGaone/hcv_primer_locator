library(ggplot2)
library(dplyr)

# Read the mutation matrix
data <- read.csv("/home/ubuntu/mutation_matrix.csv")

# Extract features (mutations) and labels
features <- data %>% select(-anti_hbc_status, -id)
labels <- data$anti_hbc_status

# Remove zero-variance columns
features <- features[, apply(features, 2, var) != 0]

# Perform PCA
pca_res <- prcomp(features, scale. = TRUE)

# Prepare data for plotting
pca_data <- as.data.frame(pca_res$x)
pca_data$status <- factor(labels, 
                          levels = c(2, 1, 0), 
                          labels = c("Positive", "Negative", "Unknown"))
pca_data$id <- data$id

# Calculate explained variance
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2) * 100

# Create PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = status, label = id)) +
  geom_point(size = 3, alpha = 0.8) +
  # stat_ellipse(aes(fill = status), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
  labs(
    title = "PCA of HBV Mutation Profiles by anti-HBc Status",
    subtitle = "Relationship of Unknown samples to Positive and Negative statuses",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "anti-HBc Status"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Positive" = "#e41a1c", "Negative" = "#377eb8", "Unknown" = "#4daf4a"))

# Save the plot
ggsave("/home/ubuntu/pca_anti_hbc.png", plot = p, width = 10, height = 8)

# Calculate Euclidean distances to centroids
centroids <- pca_data %>%
  filter(status != "Unknown") %>%
  group_by(status) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2))

unknown_samples <- pca_data %>% filter(status == "Unknown")

dist_results <- list()
for(i in 1:nrow(unknown_samples)) {
  u_pc1 <- unknown_samples$PC1[i]
  u_pc2 <- unknown_samples$PC2[i]
  
  dist_pos <- sqrt((u_pc1 - centroids$PC1[centroids$status == "Positive"])^2 + 
                   (u_pc2 - centroids$PC2[centroids$status == "Positive"])^2)
  dist_neg <- sqrt((u_pc1 - centroids$PC1[centroids$status == "Negative"])^2 + 
                   (u_pc2 - centroids$PC2[centroids$status == "Negative"])^2)
  
  dist_results[[i]] <- data.frame(
    id = unknown_samples$id[i],
    dist_to_positive = dist_pos,
    dist_to_negative = dist_neg,
    closer_to = ifelse(dist_pos < dist_neg, "Positive", "Negative")
  )
}

final_dist <- do.call(rbind, dist_results)
write.csv(final_dist, "/home/ubuntu/pca_distances.csv", row.names = FALSE)
print(final_dist)
