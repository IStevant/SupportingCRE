source(".Rprofile")

###########################################
#                                         #
#               Libraries                 #
#                                         #
###########################################
suppressPackageStartupMessages({
  library("ggplot2")
  library("cowplot")
})

###########################################
#                                         #
#               Load data                 #
#                                         #
###########################################

linkage <- read.table(snakemake@input[["linkage"]], header=TRUE)

###########################################
#                                         #
#               Plot data                 #
#                                         #
###########################################

# Link per gene

gene_split <- split(linkage, linkage$Gene)

count_correlations <- function(data) {
  pos <- sum(data$correlations > 0)
  neg <- sum(data$correlations < 0)
  return(data.frame(
    Gene = unique(data$Gene),
    Type = c(rep("positive", length(pos)), rep("negative", length(neg))),
    Count = c(pos, neg)
  ))
}

summary_links <- do.call(rbind, lapply(gene_split, count_correlations))
median_groups <- data.frame(
  median = c(
    median(summary_links$Count[summary_links$Type == "positive"]),
    median(summary_links$Count[summary_links$Type == "negative"])
  ),
  mean = c(
    mean(summary_links$Count[summary_links$Type == "positive"]),
    mean(summary_links$Count[summary_links$Type == "negative"])
  ),
  Type = c("positive", "negative") 
)

link_number_gene <- ggplot(summary_links, aes(x=Count, color=Type, fill=Type)) +
  geom_bar(data=subset(summary_links,Type=="positive"), alpha=0.5) + 
  geom_bar(data=subset(summary_links,Type=="negative"),aes(y=..count..*(-1)), alpha=0.5) + 
  scale_fill_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  scale_color_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  annotate(
    geom="text", 
    x=30, 
    y=-1000, 
    label=paste(
      "Median:", 
      round(subset(median_groups, Type=="negative")$median, digit=2),
      "\nMean:",
      round(subset(median_groups, Type=="negative")$mean, digit=2)
    ), 
    color="#2191FB"
  ) +
  annotate(
    geom="text", 
    x=30, 
    y=1000, 
    label=paste(
      "Median:", 
      round(subset(median_groups, Type=="positive")$median, digit=2),
      "\nMean:",
      round(subset(median_groups, Type=="positive")$mean, digit=2)
    ), 
    color="#EF6351"
  ) +
  scale_y_continuous(labels = \(x) scales::comma(abs(x))) +
  ggtitle("Number of linked OCRs per gene") +
  ylab("Count") +
  coord_flip() +
  theme_light() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    aspect.ratio = 0.5
  )

# Link per OCR

gene_split <- split(linkage, linkage$Peak)

count_correlations <- function(data) {
  pos <- sum(data$correlations > 0)
  neg <- sum(data$correlations < 0)
  return(data.frame(
    peak = unique(data$Peak),
    Type = c(rep("positive", length(pos)), rep("negative", length(neg))),
    Count = c(pos, neg)
  ))
}

summary_links <- do.call(rbind, lapply(gene_split, count_correlations))
median_groups <- data.frame(
  median = c(
    median(summary_links$Count[summary_links$Type == "positive"]),
    median(summary_links$Count[summary_links$Type == "negative"])
  ),
  mean = c(
    mean(summary_links$Count[summary_links$Type == "positive"]),
    mean(summary_links$Count[summary_links$Type == "negative"])
  ),
  Type = c("positive", "negative") 
)

link_number_OCR <- ggplot(summary_links, aes(x=Count, color=Type, fill=Type)) +
  geom_bar(data=subset(summary_links,Type=="positive"), alpha=0.5) + 
  geom_bar(data=subset(summary_links,Type=="negative"),aes(y=..count..*(-1)), alpha=0.5) + 
  scale_fill_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  scale_color_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  annotate(
    geom="text", 
    x=12, 
    y=-10000, 
    label=paste(
      "Median:", 
      round(subset(median_groups, Type=="negative")$median, digit=2),
      "\nMean:",
      round(subset(median_groups, Type=="negative")$mean, digit=2)
    ), 
    color="#2191FB"
  ) +
  annotate(
    geom="text", 
    x=12, 
    y=10000, 
    label=paste(
      "Median:", 
      round(subset(median_groups, Type=="positive")$median, digit=2),
      "\nMean:",
      round(subset(median_groups, Type=="positive")$mean, digit=2)
    ), 
    color="#EF6351"
  ) +
  scale_y_continuous(labels = \(x) scales::comma(abs(x))) +
  ggtitle("Number of linked genes per OCRs") +
  ylab("Count") +
  coord_flip() +
  theme_light() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    aspect.ratio = 0.5
  )

# Link distance to TSS

count_distances <- function(data) {
  pos <- data[data$correlations > 0, "Summit2TSS"]
  neg <- data[data$correlations < 0, "Summit2TSS"]
  return(data.frame(
    Type = c(rep("positive", length(pos)), rep("negative", length(neg))),
    Distance = c(pos, neg)/1000
  ))
}

summary_links <- count_distances(linkage)

link_distance <- ggplot(summary_links, aes(x=Distance, color=Type, fill=Type)) +
  geom_histogram(alpha=0.1, position="identity") +
  # geom_density(alpha=0.25) +
  scale_fill_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  scale_color_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  ggtitle("OCR distance to linked gene TSS") +
  xlab("Distance to TSS (kb)") +
  scale_y_continuous(labels = scales::comma) +
  ylab("Count") +
  theme_light() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    aspect.ratio = 0.5
  )

# Link features


count_feature <- function(data) {
  pos <- data[data$correlations > 0, "Type"]
  neg <- data[data$correlations < 0, "Type"]
  return(data.frame(
    Type = c(rep("positive", length(pos)), rep("negative", length(neg))),
    Feature = c(pos, neg)
  ))
}


summary_links <- count_feature(linkage)

link_feature <- ggplot(summary_links, aes(x=Feature, color=Type, fill=Type)) +
  # geom_histogram(aes(y=..density..), alpha=0.1, position="identity") +
  geom_bar(alpha=0.25, position = "dodge") +
  scale_fill_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  scale_color_manual(values=c(positive = "#EF6351", negative = "#2191FB")) +
  ggtitle("Location of linked OCRs") +
  xlab("Genomic features") +
  scale_y_continuous(labels = scales::comma) +
  ylab("Count") +
  theme_light() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    aspect.ratio = 0.5
  )



###########################################
#                                         #
#              Save plots                 #
#                                         #
###########################################

figure <- plot_grid(
  plotlist = list(
    link_number_gene,
    link_number_OCR,
    link_distance,
    link_feature
  ),
  labels = "AUTO",
  ncol = 2,
  align = "hv"
)

# As PDF
save_plot(
  snakemake@output[["pdf"]],
  figure,
  base_width = 25,
  base_height = 11,
  units = c("cm"),
  dpi = 300
)

# As PNG
save_plot(
  snakemake@output[["png"]],
  figure,
  base_width = 25,
  base_height = 11,
  units = c("cm"),
  dpi = 300,
  bg = "white"
)

