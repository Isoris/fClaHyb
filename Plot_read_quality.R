library(ggExtra)  # For marginal plots
library(ggplot2)
library(magick)
library(patchwork)
library(gridExtra)
library(reshape2)
library(cowplot)

setwd("C:/Users/quent/Desktop/MS_CMA/Figures/20-reads_qc")

HIFI_reads <- read.table(file = "CG_F_HIFI_None_reads_read_lengths_and_phred_scores.csv", sep = "\t", header = T)

head(HIFI_reads)

HIFI_reads$Length_.bp. <- as.numeric(HIFI_reads$Length_.bp.)
HIFI_reads$GC_content  <- as.numeric(HIFI_reads$GC_content )
HIFI_reads$Mean_PHRED_Score <- as.numeric(HIFI_reads$Mean_PHRED_Score)

HIFI_reads_sample <- HIFI_reads[sample(nrow(HIFI_reads), 100000), ]

# Main scatter plot (Mean Phred Score vs Read Length) with GC content as fill color
scatter_plot <- ggplot(HIFI_reads_sample, aes(x = Length_.bp., 
                                              y = Mean_PHRED_Score,
                                              color = GC_content)) +
  geom_point(size = 0.05) + 
  geom_density_2d(linewidth = 0.5, color = "white" ) + # Use geom_tile for the grid-like scatter plot
  scale_color_gradient(high = "#052950", low = "#0b8583", limits = c(20, 60)) +  
  #scale_color_gradient(high = "#052950", low = "#06595B") + # Color gradient from #052950 to #05595B
  labs(x = "Read length (bp)", y = "Mean PHRED Score", fill = "GC Content") +
  theme_minimal() + xlim(100, max(HIFI_reads_sample$Length_.bp.)) + scale_x_log10() + 
  ggtitle("HIFI reads dataset") +
  theme(
    #   plot.title = element_text(hjust = 0, faHI mily = "Bahnschrift", size = 20),
    axis.ticks.y = element_blank(),
    #    axis.text = element_text(family = "Bahnschrift", size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.5),
    panel.border = element_blank(),
    # Set the legend inside the plot at the upper right corner
    legend.position = c(1, 1),  # Position at top-right inside the plot
    legend.justification = c(1, 1),  # Adjust the justification to the corner
    legend.title = element_text( size=7), 
    legend.text=element_text(size=7),
    plot.title = element_text(vjust = -7)
  ) 

# Use ggMarginal to add the marginal histogram above the scatter plot
scatter_with_marginal <- ggMarginal(scatter_plot, type = "histogram", bins = 1000, fill = "#052950", 
                                    color = "#052950", alpha = 0.8)
plot(scatter_with_marginal)

ggsave(filename = "CG_F_scatter_with_marginal_HIFI.pdf", 
       plot = scatter_with_marginal, 
       device = cairo_pdf, width = 5, height = 3.5, units = "in")


ONT_reads <- read.table(file = "CG_F_nanopore_None_reads_read_lengths_and_phred_scores.csv", sep = "\t", header = T)
head(ONT_reads)

ONT_reads$Length_.bp. <- as.numeric(ONT_reads$Length_.bp.)
ONT_reads$GC_content  <- as.numeric(ONT_reads$GC_content )
ONT_reads$Mean_PHRED_Score <- as.numeric(ONT_reads$Mean_PHRED_Score)

# Randomly sample a smaller portion of the dataset for the plot
ONT_reads_sample <- ONT_reads[sample(nrow(ONT_reads), 100000), ]


# Sort reads by length in descending order
sorted_lengths <- sort(ONT_reads$Length_.bp., decreasing = TRUE)

# Calculate total length of all reads
total_length <- sum(sorted_lengths)

# Function to calculate NXX (N50, N90, etc.)
calculate_NXX <- function(sorted_lengths, threshold) {
  cumulative_lengths <- cumsum(sorted_lengths)
  NXX_index <- which(cumulative_lengths >= (threshold * total_length))[1]
  return(sorted_lengths[NXX_index])
}

# Calculate N50 (50% of total length)
N50 <- calculate_NXX(sorted_lengths, 0.5)

# Calculate N90 (90% of total length)
N1 <- calculate_NXX(sorted_lengths, 0.9)

# Max read length
max_length <- max(sorted_lengths)

cat("N50:", N50, "\n")
cat("N1:", N1, "\n")
cat("Max Read Length:", max_length, "\n")



scatter_plot_ONT <- ggplot(ONT_reads_sample, aes(x = Length_.bp., 
                                                 y = Mean_PHRED_Score,
                                                 color = GC_content)) +
  geom_point(size = 0.05, alpha = 0.5) + 
  geom_density_2d(linewidth = 0.5, color = "white" ) + 
  scale_color_gradient(high = "#603601", low = "#ba9261", limits = c(20, 60)) +  
  labs(x = "Read length (bp)", y = "Mean PHRED Score", fill = "GC Content") +
  theme_minimal() + 
  xlim(100, max(ONT_reads_sample$Length_.bp.)) +ylim(5,23)+
  scale_x_log10() + 
  ggtitle("ONT reads dataset") +
  
  # Add N50 and N1 vertical lines with geom_segment for annotations
  geom_vline(xintercept = N50, linetype = "dashed", color = "#baa891", size = 0.4) + 
  geom_vline(xintercept = N1, linetype = "dashed", color = "#baa891", size = 0.4) + 
  
  # Use geom_segment to annotate N50 with a segment
  geom_segment(aes(x = N50, xend = N50, y = 0, yend = 21), 
               linetype = "dashed", color = "black", size = 0) +
  geom_text(aes(x = N50, y = 22, label = paste("N50 =", N50)), 
            color = "black", size = 2, hjust = 0) +
  
  # Use geom_segment to annotate N1 with a segment
  geom_segment(aes(x = N1, xend = N1, y = 0, yend = 21), 
               linetype = "dashed", color = "black", size = 0) +
  geom_text(aes(x = N1, y = 22, label = paste("N1 =", N1)), 
            color = "black", size = 2, hjust = 0) +
  
  theme(
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.5),
    panel.border = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.title = element_text(size = 7), 
    legend.text = element_text(size = 7),
    plot.title = element_text(vjust = -7)
  )

# Add marginal histogram to the scatter plot
scatter_with_marginal_ONT <- ggMarginal(scatter_plot_ONT, type = "histogram", bins = 1000, 
                                        fill = "#603601", color = "#603601", alpha = 0.5)

# Save plot
ggsave(filename = "CG_F_scatter_with_marginal_ONT.pdf", 
       plot = scatter_with_marginal_ONT, 
       device = cairo_pdf, width = 5, height = 3.5, units = "in")

# Combine the two plots side by side using grid.arrange
combined <- grid.arrange(scatter_with_marginal, scatter_with_marginal_ONT, ncol = 2, nrow = 1)

plot(combined)

ggsave(filename = "CG_F_scatter_with_marginal_HIFI_ONT_v2.pdf", 
       plot = combined, 
       device = cairo_pdf, width = 10, height = 7, units = "in")

# Install cowplot and magick if not installed


# Load Genomescope image
genomescope_image <- magick::image_read("Genomescope_HIFI_Illumina.png")
# Get the dimensions of the image to calculate the new width and height after cropping
image_info <- image_info(genomescope_image)
new_width <- image_info$width  # Subtract 15 pixels from the width
new_height <- image_info$height - 20  # Subtract 15 pixels from the height

# Crop the image by removing 15 pixels from the right and 15 pixels from the top
cropped_image <- image_crop(genomescope_image, paste0(new_width, "x", new_height, "+0+15"))

# Convert the cropped image to a ggplot object
genomescope_plot <- ggdraw() + draw_image(cropped_image)

# Arrange the plots (50% for the image and 50% for the plot grid)
combined_final <- plot_grid(
  genomescope_plot,            # The external image plot
  combined,
  ncol = 1,
  nrow = 2, # 2 columns
  rel_widths = c(1, 1)   # Equal widths for the plots and image
)

# Display the combined plot with the image
print(combined_final)

ggsave(filename = "CG_F_scatter_with_marginal_HIFI_ONT_Genomesscope_v3.pdf", 
       plot = combined_final, 
       device = cairo_pdf, width =  10, height = 10, units = "in")

# Sample table data (modify this if necessary)
table_data <- data.frame(
  #  Library_Type = c("Pacific Biosciences", "Oxford Nanopore 1D ONT (NovoGene)", "HiC (Phase Genomics) DnpII RE", 
  #   "HiC (Phase Genomics) DnpII RE", "Illumina (NovoGene)", "Illumina (NovoGene)"),
  Protocol = c("HiFi","HiFi (mapped)", "ONT","ONT (mapped)", "Hi-C, PE (150) FWD", "Hi-C, PE (150) REV", "WGS, PE (150) FWD", "WGS, PE (150) REV"),
  Coverage = c(10,"~3-10 reads", 30, "~2-5 reads", 37, 37, 40, 40)
)


# Convert table data to tableGrob for plotting
table_plot <- tableGrob(table_data, rows = NULL)

# Combine Genomescope image and table below it
genomescope_with_table <- plot_grid(
  table_plot,
  genomescope_plot, 
  ncol = 2, 
  rel_heights = c(1, 1) ) # Adjust relative heights to move both up)

# Combine the plots with the new layout
combined_image_table <- plot_grid(            
  genomescope_with_table,        # First row with Genomescope and table combined
  combined,                      # Second row with HIFI and ONT combined plot
  nrow = 2,                      # Two columns: first column is "combined" plot, second is image+table
  rel_widths = c(1, 1)           # Adjust width equally for both columns
)

# Save the result
ggsave(filename = "CG_F_combined_Genomescope_Table_v3.png", 
       plot = combined_image_table, 
       device = "png", width = 8.5, height = 8.5, units = "in")


# Your existing HIFI and ONT plots
# Assuming scatter_with_marginal and scatter_with_marginal_ONT are created previously

# Hi-C data for haplotypes (as provided by you)
data <- data.frame(
  Metric = c("Sequenced Read Pairs", "No chimera found", "One or both reads unmapped", "2 alignments", "3 or more alignments", "Ligation Motif Present", "Total Unique", "Total Duplicates", "Hi-C Contacts", "Inter-chromosomal", "Intra-chromosomal", "Short Range (<20Kb)", "Long Range (>20Kb)"),
  hap0 = c(123959443, 4670875, 4670875, 102740927, 16547641, 82122529, 74765772, 27975155, 59403641, 22710062, 36693579, 17467521, 19718858),
  hap1 = c(123959443, 4616080, 4616080, 102818325, 16525038, 82122529, 74814629, 28003696, 59467557, 22740327, 36727230, 17469127, 19757703)
)

# Melt the data for easy plotting
data_long <- melt(data, id.vars = "Metric", variable.name = "Haplotype", value.name = "Value")

# Custom number formatting function
format_millions <- function(x) {
  paste0(formatC(x / 1e6, format = "f", digits = 1), " M")
}

# Create Hi-C bar plot with coord_flip
hic <- ggplot(data_long, aes(x = Metric, y = Value, fill = Haplotype)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#3b528b", "#fdEEE5")) +  
  theme_minimal() + 
  labs(title = "HiC Metrics", y = "Value", x = "") +
  scale_y_continuous(labels = format_millions) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.5),
    panel.border = element_blank(),
    legend.position = c(1, 1),  
    legend.justification = c(1, 1),
    legend.title = element_text(size = 7), 
    legend.text = element_text(size = 7),
    plot.title = element_text(vjust = 0, size = 10)
  )


# Combine ONT and HIFI plots in a vertical layout (50%)
right_side <- plot_grid(scatter_with_marginal, scatter_with_marginal_ONT, nrow = 2, rel_heights = c(1, 1))

# Combine Genomescope image and HiC plot in a vertical layout (50%)
# HiC plot should take 15-20% space, so use rel_heights
left_side <- plot_grid(genomescope_plot, hic, ncol = 1, rel_heights = c(0.8, 0.5))

# Combine left (ONT + HIFI) and right (Genomescope + HiC) side by side
final_plot <- plot_grid(left_side, right_side, ncol = 2, rel_widths = c(1,0.9))

# Display the final combined plot
print(final_plot)

ggsave(filename = "CG_F_scatter_with_marginal_HIFI_ONT_Genomesscope_hic_mm.pdf", 
       plot = final_plot, 
       device = cairo_pdf, width =  183, height = 150, units = "mm")

ggsave(filename = "CG_F_scatter_with_marginal_HIFI_ONT_Genomesscope_hic.pdf", 
       plot = final_plot, 
       device = cairo_pdf, width =  10, height = 10, units = "in")
