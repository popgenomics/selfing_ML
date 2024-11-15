#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages(library(tidyverse))

# Parse command-line arguments
args = commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
	stop("Usage: Rscript plot_selfing_sequences.R <filename.txt>")
}

# Get the input filename from the command-line argument
filename = args[1]

# Check if the input file exists
if (!file.exists(filename)) {
	stop(paste("File not found:", filename))
}

# Define the output PDF filename by replacing '.txt' with '.pdf'
output_filename = sub("\\.txt$", ".pdf", filename)

# Load and filter the data
data = tibble(read.table(filename, header = TRUE))
filtered_data = data %>%
	group_by(individual) %>%
	dplyr::filter(logLik != min(logLik))

# Find the optimal selfing rate for each individual
optimal_selfing = data %>%
	group_by(individual) %>%
	slice_max(order_by = logLik, n = 1) %>%
	ungroup() %>%
	select(individual, selfing_rate)

# Calculate the mean of the optimal selfing rates
mean_selfing_rate = mean(optimal_selfing$selfing_rate)
print(paste("Mean selfing rate:", round(mean_selfing_rate, 2)))

# Create the plot
plot = filtered_data %>%
	ggplot(aes(x = selfing_rate, y = logLik)) +
	geom_line(size = 2) +
	geom_vline(data = optimal_selfing, aes(xintercept = selfing_rate), color = "red", linetype = "dashed", size = 1) +
	facet_wrap(~ individual, scales = 'free_y') +
	ggtitle(paste(unique(data$species), '\ns (mean) = ', round(mean_selfing_rate, 2), sep = '')) +
	theme_bw(base_size = 12) +
	xlab('Selfing Rate (s)') +
	ylab('Log-Likelihood') +
	xlim(0, 1)

# Save the plot to a PDF file
ggsave(output_filename, plot, width = 10, height = 8)
print(paste("Plot saved to:", output_filename))

