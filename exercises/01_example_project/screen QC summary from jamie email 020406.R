library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(fs)
library(jsonlite)
library(purrr)


read_pycroquet_json<- function(file){
      data <- fromJSON(file)
  
  # Extract sample name
  sample_name <- data$sample_name
    total_reads <- data$total_reads
  mapped_reads <- data$mapped_to_guide_reads
  alignment_rate <- ifelse(total_reads > 0, (mapped_reads / total_reads * 100), 0)
  
  # Calculate guide coverage statistics
  total_guides <- data$total_guides
  zero_count_guides <- data$zero_count_guides
  low_count_lt_15 <- data$low_count_guides_lt_15
  mean_count <- data$mean_count_per_guide
  
  # Calculate guide coverage percentage
  guides_with_counts <- total_guides - zero_count_guides
  guide_coverage_pct <- ifelse(total_guides > 0, (guides_with_counts / total_guides * 100), 0)
  
  # Create the datagram row
  row <- data.frame(
    sample_id = sample_name,
    total_reads = total_reads,
    mapped_reads = mapped_reads,
    alignment_rate = round(alignment_rate, 2),
    unmapped_reads = data$unmapped_reads,
    guide_count = total_guides,
    guides_with_zero_counts = zero_count_guides,
    guides_with_low_counts_lt_15 = low_count_lt_15,
    mean_count_per_guide = round(mean_count, 2),
    guide_coverage_pct = round(guide_coverage_pct, 2),
    pycroquet_version = data$version,
    stringsAsFactors = FALSE
  )
  
  return(row)
  return(x)
}

# Metadata recieved from Alan, 
# Shared from a Slack message from Jamie dated 1st May 2023
# See thread in Misc general informatics
x <- read_tsv("/home/ubuntu/projects/Reproducible-Bioinformatics-Analyses/exercises/01_example_project/metadata.tsv")

files <-fs::dir_ls(here("exercises/01_example_project/pycroquet"))

qc_data <- map_dfr(files, read_pycroquet_json) 

ggplot(qc_data, aes(x = sample_id, y = mapped_reads)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Number of Mapped Reads per sample", x = "Sample ID", y = "Guide coverage percentage")

)



x
plot_list <- qc_data |> 
group_by(sample_id) |>
group_split()


for (file %in% plot_list){
  ggplot(file, aes(y = ))
}