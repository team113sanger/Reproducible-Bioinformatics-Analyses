library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(here)

## This script takes some of the data from the CDS1/CDS2 Single guide RNA 
## paper and does a mock analysis. Disclaimer it's not inspired by anyone's 
## real analysis!

## I've put together a bit of code few **purposefully** bad practices


## Your task is to fix them.
## You'll need to:
# 1) Correct any paths to relative ones (hint: use `here::here()`)
# 2) Try using pipe operator (`%>%` or |> ) instead of repeatedly assigning variables
# 3) Reorganise the code to read/write files in a more orderly structure
# 4) Capture the environment you used with renv




x <- read_tsv("/home/ubuntu/projects/Reproducible-Bioinformatics-Analyses/exercises/01_example_project/raw_count_matrix.tsv")

y<- x |> select(-c(is_unique_guide,id)) 

y <- y |> pivot_longer(cols = -c(sgrna_ids,sgrna_seqs, gene_symbol), names_to = "sample", values_to = "count")



p = ggplot(y, aes(x= sample,y=count)) + 
geom_boxplot()
ggsave(p,  filename = "/home/ubuntu/projects/Reproducible-Bioinformatics-Analyses/exercises/01_example_project/boxplot.pdf")  

df1 <- y |> 
separate(sample, into =c("cell_line", "screen"), sep = "_")

df2<- df1 |>
        separate(screen, into = c("screen", "replicate"), sep = "R")

q <- ggplot(df2, aes(x = log2(count + 1))) +
geom_density(aes(fill = `cell_line`)) +
facet_wrap(~replicate)
ggsave(q,  filename = "/home/ubuntu/projects/Reproducible-Bioinformatics-Analyses/exercises/01_example_project/density_plot.pdf")  


# Repeat the plotting for the filtered data (Guides with > 30 counts in the Plasmid sample)
guides_to_keep <- x |>  
select(-c(is_unique_guide,id)) |>
pivot_longer(cols = -c(sgrna_ids, sgrna_seqs, gene_symbol), names_to = "sample", values_to = "count") |> 
filter(sample == "Plasmid" & count >= 30) |> 
pull(sgrna_ids)
 
z <- y |> 
    filter(sgrna_ids %in% guides_to_keep) 

df1_filtered <- z |> 
separate(sample, into =c("cell_line", "screen"), sep = "_")

df2_filtered <- df1_filtered |>
        separate(screen, into = c("screen", "replicate"), sep = "R")

q2 <- ggplot(df2_filtered, aes(x = log2(count + 1))) +
geom_density(aes(fill = `cell_line`)) +
facet_wrap(~replicate)
ggsave(q2,  filename = "/home/ubuntu/projects/Reproducible-Bioinformatics-Analyses/exercises/03_example_project/density_plot_filtered.pdf")  

write_tsv(df2_filtered, "/home/ubuntu/projects/filtered_data.tsv")