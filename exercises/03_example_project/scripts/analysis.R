library(ggplot2)
library(dplyr)

# some analysis
d = read.csv("~/data/counts.txt") 
d2 = d[d$count > 10, ]  
p = ggplot(d2, aes(x=sample, y=count)) + geom_bar(stat="identity")
ggsave("~/results/plot.pdf")  


# Here we have a few analysis steps that assign the results to variables

data1 <- merge(expression_tidy, sample_metadata, by = 'sample')
data2 <- data1[data1$expression > 5, ]
data3 <- aggregate(expression ~ gene + group, data = data2, mean)
data4 <- data3[order(data3$expression, decreasing = TRUE), ]
