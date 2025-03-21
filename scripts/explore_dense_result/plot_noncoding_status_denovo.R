input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_dense_results/denovo_noncoding_status.tsv"


data <- read.table(input_file, header = TRUE, sep = "\t")


## Total number of noncoding origins
inter <- sum(data[, "intergenic"])
f0 <- sum(data[, "f.0"])
f1 <- sum(data[, "f.1"])
f2 <- sum(data[, "f.2"])
# barplot
barplot(c(inter, f0, f1, f2),
        names.arg = c("intergenic", "f+0", "f+1", "f+2"),
        col = c("#009E73", "#E69F00", "#CC79A7", "#0072B2"),
        xlab = "NC origin status", ylab = "Number of nucleotides",
        main = "For all outgroups of all de novo genes:\nlocation of each nucleotide in the noncoding (NC) region")


## Percentage of genes that have several origins
origins <- data[, c("intergenic", "f.0", "f.1", "f.2")]
# get all rows for which less than 3 cols are 0
several_origins <- nrow(origins[apply(origins, 1, function(x) sum(x == 0) < 3), ])
one_origin <- nrow(origins) - several_origins
# Camembert plot (pie chart)
pie(c(one_origin, several_origins),
    labels = c(paste("one origin (", one_origin, ")", sep = ""),
               paste("several origins (", several_origins, ")", sep = "")),
    col = c("#009E73", "#E69F00"),
    main = "For all de novo genes: number of different origins in the noncoding region")
