input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_dense_results/denovo_noncoding_status.tsv"


data <- read.table(input_file, header = TRUE, sep = "\t")
data[, c("intergenic", "f.0", "f.1", "f.2")] <- lapply(data[, c("intergenic", "f.0", "f.1", "f.2")], function(x) ifelse(x < 3, 0, x))


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
n_several_origins <- nrow(origins[apply(origins, 1, function(x) sum(x == 0) < 3), ])
n_one_origin <- nrow(origins) - several_origins
# Camembert plot (pie chart)
pie(c(n_one_origin, n_several_origins),
    labels = c(paste("one origin (", n_one_origin, ")", sep = ""),
               paste("several origins (", n_several_origins, ")", sep = "")),
    col = c("#009E73", "#E69F00"),
    main = "For all de novo genes: number of different origins in the noncoding region")


# Altframes
n_f1 <- nrow(origins[origins$f.1 > 0, ])
n_f2 <- nrow(origins[origins$f.2 > 0, ])
# Camembert plot (pie chart)
pie(c(n_f1, n_f2),
    labels = c(paste("f+1 (", n_f1, ")", sep = ""),
               paste("f+2 (", n_f2, ")", sep = "")),
    col = c("#CC79A7", "#0072B2"),
    main = "Frame of origin for all de novo genes\ncoming from altframes")


# Length distributions
len_f1 <- data[data$f.1 > 0, "f.1"]
len_f2 <- data[data$f.2 > 0, "f.2"]
# Boxplot
boxplot(len_f1, len_f2,
        names = c("f+1", "f+2"),
        col = c("#CC79A7", "#0072B2"),
        xlab = "Frame of origin", ylab = "Length of the noncoding region (nt)",
        main = "Length distribution of the noncoding region\nfor de novo genes coming from altframes")




## De novo genes with only one origin
one_origin <- origins[apply(origins, 1, function(x) sum(x == 0) == 3), ]
n_only_f1 <- nrow(one_origin[one_origin$f.1 > 0, ])
n_only_f2 <- nrow(one_origin[one_origin$f.2 > 0, ])
# Number of genes in each
pie(c(n_only_f1, n_only_f2),
    labels = c(paste("f+1 (", n_only_f1, ")", sep = ""),
               paste("f+2 (", n_only_f2, ")", sep = "")),
    col = c("#CC79A7", "#0072B2"),
    main = "Frame of origin for all de novo genes\ncoming from only one altframe")


# Length distributions
len_only_f1 <- one_origin[one_origin$f.1 > 0, "f.1"]
len_only_f2 <- one_origin[one_origin$f.2 > 0, "f.2"]
# Boxplot
boxplot(len_only_f1, len_only_f2,
        names = c("f+1", "f+2"),
        col = c("#CC79A7", "#0072B2"),
        xlab = "Frame of origin", ylab = "Length of the noncoding region (nt)",
        main = "Length distribution of the noncoding region\nfor de novo genes coming from only one altframe")


# Print all genes with several origins
f_origins <- data[, c("f.0", "f.1", "f.2")]
data_several_origins <- data[rownames(f_origins[apply(f_origins, 1, function(x) sum(x == 0) < 2),]),]
# Order alphabetically by the row "outgroup" then ascending "noncoding_match_start"
data_several_origins <- data_several_origins[order(data_several_origins$outgroup, data_several_origins$noncoding_match_start), ]
print(data_several_origins)
# Save to file
write.table(data_several_origins, file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/status_several_origins.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
