library(ggplot2)

input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequence_features.csv"
data <- read.table(input_file, header = TRUE, sep = "\t")
data_gc <- data[data$feature == "gc", ]
data_aro <- data[data$feature == "aromaticity", ]
data_inst <- data[data$feature == "instability", ]

boxplot(data_gc$cds, data_gc$trg, data_gc$denovo,
        names = c("CDS", "TRG", "De novo"),
        xlab = "Sequence type",
        ylab = "GC content (%)",
        main = "mean GC content distribution",
        col = c("#cc7f0a", "#ad4646", "#302f2f"))


boxplot(data_aro$cds, data_aro$trg, data_aro$denovo,
        names = c("CDS", "TRG", "De novo"),
        xlab = "Sequence type",
        ylab = "Aromaticity",
        main = "mean aromaticity distribution",
        col = c("#cc7f0a", "#ad4646", "#302f2f"))


boxplot(data_inst$cds, data_inst$trg, data_inst$denovo,
        names = c("CDS", "TRG", "De novo"),
        xlab = "Sequence type",
        ylab = "Instability index",
        main = "mean instability index distribution",
        col = c("#cc7f0a", "#ad4646", "#302f2f"))
