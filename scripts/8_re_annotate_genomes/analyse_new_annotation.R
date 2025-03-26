library(dplyr)

NEW_GFF <- "/home/eliott.tempez/Documents/archaea_data/complete_122/reannotated_gff/"
OLD_GFF <- "/home/eliott.tempez/Documents/archaea_data/complete_122/gff3_no_fasta/"


# List of files 
files_new <- list.files(NEW_GFF, pattern = "gff3$", full.names = TRUE)
files_old <- list.files(OLD_GFF, pattern = "gff3$", full.names = TRUE)

# Get the CDS lengths the old and new way
n_files_new <- 0
lengths_new <- c()
# For each file
for (file in files_new) {
    # Only files that don't start with GCA
    if (grepl("GCA", file)) {
        next
    }
    n_files_new <- n_files_new + 1
    data <- read.csv(file, header = FALSE, sep = "\t", comment.char = "#")
    # Get rid of the last column
    data <- data[, -ncol(data)]
    # Get the CDS
    cds <- data[data$V3 == "CDS", ]
    # Get the length of each CDS
    cds$length <- cds$V5 - cds$V4 + 1
    lengths_new <- c(lengths_new, cds$length)
}
lengths_old <- c()
n_files_old <- 0
# For each file
for (file in files_old) {
    # Only files that don't start with GCA
    if (grepl("GCA", file)) {
        next
    }
    n_files_old <- n_files_old + 1
    data <- read.csv(file, header = FALSE, sep = "\t", comment.char = "#")
    # Get the ID from the last column
    data$ID <- gsub("ID=", "", gsub(";.*", "", data$V9))
    # Get rid of the last column
    data <- data[, -9]
    # Get the CDS
    cds <- data[data$V3 == "CDS", ]
    # Remove all CDS with duplicated IDs, including all occurrences of these IDs
    duplicated_ids <- cds$ID[duplicated(cds$ID) | duplicated(cds$ID, fromLast = TRUE)]
    cds <- cds[!cds$ID %in% duplicated_ids, ]
    # Get the length of each CDS
    cds$length <- cds$V5 - cds$V4 + 1
    # Print the genes < 89
    lengths_old <- c(lengths_old, cds$length)
}


print(paste("Number of files new:", n_files_new))
print(paste("Number of files old:", n_files_old))


#hist(lengths_new[lengths_new < 300], breaks = round(length(lengths_new)/10000), main = "Length distribution of CDS\nannotated with prokka", xlab = "Length of CDS", ylab = "Number of CDS", ylim = c(0, 4000), xlim = c(0, 300))
#hist(lengths_old[lengths_old < 300], breaks = round(length(lengths_old)/10000), main = "Length distribution of CDS\nannotated with microscope", xlab = "Length of CDS", ylab = "Number of CDS", ylim = c(0, 4000), xlim = c(0, 300))

print("Old annotation summary:")
print(summary(lengths_old))
print("New annotation summary:")
print(summary(lengths_new))



############ Prokka with min length 75 ############
NEW_GFF_75 <- "/home/eliott.tempez/Documents/archaea_data/complete_122/reannotated_gff_75/"
files_new_75 <- list.files(NEW_GFF_75, pattern = "gff3$", full.names = TRUE)
for (file_new_75 in files_new_75) {
    file_new <- gsub("reannotated_gff_75", "reannotated_gff", file_new_75)
    data_new_75 <- read.csv(file_new_75, header = FALSE, sep = "\t", comment.char = "#")
    data_new <- read.csv(file_new, header = FALSE, sep = "\t", comment.char = "#")
    # Get the ID from the last column
    data_new_75$ID <- gsub("ID=", "", gsub(";.*", "", data_new_75$V9))
    data_new$ID <- gsub("ID=", "", gsub(";.*", "", data_new$V9))
    # Get rid of the last column
    data_new_75 <- data_new_75[, -9]
    data_new <- data_new[, -9]
    # Get the CDS
    cds_new_75 <- data_new_75[data_new_75$V3 == "CDS", ]
    cds_new <- data_new[data_new$V3 == "CDS", ]
    # Get the CDS we gained
    cds_gained <- anti_join(cds_new_75, cds_new, by = c("V1", "V4", "V5"))
    cds_gained$length <- cds_gained$V5 - cds_gained$V4 + 1
    print(file_new_75)
    print("New CDS gained:")
    print(cds_gained)
    print("CDS gained <90 nucl:")
    print(cds_gained[cds_gained$length < 90, ])
}
