import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

INPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/superkingdom_level.txt"

# Create empty pandas dataframe
rows = []
with open(INPUT, "r") as f:
    for line in f:
        line_split = line.strip().split("\t")
        query = line_split[0] 
        counts = line_split[1].split(";")
        for count in counts:
            kingdom, count = count.split(":")
            # Collect rows
            rows.append({"query": query, "class": kingdom, "count": int(count)})

# Create dataframe from collected rows
data = pd.DataFrame(rows, columns=["query", "class", "count"])


# Piechart of the distribution of classes
global_counts = data.groupby("class")["count"].sum().reset_index()
fig, ax = plt.subplots()
colors = plt.cm.viridis(np.linspace(0, 1, len(global_counts["class"])))
wedges, texts = ax.pie(global_counts["count"], textprops=dict(color="w"), colors=colors)
ax.legend(wedges, global_counts["class"], title="Classes", loc="center left", bbox_to_anchor=(0.5, -0.4, 0.5, 1))
ax.axis('equal')
plt.title("Distribution of superkingdoms among matches\nfor 1000 genes identified as common at the superkingdom level")


# Compute the proportion of each class for each query
data["proportion"] = data.groupby("query")["count"].apply(lambda x: x / x.sum()).reset_index(level=0, drop=True)
# Boxplot of the distribution of proportions per query
proportion_distrib = data.groupby("class")["proportion"].apply(list).reset_index()
fig, ax = plt.subplots()
box = ax.boxplot(proportion_distrib["proportion"], tick_labels=proportion_distrib["class"], patch_artist=True)

# Color each boxplot
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

plt.xlabel("Kingdom")
plt.ylabel("Proportion")
plt.title("Distribution of the proportion of matches per class\nfor 1000 genes identified as common at the superkingdom level")
#plt.show()


# Get the genes for which the Archaea matches are high in proportion
archaea_prop_90 = data[(data["class"] == "Archaea") & (data["proportion"] > 0.9)]
archaea_prop_95 = data[(data["class"] == "Archaea") & (data["proportion"] > 0.95)]
archaea_prop_99 = data[(data["class"] == "Archaea") & (data["proportion"] > 0.99)]
archaea_prop_999 = data[(data["class"] == "Archaea") & (data["proportion"] > 0.999)]
print(f"There are {len(archaea_prop_90)} genes for which the proportion of Archaea matches is > 0.9:")
print(f"There are {len(archaea_prop_95)} genes for which the proportion of Archaea matches is > 0.95:")
print(f"There are {len(archaea_prop_99)} genes for which the proportion of Archaea matches is > 0.99:")


# Plot the distribution of matches per gene for the genes with high proportion of Archaea matches
count_distrib_90 = data[data["query"].isin(archaea_prop_90["query"])].groupby("class")["count"].apply(list).reset_index()
fig, ax = plt.subplots()
box = ax.boxplot(count_distrib_90["count"], tick_labels=count_distrib_90["class"], patch_artist=True)
# Color each boxplot
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
    plt.yscale("log")
plt.xlabel("Class")
plt.ylabel("Count")
plt.title(f"Distribution of the number of matches per class\nfor {len(archaea_prop_90)} genes with > 90% of Archaea matches")


count_distrib_95 = data[data["query"].isin(archaea_prop_95["query"])].groupby("class")["count"].apply(list).reset_index()
fig, ax = plt.subplots()
box = ax.boxplot(count_distrib_95["count"], tick_labels=count_distrib_95["class"], patch_artist=True)
# Color each boxplot
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.yscale("log")
plt.xlabel("Class")
plt.ylabel("Count")
plt.title(f"Distribution of the number of matches per class\nfor {len(archaea_prop_95)} genes with > 95% of Archaea matches")


count_distrib_99 = data[data["query"].isin(archaea_prop_99["query"])].groupby("class")["count"].apply(list).reset_index()
fig, ax = plt.subplots()
box = ax.boxplot(count_distrib_99["count"], tick_labels=count_distrib_99["class"], patch_artist=True)
# Color each boxplot
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)
plt.yscale("log")
plt.xlabel("Class")
plt.ylabel("Count")
plt.title(f"Distribution of the number of matches per class\nfor {len(archaea_prop_99)} genes with > 99% of Archaea matches")
#plt.show()


# Genes for which we only have Archaea matches
archaea_only = data[data["class"] == "Archaea"]
archaea_only = archaea_only[archaea_only["proportion"] == 1]
print(f"\nThere are {len(archaea_only)} genes for which we only have Archaea matches")

# Export the gene names for which we have > 99% of Archaea matches
with open("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/archaea_prop_99_genes.txt", "w") as f:
    for gene in archaea_prop_99["query"].unique():
        f.write(f"{gene}\n")