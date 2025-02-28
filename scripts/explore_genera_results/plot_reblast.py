import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

INPUT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/superkingdom_level.txt"

data = pd.read_csv(INPUT, sep="\t", header=None)
data.columns = ["query", "class", "count"]


# Piechart of the distribution of classes
global_counts = data.groupby("class")["count"].sum().reset_index()
fig, ax = plt.subplots()
def func(pct, allvals):
    """Function to display the percentage and the count"""
    absolute = int(pct/100.*sum(allvals))
    return f"{pct:.1f}% ({absolute:d})"
colors = plt.cm.viridis(np.linspace(0, 1, len(global_counts["class"])))
wedges, texts, autotexts = ax.pie(global_counts["count"], autopct=lambda pct: func(pct, global_counts["count"]), textprops=dict(color="w"), colors=colors)
ax.legend(wedges, global_counts["class"], title="Classes", loc="center left", bbox_to_anchor=(0.5, -0.5, 0.5, 1))
ax.axis('equal')
plt.title("Distribution of superkingdoms among matches\nfor 1000 genes identified as common at the superkingdom level")
#plt.show()


# Boxplot of the distribution of counts per query
count_distrib = data.groupby("class")["count"].apply(list).reset_index()
fig, ax = plt.subplots()
box = ax.boxplot(count_distrib["count"], labels=count_distrib["class"], patch_artist=True)

# Color each boxplot
for patch, color in zip(box['boxes'], colors):
    patch.set_facecolor(color)

plt.xlabel("Class")
plt.ylabel("Count")
plt.title("Distribution of the number of matches per class\nfor 1000 genes identified as common at the superkingdom level")
plt.show()
