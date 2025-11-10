# ====================================================
# 📊 HackBio Stage 1 – Figures 1A–F
# Author: Ayşenur Akcan
# Description: Reproduction of figures using Pandas, Seaborn, and Matplotlib
# =====================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from io import StringIO

sns.set_theme(style="whitegrid", font_scale=1)

# -----------------------------------------------------
# 🔹 Figure a – Clustered Heatmap of Gene Expression
# -----------------------------------------------------
data = """gene,HBR_1,HBR_2,HBR_3,UHR_1,UHR_2,UHR_3
SULT4A1,375,343.6,339.4,3.5,6.9,2.6
MPPED1,157.8,158.4,162.6,0.7,3,2.6
PRAME,0,0,0,568.9,467.3,519.2
IGLC2,0,0,0,488.6,498,457.5
IGLC3,0,0,0,809.7,313.8,688
CDC45,2.6,1,0,155,152.5,149.9
CLDN5,77.6,88.5,67.2,1.4,2,0
PCAT14,0,0,1.2,139.8,154.4,155.1
RP5-1119A7.17,53,57.6,51.9,0,0,0
MYO18B,0,0,0,59.5,84.2,56.5
RP3-323A16.1,0,0,1.2,51.9,76.2,53.1
CACNG2,42.7,35,56.6,0,1,0
"""
df = pd.read_csv(StringIO(data), index_col=0)

sns.clustermap(
    df,
    cmap="Blues",
    linewidths=0.5,
    figsize=(6,5),
    cbar_kws={"label": "Normalized Expression"}
).fig.suptitle("a. Clustered Heatmap of Differentially Expressed Genes", y=1.05)
plt.show()

# -----------------------------------------------------
# 🔹 Figure b – Volcano Plot (Padj version)
# -----------------------------------------------------
def generate_volcano_plot(data_url, lfc_thresh=1):
    """
    Volcano plot: log2FoldChange vs -log10(Padj)
    """
    df = pd.read_csv(data_url)
    df["Padj"] = df["Padj"].replace(0, 1e-300)
    df["-log10Padj"] = -np.log10(df["Padj"])

    color_map = {"up": "green", "down": "orange", "ns": "grey"}
    label_map = {"up": "Upregulated", "down": "Downregulated", "ns": "Not Significant"}

    plt.figure(figsize=(8,6))
    for cat, color in color_map.items():
        sub = df[df["significance"] == cat]
        plt.scatter(
            sub["log2FoldChange"],
            sub["-log10Padj"],
            s=30, alpha=0.8, color=color, label=label_map[cat]
        )
    plt.axvline(x=lfc_thresh, color="black", linestyle="--", linewidth=1.2)
    plt.axvline(x=-lfc_thresh, color="black", linestyle="--", linewidth=1.2)
    plt.xlabel("log2FoldChange")
    plt.ylabel("-log10(Padj)")
    plt.title("b. Volcano Plot of Differential Gene Expression")
    plt.legend(title="Significance")
    plt.tight_layout()
    plt.show()

generate_volcano_plot(
    "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
)

# -----------------------------------------------------
# 🔹 Figure c – Scatter Plot (radius vs texture)
# -----------------------------------------------------
def plot_radius_vs_texture(data_url):
    df = pd.read_csv(data_url)
    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=df, x="radius_mean", y="texture_mean",
        hue="diagnosis", palette={"M":"C0", "B":"C1"}, alpha=0.8
    )
    plt.title("c. Texture Mean vs Radius Mean by Diagnosis")
    plt.xlabel("radius_mean")
    plt.ylabel("texture_mean")
    plt.legend(title="Diagnosis")
    plt.show()

plot_radius_vs_texture(
    "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
)

# -----------------------------------------------------
# 🔹 Figure d – Correlation Heatmap
# -----------------------------------------------------
def correlation_heatmap(data_url):
    df = pd.read_csv(data_url)
    features = [
        "radius_mean","texture_mean","perimeter_mean",
        "area_mean","smoothness_mean","compactness_mean"
    ]
    corr = df[features].corr()
    plt.figure(figsize=(8,6))
    sns.heatmap(
        corr, annot=True, cmap="Blues",
        fmt=".1f", linewidths=0.5,
        cbar_kws={"label":"Correlation Coefficient"}
    )
    plt.title("d. Correlation Heatmap of Key Mean Features")
    plt.show()

correlation_heatmap(
    "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
)

# -----------------------------------------------------
# 🔹 Figure e – Scatter (smoothness vs compactness)
# -----------------------------------------------------
def plot_smoothness_vs_compactness(data_url):
    df = pd.read_csv(data_url)
    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=df, x="smoothness_mean", y="compactness_mean",
        hue="diagnosis", palette={"M":"C0","B":"C1"}, alpha=0.8
    )
    plt.title("e. Compactness Mean vs Smoothness Mean by Diagnosis")
    plt.xlabel("smoothness_mean")
    plt.ylabel("compactness_mean")
    plt.legend(title="Diagnosis")
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.show()

plot_smoothness_vs_compactness(
    "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
)

# -----------------------------------------------------
# 🔹 Figure f – Density Plot (area distribution)
# -----------------------------------------------------
def plot_area_density(data_url):
    df = pd.read_csv(data_url)
    plt.figure(figsize=(8,6))
    sns.kdeplot(
        data=df, x="area_mean", hue="diagnosis",
        fill=True, alpha=0.25, linewidth=2,
        palette={"M":"C0","B":"C1"}
    )
    plt.title("f. Density Plot of Area Mean by Diagnosis")
    plt.xlabel("area_mean")
    plt.ylabel("Density")
    plt.legend(title="Diagnosis")
    plt.tight_layout()
    plt.show()

plot_area_density(
    "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
)
