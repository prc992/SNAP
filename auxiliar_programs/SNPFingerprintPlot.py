import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

def load_pvalue_matrix(file_path):
    """
    Reads a tab-delimited file containing a matrix of p-values.
    Converts p-values to -log10(p) for better visualization in a heatmap.
    Removes the suffix "_unique.sorted.dedup.bam" from row and column names.
    Ignores the first line if it starts with '#'.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    if lines[0].startswith("#"):
        lines = lines[1:]
    
    with open(file_path, 'w') as file:
        file.writelines(lines)
    
    df = pd.read_csv(file_path, sep="\t", index_col=0, dtype=str)
    
    # Remove specific suffixes from row and column names
    df.index = df.index.str.replace("_unique.sorted.dedup.bam", "", regex=False)
    df.columns = df.columns.str.replace("_unique.sorted.dedup.bam", "", regex=False)
    
    df.index = df.index.str.replace(".dac_filtered.dedup.unique.sorted.bam", "", regex=False)
    df.columns = df.columns.str.replace(".dac_filtered.dedup.unique.sorted.bam", "", regex=False)
    
    # Convert values to float
    df_float = df.astype(float)
    
    return df, df_float

def plot_pvalue_heatmap_sorted_dend(df_original, output_path="heatmap.jpg"):
    """
    Plots a heatmap of transformed p-values with hierarchical clustering dendrograms.
    Saves the heatmap as a .jpg file instead of displaying it.
    
    Parameters:
        df_original (pd.DataFrame): The transformed p-value matrix.
        output_path (str): File path to save the heatmap image.
    """
    row_dist = pdist(df_original, metric='correlation')
    col_dist = pdist(df_original.T, metric='correlation')
    
    row_linkage = linkage(row_dist, method='average')
    col_linkage = linkage(col_dist, method='average')
    size = max(df_original.shape[0]/2,10)
    figsize=(size, size)
    
    g = sns.clustermap(df_original, cmap="coolwarm", linewidths=0.5, annot=False, fmt=".1e",
                       cbar_kws={'label': '-log10 (p-value)'}, row_linkage=row_linkage, col_linkage=col_linkage,
                       figsize=figsize)
    
    cbar = g.ax_heatmap.collections[0].colorbar
    cbar.set_label('-log10(p-value)', rotation=270, labelpad=15)
    
    #g.ax_heatmap.set_title("P-value Similarity Heatmap with Clustering", fontsize=14, pad=20)
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate a clustered heatmap from a p-value matrix")
    parser.add_argument("-i", "--input", type=str, default="pval_out.txt", help="Input file path (default: pval_out.txt)")
    parser.add_argument("-o", "--output", type=str, default="heatmap.jpg", help="Output image file path (default: heatmap.jpg)")

    args = parser.parse_args()
    
    df_original, df_float = load_pvalue_matrix(args.input)
    plot_pvalue_heatmap_sorted_dend(df_float, args.output)

if __name__ == "__main__":
    main()