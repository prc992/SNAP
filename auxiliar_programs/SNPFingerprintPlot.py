import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from decimal import Decimal, getcontext

def load_pvalue_matrix(file_path):
    """
    Reads a tab-delimited file containing a matrix of p-values.
    Converts p-values to -log10(p) for better visualization in a heatmap.
    Removes the suffix "_unique.sorted.dedup.bam" from row and column names.
    Ignores the first line if it starts with '#'.
    """
    # Set high precision for Decimal calculations
    getcontext().prec = 100  # 100 digits precision
    
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
    
    # Convert values using high-precision Decimal arithmetic
    def precise_log10_conversion(val_str):
        """Convert p-value string to -log10(p) using maximum precision"""
        try:
            val_str = str(val_str).strip()
            
            # Use Decimal for maximum precision
            decimal_val = Decimal(val_str)
            
            if decimal_val <= 0:
                return 1000.0  # Large value for zero/negative
            
            # Calculate -log10 with high precision
            log10_val = decimal_val.log10()
            result = float(-log10_val)
            
            return result
            
        except Exception as e:
            print(f"Precision conversion error for '{val_str}': {e}")
            return 1000.0  # Fallback
    
    # Test precision on sample values
    print("Testing high-precision conversion:")
    for i in range(min(3, df.shape[0])):
        original_val = df.iloc[i, i]
        converted_val = precise_log10_conversion(original_val)
        sample_name = df.index[i]
        print(f"  {sample_name}: '{original_val}' -> {converted_val:.6f}")
    
    # Apply the precise conversion to each cell
    df_float = df.applymap(precise_log10_conversion)
    
    print(f"High-precision matrix: {df_float.shape[0]} samples, range {df_float.min().min():.3f}-{df_float.max().max():.3f}")
    print(f"Standard deviation: {df_float.std().std():.6f}")
    
    # Show sample of converted matrix with high precision
    print("\nSample of converted matrix (high precision):")
    print(df_float.iloc[:3, :3].round(3))
    
    return df, df_float

def plot_pvalue_heatmap_sorted_dend(df_original, output_path="heatmap.jpg"):
    """
    Plots a heatmap of transformed p-values with hierarchical clustering dendrograms.
    Saves the heatmap as a .jpg file instead of displaying it.
    
    Parameters:
        df_original (pd.DataFrame): The transformed p-value matrix.
        output_path (str): File path to save the heatmap image.
    """
    # Check for clustering feasibility
    if not np.isfinite(df_original.values).all():
        print("Warning: Non-finite values detected in matrix")
        return
    
    row_dist = pdist(df_original, metric='correlation')
    col_dist = pdist(df_original.T, metric='correlation')
    
    if not np.isfinite(row_dist).all() or not np.isfinite(col_dist).all():
        print("Error: Distance calculation failed, cannot perform clustering")
        return
    
    row_linkage = linkage(row_dist, method='average')
    col_linkage = linkage(col_dist, method='average')
    size = max(df_original.shape[0]/2,13)
    figsize=(size, size)
    
    g = sns.clustermap(df_original, cmap="coolwarm", linewidths=0.5, annot=False, fmt=".1e",
                       cbar_kws={'label': '-log10 (p-value)'}, row_linkage=row_linkage, col_linkage=col_linkage,
                       figsize=figsize)
    
    cbar = g.ax_heatmap.collections[0].colorbar
    cbar.set_label('-log10(p-value)', rotation=270, labelpad=15)
    
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate a clustered heatmap from a p-value matrix")
    parser.add_argument("-i", "--input", type=str, default="pval_out.txt", help="Input file path (default: pval_out.txt)")
    parser.add_argument("-o", "--output", type=str, default="heatmap.jpg", help="Output image file path (default: heatmap.jpg)")

    args = parser.parse_args()
    
    df_original, df_float = load_pvalue_matrix(args.input)
    plot_pvalue_heatmap_sorted_dend(df_float, args.output)

if __name__ == "__main__":
    main()