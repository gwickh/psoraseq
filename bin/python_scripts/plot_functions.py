#!/usr/bin/env python3
# Filename:         plot_functions.py
# Author:           Gregory Wickham
# Created:          2024-11-25
# Version:          1.0
# Date modified:    2024-11-25
# Description:      Visualise read log2FC and Jemsen-Shannon Divergence

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.special import rel_entr

def plot_log2fc(
    normalised_file_path, 
    output_path="log2fc_by_pos.png",
    smoothing_window=200
    ):
    """
    Generates plot showing log2 fold change by position from log2FC_binned_reads.csv.
    """
    
    df = pd.read_csv(normalised_file_path)     # Read the CSV file
    bins = df['Bins']    # Extract the 'bin' column for the x-axis

    plt.figure(figsize=(12, 6))  # Set figure size
    for column in df.columns[1:]:
        smoothed_values = df[column].rolling(window=smoothing_window, center=True).mean()
        plt.plot(bins, smoothed_values, label=f"{column}")

    # Add labels, title, and legend
    plt.xlabel('Position (Mb)')
    plt.ylabel('Log2 Fold Change')
    plt.title('Log2 Fold Change in Reads from Global Mean')
    plt.legend(title='Columns')
    plt.grid(True, linestyle='--', alpha=0.6)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Smoothed plot saved to {output_path}")
    plt.show()
    
def plot_JS_div(
    raw_file_path,
    output_path="jensen_shannon_divergence.png",
    ):
    """
    Normalise data, calculate JS divergence and plot data from binned_reads.csv.
    """

    # Function to calculate Jensen-Shannon divergence
    def jensen_shannon_divergence(P, Q):
        # M is the average of the two distributions
        M = 0.5 * (P + Q)
        
        # KL divergence of P from M and Q from M
        kl_p_m = np.sum(rel_entr(P, M))
        kl_q_m = np.sum(rel_entr(Q, M))
        
        # Jensen-Shannon divergence
        return 0.5 * (kl_p_m + kl_q_m)

    # Load and normalize the data
    df = pd.read_csv(raw_file_path)

    counts_data = df.drop(columns=['Bins'])
    normalized_data = counts_data/ counts_data.sum(axis=0)
    col_names = list(counts_data.columns.values)
    

    # Initialize an empty similarity matrix
    n_samples = normalized_data.shape[1]
    js_matrix = np.zeros((n_samples, n_samples))

    # Calculate the JS_D between each pair of samples
    for i in range(n_samples):
        for j in range(i + 1, n_samples):  # Calculate for upper triangle only (symmetric)
            p = normalized_data.iloc[i].values
            q = normalized_data.iloc[j].values
            js_div_value = jensen_shannon_divergence(p, q)
            js_matrix[i, j] = js_matrix[j, i] = js_div_value  # Matrix is symmetric

    # Create a DataFrame for better visualization
    js_matrix_df = pd.DataFrame(js_matrix, index=col_names, columns=col_names)

    # Generate the heatmap
    plt.figure(figsize=(10, 8))  # Adjust size to fit your data
    sns.heatmap(
        js_matrix_df, 
        annot=True, 
        cmap="coolwarm", 
        linewidths=0.5, 
        vmin=0, 
        vmax=np.log(2)
    )
    plt.title("Jensen-Shannon Divergence Between Samples")
    plt.xlabel("Sample")
    plt.ylabel("Sample")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, ha='right')
    
    # Save and show the plot
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"Jensen-Shannon divergence matrix saved to {output_path}")
    plt.show()

# Generate boxplot showing distribution of raw reads
def raw_read_boxplot(raw_file_path):
    raw_df = pd.read_csv(raw_file_path)
    raw_df = pd.read_csv("~/webber_group/Gregory_Wickham/psoraseq/bowtie2/output_test/binned_reads.csv")
    raw_df = raw_df.melt(id_vars=['Bins'], var_name='sample', value_name='read_count')
    sns.boxplot(
        raw_df,
        y = 'read_count',
        x = 'sample',
        inner = None,
        linecolor = 'black',
        color = 'white', 
        width = 0.5
    )
    sns.violinplot(
        raw_df,
        y = 'read_count',
        x = 'sample',
        inner = None,
        linecolor = 'black',
        color = 'white', 
        linewidth = 1
    )
    plt.title("Raw Count of Reads per Bin")
    plt.xlabel("Sample")
    plt.ylabel("Read Count")

# Generate boxplot showing distribution of normalised reads
def raw_read_boxplot(raw_file_path):
    normalised_df = pd.read_csv(raw_file_path)
    normalised_df = pd.read_csv("~/webber_group/Gregory_Wickham/psoraseq/bowtie2/output_test/log2FC_binned_reads.csv")
    normalised_df = normalised_df.melt(id_vars=['Bins'], var_name='sample', value_name='read_count')
    sns.boxplot(
        normalised_df,
        y = 'read_count',
        x = 'sample',
        linecolor = 'black',
        color = 'white', 
        width = 0.5
    )
    sns.violinplot(
        normalised_df,
        y = 'read_count',
        x = 'sample',
        inner = None,
        linecolor = 'black',
        color = 'white', 
        linewidth = 1
    )
    plt.title("Normalised log 2 Fold Change per Bin")
    plt.xlabel("Sample")
    plt.ylabel("log 2 Fold Change")

# Generate plots
plot_log2fc("~/webber_group/Gregory_Wickham/psoraseq/bowtie2/output_test/log2FC_binned_reads.csv")
plot_JS_div("~/webber_group/Gregory_Wickham/psoraseq/bowtie2/output_test/binned_reads.csv")

