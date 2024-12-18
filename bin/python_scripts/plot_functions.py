#!/usr/bin/env python3
# Filename:         plot_functions.py
# Author:           Gregory Wickham
# Created:          2024-11-25
# Version:          1.1
# Date modified:    2024-11-25
# Description:      Visualise read log2FC and Jensen-Shannon Divergence

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.special import rel_entr

def plot_log2fc(normalised_file_path, smoothing_window=100):
    """
    Generates plot showing log2 fold change by position from log2FC_binned_reads.csv.
    """
    df = pd.read_csv(normalised_file_path)  # Read the CSV file
    bins = df['Bins']  # Extract the 'bin' column for the x-axis

    for column in df.columns[1:]:
        # Wrap data by appending and prepending data for smoothing
        padded_data = np.concatenate(
            (df[column].iloc[-smoothing_window:].values, 
            df[column].values, 
            df[column].iloc[:smoothing_window].values)
        )
        smoothed_padded_values = pd.Series(padded_data)\
            .rolling(
                window = smoothing_window, 
                center = True)\
            .mean()
        smoothed_values = smoothed_padded_values[smoothing_window:-smoothing_window].values

        plt.plot(bins, smoothed_values, label=f"{column}")

    # Add labels, title, and legend
    plt.xlabel('Position (Mb)')
    plt.ylabel('Log2 Fold Change')
    plt.title('Log2 Fold Change in Reads from Global Mean')
    plt.legend(title='Columns')
    plt.grid(True, linestyle='--', alpha=0.6)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig("log2fc_by_pos.png")
    print("Smoothed plot saved to log2fc_by_pos.png")
    plt.show()


def plot_JS_div(raw_file_path):
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
    normalized_data = counts_data / counts_data.sum(axis=0)
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
    plt.savefig("jensen_shannon_divergence.png")
    print("Jensen-Shannon divergence matrix saved to jensen_shannon_divergence.png")
    plt.show()


def raw_read_boxplot(raw_file_path):
    """
    Generate a boxplot showing the distribution of raw reads.
    """
    raw_df = pd.read_csv(raw_file_path)
    raw_df = raw_df.melt(id_vars=['Bins'], var_name='sample', value_name='read_count')
    plt.figure(figsize=(12, 6))
    sns.boxplot(
        data = raw_df,
        y = 'read_count',
        x = 'sample',
        color = 'white',
        width = 0.5,
        linewidth = 0.8
    )
    sns.violinplot(
        data = raw_df,
        y = 'read_count',
        x = 'sample',
        color = 'white',
        linewidth = 1,
        inner = None
    )
    plt.title("Raw Count of Reads per Bin")
    plt.xlabel("Sample")
    plt.ylabel("Read Count")
    plt.tight_layout()
    plt.savefig("raw_read_boxplot.png")
    print("Boxplot saved to raw_read_boxplot.png")
    plt.show()


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Visualize read log2FC and Jensen-Shannon Divergence.")
    parser.add_argument(
        "raw_file_path", 
        help = "Path to the binned_reads.csv")
    parser.add_argument(
        "normalised_file_path", 
        help = "Path to log2FC_binned_reads.csv")
    parser.add_argument(
        "--smoothing_window", 
        type=int, 
        default=1, 
        help="Smoothing window size for log2FC plot."
        )

    # Parse arguments
    args = parser.parse_args()

    # Generate plots
    plot_log2fc(
        args.normalised_file_path, 
        smoothing_window=args.smoothing_window
    )
    plot_JS_div(args.raw_file_path)
    raw_read_boxplot(args.raw_file_path)


if __name__ == "__main__":
    main()
