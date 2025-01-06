#!/usr/bin/env python3
# Filename:         plot_functions.py
# Author:           Gregory Wickham
# Created:          2024-11-25
# Version:          1.2
# Date modified:    2025-01-06
# Description:      Visualise read binsizes (absolute and relative) and Jensen-Shannon Divergence

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.special import rel_entr

def plot_log2fc(normalised_file_path, smoothing_window=100):
    """
    Takes normalised readcounts and plots a moving average of the log2 fold change by position.
    """
    df = pd.read_csv(normalised_file_path)
    bins = df['Bins']  # Extract the 'bin' column for the x-axis

    plt.figure(figsize=(10, 6))
    for column in df.columns[1:]:
        # Wrap data by appending and prepending data for smoothing
        padded_data = np.concatenate(
            (df[column].iloc[-smoothing_window:].values, 
            df[column].values, 
            df[column].iloc[:smoothing_window].values)
        )
        padded_data = pd.DataFrame(padded_data)
        
        # Replace -inf values with NaN and interpolate to ignore them during smoothing
        padded_data = padded_data.replace(-np.inf, np.nan)
        interpolated_padded_values = padded_data.interpolate(
            method='spline', 
            order=3
        ).rolling(
            window=smoothing_window, 
            center=True
        ).mean()
        smoothed_values = interpolated_padded_values[smoothing_window:-smoothing_window].values
        plt.plot(bins, smoothed_values, label=f"{column}")
    
    # Plot aesthetics
    plt.xlabel('Position (Mb)')
    plt.ylabel('Log2 Fold Change')
    plt.title('Log2 Fold Change in Reads from Sample Mean')
    plt.legend(title='Sample', bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig("log2fc_by_pos.png")
    plt.close()

def relative_log2FC(raw_file_path, smoothing_window=100):
    """
    Takes raw readcounts and calculates readcount relative to control and plots a moving 
    average of the log2 fold change by position.
    """
    # Load and normalize the data
    df = pd.read_csv(raw_file_path)
    bins = df['Bins']
    counts_data = df.drop(columns=['Bins'])
    normalized_data = counts_data / counts_data.sum(axis=0)

    # Identify "test" sample and matching "control"
    test_columns = [col for col in normalized_data.columns if "test" in col]
    control_columns = {col.replace("test", "control"): col for col in test_columns}

    relative_data = pd.DataFrame()
    for control_col, test_col in control_columns.items():
        if control_col in normalized_data.columns:
            # Compute the ratio and store in a new column
            relative_data[test_col] = normalized_data[test_col] / normalized_data[control_col]
    relative_df = pd.concat([bins, relative_data], axis=1)
    relative_df.to_csv("relative_change.csv", index=False)
    
    if len(relative_df.columns) > 1: # Plot only if there are more than the Bins column
        plt.figure(figsize=(10, 6))
        for column in relative_df.columns[1:]:
            # Wrap data by appending and prepending data for smoothing
            padded_data = np.concatenate(
                (relative_df[column].iloc[-smoothing_window:].values, 
                relative_df[column].values, 
                relative_df[column].iloc[:smoothing_window].values)
            )
            padded_data = pd.DataFrame(padded_data)
            
            # Replace -inf values with NaN and interpolate to ignore them during smoothing
            padded_data = padded_data.replace([np.inf, -np.inf, 0, ''], np.nan)
            interpolated_padded_values = padded_data.interpolate(
                method='spline', 
                order=3
            ).rolling(
                window=smoothing_window, 
                center=True
            ).mean()
            smoothed_values = interpolated_padded_values[smoothing_window:-smoothing_window].values
            log2fc_values = np.log2(smoothed_values)
            plt.plot(bins, log2fc_values, label=f"{column}")
            
        # Plot aesthetics
        plt.xlabel('Position (Mb)')
        plt.ylabel('Relative Log2 Fold Change')
        plt.title('Change in Normalised Reads from Control')
        plt.legend(title='Sample', bbox_to_anchor=(1.04, 1), loc="upper left")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        plt.savefig("relative_change.png")
        plt.close()
    else:
        return
    
def plot_JS_div(raw_file_path):
    """
    Takes raw readcounts, normalises them and calculates the Jensen-Shannon divergence 
    between samples.
    """

    # Function to calculate Jensen-Shannon divergence as D_JS = [D_KL(P||M)/2 + D_KL(Q||M)]/2
    def jensen_shannon_divergence(P, Q):
        # Calculate the average distribution M
        M = 0.5 * (P + Q)

        # Calulate KL divergence of P from M and Q from M
        kl_p_m = np.sum(rel_entr(P, M))
        kl_q_m = np.sum(rel_entr(Q, M))

        # Return Jensen-Shannon divergence
        return 0.5 * (kl_p_m + kl_q_m)

    # Load and normalize the data
    df = pd.read_csv(raw_file_path)
    counts_data = df.drop(columns=['Bins'])
    normalized_data = counts_data / counts_data.sum(axis=0)
    col_names = list(counts_data.columns.values)

    # Initialize empty similarity matrix
    n_samples = normalized_data.shape[1]
    js_matrix = np.zeros((n_samples, n_samples))

    # Calculate the JS_D between each pair of samples
    for i in range(n_samples):
        for j in range(i + 1, n_samples):  # Calculate for upper triangle only
            p = normalized_data.iloc[:, i].values
            q = normalized_data.iloc[:, j].values
            js_div_value = jensen_shannon_divergence(p, q)
            js_matrix[i, j] = js_matrix[j, i] = js_div_value  # Matrix is symmetric

    js_matrix_df = pd.DataFrame(js_matrix, index=col_names, columns=col_names)

    # Generate the heatmap
    plt.figure(figsize=(10, 8)) 
    sns.heatmap(
        js_matrix_df,
        annot=True,
        cmap="coolwarm",
        linewidths=0.5
    )
    # Plot aesthetics
    plt.title("Jensen-Shannon Divergence Between Samples")
    plt.xlabel("Sample")
    plt.ylabel("Sample")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0, ha='right')
    plt.tight_layout()
    plt.savefig("jensen_shannon_divergence.png")
    plt.close()

def raw_read_boxplot(raw_file_path):
    """
    Generates a boxplot showing the distribution of readcount bin sizes.
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
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.title("Raw Count of Reads per Bin")
    plt.xlabel("Sample")
    plt.ylabel("Read Count")
    plt.savefig("raw_read_boxplot.png")
    plt.close()
        
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
    args = parser.parse_args()

    # Generate plots
    plot_log2fc(
        args.normalised_file_path, 
        smoothing_window=args.smoothing_window
    )
    relative_log2FC(
        args.raw_file_path, 
        smoothing_window=args.smoothing_window
    )
    plot_JS_div(args.raw_file_path)
    raw_read_boxplot(args.raw_file_path)

if __name__ == "__main__":
    main()
