"""
Filename:         plot_functions.py
Author:           Gregory Wickham
Created:          2024-11-25
Version:          1.3
Date modified:    2025-01-10
Description:      Visualise binned read abundances and plot Jensen-Shannon Divergence between samples.
"""

import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.special import rel_entr
from datetime import datetime

class DataProcessing:
    """"
    Class to process data for plotting by interpolating missing values and smoothing with moving average.
    """
    
    def __init__(self, file_path, smoothing_window=100):
        self.file_path = file_path
        self.smoothing_window = smoothing_window


    def interpolate_missing(self, df: pd.DataFrame) -> pd.DataFrame:
        """Interpolte NaN, infinity, and zero values with cubic spline."""
        df = df.replace([np.inf, -np.inf, 0, ''], np.nan)
        return df.interpolate(method='spline', order=3)
    
    
    def smooth_data(self, data: pd.Series) -> np.ndarray:
        """Helper function to pad and smooth the data with moving average."""
        padded_data = np.concatenate(
            (data.iloc[-self.smoothing_window:].values, 
             data.values, 
             data.iloc[:self.smoothing_window].values)
        )
        padded_data = pd.DataFrame(padded_data)
        padded_data = padded_data.replace(-np.inf, np.nan)
        interpolated_padded_values = padded_data.interpolate(
            method='spline', 
            order=3
            ).rolling(
                window=self.smoothing_window, 
                center=True
                ).mean()
        return interpolated_padded_values[self.smoothing_window:-self.smoothing_window].values


class PlotUtils:
    """
    Class to read CSV files and save timestamped plots.
    """
    
    @staticmethod
    def read_csv(file_path: str) -> pd.DataFrame:
        """Check file existence and read the CSV."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
        return pd.read_csv(file_path)
    
    
    @staticmethod
    def save_plot(filename: str):
        """Save the plot with a timestamp to avoid overwriting."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        plt.savefig(f"{timestamp}_{filename}")
        

class PlotFunctions:
    """
    Class to plot log2 fold change, relative log2 fold change, Jensen-Shannon Divergence and binsize distribution.
    """
    
    def __init__(self, data_processor, plot_utils, raw_file_path, normalised_file_path):
        self.data_processor = data_processor
        self.plot_utils = plot_utils
        self.raw_file_path = raw_file_path
        self.normalised_file_path = normalised_file_path
    
    def plot_log2fc(self):
        """Plot log2 fold change for normalized data."""
        df = self.plot_utils.read_csv(self.normalised_file_path)
        bins = df['Bins']

        plt.figure(figsize=(10, 6))
        for column in df.columns[1:]:
            smoothed_values = self.data_processor.smooth_data(df[column])
            plt.plot(bins, smoothed_values, label=f"{column}")
        
        plt.xlabel('Position (Mb)')
        plt.ylabel('Log2 Fold Change')
        plt.title('Log2 Fold Change in Reads from Sample Mean')
        plt.legend(title='Sample', bbox_to_anchor=(1.04, 1), loc="upper left")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.tight_layout()
        self.plot_utils.save_plot("log2fc_by_pos.png")
        plt.close()


    def relative_log2FC(self):
        """Plot relative log2 fold change between test and control samples."""
        df = self.plot_utils.read_csv(self.raw_file_path)
        bins = df['Bins']
        counts_data = df.drop(columns=['Bins'])
        normalized_data = counts_data / counts_data.sum(axis=0)

        test_columns = [col for col in normalized_data.columns if "test" in col]
        control_columns = {col.replace("test", "control"): col for col in test_columns}

        relative_data = pd.DataFrame()
        for control_col, test_col in control_columns.items():
            if control_col in normalized_data.columns:
                relative_data[test_col] = normalized_data[test_col] / normalized_data[control_col]
        relative_df = pd.concat([bins, relative_data], axis=1)
        relative_df.to_csv("relative_change.csv", index=False)
        
        if len(relative_df.columns) > 1:
            plt.figure(figsize=(10, 6))
            for column in relative_df.columns[1:]:
                interpolated_values = self.data_processor.interpolate_missing(relative_df[column])
                smoothed_values = self.data_processor.smooth_data(interpolated_values)
                log2fc_values = np.log2(smoothed_values)
                plt.plot(bins, log2fc_values, label=f"{column}")
                
            plt.xlabel('Position (Mb)')
            plt.ylabel('Relative Log2 Fold Change')
            plt.title('Change in Normalised Reads from Control')
            plt.legend(title='Sample', bbox_to_anchor=(1.04, 1), loc="upper left")
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            self.plot_utils.save_plot("relative_change.png")
            plt.close()


    def jensen_shannon_divergence(self, P: np.ndarray, Q: np.ndarray) -> float:
        """Compute Jensen-Shannon Divergence between two distributions."""
        M = 0.5 * (P + Q)
        kl_p_m = np.sum(rel_entr(P, M))
        kl_q_m = np.sum(rel_entr(Q, M))
        return 0.5 * (kl_p_m + kl_q_m)


    def plot_JS_div(self):
        """Plot the Jensen-Shannon Divergence heatmap."""
        df = self.plot_utils.read_csv(self.raw_file_path)
        counts_data = df.drop(columns=['Bins'])
        normalized_data = counts_data / counts_data.sum(axis=0)
        col_names = list(counts_data.columns.values)

        n_samples = normalized_data.shape[1]
        js_matrix = np.zeros((n_samples, n_samples))

        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                p = normalized_data.iloc[:, i].values
                q = normalized_data.iloc[:, j].values
                js_div_value = self.jensen_shannon_divergence(p, q)
                js_matrix[i, j] = js_matrix[j, i] = js_div_value

        js_matrix_df = pd.DataFrame(js_matrix, index=col_names, columns=col_names)

        plt.figure(figsize=(10, 8)) 
        sns.heatmap(
            js_matrix_df,
            annot=True,
            cmap="coolwarm",
            linewidths=0.5,
            mask = np.triu(np.ones_like(js_matrix, dtype=bool))
        )
        plt.title("Jensen-Shannon Divergence Between Samples")
        plt.xlabel("Sample")
        plt.ylabel("Sample")
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0, ha='right')
        plt.tight_layout()
        self.plot_utils.save_plot("jensen_shannon_divergence.png")
        plt.close()


    def raw_read_boxplot(self):
        """Create a boxplot for raw read counts."""
        raw_df = self.plot_utils.read_csv(self.raw_file_path)
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
        self.plot_utils.save_plot("raw_read_boxplot.png")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Visualize read log2FC and Jensen-Shannon Divergence.")
    parser.add_argument("raw_file_path", help="Path to the binned_reads.csv")
    parser.add_argument("normalised_file_path", help="Path to log2FC_binned_reads.csv")
    parser.add_argument("--smoothing_window", type=int, default=100, help="Smoothing window size for log2FC plot.")
    args = parser.parse_args()

    # Instantiate DataProcessing and PlotUtils class objects
    data_processor = DataProcessing(args.raw_file_path)
    plot_utils = PlotUtils()
    
    # Pass objects to PlotFunctions with the correct file paths
    plotter = PlotFunctions(data_processor, plot_utils, args.raw_file_path, args.normalised_file_path)
    
    plotter.plot_log2fc()
    plotter.relative_log2FC()
    plotter.plot_JS_div()
    plotter.raw_read_boxplot()

if __name__ == "__main__":
    main()
