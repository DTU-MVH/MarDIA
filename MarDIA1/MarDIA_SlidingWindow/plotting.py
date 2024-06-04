# plotting.py

#importing libraries
import os
import shutil
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Plotting settings
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.width', 1000)  # Adjust the display width (optional)


def plot_MS2_DDA(df, num_plots: int, directory_name: str):

    # If the directory exists, remove it and all its contents
    if os.path.exists(directory_name):
        shutil.rmtree(directory_name)

    # Create a directory to save the plots
    os.makedirs(directory_name)

    # Initialize variables
    plot_data = []
    plot_count = 0
    mz_value = None
    max_intensity = None
    peak_intensity_index = None
    peak_retention_time = None

    # Iterate over the rows of the dataframe
    for i, row in df.iterrows():
        if row['type'] == 'ms1':

            # If there is data to plot and we haven't reached the limit, create a plot
            if plot_data and plot_count < num_plots:
                plt.figure()
                for data in plot_data:
                    plt.plot([data[0], data[0]], [0, data[1]], label=f'mz={data[0]}')  # Plot as pillars
                plt.ylim(bottom=0)  # Make y-axis start from 0
                plt.xlabel('mz')
                plt.ylabel('Intensity')
                plt.title(f"plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}")
                plt.savefig(f'{directory_name}/plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}.png')
                plt.close()
                plot_count += 1

            # Start a new plot
            max_intensity = max(row['intensity'])  # Find the maximum peak intensity
            peak_intensity_index = row['intensity'].index(max_intensity)  # Find the index of the maximum peak intensity
            peak_retention_time = row['retention time'][peak_intensity_index]  # Get the retention time at that index
            
            plot_data = []
            mz_value = row['mz']
        
        elif row['type'] == 'ms2':
            # Add the data to the current plot
            x = row['mz']
            y = max(row['intensity'])
            plot_data.append((x, y))

    # Get the last plot
    if plot_data and plot_count < num_plots:
        plt.figure()
        for data in plot_data:
            plt.plot([data[0], data[0]], [0, data[1]], label=f'mz={data[0]}')  # Plot as pillars
        plt.ylim(bottom=0)  # Make y-axis start from 0
        plt.xlabel('mz')
        plt.ylabel('Intensity')
        plt.title(f"plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}")
        plt.savefig(f'{directory_name}/plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}.png')
        plt.close()


def plot_MS1_MS2(df, num_plots : int, directory_name: str):

    # If the directory exists, remove it and all its contents
    if os.path.exists(directory_name):
        shutil.rmtree(directory_name)

    # Create a directory to save the plots
    os.makedirs(directory_name)

    # Initialize variables
    ms1_row = None
    plot_count = 0
    mz_value = None
    max_intensity = None
    peak_intensity_index = None
    peak_retention_time = None

    # Iterate over the dataframe rows
    for row in df.itertuples():
        if row['type'] == 'ms1':

            # If there's a previous ms1 row, save the plot and start a new one
            if ms1_row is not None:
                plt.xlabel('Retention time (RT)')
                plt.ylabel('Log(Intensity)')
                plt.title(f"plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}")
                plt.savefig(f'{directory_name}/plot_{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}.png')
                plt.close()
                plot_count += 1
                if plot_count >= num_plots:
                    break
            
            #Update variables
            max_intensity = max(row['intensity'])  # Find the maximum peak intensity
            peak_intensity_index = row['intensity'].index(max_intensity)  # Find the index of the maximum peak intensity
            peak_retention_time = row['retention time'][peak_intensity_index]  # Get the retention time at that index
            mz_value = row["mz"]
            ms1_row = row
            plt.plot(ms1_row['retention time'], ms1_row['intensity'], 'r-', linewidth=2)
        elif row['type'] == 'ms2' and ms1_row is not None:
            plt.plot(row['retention time'], row['intensity'], 'b-')

    # Save the last plot
    if ms1_row is not None and plot_count < num_plots:
        plt.xlabel('Retention time (RT)')
        plt.ylabel('Log(Intensity)')
        plt.title(f"plot{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}")
        plt.savefig(f'{directory_name}/plot_{plot_count}_mz{mz_value}_peak_rt_{peak_retention_time}.png')
        plt.close()


def plot_histogram(intensities, bins=25, file_name: str="histogram1", folder_name: str="histograms"):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    plt.hist(intensities, bins)
    plt.savefig(f'{folder_name}/{file_name}.png')
    plt.close()


# Visulisation methods for dataframes
def plot_dataframe_head(df,row_number: int=10):
    print(df.head(row_number))

def plot_dataframe_tail(df,row_number: int=10):        
    print(df.tail(row_number))

def plot_dataframe_shape(df):
    print(df.shape)

def plot_dataframe_info(df):
    print(df.info())



# Plot code:

# Create DDA plots over the ms2 data
#plot.plot_MS2_DDA(df_result_correlation,num_plots=20,directory_name="DDA_MS2_plots")


#for plotting
#log transform intensity array
#df_result_correlation_log = utils.log_transform_array_in_dataframe(df_result_correlation, array_to_log="intensity")

#displaying precursor and fragments in the same plot
#plot.plot_MS1_MS2(df_result_correlation_log, num_plots=20, directory_name="plots1_folder")