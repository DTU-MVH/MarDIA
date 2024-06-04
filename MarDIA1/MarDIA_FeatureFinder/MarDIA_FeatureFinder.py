#!/usr/bin/env python3

#Created by Martin Vennick Haugbølle

#Import libraries
import pyopenms as oms
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tabulate import tabulate
from tqdm import tqdm
from scipy.stats import pearsonr
import os
import argparse

# Importing the Plotting class from plotting.py
from plotting import (plot_MS2_DDA,
                      plot_MS1_MS2,
                      plot_histogram,
                      plot_dataframe_head,
                      plot_dataframe_tail,
                      plot_dataframe_shape,
                      plot_dataframe_info)


# Importing from utils.py
from utils import (get_MS1_and_MS2_intensity_array_from_exp,
                   add_peak_retention_time,
                   save_as_mzML_for_MS2_only,
                   save_as_mzML_for_MS1_and_MS2_included,
                   get_file_path,
                   create_floating_point_range,
                   find_number_in_range, round_to_nearest,
                   log_transform_array_in_dataframe,
                   add_unique_id_col_to_dataframe,
                   explode_dataframe,
                   get_intensity_array_of_experiment,
                   log_array,
                   find_noise_level_from_intensities,
                   remove_rows_under_point_across_the_peak_limit,
                   centroid_peaks)

#importing functions from dataloader.py
from dataloader import (load_mzML_file_onto_experiment,
                        window_mower_filter,
                        threshold_mower_filter,
                        nlargest_filter,
                        convert_experiment_to_dataframe,
                        convert_experiment_to_dataframe_for_feature_detection)

#importing functions from algorithm.py
from algorithm import (split_rows_into_seperate_peaks,
                    remove_same_scan_duplicates,
                    correlate_fragments_to_precursors,
                    filter_fragment_count_pr_precursor,
                    vectorized_correlations,
                    feature_detection_ms1,
                    peak_picking_for_feature_detection_ms1)


def main(
        input_mzfile,
        mz_window_lower_bound,
        mz_window_upper_bound,
        range_start, range_end,
        peaks_pr_spectra,
        intensity_cutoff_ms1,
        intensity_cutoff_ms2,
        round_nb,
        pearson_min_score,
        points_across_the_peak_lim,
        fragment_count_pr_precursor,
        out_file):


    # Loading the mzML file onto the pyOpenMs object: experiment
    print("Loading mzML file onto experiment object:")
    exp = load_mzML_file_onto_experiment(mzML_file=input_mzfile)

    # Feature detection for MS1
    print("Feature detection for MS1:")
    features_ms1 = feature_detection_ms1(input_mzfile)

    # Centroiding the spectra
    exp = centroid_peaks(exp)

    # Finding the noise level from the bimodal distribution of intensity values.
    log_intensities_ms1 = log_array(get_intensity_array_of_experiment(exp, ms_level=1))
    log_intensities_ms2 = log_array(get_intensity_array_of_experiment(exp,ms_level=2))

    noise_level_ms1 = np.exp(find_noise_level_from_intensities(intensities=log_intensities_ms1))
    noise_level_ms2 = np.exp(find_noise_level_from_intensities(intensities=log_intensities_ms2))

    print(f"Noise level for MS1: {noise_level_ms1}")
    print(f"Noise level for MS2: {noise_level_ms2}")

    #plot_histogram(log_intensities_ms2, bins=25, file_name="ms2_histogram", folder_name="histograms")
    #plot_histogram(log_intensities_ms1, bins=25, file_name="ms1_histogram", folder_name="histograms")

    # Set noise level if not set by user, set calculated noise level
    if intensity_cutoff_ms1 == None:
        intensity_cutoff_ms1 = float(noise_level_ms1)

    if intensity_cutoff_ms2 == None:
        intensity_cutoff_ms2 = float(noise_level_ms2)

    # To be able to work with the data. 2 dataframes, from the experiment object, is created. Each row in the dataframe represents a peak in MS/MS scan and the scans are sorted to only include ms1 in one dataframe and ms2 in the other. The columns are specified here: columns=['type_spec', "ms1_nb", 'mz', 'intensity', 'retention time']. The ms1_nb is a integer assigned to all peaks to keep track of the ms1 scan that the fragments are coming from.
    print("Converting experiment to dataframe:")
    df_ms1, df_ms2 = convert_experiment_to_dataframe_for_feature_detection(
        experiment=exp,
        mz_window_lower_bound=mz_window_lower_bound,
        mz_window_upper_bound=mz_window_upper_bound,
        range=[range_start,range_end],
        peaks_pr_spectra=peaks_pr_spectra,
        intensity_cutoff_ms1=intensity_cutoff_ms1,
        intensity_cutoff_ms2=intensity_cutoff_ms2,
        round_nb=round_nb)


    # Selecting the feature for ms1
    print("Selecting the features for MS1:")
    df_ms1 = peak_picking_for_feature_detection_ms1(df_ms1, features_ms1)


    # To remove duplicate peaks in the same scan, grouping by mz values that have the same retention time and keep only the row with the maximum intensity in each group, is done.
    #df_ms1 = df_ms1.loc[df_ms1.groupby(['mz', 'retention time'])['intensity'].idxmax()]

    df_ms2 = df_ms2.loc[df_ms2.groupby(['mz', 'retention time'])['intensity'].idxmax()]


    # To have all information of all the peaks for an mz-value in one row, 
    # we group by mz and create arrays of all the intensities and retention times for that mz,
    # as an array of all the intensities with that mz-value. 
    # The intensities and retention times are of correcponding index in their arrays. 
    # Thus the intensity at index 0 was recorded at the retention time at index 0 in the retention time array.
    #df_ms1 = df_ms1.groupby('mz').agg({"ms1_nb" : list, 'intensity': list, 'retention time': list}).reset_index()

    df_ms2 = df_ms2.groupby('mz').agg({"ms1_nb" : list, 'intensity': list, 'retention time': list}).reset_index()

    # To seperate individual peaks, as all the peaks for a mz value is now in one row, 
    # we split the rows into multiple rows if the retention time between the peaks are more than 45 seconds.
    # Note: Should we use a sliding window for picking the peaks? We should definetely look into potentially improving this step
    print("Selecting the seperate peaks:")
    #df_ms1 = split_rows_into_seperate_peaks(df_ms1)

    df_ms2 = split_rows_into_seperate_peaks(df_ms2)

    # Creating a new column in the dataframes to sort the dataframes 
    # by the first retention time in the retention time array. 
    # This is important for the correlation step.
    df_ms1['sort_column'] = df_ms1['retention time'].apply(lambda x: x[0])

    # Sorting dataframe in terms of their retention time in ascending order. 
    # Important for the logic in the correlation step.
    sorted_df_ms1 = df_ms1.sort_values(by='sort_column', ascending=True)

    # Also for ms2
    df_ms2['sort_column'] = df_ms2['retention time'].apply(lambda x: x[0])

    sorted_df_ms2 = df_ms2.sort_values(by='sort_column', ascending=True)


    # All the rows in the dataframe have arrays for ms1_nb, intensity and retention time based on the mz value. To remove duplicates of the same mz value in the same scan under the sensitivity value, we find the duplicates and take the mean of the intensity and retention time arrays. This is done to ensure that the same mz value in the same scan is only represented once.
    # Note: Should we use the max intensity? (i think yes, because the max intensity is the most important for the precursor and another could just dilute the intensity.)
    # Not neseccary to do this step, as we have already removed duplicates in the same scan.
    sorted_df_ms1 = remove_same_scan_duplicates(sorted_df_ms1)

    sorted_df_ms2 = remove_same_scan_duplicates(sorted_df_ms2)

    # Removing all rows with less than 3 points across the peak
    sorted_df_ms1 = remove_rows_under_point_across_the_peak_limit(sorted_df_ms1, points_across_the_peak_lim)
    sorted_df_ms2 = remove_rows_under_point_across_the_peak_limit(sorted_df_ms2, points_across_the_peak_lim)


    print("Correlating fragments to precursors:")
    
    # To figure out which fragments are correlated to which precursors, we correlate the intensities of the fragments to the intensities of the precursors through a pearson correlation. We only keep the fragments that have a correlation of a user determined threshold or higher. We also only correlate the fragments that have a overlap in retention times with the precursor of a user defined number or more points across the peak.
    # Note: Should we use the p-value with a 0.05 cuttoff?
    df_result_correlation = vectorized_correlations(
        sorted_df_ms1,
        sorted_df_ms2,
        pearson_min_score=pearson_min_score,
        points_across_the_peak_lim=points_across_the_peak_lim)

    # Cleaning the dataframe for precursers with under 3 fragments correlated to them.
    df_result_correlation = filter_fragment_count_pr_precursor(df_result_correlation, fragment_count_pr_precursor=fragment_count_pr_precursor)

    # Adding peak retention time to the dataframe, so that all ms2 entries have the same peak retention time as the precursor they're correlated to.
    df_result_correlation = add_peak_retention_time(df_result_correlation)

    # Sorting the dataframe by peak retention time and after ms1_nb, this ensures all fragments are under the correct precursor and they are sorted by their retention time, starting with the earliest.
    df_result_correlation = df_result_correlation.sort_values(by=["peak_retention_time", "ms1_nb"], ascending=True)

    ### store as mzML ###
    save_as_mzML_for_MS1_and_MS2_included(df_result_correlation, out_file=out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process mzML file.')
    parser.add_argument('--input_file', type=str, required=True, help='Input mzML file or path to input mzML file')
    parser.add_argument('--output_file', type=str, default="output", help='Output file or path to output file')
    parser.add_argument('--mz_window_lower_bound', type=int, default=250, help='Lower bound for mz window')
    parser.add_argument('--mz_window_upper_bound', type=int, default=750, help='Upper bound for mz window')
    parser.add_argument('--range_start', type=int, default=0, help='Start of range')
    parser.add_argument('--range_end', type=int, default=-1, help='End of range')
    parser.add_argument('--peaks_pr_spectra', type=int, default=200, help='Number of peaks per spectra')
    parser.add_argument('--intensity_cutoff_ms1', type=float, default=None, help='Intensity cutoff for MS1')
    parser.add_argument('--intensity_cutoff_ms2', type=float, default=None, help='Intensity cutoff for MS2')
    parser.add_argument('--round_nb', type=int, default=3, help='Number of rounds')
    parser.add_argument('--pearson_min_score', type=float, default=0.95, help='Minimum Pearson correlation score')
    parser.add_argument('--points_across_the_peak_lim', type=int, default=3, help='Limit for points across the peak')
    parser.add_argument('--fragment_count_pr_precursor', type=int, default=3, help='Fragment count per precursor')


    args = parser.parse_args()

    # Checking and creating absolute paths for the input and output files
    if args.output_file == "output":
        output_file = get_file_path(args.output_file + "_" + os.path.basename(args.input_file) + str(args.pearson_min_score) + ".mzML")
    else:
        output_file = get_file_path(args.output_file)
    input_file = get_file_path(args.input_file)
    
    main(input_file,
         args.mz_window_lower_bound,
         args.mz_window_upper_bound,
         args.range_start,
         args.range_end,
         args.peaks_pr_spectra,
         args.intensity_cutoff_ms1,
         args.intensity_cutoff_ms2,
         args.round_nb,
         args.pearson_min_score,
         args.points_across_the_peak_lim,
         args.fragment_count_pr_precursor,
         output_file)