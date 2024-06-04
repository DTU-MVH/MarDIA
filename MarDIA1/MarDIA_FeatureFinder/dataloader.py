#!/usr/bin/env python3

#Created by Martin Vennick Haugbølle

#Import libraries
import pyopenms as oms
import numpy as np
import pandas as pd
from tqdm import tqdm
from utils import round_to_nearest
from matplotlib import pyplot as plt


def load_mzML_file_onto_experiment(mzML_file):
    #initialize experiment
    exp = oms.MSExperiment()

    #loading file
    oms.MzMLFile().load(mzML_file, exp)

    return exp

def window_mower_filter(spec, window_size=1.5, peak_count=1, movetype="slide"):
    # Apply window mowing
    window_mower_filter = oms.WindowMower()

    # Set parameters
    params = oms.Param()
    # Defines the m/z range of the sliding window
    params.setValue("windowsize", window_size, "")   #getting the 4 isotopes
    # Defines the number of highest peaks to keep in the sliding window
    params.setValue("peakcount", peak_count, "")
    # Defines the type of window movement: jump (window size steps) or slide (one peak steps)
    params.setValue("movetype", movetype, "")

    # Apply window mowing
    window_mower_filter.setParameters(params)
    window_mower_filter.filterPeakSpectrum(spec)

    return spec

def threshold_mower_filter(spec, threshold: float):
    # Apply threshold mowing
    threshold_mower_filter = oms.ThresholdMower()

    # Set parameters
    params = oms.Param()
    
    # Defines the minimum intensity value to keep a peak
    params.setValue("threshold", threshold, "")

    # Apply threshold mowing
    threshold_mower_filter.setParameters(params)
    threshold_mower_filter.filterPeakSpectrum(spec)

    return spec

def nlargest_filter(spec, n: int):
    
    # Apply N-Largest filter
    nlargest_filter = oms.NLargest()

    # Set parameters
    params = oms.Param()
    params.setValue("n", n, "")

    # Apply N-Largest filter
    nlargest_filter.setParameters(params)
    nlargest_filter.filterPeakSpectrum(spec)

    return spec

def convert_experiment_to_dataframe(experiment, mz_window_lower_bound=250, mz_window_upper_bound=750, range=list(), peaks_pr_spectra: int=200, intensity_cutoff_ms1: float=5000.0, intensity_cutoff_ms2: float=1000.0, round_nb=3):
    # Initialize lists to store the data and setting counters
    df_rows_ms1 = []
    df_rows_ms2 = []
    ms1_nb = 0                     #The ms1_nb is a integer assigned to all peaks to keep track of the ms1 scan that the fragments are coming from.
    count = 0

    # Check if the range is valid
    try:
        if len(range) > 2:
            ValueError("The range parameter only has only a lower and upper limit, thus the list should only have lenght = 2")
        if range[0] < 0:
            ValueError("The lower limit of the range cannot be below 0")
        if range[1] > (len(experiment.getSpectra())):
            ValueError("The upper limit of the range of spectra to analyze is more than the maximum amount of spectra")
    except:    
        #if no range was chosen take the full experiment as standard
        range == [0,len(experiment.getSpectra())]


    # Iterate over the spectra in the experiment
    for spec in tqdm(experiment, total=len(experiment.getSpectra())):
        count += 1
        
        # Keeping in the range
        if count < range[0]:
            continue
        if count == range[1]:
            break
        
        # Removing the noise from the spectra by applying a filter that only keeps the n highest peaks in the spectra
        spec = nlargest_filter(spec=spec, n=peaks_pr_spectra)

        # Get metadata from exp object
        type_spec = spec.getMSLevel()

        # Filtering for ms1
        if type_spec == 1:
            
            # Apply window mowing to the spectrum to remove under threshold peaks
            spec = threshold_mower_filter(spec=spec, threshold=intensity_cutoff_ms1)

            # Apply the window mowing filter to the spectrum to only get 1 peak per mz value inside the window
            spec = window_mower_filter(spec=spec, window_size=0.04, peak_count=1, movetype="slide")
            
            # Getting spectra data
            mz, intensity = spec.get_peaks()
            rt = spec.getRT()

            # Keeping track of the ms1 scan to have a unique identifier for each ms1 scan
            ms1_nb += 1

            # Iterate over the peaks in the spectrum and only keeping peaks within the mz window.
            for n, mz_val in np.ndenumerate(mz):
                # Check if the precursor is in the mz window
                if not mz_window_lower_bound <= mz_val <= mz_window_upper_bound:
                    continue
                else:
                    df_rows_ms1.append({'type_spec': type_spec,"ms1_nb" : ms1_nb, 'mz': round_to_nearest(num=mz_val), 'intensity': round(intensity[n],round_nb), 'retention time': round(rt,round_nb)})

        # Filtering for ms2
        elif type_spec == 2:
            
            # Apply window mowing to the spectrum to remove under threshold peaks
            spec = threshold_mower_filter(spec=spec, threshold=intensity_cutoff_ms2)
            
            # Apply the window mowing filter to the spectrum to only get 1 peak per mz value inside the window
            spec = window_mower_filter(spec=spec, window_size=0.04, peak_count=1, movetype="slide")

            # Getting spectra data
            mz, intensity = spec.get_peaks()
            rt = spec.getRT()

            for n, mz_val in np.ndenumerate(mz):
                # Check if the fragment is in the mz window
                if not mz_window_lower_bound <= mz_val <= mz_window_upper_bound:
                    continue
                else:
                    df_rows_ms2.append({'type_spec': type_spec, "ms1_nb" : ms1_nb, 'mz': round_to_nearest(num=mz_val), 'intensity': round(intensity[n],round_nb), 'retention time': round(rt,round_nb)})
    
    return pd.DataFrame(df_rows_ms1), pd.DataFrame(df_rows_ms2)


def convert_experiment_to_dataframe_for_feature_detection(experiment, mz_window_lower_bound=250, mz_window_upper_bound=750, range=list(), peaks_pr_spectra: int=200, intensity_cutoff_ms1: float=5000.0, intensity_cutoff_ms2: float=1000.0, round_nb=3):
    # Initialize lists to store the data and setting counters
    df_rows_ms1 = []
    df_rows_ms2 = []
    ms1_nb = 0                     #The ms1_nb is a integer assigned to all peaks to keep track of the ms1 scan that the fragments are coming from.
    count = 0

    # Check if the range is valid
    try:
        if len(range) > 2:
            ValueError("The range parameter only has only a lower and upper limit, thus the list should only have lenght = 2")
        if range[0] < 0:
            ValueError("The lower limit of the range cannot be below 0")
        if range[1] > (len(experiment.getSpectra())):
            ValueError("The upper limit of the range of spectra to analyze is more than the maximum amount of spectra")
    except:    
        #if no range was chosen take the full experiment as standard
        range == [0,len(experiment.getSpectra())]


    # Iterate over the spectra in the experiment
    for spec in tqdm(experiment, total=len(experiment.getSpectra())):
        count += 1
        
        # Keeping in the range
        if count < range[0]:
            continue
        if count == range[1]:
            break
        
        # Removing the noise from the spectra by applying a filter that only keeps the n highest peaks in the spectra
        spec = nlargest_filter(spec=spec, n=peaks_pr_spectra)

        # Get metadata from exp object
        type_spec = spec.getMSLevel()

        # Filtering for ms1
        if type_spec == 1:
            
            # Apply window mowing to the spectrum to remove under threshold peaks
            spec = threshold_mower_filter(spec=spec, threshold=intensity_cutoff_ms1)

            # Apply the window mowing filter to the spectrum to only get 1 peak per mz value inside the window
            spec = window_mower_filter(spec=spec, window_size=0.04, peak_count=1, movetype="slide")
            
            # Getting spectra data
            mz, intensity = spec.get_peaks()
            rt = spec.getRT()

            # Keeping track of the ms1 scan to have a unique identifier for each ms1 scan
            ms1_nb += 1

            # Iterate over the peaks in the spectrum and only keeping peaks within the mz window.
            for n, mz_val in np.ndenumerate(mz):
                # Check if the precursor is in the mz window
                if not mz_window_lower_bound <= mz_val <= mz_window_upper_bound:
                    continue
                else:
                    df_rows_ms1.append({'type_spec': type_spec,"ms1_nb" : ms1_nb, 'mz': mz_val, 'intensity': round(intensity[n],round_nb), 'retention time': round(rt,round_nb)})

        # Filtering for ms2
        elif type_spec == 2:
            
            # Apply window mowing to the spectrum to remove under threshold peaks
            spec = threshold_mower_filter(spec=spec, threshold=intensity_cutoff_ms2)
            
            # Apply the window mowing filter to the spectrum to only get 1 peak per mz value inside the window
            spec = window_mower_filter(spec=spec, window_size=0.04, peak_count=1, movetype="slide")

            # Getting spectra data
            mz, intensity = spec.get_peaks()
            rt = spec.getRT()

            for n, mz_val in np.ndenumerate(mz):
                # Check if the fragment is in the mz window
                if not mz_window_lower_bound <= mz_val <= mz_window_upper_bound:
                    continue
                else:
                    df_rows_ms2.append({'type_spec': type_spec, "ms1_nb" : ms1_nb, 'mz': round_to_nearest(num=mz_val), 'intensity': round(intensity[n],round_nb), 'retention time': round(rt,round_nb)})
    
    return pd.DataFrame(df_rows_ms1), pd.DataFrame(df_rows_ms2)