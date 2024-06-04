#!/usr/bin/env python3

#Created by Martin Vennick Haugbølle

#Import libraries
import pyopenms as oms
import numpy as np
import pandas as pd
from tqdm import tqdm
import os


def get_MS1_and_MS2_intensity_array_from_exp(exp_name):
    all_intensities_ms1 = np.array([])
    all_intensities_ms2 = np.array([])

    for spec in tqdm(exp_name):
        type_spec = spec.getMSLevel()
        mz, intensity = spec.get_peaks()
        if type_spec == 1:
            all_intensities_ms1 = np.concatenate((all_intensities_ms1, intensity), axis=None)
        if type_spec == 2:
            all_intensities_ms2 = np.concatenate((all_intensities_ms2, intensity), axis=None)
    
    return all_intensities_ms1, all_intensities_ms2

def add_peak_retention_time(df):
    results = []
    last_ms1_peak_retention_time = None

    for _, row in df.iterrows():
        if row['type'] == 'ms1':
            # Find the retention time corresponding to the peak intensity
            max_intensity_index = row['intensity'].index(max(row['intensity']))
            peak_retention_time = row['retention time'][max_intensity_index]
            mz = row['mz']

            # Add a new dictionary to results
            results.append({'peak_retention_time': peak_retention_time})

            # Store the peak retention time for 'ms2' entries
            last_ms1_peak_retention_time = peak_retention_time
        elif row['type'] == 'ms2' and last_ms1_peak_retention_time is not None:
            # Add a new dictionary to results with the last 'ms1' peak retention time
            results.append({'peak_retention_time': last_ms1_peak_retention_time})

    # Convert results to a DataFrame
    df_results = pd.DataFrame(results)

    return pd.concat([df.reset_index(drop=True),df_results.reset_index(drop=True)], axis=1)

def save_as_mzML_for_MS2_only(df, out_file: str):
    ### Creating mzML file ###
    # # MSSpectrum instances.
    exp_out = oms.MSExperiment()

    # #intitializing variables
    prev_rt = -1
    prev_ms1_mz = -1

    # #we can run whole loop based on the RT's as all peaks are sorted by RT
    for i, row in df.iterrows():

        if row["type"] == "ms1":
            current_ms1_mz = row["mz"]

            # Setting the precursor metadata
            p = oms.Precursor()
            p.setMZ(current_ms1_mz)
            continue

        current_rt = row["peak_retention_time"]

        if current_rt != prev_rt or current_ms1_mz != prev_ms1_mz:
            try:
                exp_out.addSpectrum(spectrum)
            except:
                pass
            
            spectrum = oms.MSSpectrum()


            # setting metadata
            spectrum.setRT(current_rt)
            spectrum.setMSLevel(int((row["type"][2:])))         #checking for mistakes as if not 2 then something is wrong with the initial data
            spectrum.setPrecursors([p])

        
        # Creating peak object and adding metadata
        peak = oms.Peak1D()
        peak.setMZ(row["mz"])
        peak.setIntensity(max(row["intensity"]))
        spectrum.push_back(peak)                        #push peak to the spectrum

        
        # updating the prev_rt and prev_ms1_mz
        prev_rt = current_rt
        prev_ms1_mz = current_ms1_mz
    
    # Store as mzML
    oms.MzMLFile().store(out_file, exp_out)

    return print(f"All done with saving the spectra to an .mzMl file which is stored at {out_file}")

def save_as_mzML_for_MS1_and_MS2_included(df, out_file: str):
    ### Creating mzML file ###
    # # MSSpectrum instances.
    exp_out = oms.MSExperiment()

    ms2_spectrums = {}

    # #intitializing variables
    prev_rt = -1
    prev_ms1_mz = -1
    prev_ms1_rt = -1

    # #we can run whole loop based on the RT's as all peaks are sorted by RT
    for i, row in df.iterrows():

        current_rt = row["peak_retention_time"]

        if row["type"] == "ms1":
            current_ms1_mz = row["mz"]
            current_ms1_rt = row["peak_retention_time"]
            
            if current_ms1_rt != prev_ms1_rt:
                try:
                    exp_out.addSpectrum(MS1_spectrum)

                    for spec in ms2_spectrums.values():
                        exp_out.addSpectrum(spec)
                    
                    ms2_spectrums = {}
                except:
                    pass

                MS1_spectrum = oms.MSSpectrum()

                # setting metadata
                MS1_spectrum.setRT(current_rt)
                MS1_spectrum.setMSLevel(int((row["type"][2:])))
            
            # Creating peak object and adding metadata
            peak = oms.Peak1D()
            peak.setMZ(row["mz"])
            peak.setIntensity(max(row["intensity"]))
            MS1_spectrum.push_back(peak)

            # Setting the precursor metadata
            p = oms.Precursor()
            p.setMZ(current_ms1_mz)
            
            prev_ms1_rt = current_ms1_rt
            continue

        if current_rt != prev_rt or current_ms1_mz != prev_ms1_mz:

            #creating spectrum
            ms2_spectrums[current_ms1_mz] = oms.MSSpectrum()

            # setting metadata
            ms2_spectrums[current_ms1_mz].setRT(current_rt)
            ms2_spectrums[current_ms1_mz].setMSLevel(int((row["type"][2:])))         #checking for mistakes as if not 2 then something is wrong with the initial data
            ms2_spectrums[current_ms1_mz].setPrecursors([p])

        
        # Creating peak object and adding metadata
        peak = oms.Peak1D()
        peak.setMZ(row["mz"])
        peak.setIntensity(max(row["intensity"]))
        ms2_spectrums[current_ms1_mz].push_back(peak)                        #push peak to the spectrum

        
        # updating the prev_rt and prev_ms1_mz
        prev_rt = current_rt
        prev_ms1_mz = current_ms1_mz
    
    # Get the final ms1 in the script
    exp_out.addSpectrum(MS1_spectrum)
    
    # Get the final ms2 in the script
    for spec in ms2_spectrums.values():
        exp_out.addSpectrum(spec)
    
    # Store as mzML
    oms.MzMLFile().store(out_file, exp_out)

    return print(f"All done with saving the spectra to an .mzMl file which is stored at {out_file}")

def get_file_path(file_path):
    if os.path.isabs(file_path):
        return file_path
    else:
        return os.path.join(os.getcwd(), file_path)

def create_floating_point_range(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

def find_number_in_range(lst, lower, upper):
    window = [num for num in lst if lower <= num < upper]
    if window:
        return min(window, key=lst.index)
    return None

def round_to_nearest(num, peak_mz_range: float=0.01):
    return round(num / peak_mz_range) * peak_mz_range

def log_transform_array_in_dataframe(df, array_to_log: str):
    df[array_to_log] = df[array_to_log].apply(lambda array: np.log(array).tolist())
    return df

def add_unique_id_col_to_dataframe(df):
    df['unique_id'] = range(1, len(df) + 1)
    return df

def explode_dataframe(df, list_columns):
        df = df.reset_index(drop=True)
        for column in list_columns:
            df = df.explode(column)
        return df

def get_intensity_array_of_experiment(experiment, ms_level):
        #make array of all intensity values from the peaks
    intensities = []

    for spec in experiment:
        if spec.getMSLevel() != ms_level:
            continue
        mz, intensity = spec.get_peaks()
        #print(intensity)
        for intense in intensity:
            try:
                float(intense)
                if intense <= 0:
                    continue
            except:
                continue
            intensities.append(intense)
    
    # Making the intensities into a numpy array
    intensities = np.array(intensities)

    return intensities

def log_array(array):
    return np.log(array)

def find_noise_level_from_intensities(intensities: np.array):

    # Compute the histogram
    hist_values, bin_edges = np.histogram(intensities, bins=25)

    # Compute the derivative of the histogram values
    hist_derivative = np.diff(hist_values)

    # Find where the derivative changes sign
    sign_changes = np.where(np.diff(np.sign(hist_derivative)) > 0)[0]

    # If there are multiple local minima, choose the first one
    if len(sign_changes) > 0:
        min_index = sign_changes[0]
    else:
        min_index = None

    # Find the corresponding intensity value
    if min_index is not None:
        threshold = bin_edges[min_index]
    else:
        threshold = 0
    
    return threshold

def remove_rows_under_point_across_the_peak_limit(df, points_across_the_peak_limit: int):
    # Remove rows where the number of points across the peak is less than the threshold
    return df[df['intensity'].apply(len) >= points_across_the_peak_limit]

def centroid_peaks(experiment):
    centroided_spectra = oms.MSExperiment()
    # input, output, chec_spectrum_type (if set, checks spectrum type and throws an exception if a centroided spectrum is passed)
    oms.PeakPickerHiRes().pickExperiment(
    experiment, centroided_spectra, True)
    return centroided_spectra