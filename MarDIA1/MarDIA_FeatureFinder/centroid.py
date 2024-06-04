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
import shutil
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
                   remove_rows_under_point_across_the_peak_limit)

#importing functions from dataloader.py
from dataloader import (load_mzML_file_onto_experiment,
                        window_mower_filter,
                        threshold_mower_filter,
                        nlargest_filter,
                        convert_experiment_to_dataframe)

#importing functions from algorithm.py
from algorithm import (split_rows_into_seperate_peaks,
                    remove_same_scan_duplicates,
                    correlate_fragments_to_precursors,
                    filter_fragment_count_pr_precursor,
                    vectorized_correlations)


# Loading the mzML file onto the pyOpenMs object: experiment
print("Loading mzML file onto experiment object:")
input_mzfile = "Sbrodae_DIA_top200.mzML"
exp = load_mzML_file_onto_experiment(mzML_file=input_mzfile)

centroided_spectra = oms.MSExperiment()
 
# input, output, chec_spectrum_type (if set, checks spectrum type and throws an exception if a centroided spectrum is passed)
oms.PeakPickerHiRes().pickExperiment(
    exp, centroided_spectra, True)

# Save the centroided spectra as an mzML file
output_mzfile = "Sbrodae_DIA_top200_centroided.mzML"
oms.MzMLFile().store(output_mzfile, centroided_spectra)

print(f"Centroided spectra saved as {output_mzfile}")
