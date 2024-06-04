# README.txt

## Script Information

**Script Name**: mzML Processing Script

**Author**: Martin Vennick Haugbølle

**Description**: This script processes mzML files to detect features, filter peaks, and correlate fragments to precursors using Pearson correlation. The final output is saved as an mzML file.

---

## Requirements

To run this script, you need the following Python libraries:

- pyopenms
- numpy
- pandas
- matplotlib
- tabulate
- tqdm
- scipy
- os
- argparse

You can install these libraries using the following command:

```bash
pip install pyopenms numpy pandas matplotlib tabulate tqdm scipy
```

---

## Usage

To use the script, follow these steps:

1. **Prepare the Input mzML File**: Ensure you have the input mzML file you wish to process.

2. **Run the Script**: Use the following command to execute the script:

```bash
python script_name.py --input_file <path_to_input_file> --output_file <path_to_output_file> [other options]
```

### Command-Line Arguments

- `--input_file` (str): Required. Path to the input mzML file.
- `--output_file` (str): Optional. Path to the output file. Default is "output".
- `--mz_window_lower_bound` (int): Optional. Lower bound for mz window. Default is 250.
- `--mz_window_upper_bound` (int): Optional. Upper bound for mz window. Default is 750.
- `--range_start` (int): Optional. Start of the range. Default is 0.
- `--range_end` (int): Optional. End of the range. Default is -1.
- `--peaks_pr_spectra` (int): Optional. Number of peaks per spectra. Default is 200.
- `--intensity_cutoff_ms1` (float): Optional. Intensity cutoff for MS1. Default is None.
- `--intensity_cutoff_ms2` (float): Optional. Intensity cutoff for MS2. Default is None.
- `--round_nb` (int): Optional. Number of rounds. Default is 3.
- `--pearson_min_score` (float): Optional. Minimum Pearson correlation score. Default is 0.95.
- `--points_across_the_peak_lim` (int): Optional. Limit for points across the peak. Default is 3.
- `--fragment_count_pr_precursor` (int): Optional. Fragment count per precursor. Default is 3.

### Example Usage

```bash
python script_name.py --input_file example.mzML --output_file output.mzML --mz_window_lower_bound 200 --mz_window_upper_bound 800 --peaks_pr_spectra 150 --pearson_min_score 0.9
```

---

## Script Details

### Main Function

The main function of the script is responsible for the following tasks:

1. **Loading the mzML File**: Loads the mzML file onto the pyOpenMS object.
2. **Feature Detection for MS1**: Detects features for MS1 level.
3. **Centroiding the Spectra**: Centroids the spectra.
4. **Finding Noise Levels**: Finds the noise levels for MS1 and MS2.
5. **Converting Experiment to DataFrame**: Converts the experiment object to DataFrame.
6. **Peak Picking for Feature Detection**: Picks peaks for feature detection in MS1.
7. **Removing Duplicate Peaks**: Removes duplicate peaks in the same scan.
8. **Grouping Peaks**: Groups peaks by mz values.
9. **Splitting Rows into Separate Peaks**: Splits rows into separate peaks.
10. **Correlating Fragments to Precursors**: Correlates fragments to precursors using Pearson correlation.
11. **Saving Results**: Saves the results as an mzML file.

### Imported Functions

- From `plotting.py`: Plotting functions for MS2 DDA, MS1 MS2, histogram, and dataframes.
- From `utils.py`: Utility functions for intensity arrays, retention times, mzML saving, logging, and more.
- From `dataloader.py`: Functions to load mzML files, apply filters, and convert experiments.
- From `algorithm.py`: Functions for peak splitting, duplicate removal, correlation, and feature detection.

### Output

The script generates an output mzML file containing the processed data, including detected features and correlated fragments.

---

## Notes

- Ensure all required libraries are installed before running the script.
- Adjust the command-line arguments as per your requirements.
- The default values for the arguments are set to typical ranges, but they can be modified based on specific needs.

---

For any issues or further assistance, please contact the author.