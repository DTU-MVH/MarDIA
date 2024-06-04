#!/usr/bin/env python3

#Created by Martin Vennick Haugbølle

#Import libraries
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr
import pyopenms as oms

def split_rows_into_seperate_peaks(df, time_between_scans=5, maximum_peak_elution_time=45):
    rows = []
    unequal_rows = []

    for _, row in tqdm(df.iterrows(), total=len(df)):
        
        # Getting the mz value of the row
        mz = row['mz']
        
        # Storing the arrays of the dataframe
        ms1_nb = np.array(row['ms1_nb'])
        intensities = np.array(row['intensity'])
        times = np.array(row['retention time'])

        # As the ms1_nb, intensities and retention times are all related to each other, we sort the indices of the times array to sort the other arrays as well.
        sorted_indices = np.argsort(times)          
        
        # Sorting the ms1_nb, intensities and times arrays from the sorted indices array to ensure correlation between them.
        ms1_nb = ms1_nb[sorted_indices]
        intensities = intensities[sorted_indices]
        times = times[sorted_indices]

        # Initialize the start variable
        start = 0
        times_ray = []

        # Iterate over the retention times and create cut-offs when the difference between two retention times are more than 45 seconds
        for end in range(0, len(times)):
            
            if times_ray == []:
                times_ray.append(end)
                continue

            if times[end] - times[times_ray[-1]] > time_between_scans or times[end] - times[times_ray[0]] > maximum_peak_elution_time:
                new_row = {'mz': mz, "ms1_nb" : list(ms1_nb[times_ray[0]:times_ray[-1]+1]), 'intensity': list(intensities[times_ray[0]:times_ray[-1]+1]), 'retention time': list(times[times_ray[0]:times_ray[-1]+1])}
            
                # If the length of the intensity array and the retention time array is not equal, we store the row in a list to check for errors later
                if len(new_row['intensity']) != len(new_row['retention time']):
                    unequal_rows.append(new_row)
                    times_ray = []
                else:
                    rows.append(new_row)
                    times_ray = []
            
            times_ray.append(end)


        # Create a new row for the last part of the arrays
        if times_ray != []:
            new_row = {'mz': mz, "ms1_nb" : list(ms1_nb[times_ray[0]:times_ray[-1]+1]), 'intensity': list(intensities[times_ray[0]:times_ray[-1]+1]), 'retention time': list(times[times_ray[0]:times_ray[-1]+1])}        
        
            if len(new_row['intensity']) != len(new_row['retention time']):
                unequal_rows.append(new_row)
            else:
                rows.append(new_row)
    
    # Check if uneq_rows is empty to ensure all intensity arrays and retention time arrays are of equal lenght.
    if unequal_rows != []:
        raise Exception("There's some rows where the number of intensity values and retention times are not equal. This indicates a lack of data or a mistake in the original data.")

    return pd.DataFrame(rows)

def remove_same_scan_duplicates(df):
    rows = []
    for index, row in df.iterrows():
        _, idx = np.unique(row['ms1_nb'], return_index=True)
        unique_ms1_nb = [row['ms1_nb'][i] for i in sorted(idx)]
        new_intensity = []
        new_retention_time = []
        for ms1_nb in unique_ms1_nb:
            indices = [i for i, x in enumerate(row['ms1_nb']) if x == ms1_nb]
            corresponding_intensity = [row['intensity'][i] for i in indices]
            corresponding_retention_time = [row['retention time'][i] for i in indices]
            new_intensity.append(np.max(corresponding_intensity))                # Taking the max intensity as it's the most important for the precursor and we avoid diluting the intensity.
            new_retention_time.append(np.mean(corresponding_retention_time))     # Always the same retention time as it's the same scan
        rows.append({'mz': row['mz'], 'ms1_nb': unique_ms1_nb, 'intensity': new_intensity, 'retention time': new_retention_time, 'sort_column': row['sort_column']})
    new_df = pd.DataFrame(rows)
    return new_df

# def prep_for_correlation(df_ms1, df_ms2):
#     # add row called min_intensity and max_intensity for both dataframes
#     df_ms1['min_intensity'] = df_ms1['intensity'].apply(min)
#     df_ms1['max_intensity'] = df_ms1['intensity'].apply(max)
#     df_ms2['min_intensity'] = df_ms2['intensity'].apply(min)
#     df_ms2['max_intensity'] = df_ms2['intensity'].apply(max)

#     # add row called lenght_intensity_array for both dataframes
#     df_ms1['lenght_intensity_array'] = df_ms1['intensity'].apply(len)
#     df_ms2['lenght_intensity_array'] = df_ms2['intensity'].apply(len)


def correlate_fragments_to_precursor_trying(df_ms1, df_ms2, pearson_min_score=0.5, time_between_ms1_scans=5.0, points_across_the_peak_lim=3):

    # Initialize a list to store the results
    results = []

    # set limits
    pearson_correlation_coeficient_lim_pos = pearson_min_score
    pearson_correlation_coeficient_lim_neg = -pearson_min_score
    ms1_nb = 0   #to keep track of the correlated fragments to precursor

    # Iterate over each row in df_ms1
    for index1, row1 in tqdm(df_ms1.iterrows(), total=len(df_ms1)):
        
        # Skip ms1 entries with intensity arrays that are too short to be correlated
        if len(row1['intensity']) < points_across_the_peak_lim:
            continue

        #update ms1_nb
        ms1_nb += 1

        # Add the ms1 entry to the results
        results.append({'mz': row1['mz'], 'intensity': row1['intensity'], 'retention time': row1['retention time'], 'type': 'ms1', "ms1_nb": ms1_nb})

        # Calculate these once to avoid repeated calculations
        min_retention_time1 = np.min(row1['retention time'])
        max_retention_time1 = np.max(row1['retention time'])

        # Iterate over each row in df_ms2
        for row2 in df_ms2.itertuples():
            
            # Check if the retention times overlap
            max_retention_time2 = np.max(row2._asdict()['retention time'])
            if min_retention_time1 > (max_retention_time2 + time_between_ms1_scans):
                continue

            #Since the dataframe is sorted according to retention time we can add logic that will not keep looking through the dataframe if the ms2 entries have higher retention times than the ms1 we are looking at.
            min_retention_time2 = np.min(row2._asdict()['retention time'])
            if (max_retention_time1 + time_between_ms1_scans) < min_retention_time2:
                break

            ###Correlation of fragments to precursors###

            # Filter ms1_nb and ms2_nb to only include matching numbers
            matching_nb = set(row1['ms1_nb']).intersection(set(row2._asdict()['ms1_nb']))
            
            # Check if the minimum amount of points match
            if len(matching_nb) < points_across_the_peak_lim:
                continue

            # Filter the intensity and retention time (can be excluded) lists based on the matching_nb
            indices1 = [i for i, x in enumerate(row1['ms1_nb']) if x in matching_nb]
            indices2 = [i for i, x in enumerate(row2._asdict()['ms1_nb']) if x in matching_nb]
            filtered_intensity1 = [row1['intensity'][i] for i in indices1]
            filtered_intensity2 = [row2._asdict()['intensity'][i] for i in indices2]
            #filtered_retention_time1 = [row1['retention time'][i] for i in indices1]
            filtered_retention_time2 = [row2._asdict()['retention time'][i] for i in indices2]

            # Calculate the Pearson correlation coefficient for the intensity arrays
            pearson_correlation_coeficient, p_value = pearsonr(filtered_intensity1, filtered_intensity2)

            # checking if the correlation is any good
            if not pearson_correlation_coeficient > pearson_correlation_coeficient_lim_pos or pearson_correlation_coeficient < pearson_correlation_coeficient_lim_neg:
                continue
            
            # Add the fragment to the results list right after its precursor
            results.append({'mz': row2._asdict()['mz'], 'intensity': filtered_intensity2, 'retention time': filtered_retention_time2, 'correlation': pearson_correlation_coeficient, 'type': 'ms2', "ms1_nb": ms1_nb})

    # Convert the results list to a DataFrame and return it
    return pd.DataFrame(results)


def vectorized_correlations(df_ms1, df_ms2, pearson_min_score=0.5, time_between_ms1_scans=5.0, points_across_the_peak_lim=3):
    

    # Initialize a list to store the results
    results = []

    # set limits
    pearson_correlation_coeficient_lim_pos = pearson_min_score
    ms1_nb = 0   #to keep track of the correlated fragments to precursor

    # Calculate the min and max retention time for each row in df_ms1 and df_ms2
    df_ms1['min_retention_time'] = df_ms1['retention time'].apply(np.min)
    df_ms1['max_retention_time'] = df_ms1['retention time'].apply(np.max)
    df_ms2['min_retention_time'] = df_ms2['retention time'].apply(np.min)
    df_ms2['max_retention_time'] = df_ms2['retention time'].apply(np.max)


    #renaming columns to work with the itertuples function
    df_ms1 = df_ms1.rename(columns={"retention time": "retention_time"})
    
    # Iterate over each row in df_ms1
    for row1 in tqdm(df_ms1.itertuples(), total=len(df_ms1)):

        #update ms1_nb
        ms1_nb += 1

        # Add the ms1 entry to the results
        results.append({'mz': row1.mz, 'intensity': row1.intensity, 'retention time': row1.retention_time, 'type': "ms1", "ms1_nb": ms1_nb})

        # Making a copy of the ms2 dataframe for applying to it without changing the original
        working_df = df_ms2.copy()
        
        # Filter out rows where the retention times do not overlap
        working_df = working_df[(row1.min_retention_time <= (working_df['max_retention_time'] + time_between_ms1_scans)) & ((row1.max_retention_time + time_between_ms1_scans) >= working_df['min_retention_time'])]

        # find matching numbers for all rows in df_ms2 using vectorization
        ms1_nb_set = set(row1.ms1_nb)
        #print("ms1_nb_set", ms1_nb_set)
        working_df['matching_nb'] = working_df['ms1_nb'].apply(lambda x: set(x).intersection(ms1_nb_set))

        # Filter rows where the matching_nb is under points_across_the_peak_lim
        working_df = working_df[(working_df['matching_nb'].apply(type) == set) & (working_df['matching_nb'].apply(len) >= points_across_the_peak_lim)]
        
        if working_df["matching_nb"].empty:
            continue

        def get_indices1(row_ms1, row_ms2):
            indices = [i for i, n in enumerate(row_ms1) if n in row_ms2["matching_nb"]]
            return indices
        
        row_1 = row1.ms1_nb
        working_df['indices1'] = working_df.apply(lambda x: get_indices1(row_1,x), axis=1)


        def get_indices2(row):
            return [i for i, x in enumerate(row["ms1_nb"]) if x in row["matching_nb"]]

        working_df['indices2'] = working_df.apply(get_indices2, axis=1)
        

        # finding filtered intensity arrays
        working_df['filtered_intensity1'] = working_df.apply(lambda x: [row1.intensity[i] for i in x['indices1']], axis=1)
        working_df['filtered_intensity2'] = working_df.apply(lambda x: [x['intensity'][i] for i in x['indices2']], axis=1)

        # finding the filtered rentention time arrays
        working_df['filtered_retention_time2'] = working_df.apply(lambda x: [x['retention time'][i] for i in x['indices2']], axis=1)

        # Calculate the Pearson correlation coefficient for the intensity arrays
        working_df['pearson_correlation_coeficient'] = working_df.apply(lambda x: pearsonr(x['filtered_intensity1'], x['filtered_intensity2'])[0], axis=1)

        # Filter out rows where the Pearson correlation coefficient is not good enough
        working_df = working_df[(working_df['pearson_correlation_coeficient'] > pearson_correlation_coeficient_lim_pos)]

        # Add the fragment to the results list right after its precursor
        for row2 in working_df.itertuples():
            results.append({'mz': row2.mz, 'intensity': row2.filtered_intensity2, 'retention time': row2.filtered_retention_time2, 'correlation': row2.pearson_correlation_coeficient, 'type': 'ms2', "ms1_nb": ms1_nb})
        
    # Convert the results list to a DataFrame and return it
    return pd.DataFrame(results)


def correlate_fragments_to_precursors(df_ms1, df_ms2, pearson_min_score=0.5, time_between_ms1_scans=5.0, points_across_the_peak_lim=3):

    # Initialize a list to store the results
    results = []

    # set limits
    pearson_correlation_coeficient_lim_pos = pearson_min_score
    pearson_correlation_coeficient_lim_neg = -pearson_min_score
    ms1_nb = 0   #to keep track of the correlated fragments to precursor

    # Iterate over each row in df_ms1
    for index1, row1 in tqdm(df_ms1.iterrows(), total=len(df_ms1)):
        
        # Skip ms1 entries with intensity arrays that are too short to be correlated
        if len(row1['intensity']) < points_across_the_peak_lim:
            continue

        #update ms1_nb
        ms1_nb += 1

        # Add the ms1 entry to the results
        results.append({'mz': row1['mz'], 'intensity': row1['intensity'], 'retention time': row1['retention time'], 'type': 'ms1', "ms1_nb": ms1_nb})

        # Iterate over each row in df_ms2
        for index2, row2 in df_ms2.iterrows():
            
            # Check if the retention times overlap
            if (np.min(row1['retention time'])) > (np.max(row2['retention time'])+time_between_ms1_scans):
                continue

            #Since the dataframe is sorted acording to retention time we can add logic that will not keep looking through the dataframe if the ms2 entries have higher retention times than the ms1 we are looking at.
            if (np.max(row1['retention time'])+time_between_ms1_scans) < np.min(row2['retention time']):
                break


            ###Correlation of fragments to precursors###

            # Filter ms1_nb and ms2_nb to only include matching numbers
            matching_nb = list(set(row1['ms1_nb']).intersection(row2['ms1_nb']))
            
            # Check if the minimum amount of points match
            if len(matching_nb) < points_across_the_peak_lim:
                continue

            # Filter the intensity and retention time (can be excluded) lists based on the matching_nb
            indices1 = [i for i, x in enumerate(row1['ms1_nb']) if x in matching_nb]
            indices2 = [i for i, x in enumerate(row2['ms1_nb']) if x in matching_nb]
            filtered_intensity1 = [row1['intensity'][i] for i in indices1]
            filtered_intensity2 = [row2['intensity'][i] for i in indices2]
            #filtered_retention_time1 = [row1['retention time'][i] for i in indices1]
            filtered_retention_time2 = [row2['retention time'][i] for i in indices2]

            # Calculate the Pearson correlation coefficient for the intensity arrays
            pearson_correlation_coeficient, p_value = pearsonr(filtered_intensity1, filtered_intensity2)

            # checking if the correlation is any good
            if not pearson_correlation_coeficient > pearson_correlation_coeficient_lim_pos or pearson_correlation_coeficient < pearson_correlation_coeficient_lim_neg:
                continue
            
            # Add the fragment to the results list right after its precursor
            results.append({'mz': row2['mz'], 'intensity': filtered_intensity2, 'retention time': filtered_retention_time2, 'correlation': pearson_correlation_coeficient, 'type': 'ms2', "ms1_nb": ms1_nb})

    # Convert the results list to a DataFrame and return it
    return pd.DataFrame(results)

def filter_fragment_count_pr_precursor(df, fragment_count_pr_precursor=3):
    # Initialize a list to store the results
    valid_entries = []

    # Iterate over the DataFrame
    i = 0
    
    while i < len(df):
        # Check if the current entry is 'ms1'
        if df.iloc[i]['type'] == 'ms1':
            # Initialize a counter for 'ms2' entries
            ms2_count = 0
            # Check the following entries
            for j in range(i+1, len(df)):
                if df.iloc[j]['type'] == 'ms2':
                    ms2_count += 1
                else:
                    break
            
            # If there are 3 or more 'ms2' entries, append them to the valid_entries DataFrame
            if ms2_count >= fragment_count_pr_precursor:
                valid_entries.extend(df.iloc[i:i+ms2_count+1].values.tolist())
            # Skip the checked entries
            i += ms2_count + 1
        else:
            # If the current entry is not 'ms1', skip it
            i += 1

    return pd.DataFrame(valid_entries, columns=df.columns)


def feature_detection_ms1(input_file: str):
    
    # Prepare data loading (save memory by only
    # loading MS1 spectra into memory)
    options = oms.PeakFileOptions()
    options.setMSLevels([1])
    fh = oms.MzMLFile()
    fh.setOptions(options)

    input_map = oms.MSExperiment()

    # Load data
    fh.load(input_file, input_map)
    input_map.updateRanges()

    ff = oms.FeatureFinder(cutoff=300, mz_tol=0.1, rt_tol=0.1)
    ff.setLogType(oms.LogType.CMD)

    # Do feature detection
    name = "centroided"
    features = oms.FeatureMap()
    seeds = oms.FeatureMap()
    params = oms.FeatureFinder().getParameters(name)
    ff.run(name, input_map, features, params, seeds)


    features.setUniqueIds()

    f0 = features[0]
    print("mz", f0.getMZ(), "RT", f0.getRT(), "intensity", f0.getIntensity())

    features_sorted = []
    for f in features:
        features_sorted.append((f.getMZ(), f.getRT(), f.getIntensity(), f.getUniqueId()))

    features_sorted = sorted(features_sorted, key = lambda element : element[1])

    return features_sorted


def peak_picking_for_feature_detection_ms1(df, features, time_between_scans=5, maximum_peak_elution_time=45):

    rows = []

    # Looping through the features build be arrays of 2d peaks
    for tuple in tqdm(features, total=len(features)):
        #print(row)
        
        # get the mz of the feature
        feature_mz = tuple[0]

        # get the retention time of the feature
        feature_rt = tuple[1]

        #print("Feature mz: ", feature_mz)
        #print("Feature rt: ", feature_rt)
        # Select the all rows in the df that are at the feature_mz with some tolerance and at feature_rt with some tolerance, tolerance are variable. The rows should be added to a list of dictionaries with the mz, ms1_nb, intensity and retention time.
        mz_tolerance = 0.5
        retention_time_tolerance = 30


        temp_rows = []
        prev_rt = None
        ms1_nbs = []
        intensities = []
        retention_times = []
        
        # rename the retention time column to retention_time
        df.rename(columns = {'retention time': 'retention_time'}, inplace = True)
        
        # sort the dataframe by retention time
        df = df.sort_values(by=['retention_time'], ascending=True)

        for row in df.itertuples():
            if abs(row.mz - feature_mz) < mz_tolerance and abs(row.retention_time - feature_rt) < retention_time_tolerance:

                # If the difference between the current and previous retention times is more than 5, append the temporary list to rows
                if prev_rt is not None and row.retention_time - prev_rt > time_between_scans:
                    rows.append({'mz': feature_mz, "ms1_nb" : ms1_nbs, 'intensity': intensities, 'retention time': retention_times})

                    print("got a hit")
                    # reset the temporary lists
                    ms1_nbs = []
                    intensities = []
                    retention_times = []


                ms1_nbs.append(row.ms1_nb)
                intensities.append(row.intensity)
                retention_times.append(row.retention_time)
                prev_rt = row.retention_time

        # Append the last group of rows
        if ms1_nbs != []:
            rows.append({'mz': feature_mz, "ms1_nb" : ms1_nbs, 'intensity': intensities, 'retention time': retention_times})

    return pd.DataFrame(rows)