# -*- coding: utf-8 -*-
"""
Created on Tue Nov  4 09:47:48 2025
2. Behavioural extraction of fiber photometry data
For sucrose conditioning paradigm - 1 reward 

This script will perform the following preprocessings steps:
- General code for loading
- determining sampling rate from the parquet files
- Exp specific: determine time locks
- Exp specific: determine behavioural outcomes
    Raw or processed such as minimal bout duration
- general code for saving

Inputs: 
    preprocessed data in parquet file
    Determine baselines
Outputs: 
    DF with colums for each time lock - containing name and indexes of time lock and baselines
    Behavioural readouts in excel or csv

Checked behavioural output of 5 mice and they are comparible with earlier scripts

Step 2 out of 4 in the photometry analysis pipeline
Meye lab, Karlijn Kooij
Version 7_11_25


@author: kkooij
"""
# enter the following parameters per analysis

# enter the directory where the raw data files are stored
# files_dir_raw = 'D://Ghrelin - D drive//Analysis scripts//New pipeline//test_data//'

files_dir_raw = 'L://adanlab//Ongoing//kkooij//ghrelin manu//Data//vgat & pitx cre//3 sucrose RPE//Pitx//Data//all data//'
#files_dir_raw = 'L://adanlab//Ongoing//kkooij//ghrelin manu//Data//vgat & pitx cre//3 sucrose RPE//Vgat//Data//'

#files_dir_raw = 'D://Ghrelin - D drive//Data//Participated trials//data//Only_w2//'
files_dir_preprocessed = files_dir_raw + '1. Preprocessed//'

# baseline determination 
baseline_sec = [-10,0] #note here which seconds the baseline should be

# behavioural specifics
nr_trials_analysed = [0,10000000]   #first tone will be deleted later as this is not good data

min_lick_delay_sec = 5    #add 5 if you want licks during cue to be counted as participated cues, of 0 if only licks during reward bins count
max_lick_delay_sec = 10   #maximum time to lick to still receive a reward
min_durations = { "first": 4, "sec": 0.01}   # Example: first DIO requires at least 0.5s licking, second requires at least 1.2s   previous 0.7 sec
min_inter_bout_interval_sec = 0.01  #  min time that considers two different licking bouts (in seconds)


TTL_1_desc = {
    "S": (10.5, "10% sucrose", 'own baseline'), #was 10.5 voor DA    - 5 voor gaba eerste 3
    "M": (12,  "30% sucrose", 'own baseline'), # was 12   voor DA    - 7 voor gaba
    "L": (14,  "3% sucrose", 'own baseline' )} # was 14   voor DA    - 9 voor gaba

#lick rate specs
time_window = [-5, 20] # seconds relative to onset
bin_size = 0.20 # 50 ms bins

# not included
# excl_lick_baseline = 0 # if 1 the code will exclude trials in which the mouse licks during baseline.
#parameters over time
#analysis_blocks = [0,20,40,60]   # or [0,20,40,60,78] with 80 trials   ONLY WORKS with minimal 60 trial (of which 58 analysed) 


###############################################################
# import packages 
import pandas as pd
import numpy as np 
import os
import time
import statistics
from scipy.signal import butter,  filtfilt
import h5py # Make sure to install the library
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Qt5Agg')   # or 'TkAgg' if Qt is not available
import matplotlib.pyplot as plt
plt.ioff() 
from matplotlib.widgets import SpanSelector, Button
import warnings
from datetime import datetime
################################################################
# functions general


def load_files_from_folder(files_dir, extention):
  # inport data from 1 folder with a certain extention 
  list_files = os.listdir(files_dir)
  files = []
  for i_files in range(len(list_files)):
      if list_files[i_files].endswith(extention):
          files.append(list_files[i_files])
          
  return files


def create_folder_with_note(base_path, description):
    """
    Create a folder if it does not exist, and write a text file inside
    noting the date and using the folder name as description.
    
    Parameters:
        base_path (str): Path where the folder should be created.
        description (str): Name of the folder to create and description for the note.
        
    Returns:
        str: Full path to the created folder.
    """
    full_path = os.path.join(base_path, description)
    
    # Create folder if it does not exist
    if not os.path.exists(full_path):
        os.makedirs(full_path)
        print(f"Folder created: {full_path}")
    else:
        print(f"Folder already exists: {full_path}")
    
    # Write note with current date and folder name as description
    date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    date = datetime.now().strftime("%Y-%m-%d")
    note_path = os.path.join(full_path, f"{description} on {date}.txt")
    with open(note_path, "w") as f:
        f.write(f"{description} on {date_str}\n")
    
    return full_path

def parse_ID(ID):
    """
    Extract the original information from an ID string created by `create_ID`.
    Handles cases where the last part (genotype) may contain a file extension.

    Args:
        ID (str): ID string in the format name_condition_week_batch_genotype[.ext]

    Returns:
        dict: Dictionary with keys 'name', 'condition', 'week', 'batch', 'genotype'
    """
    # Remove file extension if present
    ID_no_ext = os.path.splitext(ID)[0]  # removes '.parquet', '.csv', etc.
    
    parts = ID_no_ext.split("_")
    
    if len(parts) < 5:
        raise ValueError(f"ID string '{ID}' does not contain enough parts to parse.")

    # Recover the fields
    name = parts[0]
    hemiphere = parts[1]
    condition = parts[2]
    week = parts[3]
    batch = parts[4]
    genotype = "_".join(parts[5:])  # in case genotype had underscores

    return {
        "name": name,
        "hemisphere": hemiphere,
        "condition": condition,
        "week": week,
        "batch": batch,
        "genotype": genotype
    }, ID_no_ext


def determine_sampling_rate_parquet(files, files_dir):
    """
    Determine the median sampling rate across preprocessed Parquet files.
    
    Parameters:
        files (list): List of Parquet filenames.
        files_dir (str): Directory containing Parquet files.
        
    Returns:
        int: Median sampling rate (samples/sec).
    """
    sample_rates = []

    for file_name in files:
        df = pd.read_parquet(os.path.join(files_dir, file_name))
        time = df['Time(s)'].values
        rate = len(time) / (time[-1] - time[0])
        sample_rates.append(rate)

    median_rate = int(round(statistics.median(sample_rates)))
    print(f"Sampling rates across files: min={min(sample_rates):.2f}, max={max(sample_rates):.2f}, median={median_rate}")
    
    return median_rate

##################
# specific functions

def process_DIO(DIO_input, df=None, threshold=0.8):
    """
    Process DIO array(s) and return onset/offset/activated indexes.
    
    Args:
        DIO_input (np.ndarray or dict): Single DIO array or dictionary of arrays.
        df (pd.DataFrame, optional): Not used here, kept for compatibility.
        threshold (float): Voltage threshold to consider DIO "on".
    
    Returns:
        dict: keys are names (or 'DIO') and values are tuples (onsets, offsets, activated_indexes)
    """
    if isinstance(DIO_input, dict):
        return {key: DIO_get_onset_offset(arr, df=df, threshold=threshold)
                for key, arr in DIO_input.items()}
    else:
        # Single array, return as a dictionary with key 'DIO'
        return {'DIO': DIO_get_onset_offset(DIO_input, df=df, threshold=threshold)}



def get_sub_list(my_list):
    """
    Splits a sorted list of integers into consecutive sublists.
    """
    if not my_list:
        return []

    result = []
    temp = [my_list[0]]

    for prev, curr in zip(my_list, my_list[1:] + [None]):
        if curr is not None and curr == prev + 1:
            temp.append(curr)
        else:
            result.append(temp)
            if curr is not None:
                temp = [curr]

    return result


def clean_edges(index_list):
    """
    Remove suspicious first or last elements in a list of indices.
    Assumes the list is sorted.
    """
    if not index_list:
        return []

    # Remove first elements if start at 0 or consecutive anomalies at the beginning
    idx_array = np.array(index_list)
    if idx_array[0] == 0:
        diff = np.diff(idx_array)
        # Find first large jump (>1)
        jump_idx = np.argmax(diff > 1) + 1
        idx_array = idx_array[jump_idx:]

    # Remove last elements if consecutive anomalies at the end
    if idx_array[-1] == idx_array[-2] + 1:
        # If needed, additional logic could go here to remove false positives at the end
        pass

    return idx_array.tolist()


def DIO_get_onset_offset(DIO_array, df=None, threshold=0.8):
    """
    Extracts onset and offset indices from a DIO array.
    
    Args:
        DIO_array (np.ndarray): Array of DIO signal values.
        df (pd.DataFrame): Optional, not used here.
        threshold (float): Voltage threshold to consider DIO "on".
        
    Returns:
        TTL_onsets, TTL_offsets, DIO_activated_indexes
    """
    # Get all indices where DIO is above threshold
    activated_indexes = np.where(DIO_array > threshold)[0].tolist()

    # Clean edges
    activated_indexes = clean_edges(activated_indexes)

    # Split into bouts of consecutive indices
    bouts = get_sub_list(activated_indexes)

    # Extract onsets and offsets
    TTL_onsets = np.array([b[0] for b in bouts], dtype=int)
    TTL_offsets = np.array([b[-1] for b in bouts], dtype=int)

    return TTL_onsets, TTL_offsets, activated_indexes


###########################
# exp specific

def filter_by_trial_range(DIO_results, trial_onsets, nr_trials_analysed=[0, 60]):
    """
    Filter DIO onsets/offsets to only include events that fall within the index range
    defined by the start and end trial numbers.
    
    Args:
        DIO_results (dict): Dictionary with keys as DIO names and values as tuples
                            (onsets, offsets, activated_indexes)
        trial_onsets (np.ndarray): Array of trial start indices (e.g., first_onsets)
        nr_trials_analysed (list or tuple): [start_trial, end_trial] (0-based indexing)
    
    Returns:
        dict: Filtered DIO_results with the same structure.
    """
    start_trial, end_trial = nr_trials_analysed
    n_trials = len(trial_onsets)
    
    if n_trials < end_trial:
        print(f"Note: Only {n_trials} trials available, fewer than requested {end_trial}. Using all available trials.")
        end_trial = n_trials  # adjust end_trial to maximum available
    
    # Determine actual index range for filtering
    start_idx = trial_onsets[start_trial]
    end_idx = trial_onsets[end_trial - 1] if end_trial - 1 < n_trials else trial_onsets[-1]

    filtered_results = {}
    
    for name, (onsets, offsets, activated) in DIO_results.items():
        # Keep only onsets/offsets within the trial index range
        mask = (onsets >= start_idx) & (onsets <= end_idx)
        filtered_onsets = onsets[mask]
        filtered_offsets = offsets[mask]
        
        filtered_results[name] = (
            filtered_onsets,
            filtered_offsets,
            activated  # keep full activated indexes
        )
    
    return filtered_results

def concatenate_DIO_bouts(DIO_dict, sampling_rate, min_inter_bout_interval_sec=1):
    """
    Concatenate bouts for all DIO channels in a dictionary based on a minimum inter-bout interval.
    
    Args:
        DIO_dict (dict): Dictionary with keys = DIO names, values = tuples
                         (onsets, offsets, activated_indexes)
        sampling_rate (float): Sampling rate in Hz
        min_inter_bout_interval_sec (float): Merge bouts closer than this interval (seconds)
    
    Returns:
        dict: New dictionary with concatenated onsets and offsets for each DIO channel.
              The structure is the same as input: {name: (concat_onsets, concat_offsets, activated)}
    """
    
    # Convert minimum inter-bout interval from seconds to samples
    min_gap_samples = int(min_inter_bout_interval_sec * sampling_rate)
    
    # Initialize new dictionary for concatenated results
    concatenated_DIO = {}
    
    # Loop over each DIO channel
    for name, (onsets, offsets, activated) in DIO_dict.items():
        if len(onsets) == 0:
            # No events for this channel
            concatenated_DIO[name] = (np.array([], dtype=int), np.array([], dtype=int), activated)
            continue
        
        # Initialize lists with first bout
        concat_onsets = [onsets[0]]
        concat_offsets = [offsets[0]]
        
        # Loop through remaining bouts
        for i in range(1, len(onsets)):
            gap = onsets[i] - concat_offsets[-1]  # gap between previous offset and current onset
            
            if gap < min_gap_samples:
                # Bouts are too close → merge by extending the previous offset
                concat_offsets[-1] = offsets[i]
            else:
                # Bouts are far enough → keep separate
                concat_onsets.append(onsets[i])
                concat_offsets.append(offsets[i])
        
        # Save results for this channel
        concatenated_DIO[name] = (np.array(concat_onsets, dtype=int),
                                  np.array(concat_offsets, dtype=int),
                                  activated)
    
    return concatenated_DIO

    



def filter_bouts_by_duration_all(DIO_dict, sampling_rate, min_durations):
    """
    Filters bouts for all DIO channels in a dictionary based on channel-specific minimum durations.
    
    Args:
        DIO_dict (dict): {name: (onsets, offsets, activated_indexes)}
        sampling_rate (float): Sampling rate in Hz
        min_durations (dict): {name: min_duration_sec} per DIO channel.
                              If a name is not provided, default = 1 second.
                              
    Returns:
        dict: {name: (filtered_onsets, filtered_offsets, activated_indexes)}
    """
    
    filtered_dict = {}
    
    for name, (onsets, offsets, activated) in DIO_dict.items():
        # If no custom minimum provided → default 1 sec
        min_duration_sec = min_durations.get(name, 1)
        min_duration_samples = int(min_duration_sec * sampling_rate)
        
        if len(onsets) == 0:
            filtered_dict[name] = (onsets, offsets, activated)
            continue
        
        durations = offsets - onsets + 1
        valid = durations >= min_duration_samples
        
        filtered_onsets = onsets[valid]
        filtered_offsets = offsets[valid]
        
        filtered_dict[name] = (filtered_onsets, filtered_offsets, activated)
    
    return filtered_dict

# def first_and_last_lick_after_tone(tone_offsets, lick_onsets, lick_offsets, sampling_rate, min_delay_sec, max_delay_sec=10):
#     """
#     For each tone, find the first lick and the last lick of the licking bout within a time window.
    
#     Args:
#         tone_offsets (np.ndarray): Array of tone offset indices
#         lick_onsets (np.ndarray): Array of lick bout onsets (indices)
#         lick_offsets (np.ndarray): Array of lick bout offsets (indices)
#         sampling_rate (float): Sampling rate in Hz
#         max_delay_sec (float): Maximum time after tone to consider licks
    
#     Returns:
#         first_licks (np.ndarray): First lick index after each tone (np.nan if none)
#         last_licks (np.ndarray): Last lick of the bout corresponding to that first lick (np.nan if none)
#     """
#     max_delay_samples = int(max_delay_sec * sampling_rate)
#     min_delay_samples =  int(min_delay_sec * sampling_rate)
    
#     first_licks = np.full_like(tone_offsets, fill_value=np.nan, dtype=float)
#     last_licks = np.full_like(tone_offsets, fill_value=np.nan, dtype=float)
    
#     lick_idx = 0  # pointer in lick_onsets
    
#     for i, tone_off in enumerate(tone_offsets):
#         window_start = tone_off - min_delay_samples
#         window_end = tone_off + max_delay_samples
        
#         # Move pointer to first lick onset within the window
#         while lick_idx < len(lick_onsets) and lick_onsets[lick_idx] < window_start:
#             lick_idx += 1
        
#         if lick_idx >= len(lick_onsets) or lick_onsets[lick_idx] > window_end:
#             continue  # no lick in window
        
#         # First lick found
#         if lick_idx >= len(lick_onsets) or lick_onsets[lick_idx] < tone_off:
#             first_licks[i] = tone_off
#         else :
#             first_licks[i] = lick_onsets[lick_idx]
            
#         # Last lick of this bout is the offset corresponding to the first lick
#         last_licks[i] = lick_offsets[lick_idx]
        
#         # Optionally move pointer to next bout for next tone
#         lick_idx += 1
    
#     return first_licks, last_licks

def first_and_last_lick_after_tone(tone_offsets, lick_onsets, lick_offsets, sampling_rate, min_delay_sec, max_delay_sec=10):
    """
    For each tone, find the first lick and the last lick of the licking bout within a time window.
    If the first lick is before the tone, use tone offset as the first lick.
    """
    max_delay_samples = int(max_delay_sec * sampling_rate)
    min_delay_samples = int(min_delay_sec * sampling_rate)
    
    first_licks = np.full_like(tone_offsets, fill_value=np.nan, dtype=float)
    last_licks = np.full_like(tone_offsets, fill_value=np.nan, dtype=float)
    
    lick_idx = 0  # pointer in lick_onsets
    
    for i, tone_off in enumerate(tone_offsets):
        window_start = tone_off - min_delay_samples
        window_end = tone_off + max_delay_samples
        
        # Move pointer to first lick onset within the window
        while lick_idx < len(lick_onsets) and lick_onsets[lick_idx] < window_start:
            lick_idx += 1
        
        if lick_idx >= len(lick_onsets) or lick_onsets[lick_idx] > window_end:
            continue  # no lick in window
        
        # Check if first lick is before tone
        if lick_onsets[lick_idx] < tone_off:
            first_licks[i] = tone_off
            # Take last lick of the previous bout
            last_licks[i] = lick_offsets[lick_idx]
        else:
            # First lick after tone
            first_licks[i] = lick_onsets[lick_idx]
            last_licks[i] = lick_offsets[lick_idx]
        
        # Move pointer to next bout
        lick_idx += 1
    
    return first_licks, last_licks

def split_participated_trials(tone_onsets, first_licks):
    """
    Split tone trials into participated and non-participated based on licking response.
    
    Args:
        tone_onsets (np.ndarray): Tone onset (or offset) indices for each trial.
        first_licks (np.ndarray): First lick indices after each tone (NaN if none).
    
    Returns:
        participated_tones (np.ndarray): Tone onsets where licking occurred.
        non_participated_tones (np.ndarray): Tone onsets where no lick occurred.
        
        participated_trial_idx (np.ndarray): Indices of participated trials.
        non_participated_trial_idx (np.ndarray): Indices of non-participated trials.
    """
    
    # Trial index masks
    participated_mask = ~np.isnan(first_licks)
    non_participated_mask = np.isnan(first_licks)
    
    # Split tone onsets
    participated_tones = tone_onsets[participated_mask]
    non_participated_tones = tone_onsets[non_participated_mask]
    
    # Return trial indices too (useful for future alignment)
    participated_trial_idx = np.where(participated_mask)[0]
    non_participated_trial_idx = np.where(non_participated_mask)[0]
    
    return participated_tones, non_participated_tones, participated_trial_idx, non_participated_trial_idx

def first_lick_during_cue(tone_onsets, lick_onsets, sampling_rate, cue_duration_sec=5):
    """
    Detect the first lick that happens during the cue period for each tone.
    
    Args:
        tone_onsets (np.ndarray): Indices where tones start.
        lick_onsets (np.ndarray): Indices of concatenated lick onsets.
        sampling_rate (float): Sampling rate in Hz.
        cue_duration_sec (float): Duration of cue presentation in seconds (default: 5s).
        
    Returns:
        first_lick_indices (np.ndarray): Absolute indices of first lick during cue (NaN if none).
        lick_latency_sec (np.ndarray): Latency relative to tone onset in seconds (NaN if none).
    """
    
    first_lick_indices = np.full(len(tone_onsets), np.nan, dtype=float)
    lick_latency_sec   = np.full(len(tone_onsets), np.nan, dtype=float)
    
    cue_duration_samples = int(cue_duration_sec * sampling_rate)
    
    lick_pointer = 0
    
    for i, tone_i in enumerate(tone_onsets):
        cue_end = tone_i + cue_duration_samples
        
        # Move pointer to first lick that could match this tone
        while lick_pointer < len(lick_onsets) and lick_onsets[lick_pointer] < tone_i:
            lick_pointer += 1
        
        # If we ran out of licks → stop
        if lick_pointer >= len(lick_onsets):
            break
        
        # Check whether the lick occurs during the cue
        if tone_i <= lick_onsets[lick_pointer] <= cue_end:
            first_lick_indices[i] = lick_onsets[lick_pointer]
            lick_latency_sec[i]   = (lick_onsets[lick_pointer] - tone_i) / sampling_rate
    
    return first_lick_indices, lick_latency_sec


def get_licks_no_reward(DIO_dict, tone_onsets, sampling_rate, time_before_tone_sec=5, time_after_tone_sec=25):
    """
    Determine licking bouts that occur outside of cue/reward windows, using pre-concatenated bouts.
    
    Args:
        DIO_dict (dict): Dictionary with lick DIOs, e.g., {"sec": (onsets, offsets, activated)}
        tone_onsets (np.ndarray): Tone onset indices (all trials)
        sampling_rate (float): Sampling rate in Hz
        time_before_tone_sec (float): Exclude licks within this time before tone
        time_after_tone_sec (float): Exclude licks within this time after tone
    
    Returns:
        dict: {DIO_name: (lick_onsets_no_reward, bout_durations_sec)}
    """
    
    results = {}
    
    for name, (lick_onsets, lick_offsets, activated) in DIO_dict.items():
        lick_onsets = np.array(lick_onsets)
        lick_offsets = np.array(lick_offsets)
        
        # Exclude licks near tones
        mask_keep = np.ones_like(lick_onsets, dtype=bool)
        for tone_i in tone_onsets:
            mask_keep &= (lick_onsets < tone_i - int(time_before_tone_sec * sampling_rate)) | \
                         (lick_onsets > tone_i + int(time_after_tone_sec * sampling_rate))
        
        lick_onsets_no_reward = lick_onsets[mask_keep]
        lick_offsets_no_reward = lick_offsets[mask_keep]
        
        # Compute bout durations
        bout_durations_sec = (lick_offsets_no_reward - lick_onsets_no_reward) / sampling_rate
        
        results[name] = (lick_onsets_no_reward, bout_durations_sec)
    
    return results


def time_lock_omission_lick_offset(first_licks, last_licks, tone_offsets, non_participated_idx, sampling_rate):
    """
    Assign expected lick timing for omitted trials (non-participated), relative to tone offset.
    
    Args:
        first_licks (np.ndarray): First lick indices for all trials (NaN for non-participated)
        last_licks (np.ndarray): Last lick indices for all trials (NaN for non-participated)
        tone_offsets (np.ndarray): Tone offset indices for all trials
        non_participated_idx (np.ndarray): Indices of non-participated trials
        sampling_rate (float): Sampling rate in Hz
    
    Returns:
        expected_first_lick_indices (np.ndarray): Absolute indices for omitted trials
        expected_last_lick_indices (np.ndarray): Absolute indices for omitted trials
        expected_first_lick_sec (np.ndarray): Latencies relative to tone offset (seconds)
        expected_last_lick_sec (np.ndarray): Latencies relative to tone offset (seconds)
    """
    
    # Compute average latencies from participated trials
    participated_mask = ~np.isnan(first_licks)
    avg_first_latency = int(round(np.nanmean(first_licks[participated_mask] - tone_offsets[participated_mask])))
    avg_last_latency  = int(round(np.nanmean(last_licks[participated_mask]  - tone_offsets[participated_mask])))
    
    # Assign expected indices for omitted trials
    expected_first_lick_indices = tone_offsets[non_participated_idx] + avg_first_latency
    expected_last_lick_indices  = tone_offsets[non_participated_idx] + avg_last_latency
    
    # Latencies in seconds as arrays
    expected_first_lick_sec = np.full(len(non_participated_idx), avg_first_latency / sampling_rate)
    expected_last_lick_sec  = np.full(len(non_participated_idx), avg_last_latency  / sampling_rate)
    
    return (expected_first_lick_indices,
            expected_last_lick_indices,
            expected_first_lick_sec,
            expected_last_lick_sec)

def split_bouts_by_duration(
    bouts_dict,
    ttl_descriptions,
    sampling_rate,
    tolerance=0.02
):
    """
    Split bout onsets into trial types based on TTL pulse duration.

    Parameters
    ----------
    bouts_dict : dict
        {name: (onsets, offsets, activated_indexes)}
    ttl_descriptions : dict
        {'S': (duration_sec, description),
         'M': (duration_sec, description),
         'L': (duration_sec, description)}
    sampling_rate : float
        Samples per second
    tolerance : float
        Allowed deviation in seconds

    Returns
    -------
    split_dict : dict
        {name: {
            'S': [...],
            'M': [...],
            'L': [...],
            'unassigned': [...]
        }}
    """

    # Sort TTLs by expected duration (short → long)
    ttl_items = sorted(
        ttl_descriptions.items(),
        key=lambda x: x[1][0]
    )

    split_dict = {}

    for name, (onsets, offsets, activated) in bouts_dict.items():

        onsets = np.asarray(onsets)
        offsets = np.asarray(offsets)

        durations = (offsets - onsets) / sampling_rate
        result = {k: [] for k, _ in ttl_items}
        result["unassigned"] = []
        result["unassigned_durations"] = []
        

        for onset, dur in zip(onsets, durations):
            assigned = False

            for key, (expected_dur, *_) in ttl_items:
                if abs(dur - expected_dur) <= tolerance:
                    result[key].append(int(onset))
                    assigned = True
                    break

            if not assigned:
                result["unassigned"].append(int(onset))
                result["unassigned_durations"].append(int(dur))
                

        split_dict[name] = result
    print('nr unassigned ',len(result["unassigned"]))
    print(result["unassigned_durations"])
    
    return split_dict


import numpy as np

def match_rewards_to_cues(split_bouts_first, TL_first, TL_last, fs=241, delay_sec=6, tol_sec=10):
    """
    Returns reward onset and offset indices grouped by cue type (S, M, L).
    """

    delay = delay_sec * fs
    tol = tol_sec * fs

    # Flatten cues
    cue_times = []
    cue_labels = []

    for label, indices in split_bouts_first.items():
        if label == 'unassigned':
            continue
        cue_times.extend(indices)
        cue_labels.extend([label] * len(indices))

    cue_times = np.array(cue_times)
    cue_labels = np.array(cue_labels)

    # Sort cues
    order = np.argsort(cue_times)
    cue_times = cue_times[order]
    cue_labels = cue_labels[order]

    # Prepare output
    result = {
        'S': {'onset': [], 'offset': []},
        'M': {'onset': [], 'offset': []},
        'L': {'onset': [], 'offset': []}
    }

    # Loop over rewards (onset + offset together!)
    for onset, offset in zip(TL_first, TL_last):

        if np.isnan(onset) or np.isnan(offset):
            continue

        # Find preceding cue
        idxs = np.where(cue_times < onset)[0]
        if len(idxs) == 0:
            continue

        idx = idxs[-1]
        cue_time = cue_times[idx]
        cue_label = cue_labels[idx]

        diff = (onset - cue_time)

        # Check timing constraint
        if abs(diff - delay) <= tol:
            result[cue_label]['onset'].append(int(onset))
            result[cue_label]['offset'].append(int(offset))

    # Convert to numpy arrays
    for k in result:
        result[k]['onset'] = np.array(result[k]['onset'], dtype=int)
        result[k]['offset'] = np.array(result[k]['offset'], dtype=int)

    return result


################################################################
# functions for baselines 

def compute_baselines(time_lock_dict, baseline_sec, sampling_rate, reference_dict=None):
    """
    Compute baseline indices for multiple time-locks with optional reference events.
    Works for both participated and non-participated trials (NaNs are skipped).
    """
    baseline_dict = {}
    start_samples = int(round(baseline_sec[0] * sampling_rate))
    end_samples   = int(round(baseline_sec[1] * sampling_rate))
    
    for name, indices in time_lock_dict.items():
        indices = np.array(indices)
        baseline_list = []
        
        # Determine reference for this time-lock
        if reference_dict is not None and name in reference_dict and reference_dict[name] is not None:
            ref_array = np.array(reference_dict[name])
        else:
            ref_array = indices  # baseline relative to own event
        
        for i, idx in enumerate(indices):
            if np.isnan(idx):
                baseline_list.append([np.nan, np.nan])  # skip NaNs
                continue
            
            # remove NaNs from reference
            valid_refs = ref_array[~np.isnan(ref_array)]
            
            if len(valid_refs) == 0:
                baseline_list.append([np.nan, np.nan])
                continue
            
            # find closest reference
            closest_idx = np.argmin(np.abs(valid_refs - idx))
            ref_idx = valid_refs[closest_idx]
            
            start_idx = max(int(ref_idx + start_samples), 0)
            end_idx   = max(int(ref_idx + end_samples), 0)
            baseline_list.append([start_idx, end_idx])
        
        baseline_dict[name] = baseline_list
    
    return baseline_dict

######################################################
# time locks to dictionary 

def consolidate_index_data_per_mouse_rows(mouse_ID, time_lock_dict, baseline_dict):
    """
    Consolidate all time-lock indices and baseline info into multiple rows per mouse/session,
    with one row per time-lock.
    
    Args:
        mouse_ID (str): ID string for the mouse/session
        time_lock_dict (dict): Dictionary of original time-lock indices
        baseline_dict (dict): Dictionary of baseline [start,end] per time-lock
    
    Returns:
        pd.DataFrame: Each row = one time-lock, with columns:
                      mouse info, time_lock_name, indices array, baseline_start array, baseline_end array
    """
    # Parse mouse info
    mouse_info, ID = parse_ID(mouse_ID)
    
    rows = []
    
    for tl_name, indices in time_lock_dict.items():
        indices = np.array(indices)
        baseline_pairs = baseline_dict.get(tl_name, [[np.nan, np.nan]] * len(indices))
        baseline_starts = np.array([b[0] for b in baseline_pairs])
        baseline_ends   = np.array([b[1] for b in baseline_pairs])
        
        row = {
            'name': mouse_info['name'],
            'condition': mouse_info['condition'],
            'week': mouse_info['week'],
            'batch': mouse_info['batch'],
            'genotype': mouse_info['genotype'],
            'time_lock_name': tl_name,
            'time_lock_indices': indices,
            'baseline_start': baseline_starts,
            'baseline_end': baseline_ends
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return df

####################################################################
# functions behaviour


def compute_bout_durations(onsets, offsets, sampling_rate):
    """
    Compute the duration of each licking bout in seconds.
    
    Args:
        onsets (np.ndarray): Array of lick bout onsets (indices)
        offsets (np.ndarray): Array of lick bout offsets (indices)
        sampling_rate (float): Sampling rate in Hz (samples per second)
    
    Returns:
        np.ndarray: Array of bout durations in seconds, one per bout
    """
    if len(onsets) == 0:
        return np.array([], dtype=float)
    
    durations_sec = (offsets - onsets + 1) / sampling_rate
    durations_sec = np.round(durations_sec, 2)
    return durations_sec


def get_relative_lick_times(lick_idx, tone_idx, fs, window):
    """
    Compute relative lick times (in seconds) around each tone.

    Parameters
    ----------
    lick_idx : array-like
        Sample indices of detected licks.
    tone_idx : array-like
        Sample indices of tone onsets.
    fs : float
        Sampling rate (samples per second).
    window : tuple (pre, post)
        Time window around tone in seconds, e.g. (-1, 2).

    Returns
    -------
    list_of_lists : list
        For each tone, a list of relative lick times (in seconds).
    """

    lick_idx = np.asarray(lick_idx)
    tone_idx = np.asarray(tone_idx)

    pre, post = window
    pre_samples  = int(pre * fs)
    post_samples = int(post * fs)

    relative_times_per_tone = []

    for t in tone_idx:
        # lick indices inside the window around this tone
        mask = (lick_idx >= t + pre_samples) & (lick_idx <= t + post_samples)
        licks_in_window = lick_idx[mask]

        # convert to time relative to tone
        relative_times = (licks_in_window - t) / fs
        relative_times = np.round(relative_times, 2)
        relative_times_per_tone.append(relative_times.tolist())
    
    return relative_times_per_tone

def count_licks_in_window(rasterdata, window):
    start, end = window
    return [
        sum(start <= t <= end for t in trial)
        for trial in rasterdata
    ]


def extract_licks_per_bin(lick_onsets, lick_times, time_window, sampling_rate, bin_size = 0.05):
    """
    Align licks to event onsets.

    Parameters
    ----------
    lick_onsets : array-like
        Event onset times (seconds).
    lick_times : array-like
        All lick timestamps (seconds).
    time_window : tuple
        (start, end) in seconds relative to each onset.
    sampling_rate : float
        Binning frequency (Hz). Bin width = 1/sampling_rate.

    Returns
    -------
    dict
        {
            "matrix": binary lick matrix (n_trials × n_bins),
            "relative_licks": list of arrays containing exact relative lick times,
            "time": time vector corresponding to the bins
        }
    """

    start, end = time_window
    bin_width = bin_size
    n_bins = int(np.ceil((end - start) / bin_width))

    lick_times = np.asarray(lick_times) / sampling_rate
    lick_onsets = np.asarray(lick_onsets) / sampling_rate

    lick_matrix = np.zeros((len(lick_onsets), n_bins), dtype=np.uint8)
    relative_licks = []

    for i, onset in enumerate(lick_onsets):

        # Find licks in the time window
        left = np.searchsorted(lick_times, onset + start)
        right = np.searchsorted(lick_times, onset + end)

        # Relative lick times
        rel = lick_times[left:right] - onset
        relative_licks.append(rel)

        # Convert to bin indices
        bins = np.floor((rel - start) / bin_width).astype(int)

        # Prevent rare rounding issues
        bins = bins[(bins >= 0) & (bins < n_bins)]

        # One lick per bin
        lick_matrix[i, np.unique(bins)] = 1

    # Bin centers (better for plotting)
    time = start + (np.arange(n_bins) + 0.5) * bin_width

    return {
        "matrix": lick_matrix,
        "relative_licks": relative_licks,
        "time": time,
    }
####################################
# general code - import data & make folder

    
# load files
files = load_files_from_folder(files_dir_preprocessed, '.parquet')

#determine sampling rate
sampling_rate = determine_sampling_rate_parquet(files, files_dir_preprocessed)

        
#make analysed folder if this does not exists yet with current date
os.chdir (files_dir_raw)
# date = time.strftime('%x'); new_date = date[6:8]; new_date = new_date + date[0:2] ; date = new_date + date[3:5]
# files_dir_behav = files_dir_raw +  'Behaviour_extraction ' + date  + '//'
# if os.path.exists(files_dir_behav) == 0:
#     os.makedirs(files_dir_behav)
    

files_dir_behav = create_folder_with_note(files_dir_raw, '2. Behaviour_extraction')
    
all_behavior_dfs = []

# load data file for each file
for i_file in range(0,len(files)):
    print('processing ' + files[i_file] )
    ID, ID_name = parse_ID(files[i_file])
    

    # - General code for loading
    df = pd.read_parquet(files_dir_preprocessed + files[i_file] )
    
    # extract DIOS
    DIO_arrays = {}
    if "TTL_1" in df.columns:
        DIO_arrays["first"] = np.array(df["TTL_1"]).round()
    
    if "TTL_2" in df.columns:
        DIO_arrays["sec"] = np.array(df["TTL_2"]).round()
    
    
    # determine onsets and offsets
    DIO_results = process_DIO(DIO_arrays)
    first_onsets, first_offsets, first_activated = DIO_results['first']
    sec_onsets, sec_offsets, sec_activated = DIO_results['sec']
    

#########################################################
    # - Exp specific: determine time locks
    
    # exclude time locks that are boyond a cetain trial
    filtered_DIO = filter_by_trial_range(DIO_results, first_onsets, nr_trials_analysed)
    #get lick list without bout durations filter only trial number filter
    sec_onsets, sec_offsets, sec_activated = filtered_DIO['sec']


    #concaternate bouts
    concatenated_DIO = concatenate_DIO_bouts(
        DIO_dict=filtered_DIO,
        sampling_rate=sampling_rate,
        min_inter_bout_interval_sec=min_inter_bout_interval_sec )
    
    sec_onsets, sec_offsets, sec_activated = concatenated_DIO['sec']
    
    #filtered bouts  that have less then a certain duration
    filtered_DIO_durations = filter_bouts_by_duration_all(
    DIO_dict=concatenated_DIO,
    sampling_rate=sampling_rate,
    min_durations=min_durations)

    #Extract data from dictionary
    first_onsets_final, first_offsets_final, _ = filtered_DIO_durations["first"]
    sec_onsets_final, sec_offsets_final, _ = filtered_DIO_durations["sec"]
        
    
    ###################################
    #time locks

    # time lock all tones
    TL_all_tones = first_onsets_final
       
    # time lock first lick and last lick reward collection 
    expected_offsets = first_onsets_final + sampling_rate * 5
    TL_first_licks_rew, TL_last_licks_rew = first_and_last_lick_after_tone(
        tone_offsets=expected_offsets,
        lick_onsets=sec_onsets_final,
        lick_offsets=sec_offsets_final,
        sampling_rate=sampling_rate,
        min_delay_sec=min_lick_delay_sec,
        max_delay_sec=max_lick_delay_sec )
        
    #split bouts for only ttl 1
    bouts_dict_ttl1 = {"first": filtered_DIO_durations["first"]}
    split_bouts = split_bouts_by_duration(
        bouts_dict=bouts_dict_ttl1,
        ttl_descriptions=TTL_1_desc,
        sampling_rate=sampling_rate,
        tolerance=1)
    
    print('nr S bouts ', len(split_bouts['first']['S']) )
    print('nr M bouts ', len(split_bouts['first']['M']) )
    print('nr L bouts ', len(split_bouts['first']['L']) )
    #match bouts with licking bouts. 
    matched_rewards = match_rewards_to_cues(
        split_bouts['first'],
        TL_first_licks_rew,
        TL_last_licks_rew, 
        tol_sec=10)
        
    # time lock participated tones & non participated tones
    TL_participated_tones, TL_non_participated_tones, participated_idx, non_participated_idx = \
    split_participated_trials(
        tone_onsets=first_onsets_final,   # or your filtered/concatenated tones
        first_licks=TL_first_licks_rew )       # from first_and_last_lick_after_tone()
    
    
     # time lock first lick during cue
    TL_first_lick_indices_cue, beh_lick_latency_cue = first_lick_during_cue(
        tone_onsets=concatenated_DIO["first"][0],
        lick_onsets=concatenated_DIO["sec"][0],
        sampling_rate=sampling_rate,
        cue_duration_sec=5 )   
     
      # time lock omission lick   
    TL_expected_first, TL_expected_last, expected_first_sec, expected_last_sec = \
        time_lock_omission_lick_offset(
            first_licks=TL_first_licks_rew,
            last_licks=TL_last_licks_rew,
            tone_offsets=first_offsets_final,
            non_participated_idx=non_participated_idx,
            sampling_rate=sampling_rate   )
    
    
       # time lock lick no reward
    lick_no_reward_dict = get_licks_no_reward(
        DIO_dict=filtered_DIO_durations,    # already concatenated and filtered
        tone_onsets= first_onsets,  # tone onsets
        sampling_rate=sampling_rate,
        time_before_tone_sec=5,
        time_after_tone_sec=25   )
    # Extract data from dictionary
    TL_lick_no_reward, beh_lick_no_reward_sec = lick_no_reward_dict["sec"]

    ################

    #baseline calculation of the baseline_sec before each time lock

    time_lock_dict = {
        'reward onset ' + TTL_1_desc['S'][1]: matched_rewards['S']['onset'],
        'reward offset ' + TTL_1_desc['S'][1]: matched_rewards['S']['offset'],
        'reward onset ' + TTL_1_desc['M'][1]: matched_rewards['M']['onset'],
        'reward offset ' + TTL_1_desc['M'][1]: matched_rewards['M']['offset'],
        'reward onset ' + TTL_1_desc['L'][1]: matched_rewards['L']['onset'],
        'reward offset ' + TTL_1_desc['L'][1]: matched_rewards['L']['offset'],
        'first_lick_indices_cue': TL_first_lick_indices_cue,
        'rew_omission_exp_onset': TL_expected_first,
        'rew_omission_exp_offset': TL_expected_last,
        'participated_cues': TL_participated_tones,
        'all_cues': TL_all_tones,
        'lick_no_reward': TL_lick_no_reward,
        'non_participated_cues': TL_non_participated_tones}
    
    # Assign reference per time-lock (None means baseline relative to its own event)
    reference_dict = {
        'reward onset ' + TTL_1_desc['S'][1]: first_onsets,
        'reward offset ' + TTL_1_desc['S'][1]: first_onsets,
        'reward onset ' + TTL_1_desc['M'][1]: first_onsets,
        'reward offset ' + TTL_1_desc['M'][1]: first_onsets,
        'reward onset ' + TTL_1_desc['L'][1]: first_onsets,
        'reward offset ' + TTL_1_desc['L'][1]: first_onsets,
        'first_lick_indices_cue': first_onsets,
        'rew_omission_exp_onset': first_onsets,
        'rew_omission_exp_offset': first_onsets,
        # The "base_" time-locks are relative to their own event
        'participated_cues': None,
        'all_cues': None,
        'lick_no_reward': None,
        'non_participated_cues': None}
    
    baseline_dict = compute_baselines(time_lock_dict, baseline_sec, sampling_rate, reference_dict)   
    # Example access     print("Baseline for first rewarded lick, first trial:", baseline_dict['first_licks_rew'][0])

     ###############
     # code to transform it into a df
    df_timelocks = consolidate_index_data_per_mouse_rows(files[i_file], time_lock_dict, baseline_dict)
    

    
#################################################################################    
    # - Exp specific: determine behavioural outcomes
    #     Raw or processed such as minimal bout duration
    

    # rewards collected
    beh_nr_rew_collected = len(TL_participated_tones)
    beh_missed_trials = len(TL_non_participated_tones)

    # latencies to lick during cue & reward
    # beh_lick_latency_cue =  lick latency cue, determined in first_lick_during_cue
    beh_lick_latency_cue_pt = beh_lick_latency_cue
    beh_lick_latency_cue_avr = round(np.nanmean(beh_lick_latency_cue),2)
    beh_lick_latency_rew_avr = round(np.mean(expected_first_sec),2)
    beh_lick_rew_offset_avr = round(np.mean(expected_last_sec),2)
    
    # compute bout durations of general licking bout
    beh_bout_duration_licks_pt = compute_bout_durations( onsets=sec_onsets_final,  offsets=sec_offsets_final,  sampling_rate=sampling_rate )
    beh_bout_duration_licks_avr = round(np.nanmean(beh_bout_duration_licks_pt),2)
    # compute bout durations of reward licking bout
    beh_bout_duration_rew_pt = compute_bout_durations( onsets=TL_first_licks_rew,  offsets=TL_last_licks_rew,  sampling_rate=sampling_rate )
    beh_bout_duration_rew_avr = round(np.nanmean(beh_bout_duration_rew_pt),2)
    
    # licking without reward
    beh_nr_lick_bouts_no_reward = len(TL_lick_no_reward)
    beh_lick_no_reward_sec_avr = round(np.mean(beh_lick_no_reward_sec),2)
    

    # calculates time of lick relative to the tone which is saved in the rasterdata variable
    window = (-20, 30) 
    rasterdata_licks_all = get_relative_lick_times(sec_onsets, TL_all_tones, sampling_rate, window)
    rasterdata_licks_part = get_relative_lick_times(sec_onsets, TL_participated_tones, sampling_rate, window)

    # total licks and licks during specific moments cue, reward, outside
    beh_total_licks = len(sec_onsets)
    beh_licks_cue = sum(count_licks_in_window(rasterdata_licks_all, (0, 5)))
    beh_licks_rew = sum(count_licks_in_window(rasterdata_licks_all, (5, 15)))
    beh_licks_outside = sum(count_licks_in_window(rasterdata_licks_all, (-20, 0))) + sum(count_licks_in_window(rasterdata_licks_all, (15, 30)))
    
    # calculate licks per second   combined for all trials
    beh_licks_cue_per_s = round(beh_licks_cue / 5,2)
    beh_licks_rew_per_s = round(beh_licks_rew / 10,2)
    beh_licks_outside_per_s = round(beh_licks_outside / 35,2)
    
    # 3 licks graphs
    onsets_per_conc = [split_bouts['first']['S'], split_bouts['first']['M'], split_bouts['first']['L']]
    onsets_per_conc = [matched_rewards['S']['onset'],  matched_rewards['M']['onset'], matched_rewards['L']['onset']]
    concentration_labels = ['S', 'M', 'L']  # or use more descriptive names


    licks_per_concentration = {}
    for label, onsets in zip(concentration_labels, onsets_per_conc):
        licks_per_concentration[label] = extract_licks_per_bin(
            lick_onsets=onsets,
            lick_times=sec_onsets,
            time_window=time_window,
            sampling_rate=sampling_rate,   # only used for conversion below
            bin_size=bin_size,                 # 50 ms bins
        )


    # save all data in a dictionary 
    behavior_dict = {
        'nr_rew_collected': beh_nr_rew_collected,
        'missed_trials': beh_missed_trials,
        'nr reward ' + TTL_1_desc['S'][1]: len(matched_rewards['S']['onset']),
        'nr reward  ' + TTL_1_desc['M'][1]: len(matched_rewards['M']['onset']),
        'nr reward ' + TTL_1_desc['L'][1]: len(matched_rewards['L']['onset']),
        'lick_latency_cue_avr': beh_lick_latency_cue_avr,
        'lick_latency_rew_avr': beh_lick_latency_rew_avr,
        'lick_rew_offset_avr': beh_lick_rew_offset_avr,
        'bout_duration_licks_avr': beh_bout_duration_licks_avr,
        'bout_duration_rew_avr': beh_bout_duration_rew_avr,
        'nr_lick_bouts_no_reward': beh_nr_lick_bouts_no_reward,
        'lick_no_reward_sec_avr': beh_lick_no_reward_sec_avr,
       
        # lick counts
        'total_licks': beh_total_licks,
        'licks_cue': beh_licks_cue,
        'licks_cue_per_s': beh_licks_cue_per_s,
        'licks_rew': beh_licks_rew,
        'licks_rew_per_s': beh_licks_rew_per_s,
        'licks_outside': beh_licks_outside,
        'licks_outside_per_s': beh_licks_outside_per_s,
        
        # data per trial
        'bout_duration_licks_pt': beh_bout_duration_licks_pt,
        'lick_latency_cue_pt': beh_lick_latency_cue_pt,
        'rasterdata_licks_all': rasterdata_licks_all,
        'rasterdata_licks_part': rasterdata_licks_part,
        'bout_duration_rew_pt': beh_bout_duration_rew_pt,
        
        # licks in time
        'lick_rate_S': licks_per_concentration["S"]["matrix"].mean(axis=0) /bin_size,
        'lick_rate_M': licks_per_concentration["M"]["matrix"].mean(axis=0) /bin_size,
        'lick_rate_L': licks_per_concentration["L"]["matrix"].mean(axis=0) /bin_size,
        'lick_rate_time': licks_per_concentration["S"]["time"]
    }
    

    # Merge dictionaries with ID fields first
    behavior_dict = {**ID, **behavior_dict}  #** means unpacking of dictionary
    
    # store in dataframe
    behavior_df = pd.DataFrame({k: [v] for k, v in behavior_dict.items()})
    
    # Convert to DataFrame (1 row)
    all_behavior_dfs.append(behavior_df)

    
    ######################################################################################
    # - general code for saving
       
    # Save DataFrame to Parquet and csv
    df_timelocks.to_parquet(os.path.join(files_dir_behav, f"{ID_name}.parquet"))
    df_timelocks.to_csv(os.path.join(files_dir_behav, f"{ID_name}.csv"),  index=False)
    
    print(f"✅ DataFrame saved as: {files[i_file]}")



# Combine all into a single DataFrame
combined_behavior_df = pd.concat(all_behavior_dfs, ignore_index=True)

# Save to Excel, CSV, and Parquet
combined_behavior_df.to_excel(os.path.join(files_dir_behav,"combined_behavior.xlsx"), index=False)
combined_behavior_df.to_csv(os.path.join(files_dir_behav,"combined_behavior.csv"), index=False)
combined_behavior_df.to_parquet(os.path.join(files_dir_behav,"combined_behavior.parquet"), index=False)



