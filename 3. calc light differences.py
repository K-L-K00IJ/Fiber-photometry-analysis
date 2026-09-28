# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 17:41:50 2025
3. calculation of light difference of fiber photometry data
general script

This script will perform the following steps:
- calcuate Delta F/F in different senarios 
- calculate z_scores

Inputs: 
    preprocessed data in parket file
    DF with colums for each time lock - containing name and indexes of time lock and baselines

Output:
    Per mouse per timelock 1 DF of per trial results and 1 df of all average scores
    Both deltaF/F and z_scores are calcuated
    

Step 3 out of 4 in the photometry analysis pipeline
Meye lab, Karlijn Kooij
Version 5_12_25

#update 13-3-26
- added the same baseline for z-score calculations

#update 20-5-26
add forced baseline for z-score


@author: kkooij
"""

# enter the following parameters per analysis

# enter the directory where the raw data files are stored
files_dir_general =    'L://adanlab//Ongoing//apetitbon//Projects//Lep_DA_CeA//Pilot WP 01 - 19//DORIC-Pilot_GRAB-DA3h_CeA_LepRcre_VTA_opto-ChRimson//20260423//'
#files_dir_general ='L://meyelab//Ongoing//o_Project_PFC_Binge Ensembles//02_ResearchData//Fiber Photometry Analysis//Food Restricted//Second Recordings//Analysis Karlijn//'

folder_Preprocessed = files_dir_general + '1. Preprocessed//'
folder_behaviour = files_dir_general + '2. Behaviour_extraction_selected_TS//'

# window around time lock
pre_seconds = 5
post_seconds = 10

sampling_rate = 241

# when z-score is based on whole trace, do you want to force the baseline of z to be 0?
ex_force_baseline_z = 0




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



##############################################################

def load_files_from_folder(files_dir, extention):
  # inport data from 1 folder with a certain extention 
  list_files = os.listdir(files_dir)
  files = []
  for i_files in range(len(list_files)):
      if list_files[i_files].endswith(extention):
          files.append(list_files[i_files])
          
  return files

def make_folder_date(files_dir, folder_name):
    #make analysed folder if this does not exists yet with current date
    os.chdir (files_dir)
    date = time.strftime('%x'); new_date = date[6:8]; new_date = new_date + date[0:2] ; date = new_date + date[3:5]
    files_dir_behav = files_dir +  folder_name + date  + '//'
    if os.path.exists(files_dir_behav) == 0:
        os.makedirs(files_dir_behav)
        
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






####################################
# delta F/F calc

def compute_deltaf_with_windows(data_preprocessed, data_behaviour,ID, sampling_rate, pre_seconds, post_seconds):
    """
    Compute ΔF/F around time-locks with pre/post windows
    """
    
    photo = data_preprocessed["photosignal"].values
    iso   = data_preprocessed["isobestic"].values
    n = len(photo)


    # Number of samples before and after
    pre_samples  = int(pre_seconds * sampling_rate)
    post_samples = int(post_seconds * sampling_rate)
    
    # Time vector relative to lock
    window_len = pre_samples + post_samples + 1
    time_vector = np.linspace(-pre_seconds, post_seconds, window_len)
    
    all_timelocks = []
    all_timelocks_dfs = []

    for _, row in data_behaviour.iterrows():
        #Clean the single center index
        lock_index = row.get("time_lock_indices", np.nan)
        b0 = row.get("baseline_start", np.nan)
        b1 = row.get("baseline_end",   np.nan)

        # remove nans
        lock_index = lock_index[~np.isnan(lock_index)]
        b0 = b0[~np.isnan(b0)]
        b1 = b1[~np.isnan(b1)]
        
        
        # Skip if no valid lock or baseline
        if len(lock_index) == 0 or len(b0) == 0 or len(b1) == 0:
        
            nan_window = np.full(window_len, np.nan)
        
            time_lock_res = {
                "time_lock_name": row.get("time_lock_name", None),
        
                # keep same structure
                "dff_photo": [nan_window.copy()],
                "dff_iso": [nan_window.copy()],
                "time_vec": nan_window.copy(),
            }
        
            time_lock_res = {**ID, **time_lock_res}
        
            all_timelocks.append(time_lock_res)
            continue
        
        dff_photo_list = []
        dff_iso_list = []
        
        n_trials = min(len(lock_index), len(b0), len(b1))
        
        # loop though each trial
        for i in range(n_trials):
            
            start = int(lock_index[i] - pre_samples)
            end   = int(lock_index[i] + post_samples + 1)   
            base_start = int(b0[i])
            base_end   = int(b1[i])
            
         # Skip trials too close to beginning/end for either window or baseline
            if ( start < 0 or end > n or  base_start < 0 or base_end > n or base_end <= base_start  ):
                               
                continue

            
            # Extract windows from signals
            photo_seg = photo[start:end]
            iso_seg   = iso[start:end]

            # 3) Baseline handling
            photo_base = photo[int(b0[i]):int(b1[i])]
            iso_base   = iso[int(b0[i]):int(b1[i])]

            # 4) ΔF/F computation (vectorized)
            ref_photo = photo_base.mean()
            ref_iso   = iso_base.mean()

            dff_photo = ((photo_seg - ref_photo) / ref_photo) * 100
            dff_iso   = ((iso_seg - ref_iso) / ref_iso) * 100
            
            dff_photo_list.append(dff_photo)
            dff_iso_list.append(dff_iso)
            
            
            expected_len = pre_samples + post_samples + 1

            if len(dff_photo) != expected_len:
                print(
                    f"BAD WINDOW | lock={lock_index[i]} "
                    f"start={start}, end={end}, "
                    f"len={len(photo_seg)}"
                )
        
        # ---> Convert lists of arrays to 2D matrices
        # dff_photo_list = np.vstack(dff_photo_list)
        # dff_iso_list   = np.vstack(dff_iso_list)
           
        # Store time lock results
        time_lock_res = {
            "time_lock_name": row.get("time_lock_name", None),
            "dff_photo": dff_photo_list,
            "dff_iso": dff_iso_list,
            "time_vec": time_vector }
    
        # Merge dictionaries with ID fields first
        time_lock_res = {**ID, **time_lock_res}  #** means unpacking of dictionary
        
        all_timelocks.append(time_lock_res)

    # # Convert list of dicts to DataFrame
    all_timelocks_df = pd.DataFrame(all_timelocks)

    return all_timelocks_df




def compute_deltaf_exp_fit(df, data_preprocessed, data_behaviour,ID, sampling_rate, pre_seconds, post_seconds):
    """
    Compute ΔF/F around time-locks with pre/post windows and add to existing df
    """
    
    MC = data_preprocessed["MC"].values
    gfp_expfit = data_preprocessed["gfp_expfit"].values
    
    MC_dF_F = 100*(MC/gfp_expfit)
    
    new_df = extract_trials_from_signal( signal=MC_dF_F, data_behaviour=data_behaviour, ID=ID, sampling_rate=sampling_rate, pre_seconds=pre_seconds, post_seconds=post_seconds, value_name="dff_MC")
    
    if len(new_df) > 0:
        df = df.merge( new_df, on=list(ID.keys()) + ["time_lock_name"],   how="left" )
            
    return df


def extract_trials_from_signal(signal,data_behaviour, ID, sampling_rate, pre_seconds, post_seconds, value_name):
    """
    Extract trial-aligned windows from a 1D signal.

    Returns a DataFrame with one row per time_lock_name and a column
    named `value_name` containing a list of trial arrays.
    """

    pre_samples  = int(pre_seconds * sampling_rate)
    post_samples = int(post_seconds * sampling_rate)

    all_results = []

    for _, row in data_behaviour.iterrows():

        lock_index = row.get("time_lock_indices", np.nan)
        lock_index = lock_index[~np.isnan(lock_index)]

        if len(lock_index) == 0:
            continue

        trial_list = []

        for i in range(len(lock_index)):

            start = int(lock_index[i] - pre_samples)
            end   = int(lock_index[i] + post_samples + 1)

            if start < 0 or end > len(signal):
                continue

            expected_len = pre_samples + post_samples + 1
            if len(signal[start:end]) != expected_len:
                print(
                    f"BAD WINDOW | lock={lock_index[i]} "
                    f"start={start}, end={end}, "
                    f"len={len(signal[start:end])}"
                )
        
        
            trial_list.append(signal[start:end])

        res = {"time_lock_name": row.get("time_lock_name", None), value_name: trial_list }

        res = {**ID, **res}
        all_results.append(res)

    return pd.DataFrame(all_results)


# old for lists
def compute_z_scores(res_pt):
    # Create empty columns first
    res_pt['z_score_photo'] = [[] for _ in range(len(res_pt))]
    res_pt['z_score_iso']   = [[] for _ in range(len(res_pt))]

    channels = ['_photo', '_iso']

    for i_timelock in range(len(res_pt)):
        for ch in channels:
            all_trials = res_pt.loc[i_timelock, 'dff' + ch]
            z_all_trials = []
            
            for trial_array in all_trials:
                array = np.array(trial_array)
                mean_val = np.nanmean(array)
                std_val = np.nanstd(array)
                
                if std_val == 0:
                    z_trial = np.zeros_like(array)
                else:
                    z_trial = (array - mean_val) / std_val
                
                z_all_trials.append(z_trial)
            
            res_pt.at[i_timelock, 'z_score' + ch] = z_all_trials
    return res_pt

# def compute_z_scores(res_pt):
#     # Prepare output columns
#     res_pt['z_score_photo'] = [None] * len(res_pt)
#     res_pt['z_score_iso']   = [None] * len(res_pt)

#     for i in range(len(res_pt)):
#         for suffix in ["_photo", "_iso"]:
            
#             trials_np = res_pt.at[i, 'dff' + suffix]  # Already a NumPy matrix (n_trials, n_timepoints)

#             # Compute mean and std per trial (axis=1)
#             means = np.nanmean(trials_np, axis=1, keepdims=True)
#             stds  = np.nanstd(trials_np, axis=1, keepdims=True)

#             # Avoid division by zero
#             stds_safe = np.where(stds == 0, 1, stds)

#             # Vectorized z-scoring
#             z_np = (trials_np - means) / stds_safe

#             res_pt.at[i, 'z_score' + suffix] = z_np

#     return res_pt
            


def compute_z_scores_whole_trace(df, data_preprocessed, data_behaviour,ID, sampling_rate, pre_seconds, post_seconds):
    
    MC = data_preprocessed["MC"].values
    
    MC_zscored = (MC-np.mean(MC))/np.std(MC)
    
    new_df = extract_trials_from_signal( signal=MC_zscored, data_behaviour=data_behaviour, ID=ID, sampling_rate=sampling_rate, pre_seconds=pre_seconds, post_seconds=post_seconds, value_name="z_score_MC")
    
    if len(new_df) > 0:
        df = df.merge( new_df, on=list(ID.keys()) + ["time_lock_name"],   how="left" )
    else:
           # create placeholder rows
        empty_rows = []
        pre_samples  = int(pre_seconds * sampling_rate)
        post_samples = int(post_seconds * sampling_rate)
        window_len = pre_samples + post_samples + 1
        
        for _, row in data_behaviour.iterrows():
            empty_row = {
                **ID,
                "time_lock_name": row.get("time_lock_name", None),
            
                "dff_MC": [np.full(window_len, np.nan)],
                "z_score_MC": [np.full(window_len, np.nan)],
            }
        
            empty_rows.append(empty_row)
        
        empty_df = pd.DataFrame(empty_rows)
        
        # merge placeholder
        df = df.merge(
            empty_df,
            on=list(ID.keys()) + ["time_lock_name"],
            how="left"
        )
    return df, MC_zscored

def compute_z_score_whole_trace_baseline(df, data_preprocessed, data_behaviour, ID, sampling_rate, pre_seconds, post_seconds):
    """
    Compute z-scores using whole-trace std, then baseline-correct each trial.
    """
    # 1) Whole trace for MC channel
    MC = data_preprocessed["MC"].values
    mc_mean = np.mean(MC)
    mc_std  = np.std(MC)
    
    if mc_std == 0:
        MC_z = np.zeros_like(MC)
    else:
        MC_z = (MC - mc_mean) / mc_std
    
    # 2) Extract trials from z-scored signal
    df_z = extract_trials_from_signal(
        signal=MC_z,
        data_behaviour=data_behaviour,
        ID=ID,
        sampling_rate=sampling_rate,
        pre_seconds=pre_seconds,
        post_seconds=post_seconds,
        value_name="z_score_MC_whole"
    )
    
    # 3) Merge with existing df
    df = df.merge(df_z, on=list(ID.keys()) + ["time_lock_name"], how="left")
    
    # 4) Baseline correction per trial
    df['z_score_MC_baseline'] = [[] for _ in range(len(df))]

    for i in range(len(df)):
        trials = df.loc[i, "z_score_MC_whole"]
        time_lock_name = df.loc[i, "time_lock_name"]
        
        # Extract baseline arrays for this time-lock
        baseline_row = data_behaviour.loc[data_behaviour['time_lock_name'] == time_lock_name]
        
        # Assuming one row per time-lock
        b0 = np.array(baseline_row['baseline_start'].values[0])
        b1 = np.array(baseline_row['baseline_end'].values[0])
        
        z_trials = []
        n_trials = min(len(trials), len(b0), len(b1))
        
        for t in range(n_trials):
            trial_array = np.array(trials[t])
            baseline_start = int(b0[t])
            baseline_end   = int(b1[t])
            
            if baseline_end <= baseline_start:
                z_trials.append(np.zeros_like(trial_array))
                continue
            
            # Compute baseline mean from original MC trace
            baseline_mean = np.nanmean(MC[baseline_start:baseline_end])
            
            # Baseline-correct using whole-trace std
            z_corrected = (trial_array - baseline_mean) / mc_std
            z_trials.append(z_corrected)
        
        df.at[i, 'z_score_MC'] = z_trials
    
    return df

def make_df_avr(df, col_avr):
    
    # loop trough columns to average
    for i_col in range(len(col_avr)):
        for i_tl in range(len(df[col_avr[i_col]] )):
            # make average of np.matrix
            df.at[i_tl, col_avr[i_col]] = np.nanmean(df.loc[i_tl,col_avr[i_col]], axis=0)
    
    return df

def force_baseline_z(res_pt,MC_z_trace, data_behaviour):

    z_columns = ['z_score_MC']

    for col in z_columns:

        for i_row in range(len(data_behaviour)):

            b0 = data_behaviour.iloc[i_row]["baseline_start"]
            b1 = data_behaviour.iloc[i_row]["baseline_end"]

            # skip invalid rows
            if b0 is None or b1 is None:
                continue

            cur_trials = res_pt.at[i_row, col]

            # skip missing data
            if cur_trials is None:
                continue

            # prevent index mismatch
            n_trials = min(len(cur_trials), len(b0), len(b1))

            for i_trial in range(n_trials):
                cur_trial = cur_trials[i_trial]

                start = b0[i_trial]
                end = b1[i_trial]

                # skip nan indices
                if np.isnan(start) or np.isnan(end):
                    print('baseline index contains nans')
                    continue

                start = int(start)
                end = int(end)

                baseline_vals = MC_z_trace[start:end]

                # skip empty/all-nan baseline
                if len(baseline_vals) == 0:
                    print('len vals')
                    continue

                cur_baseline = np.nanmean(baseline_vals)

                # subtract baseline
                cur_trials[i_trial] = cur_trial - cur_baseline

            # assign back
            res_pt.at[i_row, col] = cur_trials

    return res_pt

####################################
# Clear all variables (Spyder / IPython)
# try:
#     from IPython import get_ipython
#     get_ipython().magic("reset -f")
# except Exception:
#     pass


# load files
files_preprocessed = load_files_from_folder(folder_Preprocessed, '.parquet')
files_behaviour = load_files_from_folder(folder_behaviour , '.parquet')

# convert behaviour list to a set for faster lookup
behaviour_set = set(files_behaviour)
    
files_dir_light = create_folder_with_note(files_dir_general, '3. Light difference')

res_pt_list = []; res_avr_list = []

# load data file for each file
for pre_file in files_preprocessed:
    print('processing ' + pre_file )
    ID, ID_name = parse_ID(pre_file)
    

    # - General code for loading
    df_preprocessed = pd.read_parquet(folder_Preprocessed + pre_file )
    df_behaviour =    pd.read_parquet(folder_behaviour + pre_file)
    

    #delta f/f calc 
    #old method with baseline based on window
    res_pt = compute_deltaf_with_windows(df_preprocessed,df_behaviour,ID, sampling_rate, pre_seconds, post_seconds)
    #new method based on exp curve
    res_pt = compute_deltaf_exp_fit(res_pt,df_preprocessed,df_behaviour,ID, sampling_rate, pre_seconds, post_seconds )
    
    
    # z-score calc based on existing delta/f
    res_pt = compute_z_scores(res_pt)
    # based on whole trace
    res_pt, MC_z_trace = compute_z_scores_whole_trace(res_pt,df_preprocessed,df_behaviour,ID, sampling_rate, pre_seconds, post_seconds)
    #res_pt = compute_z_score_whole_trace_baseline(res_pt,df_preprocessed,df_behaviour,ID, sampling_rate, pre_seconds, post_seconds)
    
    #force baseline on z-score
    if ex_force_baseline_z == 1:
        res_pt = force_baseline_z(res_pt,MC_z_trace, df_behaviour)
    
    #make all averaged
    col_avr = ["dff_photo", "dff_iso", 'dff_MC', "z_score_photo", "z_score_iso", 'z_score_MC']
    res_avr = make_df_avr(res_pt.copy(deep=True), col_avr)



    # Save to Excel, CSV, and Parquet
    files_dir_light_pt = create_folder_with_note(files_dir_light, 'per trial')
    res_pt.to_parquet(os.path.join(files_dir_light_pt,pre_file), index=False)
    files_dir_light_av = create_folder_with_note(files_dir_light, 'average')
    res_avr.to_parquet(os.path.join(files_dir_light_av,pre_file), index=False)
    
    # unpack for excel and csv
    files_dir_light_csv = create_folder_with_note(files_dir_light, 'csv')
    res_avr.to_excel(os.path.join(files_dir_light_csv,ID_name+ '.xlsx'), index=False)



# see light fiber action article for adjustments 
# automatic ssample rate detection
# 2 gaven een fout melding. waarom? include