# -*- coding: utf-8 -*-
"""
1. Preprocessing of fiber photometry data

This script will perform the following preprocessings steps:
- load data files from folder
- extract columns 
- set max time
- unsample data if needed
- bleacing correction
- low pass filter
- smoothing

Inputs: raw data in doric or csv files
Output: parket file

Step 1 out of 4 in the photometry analysis pipeline
Meye lab, Karlijn Kooij
Version 4_11_25


"""
# enter the following parameters per analysis

# enter the directory in which the csv files are located
files_dir1 = 'D://TRAP//Merged analysis//FR Data//'


#data location within file
# note as list of list. 
detector = ['1', '2']  #note here which detector was coupled to the [first mouse and second mouse] analysed. Note as ['1', '2']
DIO = [  ['1', '2'], ['3', '4']] 


name_empty_files = 'empty'  #must have same length as ind of names
doric_version = 6

ind_names = [[0,5],[19,24]]
ind_hemiphere = [[15, 16], [15, 16]] 
ind_week = [38, 44]
ind_con = [[6, 8], [25, 27]] 
ind_gen = [[15, 16], [15, 16]] 
ind_batch = 16


miceperfile = 2

# Which steps to include in the processing?
Do_max_time = 1
Do_unsample_data = 1
Do_analog_to_dio = 0   # very rare cases- for older data
Do_artefact_removal =1
Do_bleaching_correction =1
Do_low_pass_filter = 1
Do_motion_correction = 1
Do_smoothing = 1

# add specs on max time Note in minutes
max_time = 200

# analog to DIO
analog_i = []   #only coded for csv files
AI_to_DIO = ['1']

# add specifics for each step

# On unsampling data    # 1 for dec 50 and 4.001 for dec 200    1.017 for diff 237 and 241
sample_factor = 1

# Artefact removal will be done manually be selecting data that needs to be removed

# bleaching correction
#double exp curve is fitted. No steps to add

#low pass filter

#Smoothing 
smoothing_time = 0.4  # add in sec, thus 0.4 is 400 ms




#######################################
# import packages 
import pandas as pd
import numpy as np 
import os
import time
import statistics
from scipy.signal import butter,  filtfilt
from scipy.stats import linregress
import h5py # Make sure to install the library
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Qt5Agg')   # or 'TkAgg' if Qt is not available
import matplotlib.pyplot as plt
plt.ioff() 
from matplotlib.widgets import SpanSelector, Button
from datetime import datetime



###############################################################################
# functions  

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


# def load_data (dataframe_1,files_dir1 , detector, DIO, i_miceperfile, doric_version ):
#     DIO_num_first = DIO[i_miceperfile][0]
#     DIO_num_sec = DIO[i_miceperfile][1] 
#     Det_num = detector[i_miceperfile]
    
#     if dataframe_1[-6::] == '.doric':
#         with h5py.File(files_dir1+dataframe_1, 'r') as f:   
#             path = '/DataAcquisition/FPConsole/Signals/Series0001/'
#             if doric_version == 5:
#                   time_     = np.array(f[path+'AIN01xAOUT01-LockIn/Time'])
#                   DIO_first      = np.array(f[path+'/DigitalIO/DIO0'+DIO_num_first])
#                   DIO_sec      = np.array(f[path+'/DigitalIO/DIO0'+DIO_num_sec])
#                   ain1aout1 = np.array(f[path+'/AIN0' + Det_num + 'xAOUT01-LockIn/Values'])
#                   ain1aout2 = np.array(f[path+'/AIN0' + Det_num +'xAOUT02-LockIn/Values'])
#             elif doric_version == 6:
#                time_     = np.array(f[path+'LockInAOUT01/Time'])
#                DIO_first      = np.array(f[path+'/DigitalIO/DIO0'+DIO_num_first ])   
#                DIO_sec      = np.array(f[path+'/DigitalIO/DIO0'+DIO_num_sec ])    
#                ain1aout1 = np.array(f[path+ 'LockInAOUT01/AIN0' + Det_num])   #isobestic det 1
#                ain1aout2 = np.array(f[path+'LockInAOUT02/AIN0' + Det_num])     #GFP
                
               
#          #remove 20 samples in the beginning and at the end of the dios
#             if len(DIO_first) != len(time_):
#                DIO_first = DIO_first[20:-20]
#                DIO_sec = DIO_sec[20:-20]
            
#          #extra check if lengths are the same. Some one difference between DIOs and others
#             if len(DIO_first) > len(time_):
#               difference = len(DIO_first) - len(time_)
#               DIO_first = DIO_first[:-difference] 
#               DIO_sec = DIO_sec[:-difference] 
        
#             if len(DIO_first) < len(time_):
#              difference = len(time_) - len(DIO_first)
#              time_ = time_[:-difference]
#              ain1aout1 = ain1aout1[:-difference]
#              ain1aout2 = ain1aout2[:-difference]  
        
#          # Create data frame
#             df = pd.DataFrame({
#                 'Time(s)' : time_,
#                 'TTL_first'    : DIO_first,
#                 'TTL_sec'    : DIO_sec,
#                 'isobestic': ain1aout1,
#                 'photosignal': ain1aout2
#                 })   
      
#     elif dataframe_1.endswith('.csv'):
#         # Construct dynamic column names
#         col_time = 'Time(s)'
#         col_dio_first = f'DI/O-{DIO_num_first}'
#         col_dio_sec = f'DI/O-{DIO_num_sec}'
#         col_iso = f"AIn-{Det_num} - Dem (AOut-1)"
#         col_photo = f"AIn-{Det_num} - Dem (AOut-2)"
        
#         # Read CSV with dynamic columns
#         df = pd.read_csv( dataframe_1, encoding='unicode_escape', skiprows=[0], usecols=[col_time, col_dio_first, col_iso, col_dio_sec, col_photo])
#         # Rename columns to standardized names
#         df = df.rename(columns={col_time: "Time(s)",col_dio_first: "TTL_first", col_iso: "isobestic", col_dio_sec: "TTL_sec",  col_photo: "photosignal"})

#     return  df


def load_data(dataframe_1, files_dir1, detector, DIO,analog_i, i_miceperfile, doric_version):
    """
    Load Doric (or CSV) data, supporting a variable number of DIOs per mouse/file.
    
    Parameters
    ----------
    dataframe_1 : str
        File name (either .doric or .csv)
    files_dir1 : str
        Folder path containing the file
    detector : list of str
        Detector IDs per mouse/file
    DIO : list of list of str
        DIO channel numbers per mouse/file, e.g., [['1','2'], ['3','4','5']]
    i_miceperfile : int
        Index of the current mouse/file
    doric_version : int
        Doric version (5 or 6)
    
    Returns
    -------
    pd.DataFrame
        Standardized dataframe with Time(s), TTL columns, isobestic, photosignal
    """

    
    DIO_nums = DIO[i_miceperfile]  # list of DIOs for this mouse/file
    Det_num = detector[i_miceperfile]
    
    if dataframe_1.endswith('.doric'):
        with h5py.File(files_dir1 + dataframe_1, 'r') as f:
            path = '/DataAcquisition/FPConsole/Signals/Series0001/'
            
            # Time vector
            if doric_version == 5:
                time_ = np.array(f[path+'AIN01xAOUT01-LockIn/Time'])
            elif doric_version == 6:
                time_ = np.array(f[path+'LockInAOUT01/Time'])
            
            # Load DIO signals dynamically
            DIO_arrays = []
            for idx, dio_num in enumerate(DIO_nums):
                dio_path = f"{path}/DigitalIO/DIO0{dio_num}"
                DIO_arrays.append(np.array(f[dio_path]))
            
            # Load detector signals
            if doric_version == 5:
                ain1aout1 = np.array(f[path+'AIN0' + Det_num + 'xAOUT01-LockIn/Values'])
                ain1aout2 = np.array(f[path+'AIN0' + Det_num + 'xAOUT02-LockIn/Values'])
            elif doric_version == 6:
                ain1aout1 = np.array(f[path+ 'LockInAOUT01/AIN0' + Det_num])  # isobestic
                ain1aout2 = np.array(f[path+ 'LockInAOUT02/AIN0' + Det_num])  # photosignal

            # Align lengths for all DIOs
            for i in range(len(DIO_arrays)):
                if len(DIO_arrays[i]) != len(time_):
                    DIO_arrays[i] = DIO_arrays[i][20:-20]  # remove start/end buffer
                
                if len(DIO_arrays[i]) > len(time_):
                    difference = len(DIO_arrays[i]) - len(time_)
                    DIO_arrays[i] = DIO_arrays[i][:-difference]
                elif len(DIO_arrays[i]) < len(time_):
                    difference = len(time_) - len(DIO_arrays[i])
                    time_ = time_[:-difference]
                    ain1aout1 = ain1aout1[:-difference]
                    ain1aout2 = ain1aout2[:-difference]

            # Create dataframe
            data_dict = {'Time(s)': time_}
            for idx, arr in enumerate(DIO_arrays):
                data_dict[f'TTL_{idx+1}'] = arr
            data_dict['isobestic'] = ain1aout1
            data_dict['photosignal'] = ain1aout2
            
            df = pd.DataFrame(data_dict)

    elif dataframe_1.endswith('.csv'):
        # Construct dynamic DIO column names
        col_time = 'Time(s)'
        col_dios = [f'DI/O-{dio_num}' for dio_num in DIO_nums]
        col_analog_i = [f'AIn-{a_i}' for a_i in analog_i]
        col_iso = f"AIn-{Det_num} - Dem (AOut-1)"
        col_photo = f"AIn-{Det_num} - Dem (AOut-2)"
        
        # Columns to read
        use_cols = [col_time] + col_dios +col_analog_i +  [col_iso, col_photo]
        
        # Read CSV
        df = pd.read_csv(dataframe_1, encoding='unicode_escape', skiprows=[0], usecols=use_cols)

        
        # Rename DIO columns dynamically
        rename_dict = {col_time: "Time(s)", col_iso: "isobestic", col_photo: "photosignal"}
        for idx, dio_col in enumerate(col_dios):
            rename_dict[dio_col] = f'TTL_{idx+1}'
        df = df.rename(columns=rename_dict)
    
    return df

def create_ID(name, hemiphere,condition, week, batch, genotype):
    # here ID number is created with the following information
    #name_condition_week_batch_genotype
    # Sanitize inputs to avoid illegal filename characters
    
    safe_name = str(name).replace(" ", "_")
    safe_hemiphere = str(hemiphere).replace(" ", "_")
    safe_condition = str(condition).replace(" ", "_")
    safe_week = str(week).replace(" ", "_")
    safe_genotype = str(genotype).replace(" ", "_")
    safe_batch = str(batch).replace(" ", "_")
    
    ID = f"{safe_name}_{safe_hemiphere}_{safe_condition}_{safe_week}_{safe_batch}_{safe_genotype}" 
    return ID

def det_sampling_rate(files, files_dir1, detector, DIO, doric_version, name_empty_files):

    #takes the median of the samples/sec off all files.
    sample_sec_list = []
    for i_file in range(0, len(files)):
        for i_miceperfile in range(0,miceperfile):
            name = files[i_file][ind_names[i_miceperfile][0]: ind_names[i_miceperfile][1]]
            if name.find(name_empty_files) < 0:
                df = load_data(files[i_file],files_dir1, detector, DIO,analog_i, i_miceperfile , doric_version)
                sample_sec_list.append(len(df.index)/df['Time(s)'].iloc[-1])
    
    sampling_rate = int(round(statistics.median(sample_sec_list)* sample_factor))
    print('sampling rates between ' + str(round(min(sample_sec_list))) + ' and ' + str(round(max(sample_sec_list))) + ' , taken median of ' + str(sampling_rate))
    
    return sampling_rate 

def resample_dataframe(df, factor):
    """
    Upsample or downsample a DataFrame by a given factor.
      - factor > 1 : upsample (repeat rows)
      - factor < 1 : downsample (skip rows)
    """
    if factor > 1:
        # Upsample: repeat each row `factor` times
        df = df.loc[df.index.repeat(int(factor))].reset_index(drop=True)
    elif factor < 1:
        # Downsample: keep every int(1/factor)-th row
        step = int(round(1 / factor))
        df = df.iloc[::step, :].reset_index(drop=True)
    return df


def trim_edge_nans(df: pd.DataFrame, column: str) -> pd.DataFrame:
    """
    Remove rows with NaNs at the beginning and end of a DataFrame,
    based on a single reference column.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    column : str
        Column used to detect leading/trailing NaNs
    
    Returns
    -------
    pd.DataFrame
        Trimmed DataFrame
    """
    if column not in df.columns:
        raise KeyError(f"Column '{column}' not found in DataFrame")

    # Boolean mask of non-NaN values in the reference column
    not_nan = df[column].notna()

    if not not_nan.any():
        # Column is all NaNs → return empty df with same columns
        return df.iloc[0:0]

    first_valid = not_nan.idxmax()
    last_valid = not_nan[::-1].idxmax()

    return df.loc[first_valid:last_valid]

def double_exponential(t, const, amp_fast, amp_slow, tau_slow, tau_multiplier):
    '''Compute a double exponential function with constant offset.
    Parameters:
    t       : Time vector in seconds.
    const   : Amplitude of the constant offset. 
    amp_fast: Amplitude of the fast component.  
    amp_slow: Amplitude of the slow component.  
    tau_slow: Time constant of slow component in seconds.
    tau_multiplier: Time constant of fast component relative to slow. 
    '''
    tau_fast = tau_slow*tau_multiplier
    return const+amp_slow*np.exp(-t/tau_slow)+amp_fast*np.exp(-t/tau_fast)


def fit_curve_4_bleaching_correction(signal,time_sec ):

    # Fit curve to signal.
    max_sig = np.max(signal)
    inital_params = [max_sig/2, max_sig/4, max_sig/4, 3600, 0.1]
    bounds = ([0      , 0      , 0      , 600  , 0],
              [max_sig, max_sig, max_sig, 36000, 1])
    data_parms, parm_cov = curve_fit(double_exponential, time_sec, signal, 
                                      p0=inital_params, bounds=bounds, maxfev=1000)
    fitted_signal = double_exponential(time_sec, *data_parms)
    denoised_signal = signal - fitted_signal +1
    return denoised_signal, fitted_signal


def analog_to_DIO(df, analog_i, AI_to_DIO):
    # mounts analog inputs as if they are DIO
    
    col_analog_i = [f'AIn-{a_i}' for a_i in analog_i]
    
    for i in range(len(col_analog_i)):
        col = col_analog_i[i]
        
        # round values and assign back
        df[col] = df[col].round()
        
        # convert all 2s to 1s
        df[col] = df[col].replace(2, 1)
        
        # rename column
        df.rename(columns={col: f'TTL_{AI_to_DIO[i]}'}, inplace=True)
    
    return df
    
    

def artefact_removal(df, ID):
    
    signal = np.array(df['photosignal'])
    control = np.array(df['isobestic'])
    time = np.array(df['Time(s)'])
    
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.3)
    
    ax.plot(time, control, label="Isosbestic")
    ax.plot(time, signal, label="GCaMP")
    fig.suptitle(ID)
    ax.legend()
    
    # Store selections as (xmin, xmax)
    selections = []
    highlight_lines = []
    
    
    # Callback: when user selects a range
    def onselect(xmin, xmax):
        xmin, xmax = sorted([xmin, xmax])  # make sure xmin <= xmax
        selections.append((xmin, xmax))
    
        # Highlight selected range on the plot
        mask = (time >= xmin) & (time <= xmax)
        (line,) = ax.plot(time[mask], control[mask], "r", lw=2)
        highlight_lines.append(line)
    
        fig.canvas.draw()
        print(f"Added selection: {xmin:.2f} to {xmax:.2f}")
    
    # Attach span selector
    span = SpanSelector(ax, onselect, "horizontal", useblit=True)
    
    
    # Define exponential decay model
    def exp_decay(x, A, k):
        return A * np.exp(-k * x)
    
    
    def impute_v2(arr, start, end, decay_rate=0.5):
        arr = np.asarray(arr)
        arr[start:end] = np.nan
        isnan = np.isnan(arr)
        n = len(arr)
        i = 0
    
        while i < n:
            if isnan[i]:
                # Find the full nan run
                start = i
                while i < n and isnan[i]:
                    i += 1
                end = i - 1
    
                left = start - 1
                right = end + 1
    
                if left >= 0 and right < n:
                    left_val = arr[left]
                    right_val = arr[right]
                    gap_size = right - left - 1
    
                    imputed = []
                    for j in range(gap_size):
                        weight = np.exp(-decay_rate * j)
                        val = left_val * weight + right_val * (1 - weight)
                        imputed.append(val)
    
                    return np.array(imputed)
    
                else:
                    raise ValueError(f"Cannot impute range {start}–{end}: missing boundary value(s).")
    
            i += 1
    
        raise ValueError("No missing range found in array.")
    
    def impute(x, y, start, end):
        x = np.arange(len(y))
        y_missing = y.copy()
        y_missing[start:end] = np.nan
    
        # Use only non-NaN values to fit the model
        mask = ~np.isnan(y_missing)
        x_fit = x[mask]
        y_fit = y_missing[mask]
    
        # Fit the exponential model to the available data
        params, _ = curve_fit(exp_decay, x_fit, y_fit, p0=(1, 0.1))  # initial guess (A, k)
        A_fit, k_fit = params
    
        # Predict missing values
        x_missing = x[~mask]
        y_predicted = exp_decay(x_missing, A_fit, k_fit)
    
        return y_predicted
    
    
    # --- Button Callbacks ---
    def undo(event):
        if selections:
            removed = selections.pop()
            line = highlight_lines.pop()
            line.remove()
            fig.canvas.draw()
            print(f"Undid selection: {removed[0]:.2f} to {removed[1]:.2f}")
        else:
            print("No selections to undo.")
    
    
    def clear(event):
        selections.clear()
        for line in highlight_lines:
            line.remove()
        highlight_lines.clear()
        fig.canvas.draw()
        print("All selections cleared.")
    
    
    def submit(event):

        plt.close()
        print("Imputing removed ranges...")
  
        for i, (xmin, xmax) in enumerate(selections, 1):
            start = max(0, np.searchsorted(time, xmin))
            end   = min(len(time), np.searchsorted(time, xmax))
    
            if end <= start:
                print(f"Skipping selection (too small): {xmin:.2f} to {xmax:.2f}")
                continue
        
            try:
                signal_imputed = impute_v2(signal, start, end)
                control_imputed = impute_v2(control, start, end)
        
                signal[start:end] = signal_imputed
                control[start:end] = control_imputed
        
                if 'TTL_first' in df.columns:
                    df.loc[start:end, 'TTL_1'] = np.nan
                if 'TTL_sec' in df.columns:
                    df.loc[start:end, 'TTL_2'] = np.nan
        
            except ValueError as e:
                # suppress or print short warning
                print(f"Warning: {e}. Skipping this range.")
                continue
            
                
    # --- Create Buttons ---
    def make_button(position, label, callback):
        ax_btn = plt.axes(position)
        btn = Button(ax_btn, label)
        btn.on_clicked(callback)
        return btn
    
    
    # Define button layout positions (left, bottom, width, height)
    btn_indo = make_button([0.1, 0.15, 0.2, 0.075], "Undo", undo)
    btn_clear = make_button([0.4, 0.15, 0.2, 0.075], "Clear", clear)
    btn_submit = make_button([0.7, 0.15, 0.2, 0.075], "Submit", submit)
    
    # plt.legend()
    plt.show(block = True)
    plt.close('all')
            
    # save back to df
    df['photosignal'] = signal 
    df['isobestic'] = control 
    
    # Remove any remaining NaNs (usually only at the start or end)
    valid_mask = (~df['photosignal'].isna()) & (~df['isobestic'].isna())
    df = df[valid_mask].reset_index(drop=True)
    
    return df
   

####################################
# script

# inport data
list_files = os.listdir(files_dir1)
files = []
for i_files in range(len(list_files)):
    if list_files[i_files][-4:] == '.csv' or list_files[i_files][-6:] == '.doric':
        files.append(list_files[i_files])
        

# copy of script? Or another way to make mark of this analysis


#make analysed folder if this does not exists yet with current date
os.chdir (files_dir1)
files_dir2 = create_folder_with_note(files_dir1, '1. Preprocessed')

# date = time.strftime('%x'); new_date = date[6:8]; new_date = new_date + date[0:2] ; date = new_date + date[3:5]
# files_dir2 = files_dir1 +  'Preprocessed ' + date  + '//'
# if os.path.exists(files_dir2) == 0:
#     os.makedirs(files_dir2)

# determine sampling rate
sampling_rate = det_sampling_rate(files, files_dir1, detector, DIO, doric_version, name_empty_files)


# load data file for each file
for i_file in range(0,len(files)):
    cur_file = files[i_file]
    week = cur_file[ind_week[0]: ind_week[1]]
    for i_miceperfile in range(0,miceperfile):
        name = cur_file[ind_names[i_miceperfile][0]: ind_names[i_miceperfile][1]]
        hemiphere = cur_file[ind_hemiphere[i_miceperfile][0]: ind_hemiphere[i_miceperfile][1]]
        condition = cur_file[ind_con[i_miceperfile][0]: ind_con[i_miceperfile][1]]
        genotype =  cur_file[ind_gen[i_miceperfile][0]: ind_gen[i_miceperfile][1]] 
        batch = cur_file[ind_batch]
        
        
        
        if name.find(name_empty_files) < 0:
            
            #load data frame for each mouse
            df = load_data(cur_file,files_dir1, detector, DIO,analog_i, i_miceperfile, doric_version)
            ID = create_ID(name,hemiphere, condition, week, batch, genotype)
              
            # start preprocessing the data
            # remove nans in beginning and end of signal
            df = trim_edge_nans(df, "photosignal")
            df = trim_edge_nans(df, "isobestic")
            
            #Plot the start point
            files_dir_fig = files_dir2 + '//figs//'
            if os.path.exists(files_dir_fig) == 0:
                os.makedirs(files_dir_fig)
                    
            plt.plot(df['photosignal']), plt.plot(df['isobestic'])
            plt.savefig(files_dir_fig + ID+ "_raw "  + ".png", transparent = True)
            plt.close()
              
            # set max time
            if Do_max_time == 1:
                time_to_index = int(max_time * 60 * (sampling_rate / sample_factor))
                if len(df) > time_to_index:
                    df = df[0:time_to_index]
                    
            # On unsampling data
            if Do_unsample_data == 1 and sample_factor != 1:
                df = resample_dataframe(df, sample_factor)
                print(f"Data {'up' if sample_factor > 1 else 'down'}sampled by factor {sample_factor}. New length: {len(df)} rows.")
    
            # analog to DIO - only for older data
            if Do_analog_to_dio == 1:
                df = analog_to_DIO(df,analog_i, AI_to_DIO )

            # Artefact removal will be done manually be selecting data that needs to be removed
            if Do_artefact_removal == 1:
                df = artefact_removal(df, ID)
            
            # bleaching correction
            if Do_bleaching_correction == 1:
                df['photosignal'], df['gfp_expfit'] = fit_curve_4_bleaching_correction(df['photosignal'], df["Time(s)"])
                df['isobestic'], iso_expfit = fit_curve_4_bleaching_correction(df['isobestic'], df["Time(s)"])
            
            
            
            #low pass filter
            if Do_low_pass_filter == 1:
                # Lowpass filter - zero phase filtering (with filtfilt) is used to avoid distorting the signal.
                b,a = butter(2, 10, btype='low', fs=sampling_rate)
                df['photosignal'] = filtfilt(b,a, df['photosignal'])
                df['isobestic'] = filtfilt(b,a, df['isobestic'])
                #high pass filter
                # b,a = butter(2, 0.001, btype='high', fs=sampling_rate)
                # df['photosignal'] = filtfilt(b,a, df['photosignal'])
                # df['isobestic'] = filtfilt(b,a, df['isobestic'])
            
            
            if Do_motion_correction == 1:
                #We now do motion correction by finding the best linear fit of the TdTomato signal to the dLight signal 
                #and subtracting this estimated motion component from the dLight signal. We will use the data that was bleaching corrected using the
                #double exponential fit as this is less likely to remove meaningful slow variation in the signals.
                slope, intercept, r_value, p_value, std_err = linregress(x=df['photosignal'] , y=df['isobestic'])
                est_motion = intercept + slope * np.array(df['isobestic'])
                df['MC'] = np.array(df['photosignal']) - est_motion
                
            #Smoothing 
            if Do_smoothing == 1:
                smoothing_ind = round(smoothing_time * sampling_rate)               
                df['photosignal'] = (df['photosignal'].rolling(window=smoothing_ind, center=True, min_periods=1).mean())
                df['isobestic'] = (df['isobestic'].rolling(window=smoothing_ind, center=True, min_periods=1).mean())
                df['MC'] = (df['MC'].rolling(window=smoothing_ind, center=True, min_periods=1).mean())
                               
            
            plt.plot(df['photosignal']), plt.plot(df['isobestic'])
            plt.savefig(files_dir_fig + ID + "_preprocessed" + ".png", transparent = True)
            plt.close('all')
            
            if Do_smoothing == 1:
                plt.plot(df['MC'])
                plt.savefig(files_dir_fig + ID + "_MC" + ".png", transparent = True)
                plt.close('all')
            
            
            
            # save data
            # Save DataFrame in Parquet format with a unique filename
            filepath = os.path.join(files_dir2, ID + '.parquet')
            df.to_parquet(filepath, index=False)
            print(f"✅ DataFrame saved as: {ID}")

        
            
            
