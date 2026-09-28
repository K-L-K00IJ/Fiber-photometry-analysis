# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 11:17:45 2025
4. Plotting light difference of fiber photometry data
general script

This script will perform plotting functions of:
    general average score of iso and gfp per time lock per mouse
    heat map for selected time locks
    calculation of averages

Inputs: 
    Per mouse per timelock 1 DF of per trial results and 1 df of all average scores
    Both deltaF/F and z_scores are calcuated

Output:
    graphs

Step 4 out of 4 in the photometry analysis pipeline
Meye lab, Karlijn Kooij
Version 5_12_25


@author: kkooij
"""

# enter the following parameters per analysis

# enter the directory where the raw data files are stored

files_dir_general =    'L://adanlab//Ongoing//apetitbon//Projects//Lep_DA_CeA//Pilot WP 01 - 19//DORIC-Pilot_GRAB-DA3h_CeA_LepRcre_VTA_opto-ChRimson//20260430//'
#files_dir_general = 'D://GRAB DA//5C//'

folder_light_signals = files_dir_general + '3. Light difference//'

folder_exp_specific_info = files_dir_general + 'exp_specific_info//'

# plot features
time_window = [-3, 5]
x_ticks = [-2.5, 0, 2.5, 5]

# exclusion
# based on file removal folder
# experiment sepcific info in excels of exp_specific_info folder

# lay out
colors_GFP_ISO = ['green','magenta', 'black']


# labels
y_label_dff = 'ΔF/F'
y_lim_dff = [-1,1]
y_lim_z = [-0.5,0.5]
y_label_z = 'Z-score'

# file extention
file_extentions = ['.png', '.pdf']


# individual plots
ex_GFP_iso_plots=1
ex_heatmaps =1
#averages
ex_averages =1
ex_bar_plots = 1
ex_within_session = 0
# statistics
ex_statistics = 1 # only in combi with ex_bar_plots




# info heat maps
time_window_heatmap = [-2.5, 5]
x_ticks_heatmap = [-2.5, 0, 2.5, 5]
tl_for_heatmap = ["onsets_first_TTL"]
calc_for_heatmap = ['z_score_MC']

# for other plots
cols_to_plot = ['dff_photo', 'dff_iso','dff_MC', 'z_score_MC']
y_labels_to_plot = [y_label_dff,y_label_dff, y_label_dff, y_label_z]
y_lims_to_plot = [y_lim_dff,y_lim_dff,y_lim_dff, y_lim_z]
plot_seperate_conditions = 0 # add 0 for conditions in the same plot and 1 for different plots

cols_to_plot = ['z_score_MC']
y_labels_to_plot = [y_label_z]
y_lims_to_plot = [y_lim_z]


#bars
window = [-5,5, 5]  # seconds
outcomes = ['z_score_MC']
time_windows = [(-5,0), (0,5)]
y_lim_bars =  [-1, 1]

#within session
trials_to_compare = [[0,5], [15,20]]

#paired on unpaired conditions
paired_con = 'unpaired'

# note here if you again execute the script or for the first time in the console
re_execute = 0



##############################################
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


# set default plot settings

#sns.set_style("white")

plt.rcParams['font.family'] = 'Arial'
font_size = 25
plt.rcParams['axes.labelsize'] = font_size
plt.rcParams['axes.titlesize'] = font_size
plt.rcParams['xtick.labelsize'] = font_size
plt.rcParams['ytick.labelsize'] = font_size
plt.rcParams['legend.fontsize'] = 18
plt.rcParams['font.weight'] = 'regular'
plt.rcParams['axes.labelweight'] = 'regular'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6

##############################################
# general functions - to kick out
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

    full_path = full_path +  '//'
    return full_path


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


def load_files_from_folder(files_dir, extention):
  # inport data from 1 folder with a certain extention 
  list_files = os.listdir(files_dir)
  files = []
  for i_files in range(len(list_files)):
      if list_files[i_files].endswith(extention):
          files.append(list_files[i_files])
          
  return files


def load_signal_data(folder_light_signals):
    # load files
    files_list = load_files_from_folder(os.path.join(folder_light_signals,'average'), '.parquet')
    
    res_pt_list = []; res_avr_list = []
    
    for i in range(len(files_list)):
        signals_pt = pd.read_parquet(os.path.join(folder_light_signals,'per trial', files_list [i]), engine= 'pyarrow', memory_map = True, use_threads=False )
        signals_avr = pd.read_parquet(os.path.join(folder_light_signals,'average', files_list [i]), engine= 'pyarrow', memory_map = True, use_threads=False)
        
        
        # stack all together now it is per mouse
        res_pt_list.append(signals_pt)
        res_avr_list.append(signals_avr)
        
        
    # Concatenate all DataFrames at once
    df_signals_pt = pd.concat(res_pt_list, ignore_index=True)
    df_signals_avr = pd.concat(res_avr_list, ignore_index=True)
    
    return df_signals_pt, df_signals_avr


############################################################

def plot_GFP_ISO_per_tl_per_rec(df,col_to_plot, color, y_lim, tl_to_xlabel, y_label, files_dir_plots, time_window=None, x_ticks=None ):
    
    # loop through df
    for i_row in range(len(df_signals_avr)):
        fig, ax = plt.subplots()       # create a new figure
        plt.ylim(y_lim[0],y_lim[1])
        
        # get right x label
        cur_tl = df.loc[i_row, 'time_lock_name']
        x_label = tl_to_xlabel.get(cur_tl, 'Time (s)')  # default label if TL not found
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        
        for i_channel in range(len(col_to_plot)):
            cur_tl = df.loc[i_row, 'time_lock_name']
            cur_x = df.loc[i_row, 'time_vec']
            cur_y = df.loc[i_row, col_to_plot[i_channel]]
            
            if cur_y is not None:
                # Apply time window if provided
                if time_window is not None:
                    mask = (cur_x >= time_window[0]) & (cur_x <= time_window[1])
                    cur_x = cur_x[mask]
                    cur_y = cur_y[mask]
    
                ax.plot(cur_x, cur_y, color[i_channel], label = col_to_plot[i_channel])
        
        # Set x-ticks if provided
        if x_ticks is not None:
            ax.set_xticks(x_ticks)
        
        #create ID
        ID = create_ID(df.loc[i_row, 'name'], df.loc[i_row, 'hemisphere'],df.loc[i_row, 'condition'], df.loc[i_row, 'week'], df.loc[i_row, 'batch'], df.loc[i_row, 'genotype'])
       
        fig.tight_layout()   # adjusts spacing automatically
        fig.savefig(files_dir_plots + ID + '_'+  cur_tl+  '.png', transparent = True)
        
        save_legend_from_ax( ax, files_dir_plots, 'GFP ISO')
        
        plt.close(fig)                 # close the figure to free memory
        

def plot_avg_signals_per_condition(df, col_to_plot, group_cols, files_dir_plots, tl_to_xlabel,
                     y_lim=(0, 1), title_to_add='_', colors=None, time_window=None,
                     x_ticks=None, y_label = '', plot_seperate = 0, plot_extentions = ['.png']):
    """
    Plot average signals per group across multiple rows for each time lock.
    
    Parameters:
        df: pandas.DataFrame with signals
        col_to_plot: list of signal column names to plot
        group_cols: list of column names to group by (e.g., ['condition'], ['genotype', 'hemisphere'])
        files_dir_plots: directory to save figures
        tl_to_xlabel: dict mapping time_lock_name to x-label
        y_lim: tuple of (ymin, ymax)
        title_to_add: string to add to filenames
        colors: list of colors for each group (optional)
        time_window: tuple (start, end) to limit x-axis
        x_ticks: list of x-tick positions (optional)
    """
    
    # Get unique time locks
    time_locks = df['time_lock_name'].unique()
    
    for tl in time_locks:
        # Filter data for current time lock
        df_tl = df[df['time_lock_name'] == tl]
        
        # Group by the requested columns
        grouped = df_tl.groupby(group_cols)
          
        fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
        ax.set_ylim(y_lim)
        # Loop through groups
        for i, (group_name, group_df) in enumerate(grouped):
            # If group_name is not a tuple, make it a tuple for consistency
            if not isinstance(group_name, tuple):
                group_name = (group_name,)
            

            
            # Average signals across rows
            for col_idx, col in enumerate(col_to_plot):
                # Stack all signals in group and average
                filtered = group_df[group_df[col].notna()]
                cur_signals = np.stack(filtered[col].values)
                mean_signal = np.nanmean(cur_signals, axis=0)
                sem = np.nanstd(cur_signals, axis=0) / np.sqrt(cur_signals.shape[0])
                x_vec = group_df.iloc[0]['time_vec']
                
                # Apply time window if provided
                if time_window is not None:
                    mask = (x_vec >= time_window[0]) & (x_vec <= time_window[1])
                    x_vec = x_vec[mask]
                    mean_signal = mean_signal[mask]
                    sem = sem[mask]
                
                # Determine main condition for the group
                main_condition = group_name[0] if isinstance(group_name, tuple) else group_name
                
                # Try TL color first
                c = colors.get(tl, None)
                if c is None or c in [0, "0"] or (isinstance(c, float) and np.isnan(c)):
                    # Fallback to group/condition color
                    c = colors.get(main_condition, None)
                    if c is None or c in [0, "0"] or (isinstance(c, float) and np.isnan(c)):
                        # Final fallback color
                        c = None

                label_name = "_".join([str(g) for g in group_name]) + f"_{col}"
                ax.plot(x_vec, mean_signal, c=c, label=label_name)
                
                # Plot SEM as shaded area
                ax.fill_between(x_vec, mean_signal - sem, mean_signal + sem, color=c, alpha=0.3)
        
        
                # Labels
                x_label = tl_to_xlabel.get(tl, 'Time (s)')
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                
                # Optional x-ticks
                if x_ticks is not None:
                    ax.set_xticks(x_ticks)
                    
                if plot_seperate == 1:
                    
                    # Save figure
                    fig_name = f"{files_dir_plots}{tl}{title_to_add}-{main_condition}"
                    for i_exnt in range(len(plot_extentions)):
                        fig.savefig(fig_name +plot_extentions[i_exnt] , transparent=True)
                        save_legend_from_ax( ax, files_dir_plots, title_to_add, plot_extentions = plot_extentions)
                    # fig, ax = plt.subplots(constrained_layout=True)
                    # ax.set_ylim(y_lim)
                    plt.close(fig)
 
            
        if plot_seperate == 0:
            # Save figure
            fig_name = f"{files_dir_plots}{tl}{title_to_add}"
            for i_exnt in range(len(plot_extentions)):
                fig.savefig(fig_name+plot_extentions[i_exnt], transparent=True)
                save_legend_from_ax( ax, files_dir_plots, title_to_add,plot_extentions =  plot_extentions)
            plt.close(fig)
            
            
            
            
        
def plot_avg_signals_per_time_lock(
    df, 
    col_to_plot, 
    condition_filter,       # list of conditions
    time_locks_to_plot,     # list of time lock names to plot
    files_dir_plots=None,
    tl_to_xlabel=None,
    y_lim=(0, 1),
    title_to_add='_',
    colors=None,            # dict mapping time_lock_name -> color
    time_window=None,
    x_ticks=None,
    y_label=''
):
    """
    Plot multiple time locks for a single condition (or group) as line plots with SEM.
    
    Parameters:
        df: pandas.DataFrame with signals
        col_to_plot: list of signal column names to plot
        condition_filter: dict of column_name -> value to select which condition/group to plot
        time_locks_to_plot: list of time lock names to plot
        files_dir_plots: directory to save figures
        tl_to_xlabel: dict mapping time_lock_name to x-label
        y_lim: tuple of (ymin, ymax)
        title_to_add: string to add to filenames
        colors: dict mapping time_lock_name -> color
        time_window: tuple (start, end) to limit x-axis
        x_ticks: list of x-tick positions
        y_label: label for y-axis
    """
    
    # Filter dataframe by the condition(s)
    df_filtered = df.copy()
    for i_con in range(len(condition_filter)): 
        df_filtered = df_filtered[df_filtered['condition'] == condition_filter[i_con]]
    
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_ylim(y_lim)
    
    for tl in time_locks_to_plot:
        df_tl = df_filtered[df_filtered['time_lock_name'] == tl]
        if df_tl.empty:
            print(f"Warning: no data for time lock '{tl}' with filter {condition_filter}")
            continue
        
        for col in col_to_plot:
            # Stack all signals for this time lock
            cur_signals = np.stack(df_tl[col].values)
            mean_signal = np.mean(cur_signals, axis=0)
            sem = np.std(cur_signals, axis=0) / np.sqrt(cur_signals.shape[0])
            x_vec = df_tl.iloc[0]['time_vec']
            
            # Apply time window if provided
            if time_window is not None:
                mask = (x_vec >= time_window[0]) & (x_vec <= time_window[1])
                x_vec = x_vec[mask]
                mean_signal = mean_signal[mask]
                sem = sem[mask]
            
            c = colors.get(tl, None) if colors else None
            ax.plot(x_vec, mean_signal, c=c, label=f"{tl}_{col}")
            ax.fill_between(x_vec, mean_signal - sem, mean_signal + sem, color=c, alpha=0.3)
    
    # Labels
    ax.set_xlabel('Time (s)' if tl_to_xlabel is None else ', '.join([tl_to_xlabel.get(tl, tl) for tl in time_locks_to_plot]))
    ax.set_ylabel(y_label)
    
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    
    if files_dir_plots:
        fig_name = f"{files_dir_plots}time_locks_{'_'.join(condition_filter[i_con])}{title_to_add}.png"
        fig.savefig(fig_name, transparent=True)
    
    ax.legend()
    plt.close(fig)        
        
        
def save_legend_from_ax( ax, save_dir, legend_title_add="", filename_prefix="legend", fig_size=(3, 1.6), plot_extentions = ['.png'] ):
    #Creates and saves a standalone legend using the handles from the provided axis.

    handles, labels = ax.get_legend_handles_labels()

    if len(handles) == 0:
        print("Warning: No legend handles found in axis — no legend was saved.")
        return

    # Create a separate figure only for the legend
    fig_legend = plt.figure(figsize=fig_size)
    fig_legend.legend( handles,   labels,    frameon=False    )
    plt.axis("off")
    fig_legend.tight_layout()   # adjusts spacing automatically

    # Build filename
    fname = f"{filename_prefix}{legend_title_add}"
    for i_exnt in range(len(plot_extentions)):
        out_path = os.path.join(save_dir , fname+ plot_extentions[i_exnt])
        fig_legend.savefig(out_path, transparent=True, dpi=300)
    plt.close(fig_legend)

def plot_heatmap_signals(df, time_locks,columns, time_window=None,  output_dir="", cmap="coolwarm", vmax=None,  vmin=None, tl_to_xlabel=None, x_ticks = None, fig_size=(6,4),   dpi=200  ):
    #Creates heatmaps for selected time locks and selected signal columns.

    for tl in time_locks:
        df_sel = df[df['time_lock_name'] == tl]

        # reference time vector
        t_full = np.array(df_sel.iloc[0]["time_vec"])

        # apply time window
        if time_window is not None:
            t_min, t_max = time_window
            idx_time = (t_full >= t_min) & (t_full <= t_max)
            t = t_full[idx_time]
        else:
            idx_time = slice(None)
            t = t_full

        for col in columns:
            if col not in df_sel.columns:
                print(f"⚠ Column not found: {col}")
                continue

            # loop over each sample (row) in this TL group
            for row_idx, row in df_sel.iterrows():
 
                 #create ID
                 ID = create_ID(df.loc[row_idx, 'name'], df.loc[row_idx, 'hemisphere'],df.loc[row_idx, 'condition'], df.loc[row_idx, 'week'], df.loc[row_idx, 'batch'], df.loc[row_idx, 'genotype'])
                
                 # each sample's row is a 1D signal vector → make it a 2D (trials × time)
                 data_signal = df_sel.loc[row_idx][col]
                 data_signal = np.vstack(data_signal)
 
                 # if already 2D (trials × time) keep shape
                 if data_signal.ndim == 1:
                     data_mat = data_signal[idx_time][None, :]  # 1 row heatmap
                 else:
                     data_mat = data_signal[:, idx_time]
 
                 fig, ax = plt.subplots(figsize=fig_size, dpi=dpi)
 
                 im = ax.imshow( data_mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, origin="upper"  )
 
                # Use provided x_ticks if given
                 if x_ticks is not None:
                    # Find closest indices in t for the requested x_ticks
                    tick_idx = [np.argmin(np.abs(t - xt)) for xt in x_ticks]
                    ax.set_xticks(tick_idx)
                    ax.set_xticklabels([str(xt) for xt in x_ticks])
                 else:
                    # Default: 6 evenly spaced ticks
                    n_ticks = 6
                    tick_idx = np.linspace(0, len(t)-1, n_ticks, dtype=int)
                    ax.set_xticks(tick_idx)
                    ax.set_xticklabels(np.round(t[tick_idx], 1))
 
                 ax.set_ylabel("Trials")
 
                 xlabel = tl_to_xlabel.get(tl, tl)
                 ax.set_xlabel(xlabel)
 
                 # colorbar
                 cbar = fig.colorbar(im, ax=ax)
                 cbar.ax.tick_params(labelsize=9)
 
                 fname = f"HM_{col}_{tl}_{ID}.png"
                 
                 fig.tight_layout()   # adjusts spacing automatically
                 fig.savefig(os.path.join(output_dir, fname), transparent=True)
                 plt.close(fig)

def calculate_within_session(df, id_cols, trials_to_compare, colors):

    df_base = df[id_cols + ['z_score_MC']]

    def to_array(x):
        arr = np.asarray(x)

        # force object arrays (ragged) to become numeric-safe structure
        if arr.dtype == object:
            arr = np.stack([np.asarray(i) for i in x])

        return arr

    def mean_first(x):
        x = to_array(x)
        return np.nanmean(x[:trials_to_compare], axis=0)

    def mean_last(x):
        x = to_array(x)
        return np.nanmean(x[-trials_to_compare:], axis=0)

    first_df = df_base[id_cols].copy()
    first_df['period'] = f'first_{trials_to_compare}'
    first_df['z_score_MC'] = df_base['z_score_MC'].apply(mean_first)

    last_df = df_base[id_cols].copy()
    last_df['period'] = f'last_{trials_to_compare}'
    last_df['z_score_MC'] = df_base['z_score_MC'].apply(mean_last)
    
    colors_dict = {
    f'first_{trials_to_compare}': colors_first_last[0],
    f'last_{trials_to_compare}': colors_first_last[1]
    }

    return pd.concat([first_df, last_df], ignore_index=True), colors_dict
########################################################

# Clear all variables (Spyder / IPython)
# try:
#     from IPython import get_ipython
#     get_ipython().magic("reset -f")
# except Exception:
#     pass

# df_signals_pt = pd.read_parquet(folder_light_signals + 'signals_per_trial.parquet' )
# df_signals_avr = pd.read_parquet(folder_light_signals + 'signals_averaged.parquet' )



#experiment specific
tl_df = pd.read_excel( folder_exp_specific_info  + "time_lock_labels.xlsx")
tl_to_xlabel_line = dict(zip(tl_df['time_lock_name'], tl_df['x_label_avr']))
tl_to_xlabel_bar = dict(zip(tl_df['time_lock_name'], tl_df['y_label_bar']))

colors_df = pd.read_excel(folder_exp_specific_info + "colors_to_conditions.xlsx",  dtype={"condition": str} )
colors_to_conditions = dict(zip(colors_df['condition'], colors_df['color']))
tl_to_color = dict(zip(tl_df['time_lock_name'], tl_df['TL_color']))
colors_dict = {**colors_to_conditions, **tl_to_color}


conditions = list(colors_df['condition'])
if 'Label_line' in colors_df.columns:  # Check if the column exists
    groups = list(colors_df['Label_line'])
else:
    groups = conditions  # Or handle it however you like if the column doesn't exist
xlabel_bar = list(colors_df['Label_bar'])



df_signals_pt, df_signals_avr = load_signal_data(folder_light_signals)
    
   
        
files_dir_plots = create_folder_with_note(files_dir_general, '4. Plots')

# set default plot settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams['font.weight'] = 'regular'
plt.rcParams['axes.labelweight'] = 'regular'
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6



##################################################################################
# individual plots


# general GFP vs ISO plot of each time lock
if ex_GFP_iso_plots == 1:
    col_to_plot = ['dff_photo', 'dff_iso', 'dff_MC']
    files_dir_dff = create_folder_with_note(files_dir_plots, 'GFP_ISO_dff')
    plot_GFP_ISO_per_tl_per_rec(df_signals_avr,col_to_plot, colors_GFP_ISO, y_lim_dff, tl_to_xlabel_line, y_label_dff, files_dir_dff, time_window, x_ticks  )
    
    col_to_plot = ['z_score_photo', 'z_score_iso', 'z_score_MC']
    
    files_dir_z = create_folder_with_note(files_dir_plots, 'GFP_ISO_z_score')
    plot_GFP_ISO_per_tl_per_rec(df_signals_avr,col_to_plot, colors_GFP_ISO, y_lim_z, tl_to_xlabel_line, y_label_z, files_dir_z,  time_window, x_ticks  )
    


# heat map
if ex_heatmaps == 1: 
    files_dir_heat_map = create_folder_with_note(files_dir_plots, 'Heat maps')
    
    plot_heatmap_signals(
        df=df_signals_pt,
        time_locks=tl_for_heatmap,
        columns=calc_for_heatmap,
        time_window=time_window_heatmap,
        output_dir=files_dir_heat_map,
        tl_to_xlabel=tl_to_xlabel_line, 
        x_ticks = x_ticks_heatmap)


##################################################################################
# average plots

if ex_averages == 1: 

    # making average plots based on conditions
    files_dir_dff_avr = create_folder_with_note(files_dir_plots, 'Averages line')
        
    for col in cols_to_plot:
        
        plot_avg_signals_per_condition(
            df=df_signals_avr,
            col_to_plot=[col],
            group_cols=['condition'],
            files_dir_plots=files_dir_dff_avr,
            tl_to_xlabel=tl_to_xlabel_line,
            y_lim= y_lims_to_plot[cols_to_plot.index(col)],
            title_to_add = cols_to_plot[cols_to_plot.index(col)], 
            colors=colors_dict,
            time_window=time_window,
            x_ticks= x_ticks,
            y_label = y_labels_to_plot[cols_to_plot.index(col)], plot_seperate= plot_seperate_conditions, plot_extentions = file_extentions)
    
##########################################################################
if ex_within_session == 1:
    # making average plots based on conditions
    files_dir_session = create_folder_with_note(files_dir_plots, 'Within session')
    #trials_to_compare
    
    # create df where you make the trials to compare for 1 condition. Then we group option can be used as between trials
    id_cols = ['name', 'hemisphere', 'condition', 'week',
    'batch', 'genotype', 'time_lock_name', 'time_vec']
    
    nr_trials_first_last = 3
    colors_first_last = ['black', 'grey']

    within_sessions_df, colors_session_dict = calculate_within_session(df_signals_pt, id_cols, nr_trials_first_last, colors_first_last)

    for col in cols_to_plot:
        for con in conditions:
            cur_df = within_sessions_df[within_sessions_df['condition'] == con]
            plot_avg_signals_per_condition(
                df=cur_df,
                col_to_plot=[col],
                group_cols=['period'],
                files_dir_plots=files_dir_session,
                tl_to_xlabel=tl_to_xlabel_line,
                y_lim= y_lims_to_plot[cols_to_plot.index(col)],
                title_to_add = cols_to_plot[cols_to_plot.index(col)] + '_' + con , 
                colors=colors_session_dict,
                time_window=time_window,
                x_ticks= x_ticks,
                y_label = y_labels_to_plot[cols_to_plot.index(col)], plot_seperate= 0)
    
        
#######################################################################################
# bar plots

from pathlib import Path
import numpy as np
import pandas as pd
import os

def aggregate_signals_in_time_windows(
    df,
    outcome_cols,
    window,
    meta_cols,
    output_dir,
    output_excel=None,
    mode="average"  # "average", "peak", "dip"
):
    """
    Aggregate session-level time-locked signals in defined time windows.

    Parameters
    ----------
    mode : str
        "average" -> mean
        "peak"    -> max
        "dip"     -> min
    """

    # Select aggregation function
    if mode == "average":
        agg_func = np.nanmean
        default_name = "averaged_signals_windows.xlsx"
    elif mode == "peak":
        agg_func = np.nanmax
        default_name = "peak_signals_windows.xlsx"
    elif mode == "dip":
        agg_func = np.nanmin
        default_name = "dip_signals_windows.xlsx"
    else:
        raise ValueError("mode must be 'average', 'peak', or 'dip'")

    # Handle output filename
    if output_excel is None:
        output_excel = default_name

    start_sec, end_sec, step_sec = window
    time_bins = np.arange(start_sec, end_sec, step_sec)

    results = []

    for _, row in df.iterrows():

        time_vec = np.asarray(row["time_vec"])

        # Copy metadata
        base_row = {col: row[col] for col in meta_cols}

        for outcome in outcome_cols:

            signal = row[outcome]

            # Handle invalid signals
            if not isinstance(signal, (np.ndarray, list)) or len(signal) != len(time_vec):
                for t0 in time_bins:
                    t1 = t0 + step_sec
                    base_row[f"{outcome}_{t0:.1f}_{t1:.1f}_s"] = np.nan
                continue

            signal = np.asarray(signal)

            for t0 in time_bins:
                t1 = t0 + step_sec

                idx = np.where((time_vec >= t0) & (time_vec < t1))[0]

                if len(idx) == 0:
                    val = np.nan
                else:
                    val = agg_func(signal[idx])

                base_row[f"{outcome}_{t0:.1f}_{t1:.1f}_s"] = val

        results.append(base_row)

    out_df = pd.DataFrame(results)

    # Save to Excel
    os.makedirs(output_dir, exist_ok=True)
    out_df.to_excel(os.path.join(output_dir, output_excel), index=False)

    return out_df


def make_window_colname(outcome, t0, t1):
    """
    Create column names that exactly match the averaging output.
    """
    return f"{outcome}_{t0:.1f}_{t1:.1f}_s"



def plot_per_timepoint_per_condition(df, outcome, time_windows, condition_col, conditions, subject_col, colors,x_label, ylabel, output_dir, paired_con='unpaired', show_individuals=True, y_lim=None, tl_to_xlabel = None, data_type = 'avr', plot_extentions = ['.png']):

    os.makedirs(output_dir, exist_ok=True)

    time_locks = sorted(df["time_lock_name"].unique())

    for tl in time_locks:

        df_tl = df[df["time_lock_name"] == tl]
        if df_tl.empty:
            continue

        for t0, t1 in time_windows:

            col = make_window_colname(outcome, t0, t1)
            if col not in df_tl.columns:
                print(f"⚠ Column not found: {col} (TL: {tl})")
                continue

            fig, ax = plt.subplots(figsize=(3.5, 5), constrained_layout=True)

            means, sems = [], []
            per_condition_subject_vals = {}

            for cond in conditions:
                vals = (
                    df_tl[df_tl[condition_col] == cond]
                    .groupby([subject_col, 'hemisphere'])[col]
                    .mean()
                    .astype(float)
                )
                per_condition_subject_vals[cond] = vals
                means.append(np.nanmean(vals))
                sems.append(np.nanstd(vals) / np.sqrt(vals.notna().sum()))

            x = np.arange(len(conditions))

            ax.bar(x, means, yerr=sems, capsize=6, color=[colors[c] for c in conditions], edgecolor="black", zorder=1)

            if show_individuals:

                if paired_con == 'unpaired':
                    for i, cond in enumerate(conditions):
                        vals = per_condition_subject_vals[cond]
                        jitter = np.random.normal(0, 0.05, size=len(vals))
                        ax.scatter(np.full(len(vals), i) + jitter, vals, color="black", alpha=0.5, s=25, zorder=3)

                elif paired_con == 'paired':
                    common_subjects = set.intersection(*[set(v.index) for v in per_condition_subject_vals.values()])
                    for subj in common_subjects:
                        ys = [per_condition_subject_vals[c].loc[subj] for c in conditions]
                        ax.plot(x, ys, color="black", alpha=0.4, linewidth=1, zorder=2)
                        #ax.scatter(x, ys, color="black", s=25, zorder=3)

                else:
                    raise ValueError("paired_con must be 'paired' or 'unpaired'")

            ax.set_xticks(x)
            ax.set_xticklabels(x_label)
            
            ylabel_tl = tl_to_xlabel.get(tl, tl)
            ax.set_ylabel(ylabel_tl + ' (' + ylabel[0] + '$_{' + ylabel[1] + '}$)')


            if y_lim is not None:
                ax.set_ylim(y_lim)
            
            #fig.tight_layout()
            #ax.set_title(f"{tl} ({t0}–{t1} s)")

            fname = f"{tl}_{outcome}_{t0}_{t1}_{data_type}"
            for i_exnt in range(len(plot_extentions)):
                fig.savefig(os.path.join(output_dir, fname + plot_extentions[i_exnt]), dpi=300, transparent=True)
            plt.close(fig)

def plot_all_timepoints_per_condition(df, outcome, time_windows, condition_col, conditions, subject_col, colors, ylabel, output_dir, show_individuals=True, y_lim = None, tl_to_xlabel = None, data_type = 'avr', plot_extentions = ['.png']):
    os.makedirs(output_dir, exist_ok=True)

    x_labels = [f"{t0}–{t1}" for t0, t1 in time_windows]
    x = np.arange(len(time_windows))

    # ---- loop over all time locks ----
    time_locks = sorted(df["time_lock_name"].unique())

    for tl in time_locks:
        df_tl = df[df["time_lock_name"] == tl]

        if df_tl.empty:
            continue
        for cond in conditions:

            fig, ax = plt.subplots(figsize=(3.5, 5), constrained_layout=True)

            cond_df = df_tl[df_tl[condition_col] == cond]

            if cond_df.empty:
                plt.close(fig)
                continue

            means = []
            sems = []
            subject_vals = []

            # ---- compute per-window stats (1 value per subject) ----
            for t0, t1 in time_windows:
                col = make_window_colname(outcome, t0, t1)

                if col not in cond_df.columns:
                    means.append(np.nan)
                    sems.append(np.nan)
                    subject_vals.append(pd.Series(dtype=float))
                    continue

                vals = (
                    cond_df
                    .groupby(subject_col)[col]
                    .mean()
                    .astype(float) )

                subject_vals.append(vals)
                means.append(np.nanmean(vals))
                sems.append(np.nanstd(vals) / np.sqrt(vals.notna().sum()))

            ax.bar(
                x,
                means,
                yerr=sems,
                capsize=6,
                color=colors.get(cond, "grey"),
                edgecolor="black"
            )

            # ---- individual subjects ----
            if show_individuals:
                subj_df = (
                    cond_df
                    .groupby([subject_col, 'hemisphere'])
                    [[make_window_colname(outcome, t0, t1) for t0, t1 in time_windows]]
                    .mean()
                )

                for _, row in subj_df.iterrows():
                    ax.plot(
                        x,
                        row.values.astype(float),
                        color="black",
                        alpha=0.3,
                        linewidth=0.8,
                        zorder=2)

            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
            ylabel_tl = tl_to_xlabel.get(tl, tl)
            ax.set_ylabel(ylabel_tl + ' (' + ylabel[0] + '$_{' + ylabel[1] + '}$)')
            
            # ---- manual y-limits if provided ----
            if y_lim is not None:
                ax.set_ylim(y_lim)

            # ---- title includes time lock ----
            #ax.set_title(f"{cond} – {tl}")
            #fig.tight_layout()

            fname = f"{tl}_{outcome}_{cond}_{data_type}"
            for i_exnt in range(len(plot_extentions)):
                fig.savefig(os.path.join(output_dir, fname + plot_extentions[i_exnt]),  dpi=300, transparent=True )
            plt.close(fig)

def plot_bar_over_time(df, outcome, time_windows, condition_col, conditions, subject_col, colors,xlabel, ylabel, output_dir, show_individuals=True, bar_width=0.35, y_lim=None, tl_to_xlabel = None,data_type = 'avr', plot_extentions = ['.png']):
    """Bar plot of averaged signals over time windows per condition, per time lock, with optional y-limits."""
    import os, numpy as np
    os.makedirs(output_dir, exist_ok=True)
    time_locks = sorted(df["time_lock_name"].unique())

    for tl in time_locks:
        df_tl = df[df["time_lock_name"] == tl]
        if df_tl.empty: 
            continue

        data_cols = [make_window_colname(outcome, t0, t1) for (t0, t1) in time_windows]
        fig, ax = plt.subplots(figsize=(3.5, 5))
        opacity = 0.85
        n_cond = len(conditions)
        n_time = len(time_windows)
        x_base = np.arange(n_time)

        for i_cond, cond in enumerate(conditions):
            cond_df = df_tl[df_tl[condition_col] == cond]
            if cond_df.empty: 
                continue

            means = [np.nanmean(cond_df[col].astype(float)) for col in data_cols]
            sems = [np.nanstd(cond_df[col].astype(float)) / np.sqrt(cond_df[col].notna().sum()) for col in data_cols]
            x_pos = x_base + i_cond * bar_width

            ax.bar(
                x_pos,
                means,
                yerr=sems,
                width=bar_width,
                color=colors.get(cond, "grey"),
                edgecolor="black",
                capsize=6,
                alpha=opacity,
                label=cond
            )

            if show_individuals:
                for _, row in cond_df.iterrows():
                    y_vals = [row[col] for col in data_cols]
                    ax.plot(
                        x_pos,
                        y_vals,
                        color="black",
                        alpha=0.3,
                        linewidth=0.8
                    )

        ax.set_xticks(x_base + bar_width * (n_cond - 1) / 2)
        ax.set_xticklabels([f"{t0}–{t1}s" for t0, t1 in time_windows])
        ylabel_tl = tl_to_xlabel.get(tl, tl)
        ax.set_ylabel(ylabel_tl + ' (' + ylabel[0] + '$_{' + ylabel[1] + '}$)')
        #ax.set_title(tl)  # <-- title is now the time lock name

        if y_lim is not None: 
            ax.set_ylim(y_lim)

        #plt.tight_layout()
        fname = f"{tl}_{outcome}_{data_type}_all_con_all_times"
        for i_exnt in range(len(plot_extentions)):
            plt.savefig(os.path.join(output_dir, fname + plot_extentions[i_exnt]), dpi=300, transparent=True)
        plt.close(fig)
    
#make averages
# making average plots based on conditions
files_dir_dff_avr = create_folder_with_note(files_dir_plots, 'bar plots')
outcome_cols = [ 'dff_photo', 'dff_iso', 'dff_MC','z_score_photo', 'z_score_iso', 'z_score_MC']
ID_cols  =  ['name', 'hemisphere', 'condition', 'week', 'batch', 'genotype', 'time_lock_name' ]

# Average
df_avg = aggregate_signals_in_time_windows(
    df=df_signals_avr,
    outcome_cols=outcome_cols,
    window=window,
    meta_cols=ID_cols,
    output_dir=files_dir_dff_avr,
    mode="average")

# Peak (max)
df_peak = aggregate_signals_in_time_windows(
    df=df_signals_avr,
    outcome_cols=outcome_cols,
    window=window,
    meta_cols=ID_cols,
    output_dir=files_dir_dff_avr,
    mode="peak")

# Dip (min)
df_dip = aggregate_signals_in_time_windows(
    df=df_signals_avr,
    outcome_cols=outcome_cols,
    window=window,
    meta_cols=ID_cols,
    output_dir=files_dir_dff_avr,
    mode="dip")

if ex_bar_plots == 1:

    
    #df_avg = pd.read_excel( files_dir_dff_avr  + "averaged_signals_windows_reward_baseline.xlsx")
    # df_peak = pd.read_excel( files_dir_dff_avr  + "peak_signals_windows_reward_baseline.xlsx")
    # df_dip = pd.read_excel( files_dir_dff_avr  + "dip_signals_windows_reward_baseline.xlsx")
    
    
    dfs_to_plot = [df_avg, df_peak,df_dip ]
    df_labels_add = ['avr','max',' min']
    data_type = ['avr', 'peak', 'dip']

    for i_df in range(len(dfs_to_plot)):
        for outcome in outcomes:
            if outcome == 'dff':
                df_label = ['ΔF/F',df_labels_add[i_df]]
            else:
                 df_label = ['Z', df_labels_add[i_df]]

            plot_per_timepoint_per_condition(
                df=dfs_to_plot[i_df],
                outcome=outcome,
                time_windows=time_windows,
                condition_col='condition',
                conditions=conditions,
                subject_col='name',
                colors=colors_to_conditions,
                x_label = xlabel_bar,
                ylabel=df_label,
                tl_to_xlabel = tl_to_xlabel_bar,
                output_dir=files_dir_dff_avr,paired_con=paired_con, y_lim = y_lim_bars, data_type = data_type[i_df], plot_extentions=file_extentions   )
        
            plot_all_timepoints_per_condition(
                df=dfs_to_plot[i_df],
                outcome=outcome,
                time_windows=time_windows,
                condition_col='condition',
                conditions=conditions,
                subject_col='name',
                colors=colors_to_conditions,
                ylabel=df_label,
                output_dir=files_dir_dff_avr,
                y_lim = y_lim_bars,tl_to_xlabel = tl_to_xlabel_bar, data_type = data_type[i_df], plot_extentions= file_extentions)
            
            plot_bar_over_time(
                df=dfs_to_plot[i_df],
                outcome=outcome,
                time_windows=time_windows,
                condition_col='condition',
                conditions=conditions,
                subject_col='name',
                colors=colors_to_conditions,
                xlabel = xlabel_bar,
                ylabel=df_label, 
                output_dir=files_dir_dff_avr,
                show_individuals=False, y_lim = y_lim_bars, tl_to_xlabel = tl_to_xlabel_bar,data_type = data_type[i_df], plot_extentions=file_extentions
            )


# %%
def stats_time_lock(df, outcome, time_windows, condition_col, conditions,
                    subject_col, output_dir,
                    paired_con='paired',
                    test_over_time=True,
                    test_between_conditions=True,
                    data_type='avr'):
    """
    Efficient and unified statistical testing:
    - Within-condition (paired over time)
    - Between-condition (paired/unpaired)
    - Automatically computes omnibus + post-hoc
    """

    import os, pandas as pd, numpy as np, scipy.stats as stats, pingouin as pg

    os.makedirs(output_dir, exist_ok=True)
    results_all = []

    time_locks = sorted(df["time_lock_name"].unique())

    # ------------------------
    # Helper: 2-group test
    # ------------------------
    def two_group_test(v1, v2, paired=True):
        if paired:
            diff = v1 - v2
            if stats.shapiro(diff).pvalue > 0.05:
                stat, p = stats.ttest_rel(v1, v2)
                method = "Paired t-test"
            else:
                stat, p = stats.wilcoxon(v1, v2)
                method = "Wilcoxon signed-rank test"
            df_val = len(v1) - 1
        else:
            normal = stats.shapiro(v1).pvalue > 0.05 and stats.shapiro(v2).pvalue > 0.05
            levene_p = stats.levene(v1, v2).pvalue
            if normal:
                stat, p = (stats.ttest_ind(v1, v2, equal_var=levene_p>0.05) 
                           if levene_p>0.05 else stats.ttest_ind(v1, v2, equal_var=False))
                method = "Independent t-test" if levene_p>0.05 else "Welch t-test"
            else:
                stat, p = stats.mannwhitneyu(v1, v2)
                method = "Mann–Whitney U"
            df_val = len(v1) + len(v2) - 2
        return stat, p, method, df_val

    # ------------------------
    # Helper: multi-group test
    # ------------------------
    def multi_group_test(data_matrix, labels, paired=True):
        normal = all(stats.shapiro(col).pvalue > 0.05 for col in data_matrix.T)

        if paired:
            df_long = pd.DataFrame(data_matrix, columns=labels)
            df_long[subject_col] = range(len(df_long))
            df_long = df_long.melt(id_vars=subject_col, var_name="Factor", value_name="Value")
            if normal:
                aov = pg.rm_anova(dv="Value", within="Factor", subject=subject_col,
                                  data=df_long, detailed=True, correction=True)
                method = "Repeated measures ANOVA"
                posthoc = pg.pairwise_tests(dv="Value", within="Factor", subject=subject_col,
                                            data=df_long, parametric=True, padjust='holm')
                posthoc_method = "Paired t-test"
            else:
                stat, p = stats.friedmanchisquare(*data_matrix.T)
                aov = pd.DataFrame({"F":[stat], "p-unc":[p]})
                method = "Friedman test"
                posthoc = pg.pairwise_tests(dv="Value", within="Factor", subject=subject_col,
                                            data=df_long, parametric=False, padjust='holm')
                posthoc_method = "Wilcoxon"
        else:
            df_long = pd.DataFrame({ "Factor": np.repeat(labels, len(data_matrix)),
                                     "Value": data_matrix.flatten() })
            levene_p = stats.levene(*[data_matrix[:,i] for i in range(data_matrix.shape[1])]).pvalue
            if normal:
                aov = pg.anova(dv="Value", between="Factor", data=df_long, detailed=True) if levene_p>0.05 \
                      else pg.welch_anova(dv="Value", between="Factor", data=df_long, detailed=True)
                method = "One-way ANOVA" if levene_p>0.05 else "Welch ANOVA"
                posthoc = pg.pairwise_tests(dv="Value", between="Factor", data=df_long, parametric=True, padjust='holm')
                posthoc_method = "Independent t-test"
            else:
                stat, p = stats.kruskal(*[data_matrix[:,i] for i in range(data_matrix.shape[1])])
                aov = pd.DataFrame({"F":[stat], "p-unc":[p]})
                method = "Kruskal-Wallis"
                posthoc = pg.pairwise_tests(dv="Value", between="Factor", data=df_long, parametric=False, padjust='holm')
                posthoc_method = "Mann–Whitney U"

        stat = aov['F'][0]
        p = aov['p-unc'][0]
        df_val = f"{int(aov['ddof1'][0])}, {int(aov['ddof2'][0])}" if 'ddof1' in aov.columns else np.nan
        return stat, p, method, df_val, df_long, posthoc, posthoc_method

    # ------------------------
    # Helper: append posthoc
    # ------------------------
    def add_posthoc(posthoc, tl, test_name, subtest, means, method_name):
        for _, r in posthoc.iterrows():
            stat_val = r.get('T', r.get('U', r.get('W', np.nan)))
            results_all.append({
                "Time_lock": tl,
                "Condition": f"{r['A']} vs {r['B']}",
                "Test": test_name + " (posthoc)",
                "Subtest": subtest,
                "Method": method_name,
                "Mean_1": means.get(r['A'], np.nan),
                "Mean_2": means.get(r['B'], np.nan),
                "Mean_3": np.nan,
                "df": r.get('dof', np.nan),
                "Statistic": stat_val,
                "p_value": r.get('p-corr', np.nan)
            })

    # ------------------------
    # Main loop
    # ------------------------
    for tl in time_locks:
        df_tl = df[df["time_lock_name"] == tl]
        if df_tl.empty: continue

        # --- WITHIN CONDITION ---
        if test_over_time:
            for cond in conditions:
                cond_df = df_tl[df_tl[condition_col]==cond]
                if cond_df.empty: continue
                subj_data = np.array([[
                    rows[make_window_colname(outcome, t0, t1)].mean() for t0,t1 in time_windows
                ] for (_, hemi), rows in cond_df.groupby([subject_col,'hemisphere'])])

                labels = [f"{t0}-{t1}s" for t0,t1 in time_windows]
                if subj_data.shape[1]==2:
                    stat, p, method, df_val = two_group_test(subj_data[:,0], subj_data[:,1], paired=True)
                    results_all.append({
                        "Time_lock": tl, "Condition": cond, "Test":"Within condition",
                        "Subtest": f"{time_windows[0][0]} - {time_windows[0][1]} vs {time_windows[1][0]} - {time_windows[1][1]}s",
                        "Method": method, "Statistic": stat, "p_value": p,
                        "Mean_1": subj_data[:,0].mean(), "Mean_2": subj_data[:,1].mean(),
                        "Mean_3": np.nan, "df": df_val
                    })
                else:
                    stat, p, method, df_val, _, posthoc, posthoc_method = multi_group_test(subj_data, labels, paired=True)
                    results_all.append({
                        "Time_lock": tl, "Condition": cond, "Test":"Within condition",
                        "Subtest": "Time windows", "Method": method, "Statistic": stat, "p_value": p,
                        "Mean_1": subj_data[:,0].mean(), "Mean_2": subj_data[:,1].mean(),
                        "Mean_3": subj_data[:,2].mean() if subj_data.shape[1]>2 else np.nan, "df": df_val
                    })
                    means = {label: subj_data[:,i].mean() for i,label in enumerate(labels)}
                    if p < 0.05:
                        add_posthoc(posthoc, tl, "Within condition", "Time windows", means, posthoc_method)

        # --- BETWEEN CONDITIONS ---
        if test_between_conditions and len(conditions)>1:
            for t0,t1 in time_windows:
                col = make_window_colname(outcome, t0, t1)
                vals = [df_tl[df_tl[condition_col]==cond].groupby([subject_col,'hemisphere'])[col].mean().values
                        for cond in conditions]
                subj_ids = [df_tl[df_tl[condition_col]==cond].groupby([subject_col,'hemisphere'])[col].mean().index
                           for cond in conditions]
                paired = (paired_con=='paired')
                if paired:
                    common = set(subj_ids[0]).intersection(*subj_ids[1:])
                    if not common: continue
                    vals = [v[[i in common for i in idx]] for v,idx in zip(vals, subj_ids)]

                labels_cond = conditions
                if len(vals)==2:
                    stat, p, method, df_val = two_group_test(vals[0], vals[1], paired=paired)
                    results_all.append({
                        "Time_lock": tl, "Condition":" vs ".join(conditions),
                        "Test":"Between conditions", "Subtest": f"{t0}-{t1}s",
                        "Method": method, "Statistic": stat, "p_value": p,
                        "Mean_1": np.mean(vals[0]), "Mean_2": np.mean(vals[1]),
                        "Mean_3": np.nan, "df": df_val
                    })
                else:
                    data_matrix = np.array(vals).T
                    stat, p, method, df_val, _, posthoc, posthoc_method = multi_group_test(data_matrix, labels_cond, paired=paired)
                    results_all.append({
                        "Time_lock": tl, "Condition":" vs ".join(conditions),
                        "Test":"Between conditions", "Subtest": f"{t0}-{t1}s",
                        "Method": method, "Statistic": stat, "p_value": p,
                        "Mean_1": np.mean(vals[0]), "Mean_2": np.mean(vals[1]),
                        "Mean_3": np.mean(vals[2]) if len(vals)>2 else np.nan, "df": df_val
                    })
                    means = {cond: np.mean(vals[i]) for i,cond in enumerate(conditions)}
                    if p < 0.05:
                        add_posthoc(posthoc, tl, "Between conditions", f"{t0}-{t1}s", means, posthoc_method)

    # ------------------------
    # Save
    # ------------------------
    df_res = pd.DataFrame(results_all)
    for col in ['Mean_1','Mean_2','Mean_3','Statistic','p_value']:
        if col in df_res.columns: df_res[col] = df_res[col].round(3)
    cols_order = ['Time_lock','Condition','Test','Subtest','Method','Mean_1','Mean_2','Mean_3','df','Statistic','p_value']
    df_res = df_res[cols_order]
    fname = f"stats_results_{data_type}.xlsx"
    df_res.to_excel(os.path.join(output_dir,fname), index=False)
    print(f"Stats results saved to {os.path.join(output_dir,fname)}")
    return df_res

if ex_statistics == 1:
    df_avg = df_avg.dropna()
    df_peak = df_peak.dropna()
    df_dip = df_dip.dropna()
    dfs_to_plot = [df_avg, df_peak,df_dip ]
    data_type = ['avr', 'peak', 'dip']

    for i_df in range(len(dfs_to_plot)):            
        df_results = stats_time_lock(
            df=dfs_to_plot[i_df],   # <-- important fix (was df_avg before)
            outcome="z_score_MC",
            time_windows=time_windows,
            condition_col="condition",
            conditions=conditions,
            subject_col="name",
            output_dir=files_dir_plots,
            paired_con=paired_con,   # <-- NEW (key argument)
            test_over_time=True,
            test_between_conditions=True,
            data_type=data_type[i_df] )
            
