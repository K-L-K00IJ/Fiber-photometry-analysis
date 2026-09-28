# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 15:52:29 2026

Merge 2 or more parquet files per mouse/session
and plot MC signal.
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
plt.ioff()

# directories
files_dir1 ='D://GRAB DA//shock//Batch 1//1. Preprocessed//149 split//'
files_dir_fig = files_dir1 + '//Figs//'  # <-- adjust if needed

os.makedirs(files_dir_fig, exist_ok=True)

# file names grouped by mouse
file_name = [
    ['1349_0_0_0_0_0_first', '1349_0_0_0_0_0_sec']]

# loop through mice
for i_mouse in range(len(file_name)):

    dfs = []  # collect parts here

    for part in file_name[i_mouse]:
        file_path = os.path.join(files_dir1, part + '.parquet')
        df_part = pd.read_parquet(file_path)
        dfs.append(df_part)

    # merge (row-wise)
    df_merged = pd.concat(dfs, axis=0, ignore_index=True)

    # optional: sort if you have a time column
    # df_merged = df_merged.sort_values('time').reset_index(drop=True)

    # create ID from first filename (before "_first")
    ID = file_name[i_mouse][0].split('_first')[0]

    # save merged parquet
    save_path = os.path.join(files_dir1, ID + '_merged.parquet')
    df_merged.to_parquet(save_path)
    print(f"Saved merged file: {save_path}")

    # plot MC channel
    plt.figure(figsize=(10, 4))
    plt.plot(df_merged['MC'], linewidth=1)
    plt.title(f'{ID} – MC signal')
    plt.xlabel('Sample')
    plt.ylabel('MC')
    plt.tight_layout()

    fig_path = os.path.join(files_dir_fig, ID + '_MC.png')
    plt.savefig(fig_path, transparent=True, dpi=300)
    plt.close()

    print(f"Saved MC plot: {fig_path}")
        