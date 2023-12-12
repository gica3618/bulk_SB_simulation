#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:44:13 2023

@author: gianni
"""

import os
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


log_folder = '/home/gianni/ARC/P2G/bulk_simulation_script/log_files_2023-10-23_05:01:08'
log_filenames = os.listdir(log_folder)

datetime_format = '%Y-%m-%dT%H:%M:%S'
output_file = 'IndexError_timestamps.txt'

start_times = []
error_occured = []
IndexError_times = []
for filename in log_filenames:
    filepath = os.path.join(log_folder,filename)
    with open(filepath,'r') as f:
        lines = f.readlines()
        start_time_str = lines[0].split(' ')[0]
        start_time = datetime.datetime.strptime(start_time_str,datetime_format)
        start_times.append(start_time)
        IndexError_line_numbers = [i for i,line in enumerate(lines) if
                                   'Exception caught: list index out of range' in
                                   line]
        n_candidates = len(IndexError_line_numbers)
        assert n_candidates <= 1
        if n_candidates == 1:
            error_occured.append(1)
            index = IndexError_line_numbers[0]
            error_time_str = lines[index].split(' ')[0]
            IndexError_times.append(datetime.datetime.strptime(error_time_str,
                                                               datetime_format))
        else:
            error_occured.append(0)

starttimes_sort_indices = sorted(range(len(start_times)),key=lambda i: start_times[i])
start_times = [start_times[i] for i in starttimes_sort_indices]
error_occured = [error_occured[i] for i in starttimes_sort_indices]
IndexError_times = sorted(IndexError_times)


IndexError_times = sorted(IndexError_times)


with open(output_file,'w') as f:
    output_time_strings = [t.strftime(datetime_format) for t in IndexError_times]
    output_str = '\n'.join(output_time_strings)
    f.write(output_str)

fig,ax = plt.subplots()
ax.plot(start_times,error_occured)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m %H:%M'))
for label in ax.get_xticklabels(which='major'):
    label.set(rotation=30, horizontalalignment='right')
ylabels = [item.get_text() for item in ax.get_yticklabels()]
def convert_ylabels(label):
    if label == '1.0':
        return 'Failed with\nIndexError'
    elif label == '0.0':
        return 'No IndexError'
    else:
        return ''        
ylabels = [convert_ylabels(yl) for yl in ylabels]
ax.set_yticklabels(ylabels)
fig.tight_layout()