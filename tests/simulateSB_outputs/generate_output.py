#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 17:55:24 2026

@author: gianni
"""

import subprocess
import pickle
import os
import glob

commands = {"success":"simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
            "failed":"simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT-4H,2026-02-20 -C c43-3"}

for ID,command in commands.items():
    output = subprocess.run(command,shell=True,universal_newlines=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    summary_file = glob.glob("*OSS_summary.txt")
    output_folder = ID
    os.mkdir(output_folder)
    filepath = os.path.join(output_folder,f"simulateSB_output_{ID}.pkl")
    with open(filepath, 'wb') as f: 
        pickle.dump(obj=output, file=f)
    if len(summary_file) > 0:
        assert len(summary_file) == 1
        summary_file = summary_file[0]
        os.rename(summary_file,os.path.join(output_folder,summary_file))