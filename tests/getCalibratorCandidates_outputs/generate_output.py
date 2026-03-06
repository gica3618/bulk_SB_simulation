#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 14:58:14 2026

@author: gianni
"""

import subprocess

commands = {"general_output":"getCalibratorCandidates.py 2025.1.00378.S_SchedBlock0.xml -e TRANSIT-1h,2026-02-20 -c bandpass -t 1800 -r 120 -C c43-3",
            "empty_candidate_list":"getCalibratorCandidates.py 2025.1.00378.S_SchedBlock0.xml -e TRANSIT-1h,2026-02-20 -c bandpass -t 1800 -r 1 -C c43-3",
            "single_candidate":"getCalibratorCandidates.py 2025.1.00378.S_SchedBlock0.xml -e TRANSIT-1h,2026-02-20 -c bandpass -t 1800 -r 23 -C c43-3"}

for ID,command in commands.items():
    output = subprocess.run(command,shell=True,universal_newlines=True,
                            stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    filepath = f"getCalibratorCandidates_output_{ID}.txt"
    with open(filepath, 'w') as f: 
        f.write(output.stdout)
