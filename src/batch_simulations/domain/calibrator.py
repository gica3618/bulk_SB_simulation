#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:00:06 2024

@author: gianni
"""

import logging


CALIBRATOR_TYPES = {
    "Polarization",
    "Bandpass",
    "Phase",
    "Check",
    "Amplitude",
    "DGC"}
calibrator_keywords = {cal:[cal,cal.lower()] for cal in CALIBRATOR_TYPES}

def classify_calibrator(name):
    candidate_cal_types = []
    for cal,keywords in calibrator_keywords.items():
        #I don't want to read "Bandpass pointing" etc.
        #these appear sometimes, e.g. Circinus_a_04_TM1 of project 2025.1.00238.S
        if any([keyword in name for keyword in keywords]) and "pointing" not in name.lower():
            candidate_cal_types.append(cal)
    if len(candidate_cal_types) == 0:
        logging.info(f"Source '{name}' is not a calibrator")
        return None
    if len(candidate_cal_types) > 1:
        raise ValueError("unable to uniquely determine calibrator type "
                         +f"(candidate_cal_types: {candidate_cal_types})")
    return candidate_cal_types[0]


class Calibrator:

    def __init__(self,name,source_name,cal_type,is_hardcoded,coordinates=None):
        self.name = name #this is something like "Phase", or "Polarization calibrator"
        self.source_name = source_name #this can be "query", or something like "J1326-5256"
        self.cal_type = cal_type
        self.is_hardcoded = is_hardcoded
        self.coordinates = coordinates
        self.check_consistency()

    def check_consistency(self):
        if self.cal_type not in CALIBRATOR_TYPES:
            raise ValueError(f"Invalid calibrator type '{self.cal_type}'. "
                             +f"Must be one of {CALIBRATOR_TYPES}.")
        if self.is_hardcoded and self.coordinates is None:
            raise ValueError("Hardcoded calibrator must have coordinates.")

    def is_DGC(self):
        return self.cal_type == "DGC"