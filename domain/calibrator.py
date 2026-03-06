#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:00:06 2024

@author: gianni
"""

import subprocess
import logging
import datetime
from dataclasses import dataclass
from typing import Optional


CALIBRATOR_TYPES = {
    "Polarization",
    "Bandpass",
    "Phase",
    "Check",
    "Amplitude",
    "DGC",
}
calibrator_keywords = {cal:[cal,cal.lower()] for cal in CALIBRATOR_TYPES}

def classify_calibrator(name):
    candidate_cal_types = []
    for cal,keywords in calibrator_keywords.items():
        if any([keyword in name for keyword in keywords]):
            candidate_cal_types.append(cal)
    if len(candidate_cal_types) == 0:
        logging.info(f"Source '{name}' is not a calibrator")
        return None
    if len(candidate_cal_types) > 1:
        raise ValueError("unable to uniquely determine calibrator type "
                         +f"(candidate_cal_types: {candidate_cal_types})")
    return candidate_cal_types[0]

@dataclass(frozen=True)
class Calibrator:
    name: str # this is something like "Phase", or "Polarization calibrator"
    source_name: str #this can be "query", or something like "J1326-5256"
    cal_type: str
    is_hardcoded: bool
    coordinates: Optional[object]

    def __post_init__(self):
        if self.cal_type not in CALIBRATOR_TYPES:
            raise ValueError(
                f"Invalid calibrator type '{self.cal_type}'. "
                f"Must be one of {CALIBRATOR_TYPES}."
            )

        if self.is_hardcoded and self.coordinates is None:
            raise ValueError(
                "Hardcoded calibrator must have coordinates."
            )

        if not self.is_hardcoded and self.coordinates is not None:
            raise ValueError(
                "Query calibrator must not have coordinates."
            )