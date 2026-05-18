#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 11:53:10 2026

@author: gianni
"""

#TODO only read calibrator query if calibrator is not hardcoded
from batch_simulations.infrastructure.calibrator_query import parse_calibrator_candidate_data


class LogFileReader:

    calibrator_keys = {"Bandpass":"bandpass ra",
                       "Phase":"phase ra",
                       "Check":"check ra",
                       "DGC":"diffgain ra"}

    def __init__(self,filepath):
        with open(filepath,"r") as f:
            self.lines = f.readlines()

    @staticmethod
    def extract_value(line,name):
        line_split = line.split(" ")
        value_str = [s for s in line_split if f"{name}=" in s]
        if len(value_str) != 1:
            raise RuntimeError
        value_str = value_str[0]
        value_split = value_str.split("=")
        if len(value_split) != 2:
            raise RuntimeError
        return float(value_split[-1])

    @staticmethod
    def read_queried_calibrator(line):
        #2026-04-28T03:01:30 [ArrayOSS/CalibratorCatalog] |[J0429+2724] 1              |   2.4|  39.5| 0.11| 0.10|-0.70+- 0.15|-1.00|  26|20260311| 0.052+- 0.009|    6.3+-   1.1|  3.2| True|0.12(0.12)|  0.0|   47|-16978.5|    nan|-1.00|       TW|
        trimmed_line = line.split("[ArrayOSS/CalibratorCatalog] ")[1]
        trimmed_line = trimmed_line.rstrip()
        return parse_calibrator_candidate_data(trimmed_line)

    def read_calibrator_query(self,calibrator):
        key = self.calibrator_keys[calibrator]
        index = [i for i, line in enumerate(self.lines) if key in line]
        if len(index) == 0:
            raise ValueError(f"could not find calibrator query for {calibrator}")
        if len(index) > 1:
            raise RuntimeError(f"unable to uniquely read calibrator query for {calibrator}")
        index = index[0]
        search_radius = self.extract_value(line=self.lines[index], name="radius")
        min_SNR = self.extract_value(line=self.lines[index], name="SNR_min")
        while True:
            index += 1
            if "SourceName (and type IDs)" in self.lines[index]:
                break
        calibrator_index = index+2
        calibrator_candidates = []
        while "|" in self.lines[calibrator_index]:
            candidate = self.read_queried_calibrator(line=self.lines[calibrator_index])
            calibrator_candidates.append(candidate)
            calibrator_index += 1
        return search_radius,min_SNR,calibrator_candidates


class CandidateAnalysis:

    def __init__(self):
        raise NotImplementedError()