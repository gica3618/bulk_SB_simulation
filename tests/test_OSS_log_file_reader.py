#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:55:51 2026

@author: gianni
"""

from batch_simulations.infrastructure.OSS_log_file_reader import LogFileReader
from pathlib import Path
import pytest


class TestLogFileReader:

    def test_extract_value(self):
        line = ("2026-04-28T03:01:30 [ArrayOSS/CalibratorCatalog]   -> Listing "
                +"ranked candidate list... [check ra=  69.521 dec=  30.079 "
                +"radius=10.0 SNR_min=15.0")
        search_radius = LogFileReader.extract_value(line=line,name="radius")
        assert search_radius == 10
        min_SNR = LogFileReader.extract_value(line=line,name="SNR_min")
        assert min_SNR == 15

    def test_read_queried_calibrator(self):
        line = ("2026-04-28T03:01:30 [ArrayOSS/CalibratorCatalog] |[J0429+2724] 1 "
                +"             |   2.4|  39.5| 0.11| 0.10|-0.70+- 0.15|-1.00|  26|"
                +"20260311| 0.145+- 0.012|   73.5+-   5.9|  2.2| True|0.12(0.12)|"
                +"  0.0|   47|-16978.5|    nan| 0.51|         |")
        candidate = LogFileReader.read_queried_calibrator(line)
        assert candidate.source_name == "J0429+2724"
        assert candidate.SNR == 73.5
        assert candidate.SNRError == 5.9
        assert candidate.isObservable
        assert candidate.Reason == ""

    @staticmethod
    def get_log_filepath(filename):
        return Path("tests/simulateSB_outputs/OSS_log_files") / filename

    def test_read_calibrator_query(self):
        reader = LogFileReader(filepath=self.get_log_filepath(filename="HA0h_log_IRAS0416_a_09_TM2.xml_OSS.txt"))
        radius,min_SNR,candidates = reader.read_calibrator_query(calibrator="DGC")
        assert radius == 40
        assert min_SNR == 20
        assert len(candidates) == 5
        assert candidates[0].source_name == "J0423-0120"
        assert candidates[1].SNR == 28.5
        assert candidates[2].Reason == "NV,S"
        assert candidates[3].Reason == "NV,S"
        assert candidates[4].Reason == "TW"
        radius,min_SNR,candidates = reader.read_calibrator_query(calibrator="Check")
        assert radius == 10
        assert min_SNR == 15
        assert candidates[0].SNRError == 3.6
        assert candidates[1].Reason == "TW"
        #read DGC on file without diffgain query:
        reader = LogFileReader(filepath=self.get_log_filepath(filename="HA0h_log_Mouse_a_04_7M.xml_OSS.txt"))
        with pytest.raises(ValueError):
            reader.read_calibrator_query(calibrator="DGC")
        radius,min_SNR,candidates = reader.read_calibrator_query(calibrator="Bandpass")
        assert radius == 120
        assert min_SNR == 50
        assert candidates[0].SNR == 364.8
        assert candidates[-1].Reason == "TW,NV,S"
        radius,min_SNR,candidates = reader.read_calibrator_query(calibrator="Phase")
        assert radius == 10
        assert min_SNR == 15
        assert candidates[1].source_name == "J1717-3342"
        assert candidates[8].Reason == "F"