#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 23:23:40 2026

@author: gianni
"""

from batch_simulations.utils.get_calibrator_candidates_wrapper import\
        GetCalibratorCandidatesWrapper,CalibratorCandidate
import pytest
from pathlib import Path
from dataclasses import fields
import numpy as np
import datetime


class DummyProcess:
    def __init__(self, stdout):
        self.stdout = stdout


class TestGetCalibratorCandidatesWrapper:

    test_output_folder = Path("tests/getCalibratorCandidates_outputs")
    mock_xml_filepath = "something.xml"
    general_get_cal_candidates = GetCalibratorCandidatesWrapper(
                                             xml_filepath=mock_xml_filepath)

    def test_get_calibrator_candidates(self):
        int_time = 2.3
        array_config = "c43-9"
        min_cmd = self.general_get_cal_candidates.construct_command(
                 integration_time=int_time, array_config=array_config)
        assert min_cmd == ["getCalibratorCandidates.py", self.mock_xml_filepath,"-t",
                           str(int_time),"-C", array_config]
        search_radius = 40.3
        epoch = "TRANSIT"
        calibrator_type = "phase"
        spectral_spec = "something"
        src = "some source"
        no_spwavg = True
        max_cmd = self.general_get_cal_candidates.construct_command(
                      integration_time=int_time, array_config=array_config,
                      search_radius=search_radius,epoch=epoch,calibrator_type=calibrator_type,
                      spectral_spec=spectral_spec,src=src,no_spwavg=no_spwavg)
        assert max_cmd == ["getCalibratorCandidates.py", self.mock_xml_filepath, "-t", str(int_time),
                           "-C", array_config, "-r", str(search_radius), "-e",
                           epoch, "-c", calibrator_type,f"--spectralSpec={spectral_spec}",
                           f"--src={src}","--no_spwavg"]
        with pytest.raises(ValueError):
            self.general_get_cal_candidates.construct_command(
                     integration_time=int_time, array_config=array_config, calibrator_type="bandpass",
                     no_spwavg=True)

    def get_calibrator_candidates(self,filename):
        filepath = self.test_output_folder / filename
        with open(filepath, 'r') as f:
            process = DummyProcess(stdout=f.read())
        return self.general_get_cal_candidates.read_calibrator_candidates(process=process)

    @staticmethod
    def assert_calibrator_equality(calibrator,expected_calibrator):
        for field in fields(expected_calibrator):
            expected_value = getattr(expected_calibrator,field.name)
            test_value = getattr(calibrator, field.name)
            if isinstance(test_value, (int, float)):
                if np.isnan(expected_value):
                    assert np.isnan(test_value)
                    continue
            assert test_value == expected_value

    def test_calibrator_candidates_results(self):
        calibrator_candidates = self.get_calibrator_candidates(
                     filename="getCalibratorCandidates_output_general_output.txt")
        assert len(calibrator_candidates) == 26
        #some random checks:
        assert calibrator_candidates[0].source_name == "J2253+1608"
        assert calibrator_candidates[0].Reason == ""
        assert np.isnan(calibrator_candidates[0].UVmin)
        assert calibrator_candidates[1].Reason == "TW"
        assert calibrator_candidates[1].EstimatedFlux == 0.847
        assert calibrator_candidates[1].EstimatedFluxError == 0.028
        assert calibrator_candidates[1].SNR == 6.5
        assert calibrator_candidates[1].SNRError == 0.2
        assert calibrator_candidates[0].isObservable
        assert not calibrator_candidates[3].isObservable
        #test all fields of one of the calibrators:
        test_calibrator = calibrator_candidates[-2]
        expected_calibrator = CalibratorCandidate(
                                source_name="J2025+3343",
                                type_IDs="1,4,25,48,64",
                                Az=318,
                                El=15.4,
                                eRa=0.15,
                                eDec=0.12,
                                specIndex=-0.7,
                                specIndex_error=0.15,
                                reduced_chi2=-1,
                                Nobs=17,
                                LastDate=datetime.date(year=2026,month=1,day=22),
                                EstimatedFlux=0.129,
                                EstimatedFluxError=0.06,
                                SNR=0.2,
                                SNRError=0.1,
                                Sep=85.8,
                                isObservable=False,
                                fShadow=0.4,
                                fShadow_noncritical=0,
                                fRes=0,
                                dDays=29,
                                UVmax=np.nan,
                                UVmin=np.nan,
                                Score=-1,
                                Reason="TW,NV,S,L")
        #in principle I could just do test_calibrator==expected_calibrator,
        #but because of NaNs it does not work, so I need a workaround:
        self.assert_calibrator_equality(calibrator=test_calibrator,
                                        expected_calibrator=expected_calibrator)

    def test_empty_output(self):
        calibrator_candidates = self.get_calibrator_candidates(
             filename="getCalibratorCandidates_output_empty_candidate_list.txt")
        assert len(calibrator_candidates) == 0

    def test_single_candidate(self):
        calibrator_candidates = self.get_calibrator_candidates(
                   filename="getCalibratorCandidates_output_single_candidate.txt")
        assert len(calibrator_candidates) == 1
        calibrator_candidate = calibrator_candidates[0]
        expected_calibrator = CalibratorCandidate(
                                source_name="J0006-0623",
                                type_IDs="1,24,25,49,64",
                                Az=15,
                                El=72.8,
                                eRa=0.1,
                                eDec=0.1,
                                specIndex=-0.45,
                                specIndex_error=0.01,
                                reduced_chi2=3.19,
                                Nobs=153,
                                LastDate=datetime.date(year=2026,month=2,day=19),
                                EstimatedFlux=0.847,
                                EstimatedFluxError=0.028,
                                SNR=3.4,
                                SNRError=0.1,
                                Sep=21.3,
                                isObservable=True,
                                fShadow=0,
                                fShadow_noncritical=0,
                                fRes=0,
                                dDays=1,
                                UVmax=-1831,
                                UVmin=np.nan,
                                Score=-1,
                                Reason="TW")
        self.assert_calibrator_equality(expected_calibrator=expected_calibrator,
                                        calibrator=calibrator_candidate)