#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 23:57:17 2026

@author: gianni
"""

from batch_simulations.domain.simulation_result import FailReason,SimulationResult
import pickle
from pathlib import Path
import pytest


class TestFailReason:

    def test_invalid_category(self):
        with pytest.raises(ValueError):
            FailReason(error_message="bla",error_summary="blabla",
                       category="some invalid category")

    def test_constructor_from_error_message(self):
        tests = [{"error_message":"Exception: Although 1 source(s) requested, got only 0 for phase query.",
                  "expected_error_summary":"no phase",
                  "expected_category":"missing calibrator"},
                 {"error_message":"Exception: Refusing the SB execution as it will exceed the limit (2.00 hours) by 150 s",
                  "expected_error_summary":"SB exceeds 2h limit by 150 s",
                  "expected_category":"exceeds 2h limit"},
                 {"error_message":"Observation.SBExecutionMode.SBExecutionError: No visible science target [Group2]",
                  "expected_error_summary":"No visible science target [Group2]",
                  "expected_category":"unobservable"},
                 {"error_message":"Observation.SBExecutionMode.SBExecutionError: All science targets in group2 are unobservable",
                  "expected_error_summary":"All science targets in group2 are unobservable",
                  "expected_category":"unobservable"},
                 {"error_message":"Observation.SBExecutionMode.SBExecutionError: Although execution of 'BandpassCalTarget on J2253+1608 at 903.162 GHz ref=topo' is required, it is not observable. And also failed to find out another viable calibrator.",
                  "expected_error_summary":"BandpassCalTarget J2253+1608 not observable",
                  "expected_category":"unobservable"},
                 {"error_message":"Observation.SBExecutionMode.SBExecutionError: Asterix und Obelix",
                  "expected_error_summary":"Observation.SBExecutionMode.SBExecutionError: Asterix und Obelix",
                  "expected_category":"unobservable"},
                 {"error_message":"unexpected response from the source catalogue",
                  "expected_error_summary":"unexpected response from the source catalogue",
                  "expected_category":"server error"},
                 {"error_message":"socket.gaierror",
                  "expected_error_summary":"socket.gaierror",
                  "expected_category":"server error"},
                 {"error_message":"bla bla some error",
                  "expected_error_summary":"bla bla some error",
                  "expected_category":"other"},
                 {"error_message":"Exception: Specified elevation [-15.212032] is out of range.",
                  "expected_error_summary":"elevation out of range",
                  "expected_category":"unobservable"}
                 ]
        for test in tests:
            fail_reason = FailReason.from_error_message(error_message=test["error_message"])
            assert fail_reason.error_summary == test["expected_error_summary"]
            assert fail_reason.category == test["expected_category"]


class TestSimulationResult:

    output_folder = Path("tests/simulateSB_outputs")

    def get_process(self,ID):
        filepath = self.output_folder / ID / f"simulateSB_output_{ID}.pkl"
        with open(filepath, "rb") as f:
            return pickle.load(f)

    def test_wrong_xml(self):
        process = self.get_process(ID="success")
        with pytest.raises(ValueError):
            SimulationResult.from_completed_process(
               executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
               process=process,output_folder=self.output_folder/"success",
               xml_filename="abc.xml")

    def test_fail_from_summary(self):
        filepath = self.output_folder / "success_but_summaryfile_reports_failed"\
                               / "simulateSB_output_success.pkl"
        with open(filepath, "rb") as f:
            process = pickle.load(f)
        result = SimulationResult.from_completed_process(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
                       process=process,
                       output_folder=self.output_folder/"success_but_summaryfile_reports_failed",
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert not result.success
        assert result.fail_reason.error_summary == 'summary file does not report success'
        assert result.fail_reason.category == "other"

    def test_success(self):
        process = self.get_process(ID="success")
        result = SimulationResult.from_completed_process(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
                       process=process,
                       output_folder=self.output_folder/"success",
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert result.success
        assert result.fail_reason is None

    def test_failed(self):
        process = self.get_process(ID="failed")
        result = SimulationResult.from_completed_process(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT-4H,2026-02-20 -C c43-3",
                       process=process,
                       output_folder=self.output_folder/"failed",
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert not result.success
        assert result.fail_reason.error_summary == "no check"
        assert result.fail_reason.category == "missing calibrator"

    def test_get_error_message_preference(self):
        #check that preference is given to stderr
        pipe = {"stdout":"message\nerror","stderr":"othermessage\nexception"}
        assert SimulationResult.get_error_message(pipe) == "exception"

    def test_get_error_message(self):
        msg = "balerrorblabla"
        pipe =  {"stdout":f"{msg}\nmessage","stderr":"othermessage\nhallo\nrassel"}
        assert SimulationResult.get_error_message(pipe) == msg

    def test_get_get_error_message_max_n_messages(self):
        for n_messages in (SimulationResult.max_messages_to_go_back-1,
                           SimulationResult.max_messages_to_go_back,
                           SimulationResult.max_messages_to_go_back+1):
            pipe = {"stdout":"\n".join(["a"]*n_messages),
                    "stderr":"\n".join(["b"]*n_messages)}
            assert SimulationResult.get_error_message(pipe) == "a"

    def test_summary_file_reports_success(self):
        success_summary_file = self.output_folder\
                                            / "log_2025.1.01422.S_LaPequen_e_10_TM1.xml_OSS_summary.txt"
        assert SimulationResult.summary_file_reports_success(success_summary_file)
        no_success_summary_file = self.output_folder\
                                            / "no_success_log_2025.1.01422.S_LaPequen_e_10_TM1.xml_OSS_summary.txt"
        assert not SimulationResult.summary_file_reports_success(no_success_summary_file)