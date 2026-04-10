#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 23:57:17 2026

@author: gianni
"""

from batch_simulations.domain.simulation_result import SimulationResult
import pickle
from pathlib import Path
import pytest


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

    def test_server_error(self):
        assert SimulationResult.server_error(None) == False
        for fail_reason in ("bla unexpected response from the source catalogue bla",
                            "bla socket.gaierror bla"):
            assert SimulationResult.server_error(fail_reason=fail_reason)
        for fail_reason in ("no server error here","bla"):
            assert not SimulationResult.server_error(fail_reason=fail_reason)

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
        assert result.fail_reason == 'summary file does not report success'

    def test_success(self):
        process = self.get_process(ID="success")
        result = SimulationResult.from_completed_process(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
                       process=process,
                       output_folder=self.output_folder/"success",
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert result.success
        assert not result.server_error
        assert result.fail_reason is None

    def test_failed(self):
        process = self.get_process(ID="failed")
        result = SimulationResult.from_completed_process(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT-4H,2026-02-20 -C c43-3",
                       process=process,
                       output_folder=self.output_folder/"failed",
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert not result.success
        assert not result.server_error
        assert result.fail_reason == "no check"

    @staticmethod
    def test_shorten_fail_reason():
        fail_reason = "Exception: Although 1 source(s) requested, got only 0 for phase query."
        assert SimulationResult.shorten_fail_reason(fail_reason) == "no phase"
        fail_reason = "Exception: Refusing the SB execution as it will "\
                                 +"exceed the limit (2.00 hours) by 150 s"
        assert SimulationResult.shorten_fail_reason(fail_reason) == "SB exceeds 2h limit by 150 s"
        fail_reason = "blabla"
        assert SimulationResult.shorten_fail_reason(fail_reason) == fail_reason

    def test_get_fail_reason_preference(self):
        #check that preference is given to stderr
        pipe = {"stdout":"message\nerror","stderr":"othermessage\nexception"}
        assert SimulationResult.get_fail_reason(pipe) == "exception"

    def test_get_fail_reason(self):
        msg = "balerrorblabla"
        pipe =  {"stdout":f"{msg}\nmessage","stderr":"othermessage\nhallo\nrassel"}
        assert SimulationResult.get_fail_reason(pipe) == msg

    def test_get_fail_reason_max_n_messages(self):
        for n_messages in (SimulationResult.max_messages_to_go_back-1,
                           SimulationResult.max_messages_to_go_back,
                           SimulationResult.max_messages_to_go_back+1):
            pipe = {"stdout":"\n".join(["a"]*n_messages),
                    "stderr":"\n".join(["b"]*n_messages)}
            assert SimulationResult.get_fail_reason(pipe) == "a"

    def test_summary_file_reports_success(self):
        success_summary_file = self.output_folder\
                                            / "log_2025.1.01422.S_LaPequen_e_10_TM1.xml_OSS_summary.txt"
        assert SimulationResult.summary_file_reports_success(success_summary_file)
        no_success_summary_file = self.output_folder\
                                            / "no_success_log_2025.1.01422.S_LaPequen_e_10_TM1.xml_OSS_summary.txt"
        assert not SimulationResult.summary_file_reports_success(no_success_summary_file)