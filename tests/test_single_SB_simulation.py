#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:37:57 2026

@author: gianni
"""

from unittest.mock import patch
import datetime
from pathlib import Path
from types import SimpleNamespace
from astropy.coordinates import Angle
from astropy import units as u
from batch_simulations.application.single_SB_simulation import SingleSBSimulation,\
    AnalysisResult, SingleSBSimulationSummary
from batch_simulations.domain.simulation_result import FailReason
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML



class TestSingleSBSimulation:

    general_sb_simulation = SingleSBSimulation(
                               project_code="2025.1.01422.S",sb_name="LaPequen_a_10_TM1",
                               date=datetime.date(2024,1,2),array_config_12m="c43-4")

    def test_build_xml_filepath(self):
        expected_filename = "2025.1.01422.S_LaPequen_a_10_TM1.xml"
        assert self.general_sb_simulation.build_xml_filepath()\
              == Path.cwd()/expected_filename

    def test_get_array_config(self):
        class FakeSB:
            is_7m = None
            requires_TP = None
        sb = FakeSB()
        sb.is_7m = False
        assert self.general_sb_simulation.get_array_config(sb)\
                               == self.general_sb_simulation.array_config_12m
        sb.is_7m = True
        sb.requires_TP = True
        assert self.general_sb_simulation.get_array_config(sb) == "aca.cm10.pm3.cfg"
        sb.requires_TP = False
        assert self.general_sb_simulation.get_array_config(sb) == "7m"

    def test_simulate_HAs_retries_on_failure(self):
        #inspired by ChatGPT
        class Result:
            def __init__(self, success):
                self.success = success
        results_with_fail = [Result(False), Result(True)]
        results_all_success = [Result(True), Result(True)]

        class FakeSingleSBSimulation:
            xml_filepath = "/tmp/file.xml"
            default_HA_step = SingleSBSimulation.default_HA_step
            fine_HA_step = SingleSBSimulation.fine_HA_step
            date = datetime.date(1978,3,3)
            simulate_HAs = SingleSBSimulation.simulate_HAs

        class FakePlanner:
            def HA_jobs(self, xml_filepath, array_config, date, step):
                return [SimpleNamespace(HA=Angle(1*u.hour)),
                        SimpleNamespace(HA=Angle(1*u.hour)+step)]

        class FakeRunnerFail:
            def run_jobs(self, jobs):
                return results_with_fail
    
        class FakeRunnerSuccess:
            def run_jobs(self, jobs):
                return results_all_success
    
        single_sim = FakeSingleSBSimulation()
        planner = FakePlanner()

        kwargs = {"planner":planner,"array_config":"c43-5"}
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerSuccess)
        assert HAs[1]-HAs[0] == SingleSBSimulation.default_HA_step
        assert results == results_all_success
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerFail)
        assert HAs[1]-HAs[0] == SingleSBSimulation.fine_HA_step
        assert results == results_with_fail

    def test_simulate(self):
        class FakeJob:
            def __init__(self,HA):
                self.HA = HA
        class FakeSingleSBSimulation:
            default_HA_step = SingleSBSimulation.default_HA_step
            fine_HA_step = SingleSBSimulation.fine_HA_step
            date = datetime.date(2026,1,3)
            HAs = [Angle(1*u.hour),Angle(1.5*u.hour)]
            sb_name = "test_SB"
            project_code = "2025.1.00001.S"
            simulate = SingleSBSimulation.simulate
            def prepare_sb(self):
                return
            def get_array_config(self,sb):
                return "c43-3" 
            def simulate_HAs(self,planner,array_config):
                results = [f"result_{HA}" for HA in self.HAs]
                return self.HAs,results
        with patch("batch_simulations.application.single_SB_simulation.JobPlanner") as MockJobPlanner:
            MockJobPlanner.return_value = None
            fake = FakeSingleSBSimulation()
            sim_summary = fake.simulate()
            assert sim_summary.date == FakeSingleSBSimulation.date
            assert sim_summary.HAs == FakeSingleSBSimulation.HAs
            assert sim_summary.simulation_results\
                       == [f"result_{HA}" for HA in FakeSingleSBSimulation.HAs]
            assert sim_summary.unexpected_error is None
            #now the case where an unexpected error is thrown:
            error_message = "something went wrong"
            def prepare_sb_malfunction():
                raise RuntimeError(error_message)
            fake.prepare_sb = prepare_sb_malfunction
            result = fake.simulate()
            assert result.date == FakeSingleSBSimulation.date
            assert result.HAs is None
            assert result.simulation_results is None
            for key in ("error_message","full_traceback"):
                assert "RuntimeError" in result.unexpected_error[key]
                assert error_message in result.unexpected_error[key]


class TestAnalysisResult:

    def test_merge(self):
        empty_analysis = AnalysisResult()
        assert empty_analysis.inspection_reasons == []
        assert not empty_analysis.should_be_Waiting
        analysis1 = AnalysisResult(inspection_reasons=["problem"],should_be_Waiting=True)
        analysis2 = AnalysisResult(inspection_reasons=["problem","catastrophy"],
                                   should_be_Waiting=True)
        analysis3 = AnalysisResult(inspection_reasons=["minor"],should_be_Waiting=False)
        for ana in (analysis1,analysis2,analysis3):
            merged = ana.merge(empty_analysis)
            assert len(merged.inspection_reasons) > 0
            assert sorted(merged.inspection_reasons) == sorted(ana.inspection_reasons)
            assert merged.should_be_Waiting == ana.should_be_Waiting
        merged = analysis1.merge(analysis2)
        assert sorted(merged.inspection_reasons) == sorted(["problem","catastrophy"])
        assert merged.should_be_Waiting 
        merged = analysis1.merge(analysis3)
        assert sorted(merged.inspection_reasons) == sorted(["problem","minor"])
        assert merged.should_be_Waiting
        merged = analysis2.merge(analysis3)
        assert sorted(merged.inspection_reasons) == sorted(["problem","catastrophy","minor"])
        assert merged.should_be_Waiting


class TestSingleSBSimulationSummary:

    def test_HA_widths(self):
        class FakeSummary:
            HA_widths = SingleSBSimulationSummary.HA_widths
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
        fake_summary = FakeSummary()
        widths = fake_summary.HA_widths()
        assert widths == [0.5,0.5+0.75,0.75+1.75,1.75]

    def test_runnable_HA_amount(self):
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
            HA_widths = SingleSBSimulationSummary.HA_widths
            simulation_results = [SimpleNamespace(success=True),SimpleNamespace(success=False),
                                  SimpleNamespace(success=True),SimpleNamespace(success=True)]
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-2*u.hour),"max":Angle(2*u.hour)})
            runnable_HA_amount = SingleSBSimulationSummary.runnable_HA_amount
        fake_summary = FakeSummary()
        widths = fake_summary.HA_widths()
        usable_HA_range = fake_summary.runnable_HA_amount()
        expected_range = widths[0]+widths[2]
        assert usable_HA_range == expected_range

    def test_analyse_runnable_HA_range(self):
        class FakeSummary:
            def runnable_HA_amount(self):
                return 2
            analyse_runnable_HA_range = SingleSBSimulationSummary.analyse_runnable_HA_range
        fake_summary = FakeSummary()
        fake_summary.min_runnable_HA_amount = 1
        ana = fake_summary.analyse_runnable_HA_range()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary.min_runnable_HA_amount = 3.3333
        ana = fake_summary.analyse_runnable_HA_range()
        assert ana.inspection_reasons == ["runnable HA range is small"]
        assert not ana.should_be_Waiting

    def test_analyse_HA_restriction(self):
        class FakeSB:
            OT_allowed_HA = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
            def __init__(self,has_hardcoded):
                self.has_hardcoded = has_hardcoded
            def any_calibrator_hardcoded(self):
                return self.has_hardcoded
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
            simulation_results = [SimpleNamespace(success=True),SimpleNamespace(success=False),
                                  SimpleNamespace(success=True),SimpleNamespace(success=True)]
            analyse_HA_restriction = SingleSBSimulationSummary.analyse_HA_restriction
            def __init__(self,has_hardcoded):
                self.sb = FakeSB(has_hardcoded=has_hardcoded)
        fake_summary = FakeSummary(has_hardcoded=True)
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary = FakeSummary(has_hardcoded=False)
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(12*u.hour)}
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == ["unnecessarily restricted HAs: -2.0"]
        assert not ana.should_be_Waiting
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(2*u.hour)}
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == ["unnecessarily restricted HAs: -2.0, 4.0"]
        assert not ana.should_be_Waiting

    def test_HA_is_allowed(self):
        class FakeSummary:
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-2.1*u.hour),"max":Angle(3*u.hour)})
            HA_is_allowed = SingleSBSimulationSummary.HA_is_allowed
        fake_summary = FakeSummary()
        for HA in (fake_summary.sb.OT_allowed_HA["min"], fake_summary.sb.OT_allowed_HA["max"],
                   Angle(-2*u.hour),Angle(0*u.hour),Angle(1.1*u.hour)):
            assert fake_summary.HA_is_allowed(HA)
        for HA in (Angle(-2.2*u.hour),Angle(3.1*u.hour)):
            assert not fake_summary.HA_is_allowed(HA)

    def test_analyse_simulation_failures(self):
        def sim_result(success,category=None):
            fail_reason = None if success else\
                   FailReason(error_message="a",error_summary="b",category=category)
            return SimpleNamespace(success=success,fail_reason=fail_reason)
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(0.5*u.hour)]
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-4*u.hour),"max":Angle(4*u.hour)})
            analyse_simulation_failures = SingleSBSimulationSummary.analyse_simulation_failures
            HA_is_allowed = SingleSBSimulationSummary.HA_is_allowed

        def assert_correct_analysis(simulation_results,expected_reasons,expected_waiting,
                                    OT_allowed_HA=None):
            fake_summary = FakeSummary()
            if OT_allowed_HA is not None:
                fake_summary.sb.OT_allowed_HA = OT_allowed_HA
            fake_summary.simulation_results = simulation_results
            ana = fake_summary.analyse_simulation_failures()
            assert sorted(ana.inspection_reasons) == sorted(expected_reasons)
            assert ana.should_be_Waiting == expected_waiting

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error"],
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error"],
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error", "simulation failure"],
                                expected_waiting=True)
        
        #making the missing calibrator failure to occur outside of allowed HA:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error",],
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="exceeds 2h limit")],
                                expected_reasons=["simulation failure"],
                                expected_waiting=True)

        #server error outisde of allowed HA, should not trigger inspection:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="server error"),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

    def test_analyse(self):
        class FakeSummary:
            analyse = SingleSBSimulationSummary.analyse
            unexpected_error = {"error_message":"kabumm",
                                "full_traceback":"the full traceback"}
        fake_summary = FakeSummary()
        fake_summary.analyse()
        assert fake_summary.analysis_result.inspection_reasons == ["unexpected error"]
        assert fake_summary.analysis_result.should_be_Waiting

    def test_analyse_no_unexpected_error(self):
        class FakeSimResult:
            def __init__(self,success,category=None):
                self.success = success
                self.fail_reason = None if category is None\
                                    else FailReason(error_message="bla",
                                                    error_summary="bla",
                                                    category=category)
        #this SB has min allowed HA -12H and max 12H
        sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
        summary_kwargs = {"sb":sb,"config":"c43-9","date":datetime.date(1974,2,22),
                          "unexpected_error":None}

        #all success:
        HAs = [Angle(-12*u.hour),Angle(-2*u.hour),Angle(-0*u.hour),Angle(2*u.hour),
               Angle(4*u.hour),]
        simulation_results = [FakeSimResult(True)]*len(HAs)
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse_no_unexpected_error()
        assert summary.analysis_result.inspection_reasons == []
        assert not summary.analysis_result.should_be_Waiting
        
        #all failures:
        simulation_results = [FakeSimResult(False,category="missing calibrator"),]*5
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse_no_unexpected_error()
        assert sorted(summary.analysis_result.inspection_reasons)\
                  == sorted(["simulation failure","runnable HA range is small"])
        assert summary.analysis_result.should_be_Waiting

        #trigger unnecessary restrikted HA
        HAs[0] = Angle(-14*u.hour)
        simulation_results = [FakeSimResult(True)]*len(HAs)
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse_no_unexpected_error()
        assert summary.analysis_result.inspection_reasons == ['unnecessarily restricted HAs: -14.0']
        assert not summary.analysis_result.should_be_Waiting