#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:01:44 2026

@author: gianni
"""

from batch_simulations.application.batch import SingleSBSimulation,SimulationCampaign,\
        CampaignResultFormatter,AnalysisResult,SingleSBSimulationSummary
from batch_simulations.infrastructure.sb_table import SBTable
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
from batch_simulations.domain.simulation_result import FailReason
import datetime
from pathlib import Path
from astropy.coordinates import Angle
from astropy import units as u
from unittest.mock import patch
import pytest
from types import SimpleNamespace
import math
import pandas as pd
import tempfile
from batch_simulations.infrastructure.table_writer import TableWriter


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
        with patch("batch_simulations.application.batch.JobPlanner") as MockJobPlanner:
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


sb_table_folder = Path("tests/input_tables")
sb_table_filepaths = {"12m":sb_table_folder / "configuration_lookup_table_cycle_12.csv",
                      "7m":sb_table_folder / "2026-02-06_schedBlockList_7M.csv"}
sb_table_array_config_12m = "c43-5"
sb_table = SBTable(input_filepaths=sb_table_filepaths,
                   array_config_12m=sb_table_array_config_12m,SB_filter=None)

class FakeFormatterForWritingTests:
    campaign = SimpleNamespace(name="test")
    writer = TableWriter()
    p2g_columns = CampaignResultFormatter.p2g_columns
    write_master_table = CampaignResultFormatter.write_master_table
    write_table_for_P2G = CampaignResultFormatter.write_table_for_P2G
    def __init__(self):
        #write a new table for each new instance, just to prevent issues
        #if any of the test modifies the table inplace
        master_table_data = {column: ['test',]*4 for column in TableWriter.column_order}
        self.master_table = pd.DataFrame(master_table_data)


class TestCampaignResultFormatter:

    test_sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
    date = datetime.date(1623,2,12)
    HAs = [Angle(-1*u.hour),Angle(1.5*u.hour)]
    fake_simulation_results = [SimpleNamespace(success=True,fail_reason=None),
                               SimpleNamespace(success=False,
                                               fail_reason=FailReason(error_message="message",
                                                                      error_summary="no dgc",
                                                                      category="missing calibrator"))]
    fake_sb_sim_summary = SimpleNamespace(sb=test_sb,config="c43-6",
                                          date=datetime.date(1254,4,8),
                                          HAs=HAs,
                                          analysis_result=AnalysisResult(inspection_reasons=["simulation failure"],
                                                                         should_be_Waiting=True),
                                          simulation_results=fake_simulation_results,
                                          unexpected_error=None)
    fake_campaign = SimpleNamespace(sb_simulation_summaries=[fake_sb_sim_summary,]*2,
                                    sb_table=SimpleNamespace(data=sb_table.data.iloc[:2]),
                                    date=datetime.date(1258,6,11))

    def test_check_consistency(self):
        fake_campaign = SimpleNamespace(sb_table = SimpleNamespace(data=[1,2,3]),
                                        sb_simulation_summaries = ["s","d","ff"])
        CampaignResultFormatter(campaign=fake_campaign)
        fake_campaign.sb_table = SimpleNamespace(data=[1,2])
        with pytest.raises(RuntimeError):
            CampaignResultFormatter(campaign=fake_campaign)

    def test_add_general_info(self):
        row = {}
        sb_table_row = sb_table.data.iloc[0]
        class FakeFormatter:
            sb_table_keys_for_output = CampaignResultFormatter.sb_table_keys_for_output
            add_general_info = CampaignResultFormatter.add_general_info
        FakeFormatter().add_general_info(row=row,sb_table_row=sb_table_row)
        assert len(row) == len(CampaignResultFormatter.sb_table_keys_for_output)
        for key in CampaignResultFormatter.sb_table_keys_for_output:
            if isinstance(sb_table_row[key], float):
                if math.isnan(sb_table_row[key]):
                    assert math.isnan(row[key])
                    continue
            assert row[key] == sb_table_row[key]

    def test_add_note_to_aod(self):
        row = {}
        CampaignResultFormatter.add_note_to_aod(row=row,sb=self.test_sb)
        assert len(row) == 1
        assert row["note_to_AoD"] == self.test_sb.metadata["note_to_AoD"]

    def test_add_hardcoded_calibrators(self):
       row = {}
       CampaignResultFormatter.add_hardcoded_calibrators(row=row,sb=self.test_sb)
       assert len(row) == 1
       hardcoded_calibrators = [c.cal_type for c in self.test_sb.calibrators if
                                c.is_hardcoded]
       expected_entry = ",".join(hardcoded_calibrators)
       assert row["hardcoded_calibrators"] == expected_entry

    def test_add_individual_calibrators(self):
        row = {}
        CampaignResultFormatter.add_individual_calibrators(row=row,sb=self.test_sb)
        assert len(row) == len(CALIBRATOR_TYPES)
        for cal_type in CALIBRATOR_TYPES:
            if cal_type in self.test_sb.cal_types:
                assert row[cal_type] == self.test_sb.get_calibrator(cal_type).source_name
            else:
                assert row[cal_type] is None

    def test_add_HA_limits(self):
        local_test_sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
        HA_DSA = self.test_sb.get_DSA_HA_limits()
        OT_HA_limits = [{"min":Angle(-12*u.hour),"max":Angle(12*u.hour)},
                        HA_DSA,
                        {"min":HA_DSA["min"]+Angle(0.5*u.hour),"max":Angle(12*u.hour)},
                        {"min":Angle(-12*u.hour),"max":HA_DSA["min"]-Angle(0.5*u.hour)},
                        {"min":HA_DSA["min"]+Angle(0.5*u.hour),
                         "max":HA_DSA["min"]-Angle(0.5*u.hour)}]
        for i,OT_HA_lim in enumerate(OT_HA_limits):
            local_test_sb.OT_allowed_HA = OT_HA_lim
            row = {}
            CampaignResultFormatter.add_HA_limits(row=row,sb=local_test_sb)
            assert len(row) == 5
            for lim in ("min", "max"):
                assert row[f"{lim}_HA_DSA"] == HA_DSA[lim].hour
                row[f"{lim}_HA_OT"] == local_test_sb.OT_allowed_HA[lim]
            if i<=1:
                assert not row["HA_is_restricted"]
            else:
                assert row["HA_is_restricted"]

    def test_add_simulation_info(self):
        class FakeFormatter:
            campaign = SimpleNamespace(date=self.date)
            add_simulation_info = CampaignResultFormatter.add_simulation_info
            @staticmethod
            def build_per_HA_summary_string(simulation_results, HAs):
                return CampaignResultFormatter.build_per_HA_summary_string(
                               simulation_results=simulation_results,HAs=HAs)
        class FakeSBSummary:
            config = "c43-7"
            simulation_results = self.fake_simulation_results
            HAs = self.HAs
            date = self.date
            def __init__(self,unexpected_error,analysis_result):
                self.unexpected_error = unexpected_error
                self.analysis_result = analysis_result
        fake_formatter = FakeFormatter()
        row = {}
        sb_sim_summary = FakeSBSummary(unexpected_error=None,
                                       analysis_result = AnalysisResult(inspection_reasons=["simulation failure"],
                                                                        should_be_Waiting=True))
        fake_formatter.add_simulation_info(row=row,sb_sim_summary=sb_sim_summary)
        assert len(row) == 5
        assert row["simulated_config"] == FakeSBSummary.config
        assert row["simulated_date"] == str(self.date)
        assert row["simulations"] == CampaignResultFormatter.build_per_HA_summary_string(
                                      simulation_results=FakeSBSummary.simulation_results,
                                      HAs=FakeSBSummary.HAs)
        assert row["inspection reasons"] == "simulation failure"
        assert row["should be Waiting"]
        #case with unexpected error:
        unexpected_error = {"error_message":"kabumm","full_traceback":"the full traceback"}
        row = {}
        sb_sim_summary = FakeSBSummary(unexpected_error=unexpected_error,
                                       analysis_result = AnalysisResult(inspection_reasons=["unexpected error"],
                                                                        should_be_Waiting=True))
        fake_formatter.add_simulation_info(row=row,sb_sim_summary=sb_sim_summary)
        assert len(row) == 4
        assert row["simulations"] == unexpected_error["error_message"]
        assert row["traceback of unexpected error"] == unexpected_error["full_traceback"]
        assert row["inspection reasons"] == "unexpected error"
        assert row["should be Waiting"]

    def test_build_per_HA_summary_string(self):
        expected_summary = ""
        for result, HA in zip(self.fake_simulation_results,self.HAs):
            if result.success:
                expected_summary += f"{HA.hour:.3g}: success\n"
            else:
                expected_summary += f"{HA.hour:.3g}: {result.fail_reason.error_summary}\n"
        expected_summary = expected_summary[:-1]
        summary = CampaignResultFormatter.build_per_HA_summary_string(
                         simulation_results=self.fake_simulation_results, HAs=self.HAs)
        assert summary == expected_summary
        
    def test_build_row(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        row = formatter.build_row(sb_sim_summary=self.fake_sb_sim_summary,
                                  sb_table_row=sb_table.data.iloc[0])
        assert len(row) > 0
        unexpected_error = {"error_message":"whoops","full_traceback":"the details"}
        fake_sb_sim_summary_failed = SimpleNamespace(
                                        sb=None,config="c43-6",
                                        date=datetime.date(1254,4,8),HAs=None,
                                        simulation_results=None,
                                        analysis_result=AnalysisResult(inspection_reasons="unexpected error",
                                                                       should_be_Waiting=True),
                                        unexpected_error=unexpected_error)
        row = formatter.build_row(sb_sim_summary=fake_sb_sim_summary_failed,
                                  sb_table_row=sb_table.data.iloc[0])
        assert len(row) > 0

    def test_create_master_dataframe(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        formatter.create_master_dataframe()
        sorted_sb_table = self.fake_campaign.sb_table.data.sort_values("code")
        assert formatter.master_table.iloc[0]["code"] == sorted_sb_table.iloc[0].code
        assert formatter.master_table.iloc[0]["note_to_AoD"] == self.test_sb.metadata["note_to_AoD"]

    def test_write_mastertable(self):
        fake_formatter = FakeFormatterForWritingTests()
        with tempfile.TemporaryDirectory() as tmpdirname:
            for out_format in ("csv","xlsx"):
                fake_formatter.write_master_table(out_format=out_format,
                                                  output_dir=tmpdirname)
                filepath = tmpdirname / Path(f"master_table_{fake_formatter.campaign.name}.{out_format}")
                assert filepath.exists()

    def test_write_table_for_p2g(self):
        fake_formatter = FakeFormatterForWritingTests()
        with tempfile.TemporaryDirectory() as tmpdirname:
            for out_format in ("csv","xlsx"):
                # filepath = tmpdirname / Path(f"test.{extension}")
                fake_formatter.write_table_for_P2G(out_format=out_format,
                                                   output_dir=tmpdirname)
                filepath = tmpdirname / Path(f"p2g_table_{fake_formatter.campaign.name}.{out_format}")
                assert filepath.exists()
                if out_format == "csv":
                    written_data = pd.read_csv(filepath)
                elif out_format == "xlsx":
                    written_data = pd.read_excel(filepath)
                else:
                    raise RuntimeError
                assert sorted(written_data.columns) == sorted(fake_formatter.p2g_columns)


class TestCampaign:

    def test_run(self):
        class FakeSimSummary:
            def __init__(self,project_code,sb_name,array_config_12m):
                self.project_code = project_code
                self.sb_name = sb_name
                self.array_config_12m = array_config_12m
            def analyse(self):
                pass
        class FakeSingleSBSimulation:
            def __init__(self,project_code,sb_name,date,array_config_12m):
                self.project_code = project_code
                self.sb_name = sb_name
                self.date = date
                self.array_config_12m = array_config_12m
            def simulate(self):
                return FakeSimSummary(project_code=self.project_code,sb_name=self.sb_name,
                                      array_config_12m=self.array_config_12m)
        date = datetime.date(1026,2,3)
        campaign = SimulationCampaign(name="test",sb_table=sb_table, date=date)
        campaign.run(single_sb_sim_cls=FakeSingleSBSimulation)
        assert len(sb_table.data) == len(campaign.sb_simulation_summaries)
        for row,sim_summary in zip(sb_table.data.itertuples(),campaign.sb_simulation_summaries):
            assert row.code == sim_summary.project_code
            assert row.sbname == sim_summary.sb_name
            assert sb_table.array_config_12m == sim_summary.array_config_12m