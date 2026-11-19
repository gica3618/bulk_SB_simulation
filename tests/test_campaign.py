#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:01:44 2026

@author: gianni
"""

from batch_simulations.application.campaign import SimulationCampaign,\
        CampaignResultFormatter
from batch_simulations.application.single_SB_simulation import AnalysisResult
from batch_simulations.infrastructure.sb_table import SBTable
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.simulation_result import FailReason
import datetime
from pathlib import Path
from astropy.coordinates import Angle
from astropy import units as u
import pytest
from types import SimpleNamespace
import math
import pandas as pd
import tempfile
from batch_simulations.infrastructure.table_writer import TableWriter
import itertools


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
    need_inspection_master_table_selection\
                 = CampaignResultFormatter.need_inspection_master_table_selection
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
    analysis_result = AnalysisResult()
    analysis_result.inspection_reasons = set(["missing calibrator"])
    analysis_result.recommended_state = "Waiting"
    fake_sb_sim_analysis = SimpleNamespace(sb=test_sb,simulated_config="c43-6",
                                           simulated_date=datetime.date(1254,4,8),
                                           HAs=HAs,
                                           analysis_result=analysis_result,
                                           simulation_results=fake_simulation_results,
                                           unexpected_error=None)
    fake_campaign = SimpleNamespace(sb_simulation_analyses=[fake_sb_sim_analysis,]*2,
                                    sb_table=SimpleNamespace(data=sb_table.data.iloc[:2]),
                                    date=datetime.date(1258,6,11))

    def test_check_consistency(self):
        fake_campaign = SimpleNamespace(sb_table = SimpleNamespace(data=[1,2,3]),
                                        sb_simulation_analyses = ["s","d","ff"])
        CampaignResultFormatter(campaign=fake_campaign)
        fake_campaign.sb_table = SimpleNamespace(data=[1,2])
        with pytest.raises(RuntimeError):
            CampaignResultFormatter(campaign=fake_campaign)

    def test_create_master_dataframe(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        formatter.create_master_dataframe()
        sorted_sb_table = self.fake_campaign.sb_table.data.sort_values("code")
        assert formatter.master_table.iloc[0]["code"] == sorted_sb_table.iloc[0].code
        assert formatter.master_table.iloc[0]["note_to_AoD"] == self.test_sb.metadata["note_to_AoD"]

    def test_build_row(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        row = formatter.build_row(sb_sim_analysis=self.fake_sb_sim_analysis,
                                  sb_table_row=sb_table.data.iloc[0])
        assert len(row) > 0
        unexpected_error = {"error_message":"whoops","full_traceback":"the details"}
        analysis_result = AnalysisResult()
        analysis_result.inspection_reasons = set(["unexpected error"])
        analysis_result.recommended_state = "Waiting"
        fake_sb_sim_analysis_failed = SimpleNamespace(
                                        sb=None,simulated_config="c43-6",
                                        simulated_date=datetime.date(1254,4,8),
                                        HAs=None,simulation_results=None,
                                        analysis_result=analysis_result,
                                        unexpected_error=unexpected_error)
        row = formatter.build_row(sb_sim_analysis=fake_sb_sim_analysis_failed,
                                  sb_table_row=sb_table.data.iloc[0])
        assert len(row) > 0

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
        #SB without hardcoded calibrators:
       row = {}
       CampaignResultFormatter.add_hardcoded_calibrators(row=row,sb=self.test_sb)
       assert len(row) == 1
       assert row["hardcoded_calibrators"] is pd.NA
       #SB with hardcoded calibrators:
       row = {}
       sb = BuildSBFromXML.build("tests/xmls/example_bandpass_equal_dgc.xml")
       CampaignResultFormatter.add_hardcoded_calibrators(row=row,sb=sb)
       hardcoded_calibrators = [f"{c.cal_type} ({c.source_name})" for c in
                                sb.calibrators if c.is_hardcoded]
       assert hardcoded_calibrators
       expected_entry = ", ".join(hardcoded_calibrators)
       assert len(row) == 1
       assert row["hardcoded_calibrators"] == expected_entry

    def test_add_query_calibrators(self):
        row = {}
        sb = BuildSBFromXML.build("tests/xmls/example_bandpass_equal_dgc.xml")
        #fake SB with all calibrators hardcoded
        for cal in sb.calibrators:
            cal.is_hardcoded = True
        assert all([c.is_hardcoded for c in sb.calibrators])
        CampaignResultFormatter.add_query_calibrators(row=row,sb=sb)
        assert row["query_calibrators"] is pd.NA
        sb = BuildSBFromXML.build("tests/xmls/example_bandpass_equal_dgc.xml")
        CampaignResultFormatter.add_query_calibrators(row=row,sb=sb)
        query_cals = [f"{c.cal_type}" for c in sb.calibrators if not c.is_hardcoded]
        expected_entry = ", ".join(query_cals)
        assert row["query_calibrators"] == expected_entry

    def test_add_HA_limits(self):
        sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
        HA_DSA = sb.get_DSA_HA_limits()
        OT_HA_limits = [{"min":Angle(-12*u.hour),"max":Angle(12*u.hour)},
                        HA_DSA,
                        {"min":HA_DSA["min"]+Angle(0.5*u.hour),"max":Angle(12*u.hour)},
                        {"min":Angle(-12*u.hour),"max":HA_DSA["min"]-Angle(0.5*u.hour)},
                        {"min":HA_DSA["min"]+Angle(0.5*u.hour),
                         "max":HA_DSA["min"]-Angle(0.5*u.hour)}]
        recommended_HA_limits = [None,{"min":Angle(-1*u.hour),"max":Angle(1*u.hour)}]
        for OT_HA_lim,rec_HA_lim in itertools.product(OT_HA_limits,recommended_HA_limits):
            sb.OT_allowed_HA = OT_HA_lim
            analysis_result = SimpleNamespace(recommended_HA_limits=rec_HA_lim)
            sim_analysis = SimpleNamespace(sb=sb,analysis_result=analysis_result)
            row = {}
            CampaignResultFormatter.add_HA_limits(row=row,sb_sim_analysis=sim_analysis)
            assert len(row) == 7
            for lim in ("min","max"):
                assert row[f"HA_{lim}"] == sb.OT_allowed_HA[lim].hour
                assert row[f"HA_{lim}_DSA"] == HA_DSA[lim].hour
                if rec_HA_lim is None:
                    assert row[f"HA_{lim}_recommended"] is pd.NA
                else:
                    assert row[f"HA_{lim}_recommended"] == rec_HA_lim[lim].hour
            assert row["HA_is_restricted"] == sb.HA_is_restricted()

    def test_add_simulation_info(self):

        class FakeFormatter:
            add_simulation_info = CampaignResultFormatter.add_simulation_info

            @staticmethod
            def build_per_HA_summary_string(simulation_results, HAs):
                return CampaignResultFormatter.build_per_HA_summary_string(
                               simulation_results=simulation_results,HAs=HAs)

        #with unexpected error:
        row = {}
        sb_sim_analysis = SimpleNamespace(unexpected_error={"error_message":"Underlix",
                                                            "full_traceback":"beim teutates"})
        fake = FakeFormatter()
        fake.add_simulation_info(row=row,sb_sim_analysis=sb_sim_analysis)
        assert len(row) == 2
        assert row["simulations"] == sb_sim_analysis.unexpected_error["error_message"]
        assert row["traceback_of_unexpected_error"] == sb_sim_analysis.unexpected_error["full_traceback"]

        #no unexpected error:
        row = {}
        simulation_results = [SimpleNamespace(success=True),
                              SimpleNamespace(success=False,fail_reason=SimpleNamespace(error_summary="idefix"))]
        HAs = [Angle(-3*u.hour),Angle(2*u.hour)]
        sb_sim_analysis = SimpleNamespace(
                              simulated_config="c43-9",simulated_date=datetime.date(1111,12,10),
                              simulation_results=simulation_results,HAs=HAs,
                              unexpected_error=None)
        fake = FakeFormatter()
        fake.add_simulation_info(row=row,sb_sim_analysis=sb_sim_analysis)
        assert len(row) == 3
        assert row["simulated_config"] == sb_sim_analysis.simulated_config
        assert row["simulated_date"] == sb_sim_analysis.simulated_date
        expected_sim_str = "-3: success\n2: idefix"
        assert row["simulations"] == expected_sim_str

    def test_add_inspection_reasons(self):
        row = {}
        inspection_reasons = set()
        CampaignResultFormatter.add_inspection_reasons(row=row,inspection_reasons=inspection_reasons)
        assert row["inspection_reasons"] is pd.NA
        row = {}
        inspection_reasons = set(["underlix","miraculix"])
        CampaignResultFormatter.add_inspection_reasons(row=row,inspection_reasons=inspection_reasons)
        assert row["inspection_reasons"] == "; ".join(inspection_reasons)

    def test_add_recommended_state(self):
        row = {}
        recommended_state = None
        CampaignResultFormatter.add_recommended_state(row=row,recommended_state=recommended_state)
        assert row["recommended_state"] is pd.NA
        row = {}
        recommended_state = "Waiting"
        CampaignResultFormatter.add_recommended_state(row=row,recommended_state=recommended_state)
        assert row["recommended_state"] == recommended_state

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


class FakeSimAnalysis:

    def __init__(self,include_2h_fail=False,inspection_reasons=None):
        ir = set() if inspection_reasons is None else inspection_reasons
        self.analysis_result = SimpleNamespace(inspection_reasons=ir)
        self.simulation_results = [SimpleNamespace(success=True,fail_reason=None),
                                   SimpleNamespace(success=False,
                                                   fail_reason=FailReason(error_message="test",
                                                                          error_summary="no check",
                                                                          category="missing calibrator"))]
        no_phase = SimpleNamespace(success=False,
                                   fail_reason=FailReason(error_message="test",
                                                          error_summary="no phase",
                                                          category="missing calibrator"))
        #add same result twice:
        self.simulation_results += [no_phase]*2
        assert len(self.simulation_results) == 4
        if include_2h_fail:
            two_hour_fail = SimpleNamespace(
                               success=False,
                               fail_reason=FailReason(error_message="test",
                                                      error_summary="SB too long by X seconds",
                                                      category="exceeds 2h limit"))
            self.simulation_results.append(two_hour_fail)


class TestCampaignResultFormatterSummaryFile:

    def test_need_inspection_selection(self):
        class FakeFormatter:
            master_table = pd.DataFrame({"inspection_reasons":[pd.NA,"problem",pd.NA,pd.NA,"Asterix"]})
            need_inspection_master_table_selection\
                   = CampaignResultFormatter.need_inspection_master_table_selection
        fake = FakeFormatter()
        need_inspection = fake.need_inspection_master_table_selection()
        assert all(need_inspection==[False,True,False,False,True])

    def test_summarize_number_of_simulated_SBs(self):
        def get_number_of_12m_and_7m_SBs():
            return {"12m":3,"7m":2}
        data = pd.DataFrame({"sb_uid":[f"uid{i}" for i in range(5)]})
        sb_table = SimpleNamespace(data=data,
                                   get_number_of_12m_and_7m_SBs=get_number_of_12m_and_7m_SBs)
        class FakeFormatter:
            campaign = SimpleNamespace(sb_table=sb_table)
            summarize_number_of_simulated_SBs = CampaignResultFormatter.summarize_number_of_simulated_SBs
        summary = FakeFormatter().summarize_number_of_simulated_SBs()
        assert summary == "nb. of simulated SBs: 5 (12m: 3; 7m: 2)\n"

    def test_summarize_unexpected_errors(self):
        class SimAnalysis:
            def __init__(self,unexpected_error):
                self.unexpected_error = unexpected_error
        class FakeFormatter:
            sb_simulation_analyses = [SimAnalysis(None),SimAnalysis(None),
                                      SimAnalysis("test"),SimAnalysis(None)]
            campaign = SimpleNamespace(sb_simulation_analyses=sb_simulation_analyses)
            summarize_unexpected_errors = CampaignResultFormatter.summarize_unexpected_errors
        summary = FakeFormatter().summarize_unexpected_errors()
        assert summary == "1 SBs yielded an unexpected error\n"

    def test_summarize_number_of_inspections(self):
        class FakeFormatter:
            master_table = pd.DataFrame({"inspection_reasons":[pd.NA,"problem",pd.NA,pd.NA,"Asterix"]})
            need_inspection_master_table_selection\
                     = CampaignResultFormatter.need_inspection_master_table_selection
            summarize_number_of_inspections = CampaignResultFormatter.summarize_number_of_inspections
        summary = FakeFormatter().summarize_number_of_inspections()
        assert summary == "2 SBs need inspection\n"

    def test_get_count_lines(self):
        reasons = ["Asterix","Obelix","Asterix","no check","no phase","no check"]
        count_lines = CampaignResultFormatter.get_count_lines(reasons)
        assert count_lines == "Asterix: 2\n"+"no check: 2\n"+"Obelix: 1\n"\
                               +"no phase: 1\n"
        reasons = []
        assert CampaignResultFormatter.get_count_lines(reasons) == "None\n"

    def test_summarize_inspection_reasons(self):
        class FakeFormatter:
            sb_simulation_analyses = [FakeSimAnalysis(inspection_reasons=None),
                                      FakeSimAnalysis(inspection_reasons=set(["error 1", "error 2"])),
                                      FakeSimAnalysis(inspection_reasons=set(["error 1"])),
                                      FakeSimAnalysis(inspection_reasons=set(["missing calibrator"]))]
            campaign = SimpleNamespace(sb_simulation_analyses=sb_simulation_analyses)
            summarize_inspection_reasons = CampaignResultFormatter.summarize_inspection_reasons
            @staticmethod
            def get_count_lines(reasons):
                return CampaignResultFormatter.get_count_lines(reasons)
        summary = FakeFormatter().summarize_inspection_reasons()
        assert summary == ("inspection reasons:\nerror 1: 2\n" + "error 2: 1\n"
                           + "missing calibrator: 1\n")

    def test_get_fail_reasons_per_SB(self):
        sb_simulation_analyses = [FakeSimAnalysis(include_2h_fail=True),
                                  FakeSimAnalysis(include_2h_fail=False)]
        fail_reasons = CampaignResultFormatter.get_fail_reasons_per_SB(sb_simulation_analyses)
        assert len(fail_reasons) == 5
        assert fail_reasons.count("no check") == 2
        assert fail_reasons.count("no phase") == 2
        assert fail_reasons.count("exceeds 2h limit") == 1

    def test_get_fail_reasons_per_SB_unexpected_error(self):
        sb_simulation_summaries = [SimpleNamespace(simulation_results=None,
                                                   unexpected_error={"error_message":"Asterix"}),]*2
        sb_simulation_summaries.append(SimpleNamespace(simulation_results=None,
                                                       unexpected_error={"error_message":"Miraculix"}))
        fail_reasons = CampaignResultFormatter.get_fail_reasons_per_SB(sb_simulation_summaries)
        assert len(fail_reasons) == 3
        assert fail_reasons.count("Asterix") == 2
        assert fail_reasons.count("Miraculix") == 1

    def test_summarize_fail_reasons(self):
        sb_simulation_analyses = [FakeSimAnalysis(inspection_reasons=None),
                                  FakeSimAnalysis(inspection_reasons=set(["error 1", "error 2"])),
                                  FakeSimAnalysis(inspection_reasons=set(["error 1"]),
                                                  include_2h_fail=True),
                                  FakeSimAnalysis(inspection_reasons=set(["error 3"]))]
        class FakeFormatter:
            campaign = SimpleNamespace(sb_simulation_analyses=sb_simulation_analyses)
            summarize_fail_reasons = CampaignResultFormatter.summarize_fail_reasons
            @staticmethod
            def get_fail_reasons_per_SB(sb_simulation_analyses):
                return CampaignResultFormatter.get_fail_reasons_per_SB(sb_simulation_analyses)
            @staticmethod
            def get_count_lines(reasons):
                return CampaignResultFormatter.get_count_lines(reasons)
        fake = FakeFormatter()
        summary = fake.summarize_fail_reasons()
        split = summary.split("\n")
        assert split[0] == "number of SBs failing with reason X at at least one hour angle (all SBs):"
        assert "no check: 4" in split[1:3]
        assert "no phase: 4" in split[1:3]
        assert split[3] == "exceeds 2h limit: 1"
        assert split[4] == "number of SBs failing with reason X at at least one hour angle (SBs needing inspection):"
        assert "no check: 3" in split[5:7]
        assert "no phase: 3" in split[5:7]
        assert split[7] == "exceeds 2h limit: 1"

    def test_write_campaign_summar(self):
        #just checking that this runs
        class FakeFormatter:
            campaign = SimpleNamespace(name="test",date=datetime.date(1475,3,13),
                                       run_time=datetime.timedelta(seconds=145))
            summarize_number_of_simulated_SBs = lambda self: "test\n"
            summarize_unexpected_errors = lambda self: "rassel\n"
            summarize_number_of_inspections = lambda self: "Asterix und Obelix\n"
            summarize_inspection_reasons = lambda self: "okbre\n"
            summarize_fail_reasons = lambda self: "hallihallo\numgi\ntanpopopo\n"
            write_campaign_summary = CampaignResultFormatter.write_campaign_summary
        fake = FakeFormatter()
        with tempfile.TemporaryDirectory() as tmpdirname:
            fake.write_campaign_summary(output_dir=tmpdirname)
            filepath = tmpdirname / Path("campaign_test_summary.txt")
            assert filepath.exists()


class TestCampaign:

    def test_run(self):
        class FakeSimAnalysis:
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
                return FakeSimAnalysis(project_code=self.project_code,
                                       sb_name=self.sb_name,
                                       array_config_12m=self.array_config_12m)
        date = datetime.date(1026,2,3)
        campaign = SimulationCampaign(name="test",sb_table=sb_table, date=date)
        campaign.run(single_sb_sim_cls=FakeSingleSBSimulation)
        assert len(sb_table.data) == len(campaign.sb_simulation_analyses)
        for row,sim_analysis in zip(sb_table.data.itertuples(),
                                    campaign.sb_simulation_analyses):
            assert row.code == sim_analysis.project_code
            assert row.sbname == sim_analysis.sb_name
            assert sb_table.array_config_12m == sim_analysis.array_config_12m