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
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
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