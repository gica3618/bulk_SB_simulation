#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:01:44 2026

@author: gianni
"""

from batch_simulations.application.batch import SingleSBSimulation,SimulationCampaign,\
        CampaignResultFormatter
from batch_simulations.infrastructure.sb_table import SBTable
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
import datetime
from pathlib import Path
from astropy.coordinates import Angle
from astropy import units as u
from unittest.mock import patch
import pytest
from types import SimpleNamespace
import math


class TestBatch:

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

    def test_simulate_date_retries_on_failure(self):
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
            def jobs(self, xml_filepath, array_config, date, step):
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
            assert "RuntimeError" in result.unexpected_error
            assert error_message in result.unexpected_error


sb_table_folder = Path("tests/input_tables")
sb_table_filepaths = {"12m":sb_table_folder / "configuration_lookup_table_cycle_12.csv",
                      "7m":sb_table_folder / "2026-02-06_schedBlockList_7M.csv"}
sb_table_array_config_12m = "c43-5"
sb_table = SBTable(input_filepaths=sb_table_filepaths,
                   array_config_12m=sb_table_array_config_12m,SB_filter=None)

class TestCampaignResultFormatter:

    test_sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
    date = datetime.date(1623,2,12)
    fake_simulation_results = [SimpleNamespace(success=True,fail_reason=None),
                               SimpleNamespace(success=False,fail_reason="something went wrong")]
    HAs = [Angle(-1*u.hour),Angle(1.5*u.hour)]
    fake_sb_sim_summary = SimpleNamespace(sb=test_sb,config="c43-6",
                                          date=datetime.date(1254,4,8),
                                          HAs=[Angle(-1.2*u.hour),Angle(2*u.hour)],
                                          simulation_results=fake_simulation_results,
                                          unexpected_error=None)
    fake_campaign = SimpleNamespace(sb_simulation_summaries=[fake_sb_sim_summary,]*2,
                                    sb_table=sb_table.data.iloc[:2],
                                    date=datetime.date(1258,6,11))

    def test_check_consistency(self):
        fake_campaign = SimpleNamespace(sb_table = [1,2,3],
                                        sb_simulation_summaries = ["s","d","ff"])
        CampaignResultFormatter(campaign=fake_campaign)
        fake_campaign.sb_table = [1,2]
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
            def build_simulation_result_summary(simulation_results, HAs):
                return CampaignResultFormatter.build_simulation_result_summary(
                               simulation_results=simulation_results,HAs=HAs)
        class FakeSBSummary:
            config = "c43-7"
            simulation_results = self.fake_simulation_results
            HAs = self.HAs
        fake_formatter = FakeFormatter()
        row = {}
        fake_formatter.add_simulation_info(row=row,sb_sim_summary=FakeSBSummary())
        assert len(row) == 3
        assert row["simulated_config"] == FakeSBSummary.config
        assert row["simulated_date"] == str(self.date)
        assert row["simulations"] == CampaignResultFormatter.build_simulation_result_summary(
                                      simulation_results=FakeSBSummary.simulation_results,
                                      HAs=FakeSBSummary.HAs)

    def test_build_simulation_result_summary(self):
        expected_summary = ""
        for result, HA in zip(self.fake_simulation_results,self.HAs):
            if result.success:
                expected_summary += f"{HA.hour:.3g}: success\n"
            else:
                expected_summary += f"{HA.hour:.3g}: {result.fail_reason}\n"
        expected_summary = expected_summary[:-1]
        summary = CampaignResultFormatter.build_simulation_result_summary(
                         simulation_results=self.fake_simulation_results, HAs=self.HAs)
        assert summary == expected_summary
        
    def test_build_row(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        row = formatter.build_row(sb_sim_summary=self.fake_sb_sim_summary,
                                  sb_table_row=sb_table.data.iloc[0])
        assert len(row) > 0

    def test_create_master_dataframe(self):
        formatter = CampaignResultFormatter(campaign=self.fake_campaign)
        formatter.create_master_dataframe()
        assert formatter.master_table.iloc[0]["code"] == sb_table.data.iloc[0].code
        assert formatter.master_table.iloc[0]["note_to_AoD"] == self.test_sb.metadata["note_to_AoD"]


class TestCampaign:

    def test_run(self):
        class FakeSingleSBSimulation:
            def __init__(self,project_code,sb_name,date,array_config_12m):
                self.project_code = project_code
                self.sb_name = sb_name
                self.date = date
                self.array_config_12m = array_config_12m
            def simulate(self):
                return f"result_{self.project_code}_{self.sb_name}_{self.date}_{self.array_config_12m}"
        date = datetime.date(1026,2,3)
        campaign = SimulationCampaign(name="test",sb_table=sb_table, date=date)
        campaign.run(single_sb_sim_cls=FakeSingleSBSimulation)
        expected_summaries = []
        for row in sb_table.data.itertuples():
            s = f"result_{row.code}_{row.sbname}_{date}_{sb_table.array_config_12m}"
            expected_summaries.append(s)
        assert campaign.sb_simulation_summaries == expected_summaries
        