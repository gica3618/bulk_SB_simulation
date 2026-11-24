#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 18:39:38 2026

@author: gianni
"""


from batch_simulations.application.state_changes import\
                                     ProTrackStateChangeWriter as Writer
from batch_simulations.infrastructure.sb_table import SBTable
import pandas as pd
import tempfile
from pathlib import Path
import pytest
import logging
from types import SimpleNamespace


sb_table_folder = Path("tests/input_tables")
sb_table_filepaths = {"12m":sb_table_folder / "configuration_lookup_table_cycle_13_20261120.csv",
                      "7m":sb_table_folder / "2026-11-20_schedBlockList_7M.csv"}
project_uid_lookup_filepath = sb_table_folder / "2026-11-20_projectList.csv"
sb_table_array_config_12m = "c43-5"


class FakeWriter:
    get_ToO_selection = Writer.get_ToO_selection
    get_waiting_for_calibrator_selection = Writer.get_waiting_for_calibrator_selection
    get_ready_selection = Writer.get_ready_selection
    get_no_hardcoded_selection = Writer.get_no_hardcoded_selection
    get_selection_for_setting_to_WaitingForP2G = Writer.get_selection_for_setting_to_WaitingForP2G
    get_selection_for_setting_to_Ready = Writer.get_selection_for_setting_to_Ready
    merge_with_up_to_date_sb_table = Writer.merge_with_up_to_date_sb_table
    update_master_table = Writer.update_master_table
    @staticmethod
    def analyse_number_of_updates(merged):
        return Writer.analyse_number_of_updates(merged)
    @staticmethod
    def remove_missing(merged):
        return Writer.remove_missing(merged)
    def __init__(self):
        #it's important to make a new table each time a new instance is created,
        #because inside some tests, the table is modified in-place
        self.up_to_date_sb_table = SBTable(input_filepaths=sb_table_filepaths,
                                           array_config_12m=sb_table_array_config_12m,
                                           project_uid_lookup_filepath=project_uid_lookup_filepath,
                                           SB_filter=None)


class TestProTrackStateChangeWriter:

    def test_merge_with_up_to_date_sb_table(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = fake_writer.up_to_date_sb_table.data.copy()
        fake_writer.master_table = fake_writer.master_table.iloc[::-1]
        fake_writer.master_table.reset_index(drop=True,inplace=True)
        fake_writer.master_table.sb_state = "Waiting"
        fake_writer.master_table.sb_state_flag = "ForCalibrator"
        merged = fake_writer.merge_with_up_to_date_sb_table()
        assert (merged["sb_state"]==fake_writer.master_table["sb_state"]).all()
        expected_updated_states = fake_writer.up_to_date_sb_table.data["sb_state"]\
                                   .iloc[::-1].reset_index(drop=True)
        assert (merged["sb_state_updated"]==expected_updated_states).all()
        assert not fake_writer.master_table["sb_state_flag"].isna().all()
        assert (merged["sb_state_flag"]
                ==fake_writer.master_table["sb_state_flag"]).all()
        expected_updated_state_flags\
                  = fake_writer.up_to_date_sb_table.data["sb_state_flag"]\
                      .iloc[::-1].reset_index(drop=True).fillna("")
        assert (merged["sb_state_flag_updated"].fillna("")
                ==expected_updated_state_flags).all()

    def test_merge_duplicated_uids(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = fake_writer.up_to_date_sb_table.data.copy()
        fake_writer.up_to_date_sb_table.data.at[1,"sb_uid"]\
                      = fake_writer.up_to_date_sb_table.data.at[0,"sb_uid"]
        with pytest.raises(ValueError):
            fake_writer.merge_with_up_to_date_sb_table()

    def test_merge_missing(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = fake_writer.up_to_date_sb_table.data.copy()
        uid = fake_writer.master_table.sb_uid[0]
        uid_index = fake_writer.up_to_date_sb_table.data.index[
                      fake_writer.up_to_date_sb_table.data['sb_uid'] == uid]
        assert len(uid_index) == 1
        fake_writer.up_to_date_sb_table.data.drop(uid_index,inplace=True)
        merged = fake_writer.merge_with_up_to_date_sb_table()
        assert merged.sb_state_updated.isna()[0]
        assert not merged.sb_state_updated[1:].isna().any()

    def test_remove_missing(self):
        up_to_date_sb_table = SimpleNamespace(data=pd.DataFrame({"sb_state":["Ready","Waiting"],
                                                                 "sb_state_flag":[pd.NA,"ForCalibrator"],
                                                                 "sb_uid":["uid:1","uid:2"]}))
        master_table = pd.DataFrame({"sb_state":["Ready","Waiting","Ready"],
                                     "sb_state_flag":[pd.NA,"ForCalibrator",pd.NA],
                                     "sb_uid":["uid:1","uid:2","uid:3"]})
        fake_writer = FakeWriter()
        fake_writer.up_to_date_sb_table = up_to_date_sb_table
        fake_writer.master_table = master_table
        merged = fake_writer.merge_with_up_to_date_sb_table()
        assert merged[merged["sb_uid"]=="uid:3"]["sb_state_updated"].isna().all()
        missing_removed = fake_writer.remove_missing(merged)
        assert len(missing_removed) == 2
        assert "uid:3" not in missing_removed.sb_uid

    # def test_check_missing_sb_state(self):
    #     merged = pd.DataFrame({"sb_state_updated":["Ready","Waiting"],
    #                            "sb_uid":["uid:1","uid:2"]})
    #     Writer.check_missing_sb_state(merged)
    #     merged.at[0,"sb_state_updated"] = pd.NA
    #     match = re.escape("SB UIDs missing in the up-to-date table: ['uid:1']")
    #     with pytest.raises(ValueError,match=match):
    #         Writer.check_missing_sb_state(merged)

    def test_analyse_number_of_updates(self,caplog):
        merged = pd.DataFrame({"sb_state":["Ready","Ready","Waiting","Waiting","Waiting"],
                               "sb_state_flag":[pd.NA,pd.NA,"ForJAO","ForJAO","ForP2G"],
                               "sb_state_updated":["Ready","FullyObserved","Ready","Waiting","Waiting"],
                               "sb_state_flag_updated":[pd.NA,pd.NA,pd.NA,"ForJAO","ForCalibrator"]})
        with caplog.at_level(logging.INFO):
            Writer.analyse_number_of_updates(merged=merged)
        messages = [record.getMessage() for record in caplog.records]
        assert messages == ["updated 3/5 states"]

    # def test_analyse_state_updates(self):
    #     #just testing that it runs
    #     merged = pd.DataFrame({"sb_state":["Ready","Waiting"],
    #                            "sb_state_flag":[pd.NA,"ForP2G"],
    #                            "sb_state_updated":["Ready","Waiting"],
    #                            "sb_state_flag_updated":[pd.NA,"ForJAO"],
    #                            "sb_uid":["uid:1","uid:2"]})
    #     fake_writer = FakeWriter()
    #     fake_writer.analyse_state_updates(merged=merged)

    def test_update_master_table(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = fake_writer.up_to_date_sb_table.data.copy()
        fake_writer.master_table.sb_state = "Waiting"
        fake_writer.master_table.sb_state_flag = "ForCalibrator"
        fake_writer.update_master_table()
        for key in ("sb_state_updated","sb_state_flag_updated"):
            assert key not in fake_writer.master_table.columns
        for key in ("sb_state","sb_state_flag"):
            assert not fake_writer.master_table[key].isna().all()
            assert (fake_writer.master_table[key].fillna("")
                     == fake_writer.up_to_date_sb_table.data[key].fillna("")).all()

    def test_get_ToO_selection(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = pd.DataFrame({"code":["2025.1.00001.T",
                                                         "2023.1.00002.S",
                                                         "2026.1.01245.T"]})
        ToO_selection = fake_writer.get_ToO_selection()
        assert list(ToO_selection) == [True,False,True]

    def test_get_waiting_for_calibrator_selection(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = pd.DataFrame({"sb_state":["Ready",
                                                             "Waiting",
                                                             "Waiting"],
                                                 "sb_state_flag":[pd.NA,
                                                                  "ForCalibrator",
                                                                  "ForP2G"]})
        waiting_for_calibrator = fake_writer.get_waiting_for_calibrator_selection()
        assert list(waiting_for_calibrator) == [False,True,False]

    def test_get_ready_selection(self):
        fake_writer = FakeWriter()
        fake_writer.master_table = pd.DataFrame({"sb_state":["Ready",
                                                             "Waiting",
                                                             "Waiting"],
                                                 "sb_state_flag":[pd.NA,
                                                                  "ForCalibrator",
                                                                  "ForP2G"]})
        ready = fake_writer.get_ready_selection()
        assert list(ready) == [True,False,False]

    def test_get_no_hardcoded_selection(self):
        fake_writer = FakeWriter()
        hardcoded_calibrators = ["bandpass (JXY)","",pd.NA,"check (JMiraculix)"]
        fake_writer.master_table = pd.DataFrame({"hardcoded_calibrators":hardcoded_calibrators})
        no_hardcoded = fake_writer.get_no_hardcoded_selection()
        assert list(no_hardcoded) == [False,True,True,False]

    @staticmethod
    def generate_row(template,**kwargs):
        out = template.copy()
        for key,value in kwargs.items():
            out[key] = value
        return out

    def test_get_selection_for_setting_to_Ready(self):
        fake_writer = FakeWriter()
        template = {"code":"2025.1.00001.S",
                    "sb_state":"Waiting",
                    "sb_state_flag":"ForCalibrator",
                    "recommended_state":"Ready",
                    "hardcoded_calibrators":""}
        hardcoded_None = self.generate_row(template=template,hardcoded_calibrators=None)
        ToO = self.generate_row(template=template,code="2013.1.01245.T")
        already_ready = self.generate_row(template=template,sb_state="Ready",
                                          sb_state_flag="")
        waiting_JAO = self.generate_row(template=template,sb_state_flag="WaitingForJAO")
        should_be_waiting = self.generate_row(template=template,recommended_state="Waiting")
        should_be_waitingforp2g = self.generate_row(template=template,recommended_state="WaitingForP2G")
        has_hardcoded = self.generate_row(template=template,
                                          hardcoded_calibrators="Polarization (J12); Bandpass (J45)")
        combination = self.generate_row(template=template,sb_state="Ready",
                                        hardcoded_calibrators="Bandpass (J23)")
        fake_writer.master_table = pd.DataFrame([template,hardcoded_None,ToO,
                                                 already_ready,waiting_JAO,
                                                 should_be_waiting,should_be_waitingforp2g,
                                                 has_hardcoded,combination])
        selection = fake_writer.get_selection_for_setting_to_Ready()
        assert list(selection) == [True,]*2+[False]*7

    def test_get_selection_for_setting_to_WaitingForP2G_for_recommended_Waiting(self):
        #if recommendation is Waiting, then the transition should only happen
        #for SBs in Ready
        fake_writer = FakeWriter()
        template = {"code":"2025.1.00001.S",
                    "sb_state":"Ready",
                    "sb_state_flag":"",
                    "recommended_state":"Waiting"}
        ToO = self.generate_row(template=template,code="2013.1.01245.T")
        already_waiting_for_cal = self.generate_row(template=template,sb_state="Waiting",
                                                    sb_state_flag="ForCalibrator")
        already_waiting_for_p2g = self.generate_row(template=template,sb_state="Waiting",
                                                    sb_state_flag="ForP2G")
        ready_recommended = self.generate_row(template=template,recommended_state="Ready")
        waiting_JAO = self.generate_row(template=template,sb_state="Waiting",
                                        sb_state_flag="WaitingForJAO")
        combination = self.generate_row(template=template,sb_state="Waiting",
                                        sb_state_flag="ForCalibrator",
                                        recommended_state="Ready")
        fake_writer.master_table = pd.DataFrame([template,ToO,already_waiting_for_cal,
                                                 already_waiting_for_p2g,
                                                 ready_recommended,waiting_JAO,
                                                 combination])
        selection = fake_writer.get_selection_for_setting_to_WaitingForP2G()
        assert list(selection) == [True] + [False]*6

    def test_get_selection_for_setting_to_WaitingForP2G_for_recommended_WaitingForP2G(self):
        #if recommendation is WaitingForP2G, then the transition should happen from
        #both Ready and WaitingForCalibrator
        fake_writer = FakeWriter()
        template = {"code":"2025.1.00001.S",
                    "sb_state":"Ready",
                    "sb_state_flag":"",
                    "recommended_state":"WaitingForP2G"}
        waitingforcalibrator = self.generate_row(template=template,sb_state="Waiting",
                                                 sb_state_flag="ForCalibrator")
        waitingforp2g = self.generate_row(template=template,sb_state="Waiting",
                                                 sb_state_flag="ForP2G")
        ToO = self.generate_row(template=template,code="2013.1.01245.T")
        ready_recommended = self.generate_row(template=template,recommended_state="Ready")
        combination = self.generate_row(template=template,sb_state="Waiting",
                                        sb_state_flag="ForCalibrator",
                                        recommended_state="Ready")
        fake_writer.master_table = pd.DataFrame([template,waitingforcalibrator,
                                                 waitingforp2g,ToO,
                                                 ready_recommended,
                                                 combination])
        selection = fake_writer.get_selection_for_setting_to_WaitingForP2G()
        assert list(selection) == [True]*2 + [False]*4

    def test_get_table_base(self):
        class Fake:
            campaign_name = "test"
            get_table_base = Writer.get_table_base
        fake = Fake()
        state_change = {"targetState":"Waiting",
                        "targetSubstate":"ForJAO"}
        base = fake.get_table_base(state_change=state_change)
        expected_base = "set_to_WaitingForJAO_test"
        assert base == expected_base
        state_change = {"targetState":"Ready",
                        "targetSubstate":""}
        base = fake.get_table_base(state_change=state_change)
        expected_base = "set_to_Ready_test"
        assert base == expected_base

    def test_write_state_change_tables(self):
        table = pd.DataFrame({"sb_uid":["uid:1","uid:2"],
                              "sb_state":["Ready","Waiting"]})
        selection = table.sb_state == "Ready"
        state_change = {"selection":selection,
                        "targetState":"Waiting",
                        "targetSubstate":"ForJAO",
                        "comment":"test"}
        class Fake:
            campaign_name = "test"
            master_table = table
            get_table_base = Writer.get_table_base
            write_state_change_tables = Writer.write_state_change_tables
        fake = Fake()
        base = fake.get_table_base(state_change=state_change)
        with tempfile.TemporaryDirectory() as tmpdirname:
            fake.write_state_change_tables(state_change=state_change,
                                           output_dir=tmpdirname)
            human_readable_filepath = tmpdirname / Path(f"{base}.csv")
            assert human_readable_filepath.exists()
            human_readable = pd.read_csv(human_readable_filepath)
            for key in table.columns:
                assert (human_readable[key] == table[selection][key]).all()
            PT_filepath = tmpdirname / Path(f"{base}_ProTrack.csv")
            assert PT_filepath.exists()
            PT = pd.read_csv(PT_filepath,header=None)
            for index,expected_value in zip((0,1,2,3),
                                            ("uid:1",state_change["targetState"],
                                             state_change["targetSubstate"],
                                             state_change["comment"])):
                assert PT.iloc[0,index] == expected_value
            assert not PT.isna().any(axis=None)

    def test_write_tables_for_ProTrack_state_changes(self):
        #just testing if this runs
        table = pd.DataFrame({"sb_uid":["uid:1","uid:2"],
                              "sb_state":["Ready","Waiting"]})
        class Fake:
            campaign_name = "test"
            master_table = table
            get_table_base = Writer.get_table_base
            write_state_change_tables = Writer.write_state_change_tables
            write_tables_for_ProTrack_state_changes = Writer.write_tables_for_ProTrack_state_changes
            def get_selection_for_setting_to_Ready(self):
                return table.sb_state == "Waiting"
            def get_selection_for_setting_to_WaitingForP2G(self):
                return table.sb_state == "Ready"
        fake = Fake()
        with tempfile.TemporaryDirectory() as tmpdirname:
            fake.write_tables_for_ProTrack_state_changes(output_dir=tmpdirname)