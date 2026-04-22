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
import re
import numpy as np
import logging


sb_table_folder = Path("tests/input_tables")
sb_table_filepaths = {"12m":sb_table_folder / "configuration_lookup_table_cycle_12.csv",
                      "7m":sb_table_folder / "2026-02-06_schedBlockList_7M.csv"}
sb_table_array_config_12m = "c43-5"


class FakeWriter:
    analyse_state_updates = Writer.analyse_state_updates
    get_ToO_selection = Writer.get_ToO_selection
    get_selection_for_setting_to_Waiting = Writer.get_selection_for_setting_to_Waiting
    get_selection_for_setting_to_Ready = Writer.get_selection_for_setting_to_Ready
    merge_with_up_to_date_sb_table = Writer.merge_with_up_to_date_sb_table
    update_master_table = Writer.update_master_table
    @staticmethod
    def check_missing_sb_state(merged):
        return Writer.check_missing_sb_state(merged)
    @staticmethod
    def analyse_number_of_updates(merged):
        return Writer.analyse_number_of_updates(merged)
    def __init__(self):
        #it's important to make a new table each time a new instance is created,
        #because inside some tests, the table is modified inplace
        self.up_to_date_sb_table = SBTable(input_filepaths=sb_table_filepaths,
                                           array_config_12m=sb_table_array_config_12m,
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
        assert (merged["sb_state_flag"].fillna("")
                ==fake_writer.master_table["sb_state_flag"].fillna("")).all()
        expected_updated_state_flags\
                  = fake_writer.up_to_date_sb_table.data["sb_state_flag"]\
                      .iloc[::-1].reset_index(drop=True).fillna("")
        assert (merged["sb_state_flag_updated"].fillna("")
                ==expected_updated_state_flags).all()

    def test_merge_duplicated_uids(self):
        fake_writer = FakeWriter()
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

    def test_check_missing_sb_state(self):
        merged = pd.DataFrame({"sb_state_updated":["Ready","Waiting"],
                               "sb_uid":["uid:1","uid:2"]})
        Writer.check_missing_sb_state(merged)
        merged.at[0,"sb_state_updated"] = pd.NA
        match = re.escape("SB UIDs missing in the up-to-date table: ['uid:1']")
        with pytest.raises(ValueError,match=match):
            Writer.check_missing_sb_state(merged)

    def test_analyse_number_of_updates(self,caplog):
        merged = pd.DataFrame({"sb_state":["Ready","Ready","Waiting","Waiting","Waiting"],
                               "sb_state_flag":[pd.NA,pd.NA,"ForJAO","ForJAO","ForP2G"],
                               "sb_state_updated":["Ready","FullyObserved","Ready","Waiting","Waiting"],
                               "sb_state_flag_updated":[pd.NA,pd.NA,pd.NA,"ForJAO","ForCalibrator"]})
        with caplog.at_level(logging.INFO):
            Writer.analyse_number_of_updates(merged=merged)
        messages = [record.getMessage() for record in caplog.records]
        assert messages == ["updated 3/5 states"]

    def test_analyse_state_updates(self):
        #just testing that it runs
        merged = pd.DataFrame({"sb_state":["Ready","Waiting"],
                               "sb_state_flag":[pd.NA,"ForP2G"],
                               "sb_state_updated":["Ready","Waiting"],
                               "sb_state_flag_updated":[pd.NA,"ForJAO"],
                               "sb_uid":["uid:1","uid:2"]})
        fake_writer = FakeWriter()
        fake_writer.analyse_state_updates(merged=merged)

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

    @staticmethod
    def generate_row(template,**kwargs):
        out = template.copy()
        for key,value in kwargs.items():
            out[key] = value
        return out

    def test_get_selection_for_setting_to_Waiting(self):
        fake_writer = FakeWriter()
        template = {"code":"2025.1.00001.S",
                    "sb_state":"Ready",
                    "should_be_Waiting":True}
        ToO = self.generate_row(template=template,code="2013.1.01245.T")
        already_waiting = self.generate_row(template=template,sb_state="Waiting")
        should_not_be_waiting = self.generate_row(template=template,should_be_Waiting=False)
        combination = self.generate_row(template=template,sb_state="Waiting",
                                        should_be_Waiting=False)
        fake_writer.master_table = pd.DataFrame([template,ToO,already_waiting,
                                                 should_not_be_waiting,
                                                 combination])
        selection = fake_writer.get_selection_for_setting_to_Waiting()
        assert list(selection) == [True,False,False,False,False]

    def test_get_selection_for_setting_to_Ready(self):
        fake_writer = FakeWriter()
        template = {"code":"2025.1.00001.S",
                    "sb_state":"Waiting",
                    "sb_state_flag":"ForCalibrator",
                    "should_be_Waiting":False,
                    "hardcoded_calibrators":""}
        ToO = self.generate_row(template=template,code="2013.1.01245.T")
        already_ready = self.generate_row(template=template,sb_state="Ready")
        waiting_JAO = self.generate_row(template=template,sb_state_flag="WaitingForJAO")
        should_be_waiting = self.generate_row(template=template,should_be_Waiting=True)
        has_hardcoded = self.generate_row(template=template,
                                          hardcoded_calibrators="Polarization; Bandpass")
        combination = self.generate_row(template=template,sb_state="Ready",
                                        hardcoded_calibrators="Bandpass")
        fake_writer.master_table = pd.DataFrame([template,ToO,already_ready,waiting_JAO,
                                                 should_be_waiting,has_hardcoded,
                                                 combination])
        selection = fake_writer.get_selection_for_setting_to_Ready()
        assert list(selection) == [True,]+[False]*6

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
            for index,expected_value in zip((0,1,2,4),
                                            ("uid:1",state_change["targetState"],
                                             state_change["targetSubstate"],
                                             state_change["comment"])):
                assert PT.iloc[0,index] == expected_value
            assert np.isnan(PT.iloc[0,3])

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
            def get_selection_for_setting_to_Waiting(self):
                return table.sb_state == "Ready"
        fake = Fake()
        with tempfile.TemporaryDirectory() as tmpdirname:
            fake.write_tables_for_ProTrack_state_changes(output_dir=tmpdirname)