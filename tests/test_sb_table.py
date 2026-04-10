#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 19:37:21 2026

@author: gianni
"""

from batch_simulations.infrastructure.sb_table import SBTable
from pathlib import Path
import pytest
import pandas as pd


class TestSBTable:

    input_folder = Path("tests/input_tables")
    input_filepaths = {"12m":input_folder / "configuration_lookup_table_cycle_12.csv",
                       "7m":input_folder / "2026-02-06_schedBlockList_7M.csv"}
    data12m = pd.read_csv(input_filepaths["12m"])
    data7m = pd.read_csv(input_filepaths["7m"])

    def test_read_12m_array_config_number(self):
        class FakeSBTable:
            def __init__(self,array_config_12m):
                self.array_config_12m = array_config_12m
            def determine_12m_array_config_number(self):
                return SBTable.determine_12m_array_config_number(self)
        for array_config_12m,number in zip(("c43-3","c43-10"),(3,10)):
            fake_table = FakeSBTable(array_config_12m=array_config_12m)
            assert fake_table.determine_12m_array_config_number() == number
        for array_config_12m in ("c-1","c43-0","c43-11","c42-5"):
            fake_table = FakeSBTable(array_config_12m=array_config_12m)
            with pytest.raises(ValueError):
                fake_table.determine_12m_array_config_number()

    def test_read_12m_input(self):
        input_filepaths = {"12m":self.input_filepaths["12m"],"7m":None}
        table = SBTable(input_filepaths=input_filepaths,array_config_12m="c43-5",
                        SB_filter=None)
        data_selection = ~self.data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))\
                         & self.data12m["selected_c5"]
        read_data = table.read_12m_input()
        assert len(read_data) == len(self.data12m[data_selection])
        assert sorted(read_data["sb_uid"])  == sorted(self.data12m[data_selection]["sb_uid"])
        assert sorted(read_data.columns) == sorted(SBTable.input_columns)
        #explicit check using the manually (with Excel Filter) extracted uids:
        uids = pd.read_csv(self.input_folder / "c5_uids.csv")
        assert sorted(read_data["sb_uid"]) == sorted(uids["sb_uid"])

    def test_read_7m_input(self):
        input_filepaths = {"12m":None,"7m":self.input_filepaths["7m"]}
        table = SBTable(input_filepaths=input_filepaths, array_config_12m=None,
                        SB_filter=None)
        read_data = table.read_7m_input()
        assert len(read_data) == len(self.data7m)
        assert sorted(read_data.columns) == sorted(SBTable.input_columns)

    def test_read_input(self):
        table = SBTable(input_filepaths=self.input_filepaths,array_config_12m="c43-4",
                        SB_filter=None)
        data_selection_12m = self.data12m["selected_c4"]\
                & ~self.data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        read_data = table.data
        assert len(read_data) == len(self.data7m) + len(self.data12m[data_selection_12m])
        assert sorted(read_data["sb_uid"])\
               == sorted(list(self.data12m[data_selection_12m]["sb_uid"])
                         + list(self.data7m["SB UID"]))
        assert sorted(read_data.columns) == sorted(SBTable.input_columns)

    def test_filter_p2g(self):
        p2g = "gianni"
        def p2g_filter(row):
            return row["p2g_account"] == p2g
        table = SBTable(input_filepaths=self.input_filepaths,array_config_12m="c43-4",
                        SB_filter=p2g_filter)
        data_selection_12m = self.data12m["selected_c4"] & (self.data12m["p2g_account"]==p2g)\
                            &  ~self.data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        data12m_filtered = self.data12m[data_selection_12m]
        assert len(data12m_filtered) > 0
        data7m_filtered = self.data7m[self.data7m["P2G"]==p2g]
        assert len(data7m_filtered) > 0
        assert len(table.data) == len(data12m_filtered) + len(data7m_filtered)
        assert sorted(table.data["sb_uid"])\
               == sorted(list(data12m_filtered["sb_uid"]) + list(data7m_filtered["SB UID"]))

    def test_filter_project_code(self):
        code = "2025.1.00240.S"
        def project_filter(row):
            return row["code"] == code
        table = SBTable(input_filepaths=self.input_filepaths,
                       array_config_12m="c43-4",SB_filter=project_filter)
        data_selection_12m = self.data12m["selected_c4"] &  (self.data12m["code"]==code)\
                            &  ~self.data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        data12m_filtered = self.data12m[data_selection_12m]
        assert len(data12m_filtered) == 0
        data7m_filtered = self.data7m[self.data7m["Project Code"]==code]
        assert len(data7m_filtered) > 0
        assert len(table.data) == len(data12m_filtered) + len(data7m_filtered)
        assert sorted(table.data["sb_uid"])\
               == sorted(list(data12m_filtered["sb_uid"]) + list(data7m_filtered["SB UID"]))