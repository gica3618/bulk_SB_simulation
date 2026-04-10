#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:01:18 2026

@author: gianni
"""

import pandas as pd
import logging


class SBTable:

    input_columns = ["sb_uid", "code", "p2g_account", "sbname", "sb_state",
                     "sb_state_flag"]
    excluded_states_12m = ('FullyObserved','ObservingTimedOut')
    #TODO filter out VLBI and solar? Or are they anyway not in the 12m lookup
    #table and the 7M table?

    def __init__(self,input_filepaths,array_config_12m,SB_filter):
        self.input_filepaths = input_filepaths
        self.array_config_12m = array_config_12m
        self.SB_filter = SB_filter
        self.read_input()

    def determine_12m_array_config_number(self):
        #assuming that the config is something like c43-2 or c43-10
        antennas,config_number = self.array_config_12m.split('-')
        if antennas != 'c43':
            raise ValueError(f"invalid antennas string: {antennas}")
        config_number = int(config_number)
        if config_number not in range(1,11):
            raise ValueError(f"invalid config number: {config_number}")
        return config_number

    def read_input(self):
        input_data_sets = []
        if self.input_filepaths["12m"] is not None:
            input_data_sets.append(self.read_12m_input())
        if self.input_filepaths["7m"] is not None:
            input_data_sets.append(self.read_7m_input())
        data = pd.concat(input_data_sets,ignore_index=True)
        logging.info(f"in total, {len(data)} input SBs")
        if self.SB_filter is not None:
            data = data[data.apply(self.SB_filter, axis=1)]
            logging.info(f"after filtering, {len(data)} input SBs left")
        self.data = data

    def read_12m_input(self):
        filepath = self.input_filepaths["12m"]
        data = pd.read_csv(filepath)
        logging.info(f"a total of {len(data)} 12m SBs")
        array_config_number = self.determine_12m_array_config_number()
        array_selection = data[f"selected_c{array_config_number}"]
        data = data[array_selection]
        logging.info(f"after array selection, {len(data)} SBs are left")
        state_selection = ~data['sb_state'].isin(self.excluded_states_12m)
        data = data[state_selection]
        logging.info(f"after filtering out SB in {self.excluded_states_12m}, {len(data)} SBs are left")
        data = data[self.input_columns]
        return data

    def read_7m_input(self):
        filepath = self.input_filepaths["7m"]
        data = pd.read_csv(filepath)
        logging.info(f"a total of {len(data)} 7m SBs")
        column_rename_mapper = {"SB UID":"sb_uid",
                                "Project Code":"code",
                                "SB Name":"sbname",
                                "P2G":"p2g_account",
                                "State":"sb_state",
                                "SubState":"sb_state_flag",
                                }
        data.rename(mapper=column_rename_mapper,inplace=True,axis="columns")
        data = data[self.input_columns]
        return data