#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 18:39:16 2026

@author: gianni
"""

import logging
from batch_simulations.infrastructure.table_writer import TableWriter


class ProTrackStateChangeWriter:
    #This is a separate class because to write the ProTrack states change files,
    #I want the most up-to-date states of the SBs, so I will use an SB table that is
    #downloaded just before writing the files for ProTrack
    #reason is that during simulations, states of SBs might change
    def __init__(self,campaign_name,master_table,up_to_date_sb_table):
        self.campaign_name = campaign_name
        self.master_table = master_table
        self.up_to_date_sb_table = up_to_date_sb_table
        self.update_master_table()

    def merge_with_up_to_date_sb_table(self):
        logging.info("ProTrackWriter: going to update master table")
        lookup = self.up_to_date_sb_table.data[["sb_uid", "sb_state", "sb_state_flag"]]
        if lookup.sb_uid.duplicated().any():
            raise ValueError("duplicated SB UIDs")
        return self.master_table.merge(lookup,on="sb_uid",how="left",
                                       suffixes=("", "_updated"))

    @staticmethod
    def remove_missing(merged):
        missing = merged["sb_state_updated"].isna()
        logging.info(f"{missing.sum()} SBs are missing in the up-to-date table:")
        missing_uids = merged.loc[missing, "sb_uid"].tolist()
        logging.info(str(missing_uids))
        logging.info("these missing SBs will not be considered for state changes")
        return merged[~missing]

    @staticmethod
    def analyse_number_of_updates(merged):
        old_state = merged["sb_state"] + merged["sb_state_flag"].fillna("")
        updated_state = merged["sb_state_updated"] + merged["sb_state_flag_updated"].fillna("")
        was_updated = old_state != updated_state
        logging.info(f"updated {was_updated.sum()}/{len(old_state)} states")

    # def analyse_state_updates(self,merged):
    #     self.check_missing_sb_state(merged=merged)
    #     self.analyse_number_of_updates(merged=merged)

    def update_master_table(self):
        merged = self.merge_with_up_to_date_sb_table()
        merged = self.remove_missing(merged)
        self.analyse_number_of_updates(merged=merged)
        for key in ("sb_state","sb_state_flag"):
            merged[key] = merged[f"{key}_updated"]
        self.master_table = merged.drop(columns=["sb_state_updated", "sb_state_flag_updated"])

    # def get_current_state(self,sb_uid):
    #     matches = self.sb_table.data[self.sb_table.data["sb_uid"] == sb_uid]
    #     if len(matches) == 0:
    #         logging.warning(f"Could not find entry for SB UID {sb_uid} in the up-to-date SB table")
    #         return None
    #     elif len(matches) == 1:
    #         return matches.iloc[0]["sb_state"], matches.iloc[0]["sb_state_flag"]
    #     else:
    #         raise ValueError(f"Expected 0 or 1 match for SB UID {sb_uid}, found {len(matches)}")

    # def create_updated_master_table(self):
    #     #attention: what should I do if the sb_table does not contain a specific
    #     #sb anymore? at least raise a Warning, no?
    #     new_rows = []
    #     for row in self.itertuples(index=False):
    #         row_dict = row._asdict()
    #         sb_uid = row_dict["sb_uid"]
    #         current_state = self.get_current_state(sb_uid=sb_uid)
    #         if current_state is None:
    #             if self.fail_on_missing_SBs:
    #                 raise RuntimeError("SB UID {sb_uid} missing in the up-to-date table, exiting")
    #             else:
    #                 logging.warning("SB UID {sb_uid} missing in the up-to-date table, will ignore SB")
    #         else:
    #             row_dict["sb_state"] = current_state[0]
    #             row_dict["sb_state_flag"] = current_state[1]
    #             new_rows.append(row_dict)
    #     self.updated_master_table = pd.DataFrame(new_rows)

    def get_ToO_selection(self):
        return self.master_table.code.str.endswith('.T')

    def get_selection_for_setting_to_Waiting(self):
        is_ToO = self.get_ToO_selection()
        #careful, need to put parenthesis:
        return ((~is_ToO) & self.master_table["should_be_Waiting"]
                & (self.master_table["sb_state"] == "Ready"))

    def get_selection_for_setting_to_Ready(self):
        #remember to only set stuff to Ready that does not have
        #hardcoded calibrators
        is_ToO = self.get_ToO_selection()
        return (
                (~is_ToO)
                & (self.master_table["sb_state"] == "Waiting")
                & (self.master_table["sb_state_flag"] == "ForCalibrator")
                & (~self.master_table["should_be_Waiting"])
                #depending on how the master_table was read, empty string might
                #be converted to NA (e.g. with read_csv), so need to cover that:
                & (self.master_table["hardcoded_calibrators"] == "")
                  |(self.master_table["hardcoded_calibrators"].isna())
                )

    def get_table_base(self,state_change):
        return f"set_to_{state_change['targetState']+state_change['targetSubstate']}_{self.campaign_name}"

    def write_state_change_tables(self,state_change,output_dir="."):
        sc = state_change
        base = self.get_table_base(state_change=state_change)
        human_readable = self.master_table[sc["selection"]]
        writer = TableWriter()
        writer.write_table_to_disk(dataframe=human_readable,filename=f"{base}.csv",
                                   output_dir=output_dir)
        writer.write_ProTrack_csv(table=self.master_table,
                                  selection=sc["selection"],
                                  targetState=sc["targetState"],
                                  targetSubstate=sc["targetSubstate"],
                                  comment=sc["comment"],
                                  filename=f"{base}_ProTrack.csv",
                                  output_dir=output_dir)

    def write_tables_for_ProTrack_state_changes(self,output_dir="."):
        state_changes = [{"selection":self.get_selection_for_setting_to_Ready(),
                          "targetState":"Ready",
                          "targetSubstate":"",
                          "comment":"batch simulations: no issues, setting to Ready"},
                         {"selection":self.get_selection_for_setting_to_Waiting(),
                          "targetState":"Waiting",
                          "targetSubstate":"ForP2G",
                          "comment":"batch simulations: simulation failure"}]
        for sc in state_changes:
            self.write_state_change_tables(state_change=sc,output_dir=output_dir)