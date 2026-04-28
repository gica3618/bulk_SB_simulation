#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:54:41 2026

@author: gianni
"""

import pandas as pd
from pathlib import Path
import logging
import datetime


class TableWriter:

    column_order = ["code",
                    "sbname",
                    "sb_uid",
                    "p2g_account",
                    "sb_state",
                    "sb_state_flag",
                    "simulated_config",
                    "simulated_date",
                    "simulations",
                    "inspection_reasons",
                    "should_be_Waiting",
                    "traceback_of_unexpected_error",
                    "note_to_AoD",
                    "hardcoded_calibrators",
                    "Phase",
                    "Polarization",
                    "Bandpass",
                    "Amplitude",
                    "Check",
                    "DGC",
                    "min_HA_OT",
                    "max_HA_OT",
                    "min_HA_DSA",
                    "max_HA_DSA",
                    "HA_is_restricted"]
    column_order_lookup = {value:i for i,value in enumerate(column_order)}

    @staticmethod
    def export_dataframe_to_excel_autofit(df,filepath,freeze_panes=None,
                                          sheet_name="Sheet1",
                                          index=False,min_width=8,
                                          max_width=80,extra_padding=2,
                                          adjust_row_height=True):
        #written by ChatGPT
        with pd.ExcelWriter(filepath, engine="xlsxwriter") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=index)
            workbook = writer.book
            worksheet = writer.sheets[sheet_name]    
            wrap_format = workbook.add_format({"text_wrap": True})
            # Iterate over columns
            for col_idx, col_name in enumerate(df.columns):
                series = df[col_name].astype(str)
                max_len = len(str(col_name))  # include header
                for val in series:
                    if val:
                        lines = val.split("\n")
                        longest_line = max(len(line) for line in lines)
                        max_len = max(max_len, longest_line)
                # Clamp width
                width = min(max(max_len + extra_padding, min_width), max_width)    
                worksheet.set_column(first_col=col_idx, last_col=col_idx,
                                     width=width, cell_format=wrap_format)
            # Optionally adjust row heights for multiline text
            if adjust_row_height:
                for row_idx in range(len(df)):
                    max_lines = 1
                    for col in df.columns:
                        val = str(df.iloc[row_idx][col])
                        lines = val.count("\n") + 1
                        max_lines = max(max_lines, lines)
                    # Approximate row height (15 is default Excel row height)
                    #later changed to 14, seems to work better
                    worksheet.set_row(row=row_idx + 1, height=14 * max_lines)

    def sort_columns(self,columns):
        return sorted(columns, key=lambda x: self.column_order_lookup[x])

    def write_table_to_disk(self,dataframe,filename,output_dir="."):
        filepath = Path(output_dir) / filename
        if filepath.exists():
            logging.info(f"{filepath} already exists, I will make a backup")
            timestamp = str(datetime.datetime.now()).replace(" ","_")
            backup_filepath = filepath.with_stem(f"{filepath.stem}_backup_{timestamp}")
            filepath.rename(backup_filepath)
            logging.info(f"made a backup at {backup_filepath}")
        columns = self.sort_columns(dataframe.columns)
        out = dataframe[columns]
        if filepath.suffix == ".csv":
            out.to_csv(filepath,index=False)
        elif filepath.suffix == ".xlsx":
            freeze_panes = (1,2)
            self.export_dataframe_to_excel_autofit(
                     df=out, filepath=filepath,freeze_panes=freeze_panes)
        else:
            raise ValueError(f"invalid output file extension '{filepath.suffix}'")

    @staticmethod
    def write_ProTrack_csv(table,selection,targetState,targetSubstate,comment,
                           filename,output_dir="."):
        out = pd.DataFrame({"sb_uid":table[selection]["sb_uid"]})
        out["targetState"] = targetState
        out["targetSubstate"] = targetSubstate
        out["timestamp"] = ''
        out["comment"] = comment
        filepath = Path(output_dir) / filename
        out.to_csv(filepath,index=False,header=False)