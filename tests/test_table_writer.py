#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:56:05 2026

@author: gianni
"""

from batch_simulations.infrastructure.table_writer import TableWriter
import pandas as pd
import tempfile
from pathlib import Path
import pytest
import glob


class TestTableWriter:

    def test_sort_columns(self):
        writer = TableWriter()
        input_columns = ["sb_uid","code","min_HA_OT","Phase"]
        expected_output_columns = ["code","sb_uid","Phase","min_HA_OT"]
        output_columns = writer.sort_columns(input_columns)
        assert output_columns == expected_output_columns

    def test_write_table_to_disk(self):
        data = {'code': ['2025.1.00001.S', '2023.1.00123.S'],
                'p2g_account': ["ronaldo", "kaka"]}
        df = pd.DataFrame(data)
        writer = TableWriter()
        with tempfile.TemporaryDirectory() as tmpdirname:
            for extension in ("csv","xlsx"):
                filename = Path(f"test.{extension}")
                writer.write_table_to_disk(dataframe=df,filename=filename,
                                           output_dir=tmpdirname)
                assert (tmpdirname/filename).exists()
            filename = Path("test.pdf")
            with pytest.raises(ValueError):
                writer.write_table_to_disk(dataframe=df,filename=filename,
                                           output_dir=tmpdirname)
        #now test the overwrite protection:
        with tempfile.TemporaryDirectory() as tmpdirname:
            for extension in ("csv","xlsx"):
                filename = Path(f"test.{extension}")
                filepath = tmpdirname / filename
                filepath.touch()
                writer.write_table_to_disk(dataframe=df,filename=filename,
                                           output_dir=tmpdirname)
                backup_filepath = filepath.with_stem(f"{filename.stem}_backup*")
                assert len(glob.glob(str(backup_filepath))) == 1

    def test_write_protrack_csv(self):
        writer = TableWriter()
        with tempfile.TemporaryDirectory() as tmpdirname:
            filename = Path("test.csv")
            master_table_data = {column: ['test',]*4 for column in
                                 TableWriter.column_order}
            table = pd.DataFrame(master_table_data)
            selection = [True,]*len(table)
            selection[0] = False
            targetState = "Waiting"
            targetSubstate = "ForJAO"
            comment = "Asterix"
            writer.write_ProTrack_csv(
                  table=table,selection=selection,targetState=targetState,
                  targetSubstate=targetSubstate,comment=comment,
                  filename=filename,output_dir=tmpdirname)
            filepath = tmpdirname / filename
            assert filepath.exists()
            written = pd.read_csv(filepath,header=None)
            assert len(written.columns) == 5
            #I set one selection to False, so I know the expected length:
            assert len(written) == len(table)-1
            assert all(written.iloc[:,1]==targetState)
            assert all(written.iloc[:,2]==targetSubstate)
            assert all(written.iloc[:,3].isna()) #timestamp
            assert all(written.iloc[:,4]==comment)