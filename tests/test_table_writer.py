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
        input_columns = ["sb_uid","code","HA_min","query_calibrators"]
        expected_output_columns = ["code","sb_uid","query_calibrators","HA_min"]
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