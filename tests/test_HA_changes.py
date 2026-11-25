#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:57:45 2026

@author: gianni
"""

import pandas as pd
import pytest
import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET
from batch_simulations.application.HA_changes import HAChangeWriter


master_table_data = [#general case
                     {"HA_min":-12, "HA_max":12,
                      "HA_min_recommended":-2,"HA_max_recommended":2.25,
                      "recommended_state":"Ready"},
                     #no recommendation
                     {"HA_min":-12, "HA_max":12,
                      "HA_min_recommended":pd.NA,"HA_max_recommended":pd.NA,
                      "recommended_state":"Ready"},
                     #no max recommendation
                     {"HA_min":-12, "HA_max":12,
                      "HA_min_recommended":-2,"HA_max_recommended":pd.NA,
                      "recommended_state":"Ready"},
                     #no min recommendation
                     {"HA_min":-12, "HA_max":12,
                      "HA_min_recommended":pd.NA,"HA_max_recommended":1.5,
                      "recommended_state":"Ready"},
                     #no update needed
                     {"HA_min":-2, "HA_max":2.25,
                      "HA_min_recommended":-2,"HA_max_recommended":2.25,
                      "recommended_state":"Ready"},
                     #no min update needed
                     {"HA_min":-2, "HA_max":12,
                      "HA_min_recommended":-2,"HA_max_recommended":2.25,
                      "recommended_state":"Ready"},
                     #no max update needed
                     {"HA_min":-12, "HA_max":2.25,
                      "HA_min_recommended":-2,"HA_max_recommended":2.25,
                      "recommended_state":"Ready"},
                     #not recommended for Ready
                     {"HA_min":-12, "HA_max":2.25,
                      "HA_min_recommended":-2,"HA_max_recommended":2.25,
                      "recommended_state":"Waiting"},
                     ]
master_table = pd.DataFrame(master_table_data)
project_uid = "projectuid"
master_table["project_uid"] = project_uid
sb_uid = "sbuid"
master_table["sb_uid"] = sb_uid
campaign_name = "test"

writer = HAChangeWriter(campaign_name=campaign_name,master_table=master_table)
expected_selection = [True,False,True,True,False,True,True,False]

def test_get_selection_for_HA_updates():
    selection = writer.get_selection_for_HA_updates()
    assert list(selection) == expected_selection

def get_unique(element,key):
    out = element.findall(key)
    assert len(out) == 1
    return out[0]

def test_add_HA_limit_changes():
    kwargs = {"project_uid":"project12","sb_uid":"sb23"}
    with pytest.raises(RuntimeError):
        writer.add_HA_limit_changes(**kwargs,changes={})
    writer.initialise_xml()
    change_examples = [{"min":2.3,"max":4},
                       {"min":-1},
                       {"max":0.3}]
    for i,changes in enumerate(change_examples):
        writer.add_HA_limit_changes(**kwargs,changes=changes)
        edit_sb = writer.xml_root.findall("EditSB")
        assert len(edit_sb) == i+1
        edit_sb = edit_sb[i]
        assert edit_sb.get('projectuid') == kwargs["project_uid"]
        sb = get_unique(element=edit_sb, key="sb")
        assert sb.get("uid") == kwargs["sb_uid"]
        for lim in ("min","max"):
            limha = sb.findall(f"{lim}ha")
            if lim in changes:
                assert len(limha) == 1
                limha = limha[0]
                assert limha.get("angle") == str(changes[lim])
                assert limha.get("unit") == "h"
            else:
                assert len(limha) == 0

def check_xml(root):
    edit_sb = root.findall("EditSB")
    assert len(edit_sb) == 5
    for i,(e_sb,row) in enumerate(zip(edit_sb,master_table[expected_selection].itertuples())):
        assert e_sb.get("projectuid") == project_uid
        sb = get_unique(e_sb, "sb")
        assert sb.get("uid") == sb_uid
        if i == 0:
            minha = get_unique(sb, "minha")
            maxha = get_unique(sb, "maxha")
            assert minha.get("angle") == "-2"
            assert maxha.get("angle") == "2.25"
        elif i == 1:
            minha = get_unique(sb, "minha")
            assert minha.get("angle") == "-2"
            assert sb.find("maxha") is None
        elif i == 2:
            maxha = get_unique(sb, "maxha")
            assert sb.find("minha") is None
            assert maxha.get("angle") == "1.5"
        elif i == 3:
            maxha = get_unique(sb, "maxha")
            assert sb.find("minha") is None
            assert maxha.get("angle") == "2.25"
        elif i == 4:
            minha = get_unique(sb, "minha")
            assert minha.get("angle") == "-2"
            assert sb.find("maxha") is None

def test_apply_HA_changes():
    writer.initialise_xml()
    selection = writer.get_selection_for_HA_updates()
    writer.apply_HA_changes(selection=selection)
    check_xml(root=writer.xml_root)

def test_write_xml():
    writer.initialise_xml()
    selection = writer.get_selection_for_HA_updates()
    writer.apply_HA_changes(selection=selection)
    with tempfile.TemporaryDirectory() as tmpdir:
        writer.write_xml(output_dir=tmpdir)
        filename = Path(writer.get_output_base() + ".xml")
        filepath = tmpdir / filename
        tree = ET.parse(filepath)
        root = tree.getroot()
        check_xml(root=root)

def test_write_human_readable_table():
    selection = writer.get_selection_for_HA_updates()
    with tempfile.TemporaryDirectory() as tmpdir:
        writer.write_human_readable_table(selection=selection,output_dir=tmpdir)
        filename = Path(writer.get_output_base() + ".csv")
        filepath = tmpdir / filename
        written_table = pd.read_csv(filepath)
        written_table = written_table.replace({pd.NA: float("nan")}).reset_index(drop=True)
        expected_table = master_table[expected_selection].replace({pd.NA: float("nan")}).reset_index(drop=True)
        pd.testing.assert_frame_equal(
            written_table,
            expected_table,
            check_like=True,
            check_dtype=False)

def test_write_HA_update_outputs():
    #just checking if it runs
    with tempfile.TemporaryDirectory() as tmpdir:
        writer.write_HA_update_outputs(output_dir=tmpdir)