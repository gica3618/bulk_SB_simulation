#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:51:28 2026

@author: gianni
"""

import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path
import logging
import numpy as np
from batch_simulations.infrastructure.table_writer import TableWriter


class HAChangeWriter:

    def __init__(self,campaign_name,master_table):
        self.campaign_name = campaign_name
        self.master_table = master_table
        self.initialise_xml()

    def initialise_xml(self):
        self.xml_root = ET.Element("BatchActions")

    @staticmethod
    def HA_equal(HA1,HA2):
        return np.isclose(HA1,HA2,atol=0.01,rtol=1e-3)

    def get_selection_for_HA_updates(self):
        table = self.master_table
        needs_update = {}
        for lim in ("min","max"):
            #pd.NA cannot be compared in boolean way, so have to convert:
            HA_lim_recommended = table[f"HA_{lim}_recommended"].to_numpy(dtype=float, na_value=np.nan)
            HA_lim = table[f"HA_{lim}"].to_numpy(dtype=float, na_value=np.nan)
            is_equal = self.HA_equal(HA_lim, HA_lim_recommended)
            is_different = np.logical_not(is_equal)
            needs_update[lim] = (table[f"HA_{lim}_recommended"].notna()
                                & is_different)
        ready_recommended = table["recommended_state"] == "Ready"
        return ready_recommended & (needs_update["min"] | needs_update["max"])

    def add_HA_limit_changes(self,project_uid,sb_uid,changes):
        if not changes:
            raise RuntimeError("empty changes")
        edit_sb = ET.SubElement(self.xml_root, "EditSB", projectuid=project_uid)
        sb = ET.SubElement(edit_sb, "sb", uid=sb_uid)
        for lim, value in changes.items():
            ET.SubElement(sb, f"{lim}ha", angle=str(value), unit="h")

    def apply_HA_changes(self,selection):
        #re-initialise xml:
        self.initialise_xml()
        for row in self.master_table[selection].itertuples():
            changes = {}
            for lim in ("min", "max"):
                recommended_HA_lim = getattr(row, f"HA_{lim}_recommended")
                current_HA_lim = getattr(row, f"HA_{lim}")
                if (pd.notna(recommended_HA_lim)
                         and not self.HA_equal(current_HA_lim, recommended_HA_lim)):
                    changes[lim] = recommended_HA_lim
            if changes:
                self.add_HA_limit_changes(project_uid=row.project_uid,
                                          sb_uid=row.sb_uid,changes=changes)

    def get_output_base(self):
        return f"HA_updates_{self.campaign_name}"

    def write_xml(self,output_dir):
        tree = ET.ElementTree(self.xml_root)
        ET.indent(tree, space="  ", level=0)
        filename = self.get_output_base() + ".xml"
        filepath = Path(output_dir) / filename
        tree.write(filepath, encoding="UTF-8", xml_declaration=True)

    def write_human_readable_table(self,selection,output_dir):
        human_readable = self.master_table[selection]
        writer = TableWriter()
        filename = self.get_output_base() + ".csv"
        writer.write_table_to_disk(dataframe=human_readable,filename=filename,
                                   output_dir=output_dir)

    def write_HA_update_outputs(self,output_dir="."):
        logging.info("going to write HA change tables")
        selection = self.get_selection_for_HA_updates()
        self.apply_HA_changes(selection=selection)
        self.write_xml(output_dir=output_dir)
        self.write_human_readable_table(selection=selection,output_dir=output_dir)