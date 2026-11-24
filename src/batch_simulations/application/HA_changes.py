#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:51:28 2026

@author: gianni
"""

import xml.etree.ElementTree as ET


class HAChangeWriter:

    def __init__(self,campaign_name,master_table):
        self.campaign_name = campaign_name
        self.master_table = master_table
        self.xml_root = ET.Element("BatchActions")

    def get_HA_change_selection(self):
        HA_recommendations_exist = (self.master_table["HA_min_recommended"].notna()
                                    | self.master_table["HA_max_recommended"].notna())
        return ((self.master_table["recommended_state"]=="Ready")
                & HA_recommendations_exist)

    def add_HA_change(self,project_uid,sb_uid,new_HA_limits):
        edit_sb = ET.SubElement(self.xml_root, "EditSB", projectuid=project_uid)
        sb = ET.SubElement(edit_sb, "sb", uid=sb_uid)
        for lim,lim_value in new_HA_limits.items():
            if lim_value is not None:
                ET.SubElement(sb, f"{lim}ha", angle=str(new_HA_limits[lim]), unit="h")

    def write_xml(self):
        tree = ET.ElementTree(self.root)
        ET.indent(tree, space="  ", level=0)
        tree.write(filename, encoding="UTF-8", xml_declaration=True)