#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:57:53 2026

@author: gianni
"""

import os
import logging
import subprocess
from pathlib import Path


class XMLDownloader:

    def __init__(self,project_code,sb_name,filepath):
        self.project_code = project_code
        self.sb_name = sb_name
        self.filepath = filepath

    def build_getsb_command(self):
        #I took this code directly from simulateSB.py and simplified it
        scriptGetSB = "/groups/science/scripts/P2G/getsb/getsb.py"
        if not os.path.isfile(scriptGetSB):
            scriptGetSB = "/users/ahirota/AIV/science/scripts/P2G/getsb/getsb.py"
        try:
            import cx_Oracle
            serverName = None
        except:
            logging.info("cx_Oracle is not available, and thus will run getsb.py"+
                         " on red-osf")
            serverName = "red-osf.osf.alma.cl"
        if serverName:
            userName = os.getenv("USER")
            cmd = ["ssh"]
            cmd.append(f"{userName}@{serverName}")
            scriptName = "PYTHONPATH=/users/ahirota/local/lib64/python2.6/"\
                          +"site-packages/cx_Oracle-5.2.1-py2.6-linux-x86_64.egg"\
                          +f":$PYTHONPATH {scriptGetSB}"
            cmd.append(f"{scriptName} -p '{self.project_code}' -s '{self.sb_name}'")
            cmd.append("-S ora.sco.alma.cl:1521/ONLINE.SCO.CL")
        else:
            cmd = [scriptGetSB]
            cmd.extend(["-p", self.project_code, "-s", self.sb_name])
            cmd.extend(["-S", "ora.sco.alma.cl:1521/ONLINE.SCO.CL"])
        return cmd

    def download_xml(self):
        cmd = self.build_getsb_command()
        logging.info("# Retrieving SB xml with the following command [%s]"\
                  % (" ".join(cmd)))
        process = subprocess.run(cmd,capture_output=True, text=True, check=True)
        xml_str = process.stdout
        Path(self.filepath).write_text(xml_str)