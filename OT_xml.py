#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 08:57:35 2024

@author: gianni
"""

import xml.etree.ElementTree as ET
import os
import logging
import subprocess
import sys
from astropy import units as u
from astropy.coordinates import Angle,SkyCoord
import calibrator


def download_xml(project_code,SB,filepath):
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
        cmd.append(f"{scriptName} -p '{project_code}' -s '{SB}'")
        cmd.append("-S ora.sco.alma.cl:1521/ONLINE.SCO.CL")
    else:
        cmd = [scriptGetSB]
        cmd.extend(["-p", project_code, "-s", SB])
        cmd.extend(["-S", "ora.sco.alma.cl:1521/ONLINE.SCO.CL"])
    logging.info("# Retrieving SB xml with the following command [%s]"\
                  % (" ".join(cmd)))
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,env=None)
    xml_str, _ = p.communicate()
    if hasattr(sys.stdout,"detach"):
        # Py3
        xml_str = xml_str.decode("utf-8")
    with open(filepath, "w") as f:
        f.write(xml_str)


class OT_XML_Reader():
    namespaces = {'sbl':'Alma/ObsPrep/SchedBlock',
                  'prj':"Alma/ObsPrep/ObsProject",
                  'val':"Alma/ValueTypes"}
    long_lat_keys = {'longitude':'ra','latitude':'dec'}

    def __init__(self,filepath):
        tree = ET.parse(filepath)
        self.root = tree.getroot()

    def find_unique_element(self,tag):
        elements = self.root.findall(tag,namespaces=self.namespaces)
        n_elements = len(elements)
        assert n_elements == 1, f'found {n_elements} matching elements for {tag}'
        return elements[0]

    def read_allowed_HA(self):
        allowed_HA = {}
        for key in ('minAllowedHA','maxAllowedHA'):
            tag = f'sbl:Preconditions/prj:{key}'
            element = self.find_unique_element(tag)
            HA = float(element.text)
            xml_unit = element.attrib['unit']
            if xml_unit == 'deg':
                unit = u.deg
            elif xml_unit == 'h':
                unit = u.hour
            else:
                raise RuntimeError(f'unknown unit {unit} for {key}')
            HA = Angle(HA,unit=unit)
            output_key = key[:3]
            allowed_HA[output_key] = HA
        assert allowed_HA['min'] < allowed_HA['max']
        return allowed_HA

    def read_RequiresTPAntenna(self):
        text = self.root.findtext('sbl:SchedulingConstraints/sbl:sbRequiresTPAntennas',
                                  namespaces=self.namespaces)
        if text is None:
            return None
        if text == 'true':
            return True
        elif text == 'false':
            return False
        else:
            raise RuntimeError(f'unknown xml content: {text}')

    def read_coordinates(self,coord_element):
        coord = {}
        for xml_key,output_key in self.long_lat_keys.items():
            element = coord_element.find(f'val:{xml_key}',namespaces=self.namespaces)
            assert element.attrib['unit'] == 'deg'
            coord[output_key] = float(element.text)
        return SkyCoord(ra=coord['ra']*u.deg,dec=coord['dec']*u.deg)

    def get_representative_coordinates(self):
        tag = 'sbl:SchedulingConstraints/sbl:representativeCoordinates'
        coord_element = self.find_unique_element(tag)
        return self.read_coordinates(coord_element=coord_element)

    def read_modeName(self):
        return self.root.findtext('sbl:modeName',namespaces=self.namespaces)

    def get_nominal_config(self):
        return self.root.findtext('sbl:SchedulingConstraints/sbl:nominalConfiguration',
                                  namespaces=self.namespaces)

    def read_is_query(self,is_query_text):
        if is_query_text == 'false':
            return False
        elif is_query_text == 'true':
            return True
        else:
            raise RuntimeError(f'unknown xml value for isQuery: {is_query_text}')

    def read_calibrators(self):
        calibrators = []
        field_sources = self.root.findall("sbl:FieldSource",namespaces=self.namespaces)
        for field_source in field_sources:
            name = field_source.findtext('sbl:name',namespaces=self.namespaces)
            candidate_cal_types = []
            for cal_type,keywords in calibrator.Calibrator.calibrator_keywords.items():
                if any([keyword in name for keyword in keywords]):
                    candidate_cal_types.append(cal_type)
            if len(candidate_cal_types) == 0:
                logging.info(f'Source name "{name}" is not a calibrator')
                continue
            else:
                assert len(candidate_cal_types) == 1
                is_query_text = field_source.findtext(
                                      'sbl:isQuery',namespaces=self.namespaces)
                is_query = self.read_is_query(is_query_text=is_query_text)
                source_name = field_source.findtext('sbl:sourceName',
                                                    namespaces=self.namespaces)
                coord_element = field_source.find('sbl:sourceCoordinates',
                                                  namespaces=self.namespaces)
                coordinates = self.read_coordinates(coord_element=coord_element)
                cal = calibrator.Calibrator(
                             name=name,source_name=source_name,cal_type=cal_type,
                             is_query=is_query,coordinates=coordinates)
                calibrators.append(cal)
        return calibrators

    def get_NotetoAoD(self):
        return self.root.findtext('prj:note',namespaces=self.namespaces)