#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 08:57:35 2024

@author: gianni
"""

import xml.etree.ElementTree as ET
import logging
import sys
from astropy import units as u
from astropy.coordinates import Angle,SkyCoord
import calibrator
from sched_block import SB


class OT_XML():
    namespaces = {'sbl':'Alma/ObsPrep/SchedBlock',
                  'prj':"Alma/ObsPrep/ObsProject",
                  'val':"Alma/ValueTypes"}
    long_lat_keys = {'longitude':'ra','latitude':'dec'}
    unit_map = {"deg": u.deg,
                "h": u.hour,
                "rad": u.rad}

    def __init__(self,filepath):
        tree = ET.parse(filepath)
        self.root = tree.getroot()

    def find_unique_element(self,tag):
        elements = self.root.findall(tag,namespaces=self.namespaces)
        n_elements = len(elements)
        if n_elements != 1:
            raise ValueError(f'found {n_elements} matching elements for {tag}')
        return elements[0]

    def read_allowed_HA(self):
        allowed_HA = {}
        for key in ('minAllowedHA','maxAllowedHA'):
            tag = f'sbl:Preconditions/prj:{key}'
            element = self.find_unique_element(tag)
            HA = float(element.text)
            xml_unit = element.attrib['unit']
            try:
                unit = self.unit_map[xml_unit]
            except KeyError:
                raise RuntimeError(f'unknown unit {xml_unit} for {key}')
            HA = Angle(HA,unit=unit)
            output_key = key[:3]
            allowed_HA[output_key] = HA
        if allowed_HA['min'] >= allowed_HA['max']:
            raise ValueError(f"HA limits defined in OT do not make sense ({allowed_HA})")
        return allowed_HA

    def read_RequiresTPAntenna(self):
        text = self.root.findtext('sbl:SchedulingConstraints/sbl:sbRequiresTPAntennas',
                                  namespaces=self.namespaces)
        if text is None: #not found
            return None
        if text == 'true':
            return True
        elif text == 'false':
            return False
        else:
            raise RuntimeError(f'unknown xml content for sbRequiresTPAntennas: {text}')

    def read_coordinates(self,coord_element):
        coord = {}
        unit_names = []
        for xml_key,output_key in self.long_lat_keys.items():
            element = coord_element.find(f'val:{xml_key}',namespaces=self.namespaces)
            unit_names.append(element.attrib['unit'])
            coord[output_key] = float(element.text)
        if len(set(unit_names)) != 1:
            raise ValueError(f'ra and dec have different units: {unit_names}')
        unit_name = unit_names[0]
        try:
            unit = self.unit_map[unit_name]
        except KeyError:
            raise RuntimeError(f'unknown unit {unit_name}')
        return SkyCoord(ra=coord['ra']*unit,dec=coord['dec']*unit)

    def get_representative_coordinates(self):
        tag = 'sbl:SchedulingConstraints/sbl:representativeCoordinates'
        coord_element = self.find_unique_element(tag)
        return self.read_coordinates(coord_element=coord_element)

    def read_modeName(self):
        return self.find_unique_element('sbl:modeName').text

    def get_nominal_configurations(self):
         configs = self.root.findall('sbl:SchedulingConstraints/sbl:nominalConfiguration',
                                     namespaces=self.namespaces)
         return [c.text for c in configs]

    @staticmethod
    def read_is_query(is_query_text):
        if is_query_text == 'false':
            return False
        elif is_query_text == 'true':
            return True
        else:
            raise RuntimeError(f'unknown xml value for isQuery: {is_query_text}')

    def get_ordered_target_part_ids(self):
        observing_groups = self.root.findall("sbl:ObservingGroup",
                                             namespaces=self.namespaces)
        if len(observing_groups) != 2: #Calibrators and Science
            raise ValueError("expected exactly 2 observing groups")
        ordered_target_partIDs = []
        for obs_group in observing_groups:
            obs_group_name = obs_group.findtext("sbl:name",namespaces=self.namespaces)
            ordered_targets = obs_group.findall("sbl:OrderedTarget",namespaces=self.namespaces)
            logging.info(f"Obs Group {obs_group_name} contains {len(ordered_targets)} ordered targets")
            for ordered_target in ordered_targets:
                target_ref = ordered_target.find("sbl:TargetRef",namespaces=self.namespaces)
                ordered_target_partIDs.append(target_ref.get("partId"))
        return ordered_target_partIDs

    def get_field_source_part_ids(self,ordered_target_partIDs):
        #find partIDs of the corresponding field sources:
        field_source_ref_partIDs = []
        for ordered_target_partID in ordered_target_partIDs:
            target = self.find_unique_element(f"sbl:Target[@entityPartId='{ordered_target_partID}']")
            field_source_ref = target.find("sbl:FieldSourceRef",namespaces=self.namespaces)
            field_source_ref_partIDs.append(field_source_ref.get("partId"))
        #some calibrators (e.g. Pol Cal) appear in both observing groups, so
        #I take the "set"
        return set(field_source_ref_partIDs)

    def build_calibrator_from_field_source(self,field_source_partID):
        field_source = self.find_unique_element(
                           f"sbl:FieldSource[@entityPartId='{field_source_partID}']")
        name = field_source.findtext("sbl:name",namespaces=self.namespaces)
        cal_type = calibrator.classify_calibrator(name=name)
        if cal_type is None:
            return
        source_name = field_source.findtext("sbl:sourceName",namespaces=self.namespaces)
        logging.info(f"found calibrator of type '{cal_type}' ({name}, {source_name})")
        is_query_text = field_source.findtext(
                                       'sbl:isQuery',namespaces=self.namespaces)
        is_query = self.read_is_query(is_query_text=is_query_text)
        is_hardcoded = not is_query
        coord_element = field_source.find('sbl:sourceCoordinates',
                                          namespaces=self.namespaces)
        coordinates = self.read_coordinates(coord_element=coord_element)
        if is_query:
            if coordinates.ra.deg != 0 or coordinates.dec.deg != 0:
                raise ValueError("expected query calibrator of have coordinates RA=0, DEC=0,"
                                 +f"but found RA={coordinates.ra.deg} deg,"
                                 +f" DEC={coordinates.dec.deg} deg")
            coordinates = None
        return calibrator.Calibrator(
                             name=name,source_name=source_name,cal_type=cal_type,
                             is_hardcoded=is_hardcoded,coordinates=coordinates)

    def read_calibrators_from_observing_groups(self):
        #sometimes there are e.g. several PolCals. therefore I use the observing
        #groups to get all calibrators that will actually be observed
        #NOTE: if there are several query calibrators of the same type (e.g. two
        #phase calibrators in query), then the current logic still just adds one calibrator,
        #because the two ordered targets corresponding to the two phase calibrators
        #point to the same field source
        ordered_target_part_ids = self.get_ordered_target_part_ids()
        field_source_part_ids = self.get_field_source_part_ids(
                                  ordered_target_partIDs=ordered_target_part_ids)
        calibrators = []
        for partID in field_source_part_ids:
            cal = self.build_calibrator_from_field_source(field_source_partID=partID)
            if cal:
                calibrators.append(cal)
        return calibrators

    def get_NotetoAoD(self):
        return self.find_unique_element('prj:note').text


class BuildSBFromXML:

    @staticmethod
    def build(xml_filepath):
        xml = OT_XML(xml_filepath)
        metadata = {"note_to_AoD":xml.get_NotetoAoD()}
        sb =  SB(
            metadata=metadata,
            calibrators=xml.read_calibrators_from_observing_groups(),
            mode_name=xml.read_modeName(),
            nominal_configs=xml.get_nominal_configurations(),
            rep_coord=xml.get_representative_coordinates(),
            OT_allowed_HA = xml.read_allowed_HA(),
            requires_TP=xml.read_RequiresTPAntenna())
        sb.consistency_checks()
        return sb


if __name__ == "__main__":
    logging_level = logging.INFO
    #logging_level = logging.ERROR
    logging.basicConfig(format='%(levelname)s: %(message)s',level=logging_level,
                        stream=sys.stdout)

    test_xml = OT_XML(filepath="tests/test_xmls/example_several_polcal_in_obsgroups.xml")
    test_xml.read_calibrators_from_observing_groups()
    
    test_xml = OT_XML(filepath="tests/test_xmls/two_phase_calibrators.xml")
    test_xml.read_calibrators_from_observing_groups()