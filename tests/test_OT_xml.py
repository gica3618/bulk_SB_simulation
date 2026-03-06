#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 16:29:20 2026

@author: gianni
"""

import sys
sys.path.append("..")
import OT_xml
import os
import pytest
import numpy as np


xml_folder = 'tests/test_xmls'

def get_xml(filename):
    filepath = os.path.join(xml_folder,filename)
    return OT_xml.OT_XML(filepath=filepath)


def test_find_unique_element():
    xml = get_xml(filename="example_polarisation_2023.1.00013.S.xml")
    project_name = xml.find_unique_element(tag="prj:name")
    assert project_name.text == "Orion-KL_a_05_TM1"
    with pytest.raises(ValueError):
        #non-existant element:
        xml.find_unique_element(tag="prj:namehallihallo")
    with pytest.raises(ValueError):
        #non-unique element
        xml.find_unique_element(tag="sbl:FieldSource")

def test_read_allowed_HA():
    #first file gives in units of deg, second in h
    for filename in ("example_polarisation_2023.1.00013.S.xml",
                      "example_polarisation_2023.1.00013.S_allowedHA_hours.xml"):
        xml = get_xml(filename=filename)
        allowed_HA = xml.read_allowed_HA()
        assert np.isclose(allowed_HA["min"].deg, -180, atol=0, rtol=1e-8)
        assert np.isclose(allowed_HA["max"].deg, 180, atol=0, rtol=1e-8)
    xml = get_xml(filename="example_polarisation_2023.1.00013.S_invalid_allowedHA_units.xml")
    with pytest.raises(RuntimeError):
        xml.read_allowed_HA()

def test_requires_TP():
    xml = get_xml(filename="example_cycle10_7m_2023.1.01099.S.xml")
    assert xml.read_RequiresTPAntenna()
    xml = get_xml(filename="example_polarisation_2023.1.00013.S.xml")
    assert not xml.read_RequiresTPAntenna()
    xml = get_xml(filename="example_VLBI_2022.1.01268.V.xml")
    assert xml.read_RequiresTPAntenna() is None
    xml = get_xml(filename="example_cycle10_7m_2023.1.01099.S_invalid_requiresTP.xml")
    with pytest.raises(RuntimeError):
        xml.read_RequiresTPAntenna()

def test_read_coordinates():
    filenames = {"deg":"example_polarisation_2023.1.00013.S.xml",
                 "rad":"example_polarisation_2023.1.00013.S_fieldsourcecoord_rad.xml",
                 "invalid":"example_polarisation_2023.1.00013.S_fieldsourcecoord_invalid.xml",
                 "diff units ra dec":"example_polarisation_2023.1.00013.S_fieldsourcecoord_radecdiff.xml"}
    for ID,filename in filenames.items():
        xml = get_xml(filename=filename)
        field_source = xml.find_unique_element("sbl:FieldSource[@entityPartId='X1535245412']")
        coord_element = field_source.find("sbl:sourceCoordinates",namespaces=xml.namespaces)
        if ID == "deg":
            coord = xml.read_coordinates(coord_element)
            assert coord.ra.deg == 65.81583635
            assert coord.dec.deg == -1.34251820278
        elif ID == "rad":
            coord = xml.read_coordinates(coord_element)
            assert np.isclose(coord.ra.rad,0.32,atol=0,rtol=1e-8)
            assert np.isclose(coord.dec.rad,-0.2,atol=0,rtol=1e-8)
        elif ID == "invalid":
            with pytest.raises(RuntimeError):
                xml.read_coordinates(coord_element)
        elif ID == "diff units ra dec":
            with pytest.raises(ValueError):
                xml.read_coordinates(coord_element)
        else:
            raise RuntimeError

def test_representative_coord():
    xml = get_xml("example_polarisation_2023.1.00013.S.xml")
    rep_coord = xml.get_representative_coordinates()
    assert rep_coord.ra.deg == 83.80885625
    assert rep_coord.dec.deg == -5.376798611111111

def test_mode_name():
    xml = get_xml("example_VLBI_2022.1.01268.V.xml")
    assert xml.read_modeName() == "Standard VLBI"
    xml = get_xml("example_polarisation_2023.1.00013.S.xml")
    assert xml.read_modeName() == "Polarization Interferometry"
    xml = get_xml("example_cycle10_7m_2023.1.01099.S.xml")
    assert xml.read_modeName() == "Standard Interferometry"

def test_nominal_config():
    xml = get_xml("example_cycle10_7m_2023.1.01099.S.xml")
    assert xml.get_nominal_configurations() == ["7M",]
    xml = get_xml("example_polarisation_2023.1.00013.S.xml")
    assert xml.get_nominal_configurations() == ["C43-7","C43-8"]

def test_read_is_query():
    assert OT_xml.OT_XML.read_is_query("true")
    assert not OT_xml.OT_XML.read_is_query("false")
    for unknown_text in ("True","False","rassel",""):
        with pytest.raises(RuntimeError):
            OT_xml.OT_XML.read_is_query(unknown_text)

def get_calibrator(calibrators,cal_type):
    candidates = [c for c in calibrators if c.cal_type == cal_type]
    assert len(candidates) == 1
    return candidates[0]

def test_read_calibrators():
    #simple case with all queries:
    xml = get_xml("example_cycle10_7m_2023.1.01099.S.xml")
    calibrators = xml.read_calibrators_from_observing_groups()
    assert len(calibrators) == 2
    assert sorted([c.name for c in calibrators]) == ["Bandpass","Phase"]
    for c in calibrators:
        assert c.source_name == "query"
        assert not c.is_hardcoded
        assert c.coordinates is None
    #case where bandpass and dgc hardcoded to the same source
    xml = get_xml("example_bandpass_equal_dgc.xml")
    calibrators = xml.read_calibrators_from_observing_groups()
    assert len(calibrators) == 4
    dgc = get_calibrator(calibrators=calibrators, cal_type="DGC")
    bandpass = get_calibrator(calibrators=calibrators, cal_type="Bandpass")
    assert dgc.source_name == bandpass.source_name == "J2253+1608"
    #cases with several pol cals:
    cases = {"2PolCalObsGroup":{"xml":get_xml("example_several_polcal_in_obsgroups.xml"),
                                "npolcal":2,"ncal":5,"PolCals": ["J1326-5256","J1650-5044"]},
             "1PolCalObsGroup":{"xml":get_xml("example_several_polcal_only_one_in_obsgroups.xml"),
                                "npolcal":1,"ncal":4,"PolCals": ["J1650-5044",]}
             }
    for ID,case in cases.items():
        calibrators = case["xml"].read_calibrators_from_observing_groups()
        assert len(calibrators) == case["ncal"]
        query_calibrators = [cal for cal in calibrators if not cal.is_hardcoded]
        assert len(query_calibrators) == case["ncal"]-case["npolcal"]
        query_names = [c.name for c in query_calibrators]
        assert sorted(query_names) == ["Bandpass","Check","Phase"]
        for qcal in query_calibrators:
            assert qcal.source_name == "query"
            assert not qcal.is_hardcoded
            assert qcal.coordinates is None
        hardcoded_calibrators = [cal for cal in calibrators if cal.is_hardcoded]
        assert len(hardcoded_calibrators) == case["npolcal"]
        for hcal in hardcoded_calibrators:
            assert hcal.name == "Polarization calibrator"
        hardcoded_source_names = [c.source_name for c in hardcoded_calibrators]
        assert sorted(hardcoded_source_names) == sorted(case["PolCals"])
        J1650 = [c for c in hardcoded_calibrators if c.source_name == "J1650-5044"][0]
        assert J1650.coordinates.ra.deg == 252.569279575
        assert J1650.coordinates.dec.deg == -50.7467252472
        if case["npolcal"] == 2:
            J1326 = [c for c in hardcoded_calibrators if c.source_name == "J1326-5256"][0]
            assert J1326.coordinates.ra.deg == 201.705121313
            assert J1326.coordinates.dec.deg == -52.939898125

def test_note_to_aod():
    xml = get_xml("example_cycle10_7m_2023.1.01099.S.xml")
    assert xml.get_NotetoAoD() is None
    xml = get_xml("example_NoteToAoD_2023.1.00578.S.xml")
    assert xml.get_NotetoAoD() == "test Note to AoD"