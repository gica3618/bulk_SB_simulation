#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 16:29:20 2026

@author: gianni
"""


from pathlib import Path
import pytest
import numpy as np
from batch_simulations.infrastructure.OT_xml import OT_XML,BuildSBFromXML
from astropy.coordinates import SkyCoord
from astropy import units as u


xml_folder = Path('tests/xmls')

def get_xml(filename):
    filepath = xml_folder / filename
    return OT_XML(filepath=filepath)


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
    assert OT_XML.read_is_query("true")
    assert not OT_XML.read_is_query("false")
    for unknown_text in ("True","False","rassel",""):
        with pytest.raises(RuntimeError):
            OT_XML.read_is_query(unknown_text)

def test_get_part_ids():
    xml = get_xml("example_polarisation_2023.1.00013.S.xml")
    target_part_IDs = xml.get_ordered_target_part_ids()
    assert target_part_IDs == ["X1307262861","X243929360","X793454005","X784838089",
                               "X702136351","X822623979","X1531303282","X1782759204",
                               "X22837845"]
    field_source_part_ids = xml.get_field_source_part_ids(
                                                ordered_target_partIDs=target_part_IDs)
    assert field_source_part_ids == set(["X1503178100","X1535245412","X1646156811",
                                         "X1201643612","X1016110489","X1273708068",
                                         "X77481490","X603351316","X1201643612"])
    xml = get_xml("2025.1.01279.S_general_SB.xml")
    target_part_IDs = xml.get_ordered_target_part_ids()
    assert target_part_IDs == ["X650245099","X1212077582","X1386938427","X237953864",
                               "X103911417","X2147388034",]
    field_source_part_ids = xml.get_field_source_part_ids(
                                                ordered_target_partIDs=target_part_IDs)
    assert field_source_part_ids == set(["X2119741686","X504760990","X569987043",
                                         "X1440250739","X1619819355","X720613460"])

def test_build_calibrator_from_field_source():
    xml = get_xml("2025.1.01279.S_general_SB.xml")
    cal = xml.build_calibrator_from_field_source(field_source_partID="X2119741686")
    assert cal is None
    cal = xml.build_calibrator_from_field_source(field_source_partID="X504760990")
    assert cal.name == "Bandpass"
    assert cal.source_name == "query"
    assert cal.cal_type == "Bandpass"
    assert cal.is_hardcoded == False
    xml = get_xml("example_polarisation_2023.1.00013.S.xml")
    cal = xml.build_calibrator_from_field_source(field_source_partID="X1201643612")
    assert cal.name == "Polarization calibrator"
    assert cal.source_name == "J0522-3627"
    assert cal.cal_type == "Polarization"
    assert cal.is_hardcoded == True
    assert cal.coordinates == SkyCoord(ra=80.7416026833*u.deg,dec=-36.4585697806*u.deg)


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

def test_sb_builder():
    filepath = xml_folder / "2025.1.01279.S_general_SB.xml"
    BuildSBFromXML.build(filepath)