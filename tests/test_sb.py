#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 15:13:23 2026

@author: gianni
"""

from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.dsa_ha_policy import DSAHourAnglePolicy
import pytest
from pathlib import Path
from scipy import constants


class TestSB:

    @staticmethod
    def generate_sb(filename):
        xml_folder = Path('tests/xmls')
        filepath = xml_folder / filename
        return BuildSBFromXML.build(filepath)

    @staticmethod
    def assert_is_not_XX(sb,exluded):
        for is_XX in ("is_7m","is_VLBI","is_Polarisation","is_B2B","is_Solar"):
            if is_XX in exluded:
                continue
            assert not getattr(sb, is_XX)

    def test_is_polarisation(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        assert sb.is_Polarisation
        self.assert_is_not_XX(sb=sb, exluded=["is_Polarisation",])

    def test_is_vlbi(self):
        sb = self.generate_sb("example_VLBI_2022.1.01268.V.xml")
        assert sb.is_VLBI
        self.assert_is_not_XX(sb=sb, exluded=["is_VLBI",])

    def test_solar(self):
        sb = self.generate_sb("example_solar_2022.1.01544.S.xml")
        assert sb.is_Solar
        self.assert_is_not_XX(sb=sb, exluded=["is_Solar",])

    def test_is_b2b(self):
        sb = self.generate_sb("2025.1.01389.S_B2B.xml")
        assert sb.is_B2B
        self.assert_is_not_XX(sb=sb, exluded=["is_B2B",])

    def test_is_7m(self):
        sb = self.generate_sb("example_NoteToAoD_2023.1.00578.S.xml")
        assert sb.is_7m
        self.assert_is_not_XX(sb=sb, exluded=["is_7m",])

    def test_general_init(self):
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        self.assert_is_not_XX(sb=sb,exluded=[])

    def test_consistency_checks(self):
        self.generate_sb("2025.1.01279.S_general_SB.xml")
        self.generate_sb("example_several_polcal_only_one_in_obsgroups.xml")
        with pytest.raises(ValueError):
            self.generate_sb("example_several_polcal_in_obsgroups.xml")
        with pytest.raises(ValueError):
            self.generate_sb("Polarisation_without_PolCal.xml")
        with pytest.raises(ValueError):
            self.generate_sb("Polarisation_PolCal_not_hardcoded.xml")
        with pytest.raises(ValueError):
            self.generate_sb("NonPolarisation_with_PolCal.xml")

    def test_has_no_PolCal(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        assert not sb.has_no_PolCal()
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        assert sb.has_no_PolCal()

    def test_has_at_least_one_pol_cal(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        assert sb.has_at_least_one_PolCal()
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        assert not sb.has_at_least_one_PolCal()

    def test_get_calibrator(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        pol_cal = sb.get_calibrator("Polarization")
        assert pol_cal.name == "Polarization calibrator"
        assert pol_cal.source_name == "J0522-3627"
        bandpass = sb.get_calibrator("Bandpass")
        assert bandpass.name == "Bandpass"
        assert bandpass.source_name == "J0538-4405"
        phase = sb.get_calibrator("Phase")
        assert phase.name == "Phase"
        assert phase.source_name == "query"
        #attempting to get PolCal from non-polarization SB, should fail:
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        with pytest.raises(RuntimeError):
            sb.get_calibrator("Polarization")

    def test_get_PolCal(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        pol_cal = sb.get_PolCal()
        expected_pol_cal = sb.get_calibrator("Polarization")
        for field,value in vars(pol_cal).items(): 
            assert value == getattr(expected_pol_cal,field)
        #attempting to get PolCal from non-polarization SB, should fail:
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        with pytest.raises(RuntimeError):
            sb.get_PolCal()

    def test_PolCal_is_hardcoded(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        assert sb.PolCal_is_hardcoded()
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        with pytest.raises(RuntimeError):
            sb.PolCal_is_hardcoded()

    def test_any_calibrator_is_harcoded(self):
        sb = self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        assert sb.any_calibrator_hardcoded()
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        assert not sb.any_calibrator_hardcoded()

    def test_DSA_HA_limits(self):
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        assert sb.get_DSA_HA_limits() == DSAHourAnglePolicy.compute_min_max_HA(sb)

    def test_add_metadata(self):
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        key,value = "test_key",12.34
        sb.add_metadata(key=key,value=value)
        assert sb.metadata[key] == value

    def test_single_execution_time(self):
        sb = self.generate_sb("single_execution_SB.xml")
        assert sb.single_execution_time() == 34.793*constants.minute
        sb = self.generate_sb("eight_executions_SB.xml")
        assert sb.single_execution_time() == 10.24982222*constants.hour / 8
        sb = self.generate_sb("Polarisation_1session.xml")
        assert sb.single_execution_time() == 3.022894444*constants.hour
        sb = self.generate_sb("Polarisation_several_sessions.xml")
        assert sb.single_execution_time() == 3.993038889*constants.hour