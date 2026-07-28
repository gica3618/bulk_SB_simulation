#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:09:03 2026

@author: gianni
"""

from batch_simulations.domain.dsa_ha_policy import DSAHourAnglePolicy
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from pathlib import Path
import numpy as np
import pytest
from astropy.coordinates import SkyCoord,Angle
from astropy import units as u
from scipy import constants


class TestDSAPolicy:

    @staticmethod
    def generate_sb(filename):
        xml_folder = Path('tests/xmls')
        filepath = xml_folder / filename
        return BuildSBFromXML.build(filepath)

    def test_general_north(self):
        sb = self.generate_sb("2025.1.01279.S_general_SB.xml")
        assert sb.rep_coord.dec.deg == 42.19541666666667
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert HA_limits["min"].hour == -3
        assert HA_limits["max"].hour == 2

    def test_general_south(self):
        sb = self.generate_sb("general_DEC-15.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -15.39436944
        assert HA_limits["min"].hour == -4
        assert HA_limits["max"].hour == 3

    def test_leading_pol_cal(self):
        #example with leading PolCal:
        sb =self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -5.376798611111111
        assert sb.rep_coord.ra.deg == 83.80885625
        assert sb.get_PolCal().coordinates.ra.deg == 80.7416026833
        #this is a case where the min needs to be adjusted, since without PolCal,
        #it would be -4
        assert HA_limits["min"].hour == -3 - (83.80885625-80.7416026833)/360*24
        assert HA_limits["max"].hour == 3

    def test_leading_pol_cal_RA0(self):
        #example with leading PolCal around RA=0
        sb = self.generate_sb("Polarisation_trailing_PolCal_around_RA0deg.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 1
        assert sb.get_PolCal().coordinates.ra.deg == 359.471941933
        delta_ra = (360+1)-359.471941933
        assert np.isclose(HA_limits["min"].hour, -3 - delta_ra/360*24, atol=0, rtol=1e-8)
        assert HA_limits["max"].hour == 3

    def test_trailing_pol_cal(self):
        #example with trailing PolCal
        sb = self.generate_sb("Polarisation_leading_PolCal.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 165.0
        assert sb.get_PolCal().coordinates.ra.deg == 165.814104429
        assert HA_limits["min"].hour == -3 + (165.814104429-165)/360*24
        assert HA_limits["max"].hour == 3

    def test_trailing_pol_cal_RA0(self):
        #example with trailing PolCal around RA=0
        sb = self.generate_sb("Polarisatoin_trailing_PolCal_around_RA0.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 355.0
        assert sb.get_PolCal().coordinates.ra.deg == 1.14856453333
        delta_ra = (360+1.14856453333) - 355
        assert HA_limits["min"].hour == -3 + delta_ra/360*24
        assert HA_limits["max"].hour == 3

    def test_leading_pol_cal_high_dec(self):
        #leading PolCal, no adjustement needed (anyway min is -3 because of high DEC)
        sb = self.generate_sb("Polarisation_leading_PolCal_highDEC.xml")
        HA_limits = DSAHourAnglePolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == 0.5238108333333333
        assert sb.rep_coord.ra.deg == 150
        assert sb.get_PolCal().coordinates.ra.deg == 139.933497596
        assert HA_limits["min"].hour == -3
        assert HA_limits["max"].hour == 2

    def test_pol_cal_leading_too_much(self):
        #PolCal leading too much
        with pytest.raises(ValueError):
            sb = self.generate_sb("Polarisation_PolCal_leads_too_much.xml")
            DSAHourAnglePolicy.compute_min_max_HA(sb)

    def test_pol_cal_trailing_too_much(self):
        #PolCal trailing too much
        with pytest.raises(ValueError):
            sb = self.generate_sb("Polarisation_PolCal_trails_too_much.xml")
            DSAHourAnglePolicy.compute_min_max_HA(sb)

    def test_get_elevation_at_ALMA_site(self):
        #just checking if it runs
        coord = SkyCoord('05h47m17.0876901s', '-51d03m59.441135s', frame='icrs')
        DSAHourAnglePolicy.get_elevation_at_ALMA_site(DEC=coord.dec,HA=Angle(-3*u.hour))

    def test_target_elevation_limits(self):
        execution_time = 1*constants.hour
        test_cases = [#above 20 deg for HA between -1.51h and 1.51h:
                      {"coord":SkyCoord("00h00m0.0s","43d46m56.772s"),
                       "HA_ok":[-1.4,0.3],
                       "HA_bad":[-1.6,1.6,1],"DGC":False},
                      #above 40 deg for HA between -2.15h and 2.15h:
                      {"coord":SkyCoord("00h00m0.0s","15d46m56.772s"),
                       "HA_ok":[-2,0,1],
                       "HA_bad":[-3,1.5,2.3],"DGC":True},
                      #above 88 deg for HA between -0.13h and 0.13h:
                      {"coord":SkyCoord("00h00m0.0s","-23d46m56.772s"),
                       "HA_ok":[-2,0.2],
                       "HA_bad":[-0.11,0.11,0,-1,-1.05],"DGC":False},
                      #above 85 deg for HA between -0.36h and 0.36h
                      {"coord":SkyCoord("00h00m0.0s","-23d46m56.772s"),
                       "HA_ok":[-2,0.5],
                       "HA_bad":[-0.35,0.35,-1.1],"DGC":True}]
        for test_case in test_cases:
            for HA in test_case["HA_ok"]:
                assert not DSAHourAnglePolicy.outside_elevation_limits(
                              DEC=test_case["coord"].dec,start_HA=Angle(HA*u.hour),
                              execution_time=execution_time,target_is_DGC=test_case["DGC"])
            for HA in test_case["HA_bad"]:
                assert DSAHourAnglePolicy.outside_elevation_limits(
                              DEC=test_case["coord"].dec,start_HA=Angle(HA*u.hour),
                              execution_time=execution_time,target_is_DGC=test_case["DGC"])