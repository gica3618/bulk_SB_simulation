#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 14:09:03 2026

@author: gianni
"""

from batch_simulations.domain.dsa_policy import DSAPolicy
from batch_simulations.domain.calibrator import Calibrator
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
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert HA_limits["min"].hour == -3
        assert HA_limits["max"].hour == 2

    def test_general_south(self):
        sb = self.generate_sb("general_DEC-15.xml")
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -15.39436944
        assert HA_limits["min"].hour == -4
        assert HA_limits["max"].hour == 3

    def test_leading_pol_cal(self):
        #example with leading PolCal:
        sb =self.generate_sb("example_polarisation_2023.1.00013.S.xml")
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
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
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 1
        assert sb.get_PolCal().coordinates.ra.deg == 359.471941933
        delta_ra = (360+1)-359.471941933
        assert np.isclose(HA_limits["min"].hour, -3 - delta_ra/360*24, atol=0, rtol=1e-8)
        assert HA_limits["max"].hour == 3

    def test_trailing_pol_cal(self):
        #example with trailing PolCal
        sb = self.generate_sb("Polarisation_leading_PolCal.xml")
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 165.0
        assert sb.get_PolCal().coordinates.ra.deg == 165.814104429
        assert HA_limits["min"].hour == -3 + (165.814104429-165)/360*24
        assert HA_limits["max"].hour == 3

    def test_trailing_pol_cal_RA0(self):
        #example with trailing PolCal around RA=0
        sb = self.generate_sb("Polarisatoin_trailing_PolCal_around_RA0.xml")
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == -59.52381111
        assert sb.rep_coord.ra.deg == 355.0
        assert sb.get_PolCal().coordinates.ra.deg == 1.14856453333
        delta_ra = (360+1.14856453333) - 355
        assert HA_limits["min"].hour == -3 + delta_ra/360*24
        assert HA_limits["max"].hour == 3

    def test_leading_pol_cal_high_dec(self):
        #leading PolCal, no adjustement needed (anyway min is -3 because of high DEC)
        sb = self.generate_sb("Polarisation_leading_PolCal_highDEC.xml")
        HA_limits = DSAPolicy.compute_min_max_HA(sb)
        assert sb.rep_coord.dec.deg == 0.5238108333333333
        assert sb.rep_coord.ra.deg == 150
        assert sb.get_PolCal().coordinates.ra.deg == 139.933497596
        assert HA_limits["min"].hour == -3
        assert HA_limits["max"].hour == 2

    def test_pol_cal_leading_too_much(self):
        #PolCal leading too much
        with pytest.raises(ValueError):
            sb = self.generate_sb("Polarisation_PolCal_leads_too_much.xml")
            DSAPolicy.compute_min_max_HA(sb)

    def test_pol_cal_trailing_too_much(self):
        #PolCal trailing too much
        with pytest.raises(ValueError):
            sb = self.generate_sb("Polarisation_PolCal_trails_too_much.xml")
            DSAPolicy.compute_min_max_HA(sb)

    def test_get_elevation_at_ALMA_site(self):
        #just checking if it runs
        coord = SkyCoord('05h47m17.0876901s', '-51d03m59.441135s', frame='icrs')
        DSAPolicy.get_elevation_at_ALMA_site(DEC=coord.dec,HA=Angle(-3*u.hour))

    # def test_elevation_ok(self):
    #     raise NotImplementedError
    #     test_cases = [#above 20 deg for HA between -1.51h and 1.51h:
    #                   {"coord":SkyCoord("00h00m0.0s","43d46m56.772s"),
    #                    "HA_ok":[-1.4,0.3],
    #                    "HA_bad":[-1.6,1.6,1],"DGC":False},
    #                   #above 40 deg for HA between -2.15h and 2.15h:
    #                   {"coord":SkyCoord("00h00m0.0s","15d46m56.772s"),
    #                    "HA_ok":[-2,0,1],
    #                    "HA_bad":[-3,1.5,2.3],"DGC":True},
    #                   #above 88 deg for HA between -0.13h and 0.13h:
    #                   {"coord":SkyCoord("00h00m0.0s","-23d46m56.772s"),
    #                    "HA_ok":[-2,0.2],
    #                    "HA_bad":[-0.11,0.11,0,-1,-1.05],"DGC":False},
    #                   #above 85 deg for HA between -0.36h and 0.36h
    #                   {"coord":SkyCoord("00h00m0.0s","-23d46m56.772s"),
    #                    "HA_ok":[-2,0.5],
    #                    "HA_bad":[-0.35,0.35,-1.1],"DGC":True}]
    #     for test_case in test_cases:
    #         for HA in test_case["HA_ok"]:
    #             assert DSAPolicy.elevation_ok(
    #                           DEC=test_case["coord"].dec,start_HA=Angle(HA*u.hour),
    #                           is_DGC=test_case["DGC"])
    #         for HA in test_case["HA_bad"]:
    #             assert not DSAPolicy.elevation_ok(
    #                           DEC=test_case["coord"].dec,start_HA=Angle(HA*u.hour),
    #                           is_DGC=test_case["DGC"])

    def test_targets_for_elevation_check(self):
        class FakeSB:
            science_targets = [{"name":"miraculix",
                                "coordinates":SkyCoord('05h47m17.08s', '-51d03m59.44s')},
                               {"name":"weissnix",
                                "coordinates":SkyCoord('10h47m17.08s', '-71d03m59.44s')}]
            calibrators = [Calibrator(name="Phase",source_name="J1326-5256",
                                      cal_type="Phase",is_hardcoded=True,
                                      coordinates=SkyCoord('05h00m00s', '-53d00m00s')),
                           Calibrator(name="Bandpass calibrator",source_name="query",
                                      cal_type="Bandpass",is_hardcoded=False,
                                      coordinates=None),
                           Calibrator(name="DGC",source_name="asterix",
                                      cal_type="DGC",is_hardcoded=True,
                                      coordinates=SkyCoord('08h00m00s', '-20d00m00s'))]
        sb = FakeSB()
        targets = DSAPolicy.targets_for_elevation_check(sb)
        assert len(targets) == 4
        for i in (0,1):
            assert targets[i] == {"name":FakeSB.science_targets[i]["name"],
                                  "DEC":FakeSB.science_targets[i]["coordinates"].dec,"DGC":False}
        assert targets[2] == {"name":"Phase","DEC":FakeSB.calibrators[0].coordinates.dec,"DGC":False}
        assert targets[3] == {"name":"DGC","DEC":FakeSB.calibrators[2].coordinates.dec,"DGC":True}   

    # def test_targets_beyond_elevation_limits(self):
    #     raise NotImplementedError()  
    #     targets = [{"name":"dgc miraculix","DEC":SkyCoord('00h00m00s', '-10d00m00s').dec,
    #                 "DGC":True},#this should not be limited by elevation
    #                {"name":"J12345678","DEC":SkyCoord('00h00m00s', '80d00m00s').dec,
    #                 "DGC":False},#this one is unobservable by ALMA
    #                {"name":"majestix","DEC":SkyCoord("00h00m0.0s","-23d46m56.772s").dec,
    #                 "DGC":True},#above 85 deg for HA between -0.36h and 0.36h
    #                {"name":"majestix","DEC":SkyCoord("00h00m0.0s","-23d46m56.772s").dec,
    #                 "DGC":False}#above 85 deg for HA between -0.36h and 0.36h
    #                ]
    #     output = SingleSBSimulation.get_targets_beyond_elevation_limits(
    #                       start_HA=Angle(-0.8*u.hour),execution_time=0.5*constants.hour,
    #                       targets=targets)
    #     assert len(output) == 2
    #     assert output[0] == "J12345678"
    #     assert output[1] == "majestix"