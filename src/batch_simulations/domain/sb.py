#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:37:27 2026

@author: gianni
"""

from collections import Counter
import logging
from batch_simulations.domain.dsa_ha_policy import DSAHourAnglePolicy


class SB:

    def __init__(self, calibrators, mode_name, nominal_configs,
                 rep_coord, OT_allowed_HA, requires_TP, metadata=None):
        self.metadata = metadata or {}
        self.calibrators = calibrators
        self.cal_types = [c.cal_type for c in self.calibrators]
        self.mode_name = mode_name
        self.is_Polarisation = "Polarization" in mode_name
        #TODO exclude VLBI and Solar from simulations?
        self.is_VLBI = 'VLBI' in mode_name
        self.is_Solar = 'Solar' in mode_name
        self.is_B2B = mode_name == 'BandToBand Interferometry'
        self.nominal_configs = nominal_configs
        self.is_7m = "7M" in self.nominal_configs
        self.rep_coord = rep_coord
        self.OT_allowed_HA = OT_allowed_HA
        self.requires_TP = requires_TP
        self.consistency_checks()

    def consistency_checks(self):
        cal_type_counts = Counter(self.cal_types)
        for cal_type,count in cal_type_counts.items():
            if count > 1:
                #note that this test cannot detect if there are several query
                #calibrators of the same type, because those calibrators all
                #point to the same field source in the xml
                raise ValueError(f"multiple calibrators of type {cal_type}")
        if self.is_Polarisation and self.has_no_PolCal():
            raise ValueError("Polarization SB without PolCal")
        if self.is_Polarisation and (not self.PolCal_is_hardcoded()):
            raise ValueError("PolCal not hardcoded")
        if (not self.is_Polarisation) and self.has_at_least_one_PolCal():
            raise ValueError("SB is not Polarisation, but PolCal is present")
        logging.info("SB passed consistency checks")

    def has_no_PolCal(self):
        return ("Polarization" not in self.cal_types)

    def has_at_least_one_PolCal(self):
        return not self.has_no_PolCal()

    def get_calibrator(self,cal_type):
        cal = [cal for cal in self.calibrators if cal.cal_type==cal_type]
        if len(cal) != 1:
            raise RuntimeError(f"cannot uniquely determine calibrator of type {cal_type}")
        return cal[0]

    def get_PolCal(self):
        return self.get_calibrator(cal_type="Polarization")

    def PolCal_is_hardcoded(self):
        return self.get_PolCal().is_hardcoded

    def any_calibrator_hardcoded(self):
        return any([c.is_hardcoded for c in self.calibrators])

    def get_DSA_HA_limits(self):
        return DSAHourAnglePolicy.compute(sb=self)

    def add_metadata(self,key,value):
        self.metadata[key] = value