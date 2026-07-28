#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:37:27 2026

@author: gianni
"""

from collections import Counter
import logging
from scipy import constants
from batch_simulations.domain.dsa_ha_policy import DSAHourAnglePolicy


class SB:

    def __init__(self, calibrators, science_targets, mode_name, nominal_configs,
                 rep_coord, OT_allowed_HA, requires_TP, total_execution_time,
                 number_of_executions,metadata=None):
        self.metadata = metadata or {}
        self.calibrators = calibrators
        self.science_targets = science_targets
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
        self.total_execution_time = total_execution_time
        self.number_of_executions = number_of_executions
        self.consistency_checks()

    def consistency_checks(self):
        cal_type_counts = Counter(self.cal_types)
        for cal_type,count in cal_type_counts.items():
            if count > 1:
                #note that this test cannot detect if there are several query
                #calibrators of the same type, because those calibrators all
                #point to the same field source in the xml
                raise ValueError(f"multiple calibrators of type {cal_type}")
        if len(self.science_targets) < 1:
            raise ValueError("at least one science target expected")
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
        return DSAHourAnglePolicy.compute_min_max_HA(sb=self)

    def add_metadata(self,key,value):
        self.metadata[key] = value

    def single_execution_time(self):
        single_exec = self.total_execution_time/self.number_of_executions
        if self.is_Polarisation:
            #polarisation execution (session) should be at least 3 hours
            min_number_of_executions = ((3*constants.hour) // single_exec) + 1
            session_duration = single_exec*min_number_of_executions
            if session_duration < 3*constants.hour:
                raise ValueError("session duration less than 3 hours")
            return session_duration
        else:
            return single_exec