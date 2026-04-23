#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:41:43 2026

@author: gianni
"""

from batch_simulations.domain.simulation_job import SimulationJob


class JobPlanner:

    def __init__(self, sb, default_HA_step, fine_HA_step):
        self.sb = sb
        self.default_HA_step = default_HA_step
        self.fine_HA_step = fine_HA_step
        self.DSA_HA = sb.get_DSA_HA_limits()

    def compute_HAs(self, step):
        if step.hour <= 0:
            raise ValueError(f"invalid HA step: {step.hour} H")
        #min_HA considered by DSA might be affected by PolCal, but max_HA is
        #always an integer, so I start with max_HA
        if not self.DSA_HA["max"].hour.is_integer():
            raise RuntimeError("expected max DSA HA to be an integer, but it is"
                               +f"{self.DSA_HA['max'].hour}h")
        HAs = []
        HA = self.DSA_HA["max"].copy() #safer to make a copy here, to avoid modifying self.DSA_HA
        while HA >= self.DSA_HA["min"]:
            HAs.append(HA)
            #don't do "HA -= step" here, as this leads to modification of all
            #HAs already in the list:
            HA = HA - step
        return sorted(HAs)

    def HA_jobs(self, xml_filepath, array_config, date, step):
        jobs = []
        HAs = self.compute_HAs(step)
        for HA in HAs:
            job = SimulationJob(xml_filepath=xml_filepath,array_config=array_config,
                                HA=HA,date=date)
            jobs.append(job)
        return HAs,jobs

    def HA_jobs_default_HA_step(self,xml_filepath,array_config,date):
        return self.HA_jobs(xml_filepath=xml_filepath, array_config=array_config,
                            date=date, step=self.default_HA_step)

    def HA_jobs_fine_HA_step(self,xml_filepath,array_config,date):
        return self.jobs(xml_filepath=xml_filepath, array_config=array_config,
                         date=date, step=self.fine_HA_step)