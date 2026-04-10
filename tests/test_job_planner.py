#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 15:38:52 2026

@author: gianni
"""

from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.job_planner import JobPlanner
from astropy import units as u
from astropy.coordinates import Angle
import numpy as np
import pytest
from pathlib import Path
import datetime


class TestJobPlanner:

    @staticmethod
    def construct_sb(filename):
        filepath = Path('tests/xmls') / filename
        return BuildSBFromXML.build(xml_filepath=filepath)

    def generate_general_job(self):
        sb = self.construct_sb("2025.1.01279.S_general_SB.xml")
        return JobPlanner(sb=sb, default_HA_step=Angle(1*u.hour),
                          fine_HA_step=Angle(0.25*u.hour))

    def test_compute_HAs_general(self):
        job_planner = self.generate_general_job()
        expected_HA_limits = job_planner.sb.get_DSA_HA_limits()
        #test step where both DSAmin and DSAmax are exactly included:
        for step in (Angle(1*u.hour),Angle(0.25*u.hour)):
            sim_HA = job_planner.compute_HAs(step=step)
            assert np.all([HA.hour for HA in sim_HA] == np.arange(expected_HA_limits["min"].hour,
                                                                  expected_HA_limits["max"].hour+0.0001,
                                                                  step.hour))
        #test the case where DSAmin is not included
        step = Angle(0.7486*u.hour)
        sim_HA = np.array([HA.hour for HA in job_planner.compute_HAs(step=step)])
        assert sim_HA[-1] == expected_HA_limits["max"].hour
        diff = np.diff(sim_HA)
        assert np.allclose(diff[0],diff,atol=0,rtol=1e-6)
        assert sim_HA[0] > expected_HA_limits["min"].hour
        assert sim_HA[0]-step.hour < expected_HA_limits["min"].hour

    def test_compute_HAs_explicitly(self):
        job_planner = self.generate_general_job()
        job_planner.DSA_HA["min"] = Angle(-3.2*u.hour)
        job_planner.DSA_HA["max"] = Angle(3*u.hour)
        sim_HA = job_planner.compute_HAs(step=Angle(0.5*u.hour))
        assert np.all([HA.hour for HA in sim_HA] == np.arange(-3,3.1,0.5))
    
    def test_negative_step(self):
        job_planner = self.generate_general_job()
        with pytest.raises(ValueError):
            job_planner.compute_HAs(step=Angle(-1.2*u.hour))
    
    def test_noninteger_max_HA_DSA(self):
        job_planner = self.generate_general_job()
        job_planner.DSA_HA["max"] = Angle(2.1*u.hour)
        with pytest.raises(RuntimeError):
            job_planner.compute_HAs(step=Angle(1*u.hour))

    def test_HA_jobs(self):
        job_planner = self.generate_general_job()
        xml_filepath = 'tests/xmls/2025.1.01279.S_general_SB.xml'
        array_config  = "c43-6"
        date = datetime.date(year=1912,month=10,day=5)
        step = Angle(0.3*u.hour)
        jobs = job_planner.HA_jobs(xml_filepath=xml_filepath, array_config=array_config,
                                   date=date, step=step)
        job_HAs = job_planner.compute_HAs(step=step)
        assert len(jobs) == len(job_HAs)
        for i,HA in enumerate(job_HAs):
            assert jobs[i].HA == HA