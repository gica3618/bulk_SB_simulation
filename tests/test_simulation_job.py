#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 13:50:27 2026

@author: gianni
"""


from batch_simulations.domain.simulation_job import SimulationJob
from astropy.coordinates import Angle
from astropy import units as u
import datetime


def test_epoch():
    job = SimulationJob(xml_filepath="a.xml",array_config="c43-3",
                        HA=Angle(1.2*u.hour),date=None)
    assert job.epoch() == "TRANSIT+1.2h"
    job.HA = Angle(-1.234567*u.hour)
    assert job.epoch() == "TRANSIT-1.23h"
    job.date = datetime.date(year=2025,month=3,day=23)
    assert job.epoch() == "TRANSIT-1.23h,2025-03-23"
    job.HA = Angle(0*u.hour)
    assert job.epoch() == "TRANSIT,2025-03-23"

def test_command():
    
    job = SimulationJob(xml_filepath="a.xml",array_config="c43-3",
                        HA=Angle(1.2*u.hour),
                        date=datetime.date(year=1902,month=2,day=2))
    assert job.command() == ["simulateSB.py", "a.xml", "TRANSIT+1.2h,1902-02-02",
                             "-C","c43-3"]