#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 11:07:41 2026

@author: gianni
"""

from batch_simulations.domain.calibrator import Calibrator,classify_calibrator,\
       CALIBRATOR_TYPES
import pytest
from astropy.coordinates import SkyCoord
from astropy import units as u


def test_calibrator():
    coord = SkyCoord(ra=10*u.deg,dec=-2*u.deg)
    Calibrator(name="test", source_name="query", cal_type="Bandpass",
               is_hardcoded=False,coordinates=None)
    Calibrator(name="test", source_name="J123", cal_type="DGC",
               is_hardcoded=True,coordinates=coord)
    with pytest.raises(ValueError):
        Calibrator(name="test", source_name="query", cal_type="Bandpass123",
                   is_hardcoded=False,coordinates=None)
    with pytest.raises(ValueError):
        Calibrator(name="test", source_name="query", cal_type="Bandpass",
                   is_hardcoded=False,coordinates=coord)
    with pytest.raises(ValueError):
        Calibrator(name="test", source_name="Jirgendwas", cal_type="Bandpass",
                   is_hardcoded=True,coordinates=None)

def test_classify_calibrator():
    for cal_type in CALIBRATOR_TYPES:
        assert classify_calibrator(cal_type) == cal_type
        assert classify_calibrator(cal_type.lower()) == cal_type
        assert classify_calibrator(f"{cal_type} calibrator") == cal_type
    assert classify_calibrator("J123") is None
    with pytest.raises(ValueError):
        classify_calibrator("Bandpass Phase calibrator hahahah")