#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 17:41:07 2026

@author: gianni
"""

from astropy.coordinates import Angle
from astropy import units as u
import batch_simulations.utils.angle_comparisons as angle_comp


b = Angle(1*u.hour)

def test_leq():
    leq = (b-Angle(1*u.hour),b,Angle(b.rad*u.rad))
    for a in leq:
        assert angle_comp.leq(a,b)
    greater = (b+Angle(0.00001*u.hour),b+Angle(2*u.hour))
    for a in greater:
        assert not angle_comp.leq(a,b)

def test_angle_geq():
    geq = (b,Angle(b.rad*u.rad),b+Angle(1*u.hour))
    for a in geq:
        assert angle_comp.geq(a,b)
    less = (b-Angle(1*u.hour),b-Angle(0.00001*u.hour))
    for a in less:
        assert not angle_comp.geq(a,b)

def test_angle_less_than():
    le = (b-Angle(1*u.hour),b-Angle(0.00001*u.hour))
    for a in le:
        assert angle_comp.less(a,b)
    geq = (b+Angle(1*u.hour),b,Angle(b.rad*u.rad))
    for a in geq:
        assert not angle_comp.less(a,b)

def test_angle_greater_than():
    ge = (b+Angle(1*u.hour),b+Angle(0.000001*u.hour))
    for a in ge:
        assert angle_comp.greater(a,b)
    leq = (b-Angle(1*u.hour),b,Angle(b.rad*u.rad))
    for a in leq:
        assert not angle_comp.greater(a,b)

def test_angle_equal():
    equal = (b,Angle(b.rad*u.rad),Angle(b.deg*u.deg))
    for a in equal:
        assert angle_comp.equal(a,b)
    not_equal = (b+Angle(1*u.hour),b+Angle(0.0001*u.hour))
    for a in not_equal:
        assert not angle_comp.equal(a,b)