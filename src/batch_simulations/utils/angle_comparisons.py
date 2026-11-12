#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:35:18 2026

@author: gianni
"""

#for comparing astropy Angle objects due to machine precision, we need these methods
#for example, when doing HA=Angle(-15,unit=u.deg), then, when doing HA.hour,
#we get np.float64(-1.0000000000000002)

from astropy import units as u
import numpy as np

tol = 1e-10 * u.rad

def leq(a, b):
    return a <= b + tol

def geq(a, b):
    return a >= b - tol

def less(a, b):
    return a < b - tol

def greater(a, b):
    return a > b + tol

def equal(a,b):
    return np.isclose(a.rad,b.rad,rtol=0,atol=tol.to(u.rad).value)