#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 23:12:34 2026

@author: gianni
"""

from batch_simulations.utils.retry import retry
import time

def test_retry():
    out = retry(operation=lambda:3, max_trials=5, sleep_time=5,
                retry_condition=lambda x: False)
    assert out == 3

def test_retry_max_trials():
    sleep_time = 1
    max_trials = 5
    start = time.time()
    out = retry(operation=lambda:3, max_trials=max_trials, sleep_time=sleep_time,
                retry_condition=lambda x: True)
    end = time.time()
    assert out == 3
    assert end-start >= max_trials*sleep_time