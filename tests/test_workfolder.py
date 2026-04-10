#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 22:22:27 2026

@author: gianni
"""

from batch_simulations.utils import workfolder
import os


def test_work_folder():
    with workfolder.WorkFolder() as wf:
        assert os.path.isdir(wf)