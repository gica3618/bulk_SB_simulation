#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:02:52 2026

@author: gianni
"""


import logging
import tempfile
import shutil


class WorkFolder:

    def __enter__(self):
        self.path = tempfile.mkdtemp(prefix="sb_sim_")
        logging.info(f"created work folder {self.path}")
        return self.path

    def __exit__(self, exc_type, exc, tb):
        logging.info(f"removing work folder {self.path}")
        shutil.rmtree(self.path)