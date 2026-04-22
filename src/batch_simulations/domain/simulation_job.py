#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 16:09:10 2026

@author: gianni
"""

from dataclasses import dataclass
from astropy.coordinates import Angle
from datetime import date


@dataclass
class SimulationJob:
    xml_filepath: str
    array_config: str
    HA: Angle
    date: date | None

    def epoch(self):
        epoch = "TRANSIT"
        if self.HA.hour != 0:
            epoch += f"{self.HA.hour:+.3g}h"
        if self.date:
            epoch += f",{self.date.isoformat()}"
        return epoch

    def command(self):
        return ["simulateSB.py", str(self.xml_filepath), self.epoch(), "-C",
                self.array_config]