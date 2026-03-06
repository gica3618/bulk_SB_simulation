#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:53:28 2026

@author: gianni
"""

import os
import logging
from job_planner import JobPlanner
from simulation_runner import SimulationRunner
import traceback
from infrastructure.xml_repository import SBFactory,XMLDownloader
from astropy.coordinates import Angle
from astropy import units as u
from dataclasses import dataclass


@dataclass
class SingleSBSimulationResult:
    dates: list
    HAs: list | None
    results: list | None
    unexpected_error: None | str


class SingleSBSimulation:

    def __init__(self,project_code,sb_name,obs_dates,array_config_12m):
        self.project_code = project_code
        self.sb_name = sb_name
        self.obs_dates = obs_dates
        self.array_config_12m = array_config_12m
        self.xml_filepath = os.path.join(
                              os.getcwd(),
                              f"{self.project_code}_{self.sb_name}.xml".replace(" ", "_"))

    def simulate(self):
        try:
            XMLDownloader.download_xml(project_code=self.project_code,
                                       sb_name=self.sb_name,
                                       filepath=self.xml_filepath)
            sb = SBFactory.from_xml(self.xml_filepath)
            array_config = self.get_array_config(sb)
            job_planner = JobPlanner(sb,default_step=self.default_HA_step,
                                     fine_step=self.fine_HA_step)
            runner = SimulationRunner()
            all_results = []
            simulated_HAs = []
            logging.info(f"Simulating SB {self.sb_name} ({self.project_code})")
            for date in self.obs_dates:
                job_planner_kwargs = {"xml_filepath":self.xml_filepath,
                                      "array_config":array_config,"date":date}
                jobs = job_planner.jobs(**job_planner_kwargs,
                                        step=self.default_HA_step)
                results = runner.run_jobs(jobs)
                if not all(r.success for r in results):
                    logging.info("At least one simulation failed — retrying with finer HA grid")
                    jobs = job_planner.jobs(**job_planner_kwargs,
                                            step=self.fine_HA_step)
                    results = runner.run_jobs(jobs)
                simulated_HAs.append([j.HA for j in jobs])
                all_results.append(results)
            return SingleSBSimulationResult(dates=self.obs_dates, HAs=simulated_HAs,
                                            results=all_results, unexpected_error=None)
        except Exception as error:
            logging.error("unexpected error occured during consideration of "
                         +f" {self.sb_name} ({self.project_code}). Traceback:")
            logging.error(traceback.format_exc())
            return SingleSBSimulationResult(dates=self.obs_dates, HAs=None,
                                            results=None, unexpected_error=str(error))


class SimulationCampaign:

    default_HA_step = Angle(1*u.hour)
    fine_HA_step = Angle(0.25*u.hour)

    def __init__(self, sb_table, obs_dates, array_config_12m):
        self.sb_table = sb_table
        self.obs_dates = obs_dates
        self.array_config_12m = array_config_12m

    def run(self):
        results = []
        for row in self.sb_table.itertuples():
            single_sb_sim = SingleSBSimulation(
                                project_code=row.code, sb_name=row.sbname,
                                obs_dates=self.obs_dates,
                                array_config_12m=self.array_config_12m)
            results.append(single_sb_sim.simulate())
        return results