#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:43:14 2026

@author: gianni
"""

import os
from utils.workfolder import WorkFolder
from infrastructure.simulator import SBSimulator
from domain.simulation_result import SimulationResult
from utils.retry import retry


class SimulationRunner:

    max_trials=10
    sleep_time=60

    def run(self,job):
        with WorkFolder() as work_folder:
            simulator = SBSimulator(work_folder=work_folder)
            def operation():
                process, command = simulator.simulate(job)
                return SimulationResult.from_completed_process(
                                       executed_command=" ".join(command),
                                       process=process,output_folder=work_folder,
                                       xml_filename=os.path.basename(job.xml_filepath))
            return retry(operation=operation,max_trials=self.max_trials,
                         sleep_time=self.sleep_time,
                         retry_condition=lambda sim: sim.server_error)

    def run_jobs(self,jobs):
        return [self.run(job) for job in jobs] 