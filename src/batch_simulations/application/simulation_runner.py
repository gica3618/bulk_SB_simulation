#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:43:14 2026

@author: gianni
"""

from pathlib import Path
from batch_simulations.utils.workfolder import WorkFolder
from batch_simulations.infrastructure.simulator import Simulator
from batch_simulations.domain.simulation_result import SimulationResult
from batch_simulations.utils.retry import retry


class SimulationRunner:

    max_trials=10
    sleep_time=60

    def run(self,job):
        with WorkFolder() as work_folder:
            simulator = Simulator(work_folder=work_folder)
            def operation():
                process, command = simulator.simulate(job)
                return SimulationResult.from_completed_process(
                                       executed_command=" ".join(command),
                                       process=process,output_folder=work_folder,
                                       xml_filename=Path(job.xml_filepath).name)
            return retry(operation=operation,max_trials=self.max_trials,
                         sleep_time=self.sleep_time,
                         retry_condition=lambda sim_result: sim_result.server_error)

    def run_jobs(self,jobs):
        return [self.run(job) for job in jobs] 