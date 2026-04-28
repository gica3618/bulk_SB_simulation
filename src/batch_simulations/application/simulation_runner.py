#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:43:14 2026

@author: gianni
"""

from pathlib import Path
import logging
from batch_simulations.utils.workfolder import WorkFolder
from batch_simulations.infrastructure.simulator import Simulator
from batch_simulations.domain.simulation_result import SimulationResult
from batch_simulations.utils.retry import retry


class SimulationRunner:

    max_trials = 10
    sleep_time = 60

    def run(self,job):
        with WorkFolder() as work_folder:
            simulator = Simulator(work_folder=work_folder)
            def operation():
                process, command = simulator.simulate(job)
                return SimulationResult.from_completed_process(
                                       executed_command=" ".join(command),
                                       process=process,output_folder=work_folder,
                                       xml_filename=Path(job.xml_filepath).name)
            def retry_condition(sim_result):
                if sim_result.success:
                    logging.info("retry condition not satisfied (simulation succeeded)")
                    return False
                if sim_result.fail_reason.category == "server error":
                    logging.info("retry condition satisfied (server error)")
                    return True
                else:
                    logging.info("simulation failed, but retry condition not satisfied")
                    return False
            return retry(operation=operation,max_trials=self.max_trials,
                         sleep_time=self.sleep_time,
                         retry_condition=retry_condition)

    def run_jobs(self,jobs):
        return [self.run(job) for job in jobs] 