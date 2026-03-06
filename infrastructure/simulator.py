#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:56:21 2026

@author: gianni
"""

import logging
import shutil
import subprocess


class Simulator:

    def __init__(self,work_folder):
        self.work_folder = work_folder

    def prepare_config_file(self,array_config):
        if array_config.endswith(".cfg"):
            logging.info(f'will use cfg file {array_config}, copying it'
                         +' to the work folder')
            shutil.copy(src=array_config,dst=self.work_folder)

    def simulate(self,job):
        command = job.command()
        self.prepare_config_file(array_config=job.array_config)
        logging.info(f"going to execute the following command:\n{' '.join(command)}")
        process = subprocess.run(command,cwd=self.work_folder,text=True,
                                 capture_output=True)
        pipe = {'stdout':process.stdout,'stderr':process.stderr}
        loggers = {'stdout':logging.info,'stderr':logging.error}
        for key,logger in loggers.items():
            logger(f'{key} of simulateSB.py:')
            logger(pipe[key])
        return process