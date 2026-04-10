#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:54:55 2026

@author: gianni
"""

from batch_simulations.infrastructure.simulator import Simulator
import tempfile
import os
from pathlib import Path


def test_prepare_config_file_std_config():
    with tempfile.TemporaryDirectory() as tmpdir:
        simulator = Simulator(work_folder=tmpdir)
        simulator.prepare_config_file_if_necessary("c43-3")
        assert not os.listdir(tmpdir)

def test_prepare_config_file_with_config_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        simulator = Simulator(work_folder=tmpdir_path)
        array_config = "some_config_file.cfg"
        p = Path(array_config)
        p.touch()
        simulator.prepare_config_file_if_necessary(array_config)
        assert (tmpdir_path / array_config).is_file()
        #delete the file:
        p.unlink()
        assert not p.exists()

class DummyProcess:
    def __init__(self, stdout, stderr):
        self.stdout = stdout
        self.stderr = stderr
        

def test_send_to_logger():
    with tempfile.TemporaryDirectory() as tmpdir:
        simulator = Simulator(work_folder=tmpdir)
        process = DummyProcess(stdout="some messages\nand some more",
                               stderr="ERROR!\nsomething went wrong")
        simulator.send_to_logger(process)


class DummyJob():
    def __init__(self, array_config):
        self.array_config = array_config

    def command(self):
        return ["echo",'Hello','World']

def test_simulate():
    with tempfile.TemporaryDirectory() as tmpdir:
        simulator = Simulator(work_folder=tmpdir)
        job = DummyJob(array_config="c43-3")
        process,command = simulator.simulate(job)
        assert process.stdout == "Hello World\n"
        assert process.stderr == ""
        assert command == job.command()