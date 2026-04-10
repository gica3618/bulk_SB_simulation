#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:13:44 2026

@author: gianni
"""

from batch_simulations.application.simulation_runner import SimulationRunner
from unittest.mock import Mock, patch
from pathlib import Path


#following test was written by ChatGPT

def test_simulation_runner_calls_retry():
    runner = SimulationRunner()

    fake_work_folder = "/tmp/work"
    fake_job = Mock()
    fake_job.xml_filepath = Path("/path/to/file.xml")
    fake_process = Mock()
    fake_command = ["run", "something"]
    fake_sim_result = Mock()
    fake_sim_result.server_error = False
    #note that I need to patch e.g. batch_simulations.application.simulation_runner.Simulator,
    #not batch_simulations.infrastructure.simulator! This is because Simulator
    #is imported in simulation_runner
    with patch("batch_simulations.application.simulation_runner.WorkFolder") as MockWorkFolder, \
         patch("batch_simulations.application.simulation_runner.Simulator") as MockSimulator,\
         patch("batch_simulations.application.simulation_runner.retry") as mock_retry, \
         patch("batch_simulations.application.simulation_runner.SimulationResult") as MockSimResult:

        # WorkFolder context manager
        MockWorkFolder.return_value.__enter__.return_value = fake_work_folder

        # Simulator behavior
        mock_simulator_instance = MockSimulator.return_value
        mock_simulator_instance.simulate.return_value = (fake_process, fake_command)

        # SimulationResult behavior
        MockSimResult.from_completed_process.return_value = fake_sim_result

        # retry just calls the operation once
        def fake_retry(operation, max_trials, sleep_time, retry_condition):
            result = operation()
            # verify retry_condition works
            assert retry_condition(result) is False
            return result
        mock_retry.side_effect = fake_retry

        result = runner.run(fake_job)
        MockSimulator.assert_called_once_with(work_folder=fake_work_folder)
        MockSimResult.from_completed_process.assert_called_once_with(
            executed_command=" ".join(fake_command),
            process=fake_process,
            output_folder=fake_work_folder,
            xml_filename=fake_job.xml_filepath.name)
        assert result == fake_sim_result