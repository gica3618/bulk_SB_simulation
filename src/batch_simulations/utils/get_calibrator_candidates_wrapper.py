#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 18:38:28 2026

@author: gianni
"""

import subprocess
import logging
from batch_simulations.infrastructure.calibrator_query import\
                                                parse_calibrator_candidate_data


class GetCalibratorCandidatesWrapper:
    
    def __init__(self,xml_filepath):
        self.xml_filepath = xml_filepath
    
    def construct_command(self,integration_time,array_config,search_radius=None,
                          epoch=None,calibrator_type=None,spectral_spec=None,
                          src=None,no_spwavg=False,maxAge=None):
        command = ["getCalibratorCandidates.py", self.xml_filepath, "-t",
                   str(integration_time),"-C", array_config]
        if search_radius is not None:
            command.extend(["-r", str(search_radius)])
        if epoch is not None:
            command.extend(["-e", epoch])
        if calibrator_type is not None:
            command.extend(["-c", calibrator_type])
        if spectral_spec is not None:
            command.append(f'--spectralSpec={spectral_spec}')
        if src is not None:
            command.append(f'--src={src}')
        if no_spwavg:
            #TODO verify that only phase and check use no_spwavg, see email 
            #conversation with Akihiko
            if calibrator_type not in ('phase', 'check'):
                raise ValueError(
                    f"no_spwavg has no effect for calibrator_type={calibrator_type}"
                )
            command.append('--no_spwavg')
        if maxAge is not None:
            command.append(f"--maxAge={maxAge}")
        return command

    def run(self,integration_time,array_config,search_radius=None,epoch=None,
            calibrator_type=None,spectral_spec=None,src=None,no_spwavg=False):
        command = self.construct_command(
                        integration_time=integration_time,array_config=array_config,
                        search_radius=search_radius,epoch=epoch,
                        calibrator_type=calibrator_type,spectral_spec=spectral_spec,
                        src=src,no_spwavg=no_spwavg)
        process = subprocess.run(command, capture_output=True, text=True, timeout=300)
        if process.returncode != 0:
            logging.error('getCalibratorCandidates.py crashed')
            logging.error(process.stderr)
            raise RuntimeError(f"getCalibatorCandidates.py crashed\n{process.stderr}")
        return self.read_calibrator_candidates(process=process)

    @staticmethod
    def read_calibrator_candidates(process):
        stdout = process.stdout.splitlines()
        calibrator_candidates = []
        first = None
        last = None
        for i,line in enumerate(stdout):
            if '-> Listing ranked candidate list...' in line:
                first = i+4
            if '6th col:' in line:
                last = i-2
            # if 'sources passed the selection criteria' in line:
            #     #example: "3 sources passed the selection criteria: [J0529-0519, J0532-0307, J0541-0541]"
            #      source_names = line.split('[')[1].replace(']','')
            #      source_names = source_names.split(',')
            #      source_names = [sn.strip() for sn in source_names]
        if first is None or last is None:
            raise RuntimeError("Could not locate candidate table in output")
        for line in stdout[first:last+1]:
            calibrator = parse_calibrator_candidate_data(line=line)
            calibrator_candidates.append(calibrator)
        return calibrator_candidates