#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 12:46:03 2023

@author: gianni
"""

import bulk_SB_simulation
from datetime import date
import sys
import logging

logging_level = logging.INFO
#logging_level = logging.ERROR
obs_dates = [date.today(),date(year=2023,month=9,day=10)]


logging.basicConfig(format='%(levelname)s: %(message)s',level=logging_level,
                    stream=sys.stdout)

# sim12m = bulk_SB_simulation.BulkSimulation12m(
#                           SB_list='lookup_table_12mSBs_20230828.csv',
#                           support_arcs=['EA',],array_config='c43-9')
# sim12m.run_simulations(obs_dates=obs_dates)
# sim12m.write_failed_simulations_to_csvfile(filepath='failed_simulations_12m.csv')

# def custom_filter_7m(data):
#     return data['Project Code'][:4] != '2023'
# sim7m = bulk_SB_simulation.BulkSimulation7m(
#                    SB_list='SBs_7m_2023-08-28.csv',custom_SB_filter=custom_filter_7m)
# sim7m.run_simulations(obs_dates=obs_dates)
# sim7m.write_failed_simulations_to_csvfile(filepath='failed_simulations_7m.csv')

#search for 7m SBs that require TP antennas:
sim7m = bulk_SB_simulation.BulkSimulation7m(
                    SB_list='SBs_7m_2023-09-15.csv',
                    array_configs=['aca.cm10.pm3.cfg','7m'],HAs=['TRANSIT-1h',])
sim7m.run_simulations(obs_dates=[date(year=2023,month=10,day=1),])
sim7m.write_failed_simulations_to_csvfile(
               filepath='failed_simulations_7m_TRANSIT-1h.csv')