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


logging.basicConfig(format='%(levelname)s: %(message)s',level=logging_level,
                    stream=sys.stdout)

#12m SBs
# obs_dates = [date.today(),date(year=2023,month=10,day=20)]
# sim12m = bulk_SB_simulation.BulkSimulation12m(
#                           SB_list='12m_SBs_2023-10-05.csv',
#                           support_arcs=None,array_config='c43-8')
# sim12m.run_simulations(obs_dates=obs_dates)
# sim12m.write_failed_simulations_to_csvfile(filepath='failed_simulations_12m.csv')
# sim12m.print_statistics()

#7m SBs
obs_dates = [date.today(),]
sim7m = bulk_SB_simulation.BulkSimulation7m(SB_list='7m_SBs_2023-10-05.csv')
sim7m.run_simulations(obs_dates=obs_dates)
sim7m.write_failed_simulations_to_csvfile(filepath='failed_simulations_7m.csv')
sim7m.print_statistics()

#search for 7m SBs that require TP antennas:
# sim7m = bulk_SB_simulation.BulkSimulation7m(
#                     SB_list='SBs_7m_2023-09-15.csv',
#                     array_configs=['aca.cm10.pm3.cfg','7m'],HAs=[-1,])
# sim7m.run_simulations(obs_dates=[date.today(),])
# sim7m.write_failed_simulations_to_csvfile(
#                filepath='failed_simulations_7m_TRANSIT-1h.csv')

#simulate specific 7m project
# obs_dates = [date.today(),]
# array_configs = ['aca.cm10.pm3.cfg',]
# project_code = '2023.1.00578.S'
# def custom_filter_7m(data):
#     return data['Project Code'] == project_code
# sim7m = bulk_SB_simulation.BulkSimulation7m(
#                     SB_list='2023-10-03_7m_schedBlockList.csv',
#                     custom_SB_filter=custom_filter_7m,array_configs=array_configs)
# sim7m.run_simulations(obs_dates=obs_dates)
# sim7m.write_failed_simulations_to_csvfile(
#                       filepath=f'failed_simulations_7m_{project_code}.csv')