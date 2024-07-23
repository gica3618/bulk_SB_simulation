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
# def project_selector(data):
#     return data['code'] == '2018.1.01243.S'
array_config = 'c43-5'
#obs_dates = [date.today(),date(year=2024,month=6,day=23)]
obs_dates = [date(year=2024,month=7,day=28),date(year=2024,month=8,day=17)]
sim12m = bulk_SB_simulation.BulkSimulation12m(
                          SB_list='12m_SBs_2024-07-19.csv',
                          support_arcs=None,array_config=array_config,
                          custom_SB_filter=None)
sim12m.run_simulations(obs_dates=obs_dates)
sim12m.write_failed_and_skipped_simulations_to_csvfile(
                  filepath=f'failed_and_skipped_simulations_12m_{array_config}.csv')
sim12m.write_statistics(f'statistics_{array_config}.txt')

'''
#7m SBs
def sb_selector(data):
    return data['SB Name'] == 'HD_16329_a_09_7M'
obs_dates = [date.today(),]
sim7m = bulk_SB_simulation.BulkSimulation7m(SB_list='7M_SBs_2024-06-12.csv',
                                            custom_SB_filter=sb_selector)
sim7m.run_simulations(obs_dates=obs_dates)
sim7m.write_failed_and_skipped_simulations_to_csvfile(
                              filepath='failed_and_skipped_simulations_7m.csv')
sim7m.write_statistics('statistics_7m.txt')

'''
'''
#check which 7m SBs require TP antennas:
# def project_selector(data):
#     return data['Project Code'] == '2023.1.01358.S'
check_TP = bulk_SB_simulation.CheckTP_for_7m_SBs(SB_list='7M_SBs_2024-06-12.csv',
                                                  obs_date=date.today())
check_TP.check_TP(check_results_filepath='7m_needs_TP.csv',
                  statistics_filepath='7mTP_general_statistics.txt')
check_TP.write_needsTP_statistics(filepath='needsTP_statistics.txt')
'''

#simulate specific 7m project
# obs_dates = [date.today(),]
# array_configs = ['aca.cm10.pm3.cfg',]
# project_code = '2023.1.00578.S'
# def custom_filter_7m(data):
#     return data['Project Code'] == project_code
# sim7m = bulk_SB_simulation.BulkSimulation7m(
#                     SB_list='7m_SBs_2023-10-16.csv',
#                     custom_SB_filter=custom_filter_7m,array_configs=array_configs)
# sim7m.run_simulations(obs_dates=obs_dates)
# sim7m.write_failed_and_skipped_simulations_to_csvfile(
#                       filepath=f'failed_simulations_7m_{project_code}.csv')
# sim7m.write_statistics('statistics.txt')

#testing some individual SBs
# def custom_filter(data):
#     return data['Project Code'] in ('2022.1.00210.S','2023.1.00570.S')
# obs_dates = [date.today(),date(year=2023,month=10,day=20)]
# sim = bulk_SB_simulation.BulkSimulation12m(
#                           SB_list='12m_SBs_2023-10-05.csv',
#                           support_arcs=None,array_config='c43-8',
#                           custom_SB_filter=custom_filter)
# sim.run_simulations(obs_dates=obs_dates)
