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
import os

logging_level = logging.INFO
#logging_level = logging.ERROR

logfile = "batch_simulations.log"
if os.path.isfile(logfile):
    os.remove(logfile)
handlers = [logging.FileHandler(logfile),logging.StreamHandler(sys.stdout)]
logging.basicConfig(format='%(levelname)s %(asctime)s: %(message)s',level=logging_level,
                    handlers=handlers)


'''
#12m SBs
# def custom_SB_filter(data):
#     return data['sb_state'] == 'Waiting' and\
#                              data['sb_state_flag'] == 'ForCalibrator'
# def custom_SB_filter(data):
#     return data['sbname'] == "G24.01_a_06_TM1"

array_config = 'c43-3'
obs_dates = [date(year=2026,month=1,day=10),date(year=2026,month=1,day=31)]
sim12m = bulk_SB_simulation.BulkSimulation12m(
                          SB_list='configuration_lookup_table_cycle_12_20250104.csv',
                          support_arcs=None,array_config=array_config,
                          custom_SB_filter=None)
sim12m.run_simulations(obs_dates=obs_dates)
sim12m.write_failed_and_skipped_simulations_to_csvfile(
                  filepath=f'failed_and_skipped_simulations_12m_{array_config}.csv')
sim12m.write_statistics(f'statistics_{array_config}.txt')
'''

#7m SBs
# def sb_selector(data):
#     # return data['State'] == 'Waiting' and data['SubState'] == 'ForCalibrator'
#     project_list = ['2023.A.00038.T','2023.A.00038.T',"2023.A.00038.T",
#                     "2023.A.00038.T","2023.A.00038.T","2023.A.00038.T","2024.1.00137.T",
#                     "2024.1.00137.T","2024.1.00137.T","2024.1.00137.T","2024.1.00137.T",
#                     "2024.1.00477.T","2024.1.00477.T","2024.1.00477.T","2024.1.00477.T",
#                     "2024.1.00477.T","2024.1.00477.T","2024.1.00477.T","2024.1.00477.T",
#                     "2024.1.00477.T","2024.A.00044.T","2024.A.00044.T","2024.A.00044.T",
#                     "2024.A.00044.T","2024.A.00044.T","2024.A.00044.T","2024.A.00044.T",
#                     "2024.A.00044.T","2024.A.00044.T","2024.A.00044.T","2024.A.00044.T",
#                     "2025.1.00552.S","2025.1.00762.S","2025.1.00891.T","2025.1.00891.T",
#                     "2025.1.00891.T","2025.1.00891.T","2025.1.00891.T","2025.1.00915.S",
#                     "2025.1.00915.S","2025.1.00915.S","2025.1.00915.S","2025.1.00915.S",
#                     "2025.1.00915.S","2025.1.00915.S","2025.1.00915.S","2025.1.01039.S",
#                     "2025.1.01084.S","2025.1.01539.S"]
#     sb_list = ["T_CrB_e_03_7M","T_CrB_e_06_7M","T_CrB_f_03_7M","T_CrB_f_06_7M",
#                "T_CrB_g_03_7M","T_CrB_g_06_7M","Tsuchins_f_06_7M","Tsuchins_g_06_7M",
#                "Tsuchins_h_06_7M","Tsuchins_i_06_7M","Tsuchins_j_06_7M","ToO_Come_a_07_7M",
#                "ToO_Come_b_06_7M","ToO_Come_b_07_7M","ToO_Come_c_06_7M","ToO_Come_d_06_7M",
#                "ToO_Come_e_06_7M","ToO_Come_f_06_7M","ToO_Come_g_06_7M","ToO_Come_h_06_7M",
#                "V462_Lup_f_01_7M","V462_Lup_f_04_7M","V462_Lup_f_07_7M","V462_Lup_g_01_7M",
#                "V462_Lup_g_04_7M","V462_Lup_g_07_7M","V462_Lup_h_01_7M","V462_Lup_h_04_7M",
#                "V462_Lup_h_07_7M","V462_Lup_i_01_7M","V462_Lup_i_04_7M","UGC_0152_a_06_7M",
#                "T_Tau_a_09_7M","T_CrB_e_06_7M","T_CrB_f_03_7M","T_CrB_f_06_7M",
#                "T_CrB_g_03_7M","T_CrB_g_06_7M","UGC02229_a_06_7M","UGC02229_a_07_7M",
#                "UGC04145_a_06_7M","UGC04197_a_07_7M","UGC08322_a_06_7M","UGC08322_a_07_7M",
#                "UGC12518_a_06_7M","UGC12518_a_07_7M","UGC_0948_a_06_7M","UM_384_a_06_7M",
#                "SMC_area_a_08_7M"]
#     assert len(project_list) == len(sb_list)
#     if data['SB Name'] in sb_list:
#         for i,sb in enumerate(sb_list):
#             if sb == data['SB Name']:
#                 if data["Project Code"] == project_list[i]:
#                     return True
#     return False

obs_dates = [date.today(),]
sim7m = bulk_SB_simulation.BulkSimulation7m(SB_list='2026-01-09_schedBlockList.csv',
                                            custom_SB_filter=None)
sim7m.run_simulations(obs_dates=obs_dates)
sim7m.write_failed_and_skipped_simulations_to_csvfile(
                              filepath='failed_and_skipped_simulations_7m.csv')
sim7m.write_statistics('statistics_7m.txt')

'''
#check which 7m SBs require TP antennas:
# def project_selector(data):
#     return data['Project Code'] == '2023.1.01358.S'
check_TP = bulk_SB_simulation.CheckTP_for_7m_SBs(SB_list='7m_2025-08-28.csv',
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
