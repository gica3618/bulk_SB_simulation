#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 09:05:57 2023

@author: gianni
"""


    # def SB_xml_summary(self):
    #     #TODO needs a test
    #     summary = {}
    #     summary["note_to_AoD"] = self.xml.get_NotetoAoD()
    #     summary["hardcoded_calibrator(s)"] = any([c.is_hardcoded for c in
    #                                              self.calibrators])
    #     for cal_type in calibrator.calibrator_types:
    #         if cal_type not in self.cal_types:
    #             summary[cal_type] = None
    #         else:
    #             cal = self.get_calibrator(cal_type=cal_type)
    #             summary[cal_type] = cal.source_name
    #     OT_allowed_HA = self.xml.read_allowed_HA()
    #     for limit in ("min","max"):
    #         summary[f"{limit}_HA_OT"] = round(OT_allowed_HA[limit].hour,6)
    #         summary[f"{limit}_HA_DSA"] = round(getattr(self,f"{limit}_HA_DSA").hour,6)
    #     return summary


    # def output_for_p2g(self):
    #     output = []
    #     for executed_sim in self.executed_simulations:
    #         out_data = {key:executed_sim.sb_metadata[key] for key in
    #                     ("sb_uid","code","p2g_account","sb_state","sb_state_flag")}
    #         out_data["sb_error"] = executed_sim.sb_error
    #         for obs_date in executed_sim.obs_dates:
    #             date_key = obs_date.isoformat()
    #             result_summary_str = ""
    #             sim_results = executed_sim.simulation_results[date_key]
    #             for HA,sim_result in zip(executed_sim.HAs[date_key],sim_results):
    #                 if sim_result.success:
    #                     result_summary_str += "success\n"
    #                 else:
    #                     result_summary_str += f"{HA.hour:.3g}h: {sim_result.fail_reason}\n"
    #             out_data[f"results {date_key}"] = result_summary_str
    #             executed_commands = [sim_result.executed_command for
    #                                  sim_result in sim_results]
    #             out_data[f"simulations {date_key}"] = "\n".join(executed_commands)
    #         output.append(out_data)


if __name__ == '__main__':
    logging_level = logging.INFO
    #logging_level = logging.ERROR
    logging.basicConfig(format='%(levelname)s: %(message)s',level=logging_level,
                        stream=sys.stdout)

    # test_SB = SB(project_code='2023.1.01430.S',name='MMS_1_a_08_TM2',
    #              array_config='c43-3')
    # test_SB.run_simulations()
    
    input_filepaths = {"12m":"tests/test_input_tables/configuration_lookup_table_cycle_12.csv",
                       "7m":"tests/test_input_tables/2026-02-06_schedBlockList_7M.csv"}
    for array_config_12m,number in zip(("c43-3","c43-10"),(3,10)):
        batch_sim = BatchSimulations(
                       input_filepaths=input_filepaths,array_config_12m=array_config_12m)