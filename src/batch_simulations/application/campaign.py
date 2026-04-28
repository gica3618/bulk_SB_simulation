#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:53:28 2026

@author: gianni
"""


import logging
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
from batch_simulations.infrastructure.table_writer import TableWriter
from batch_simulations.application.single_SB_simulation import SingleSBSimulation
import pandas as pd
import datetime
from pathlib import Path
from collections import Counter


class SimulationCampaign:

    def __init__(self, name, sb_table, date):
        self.name = name
        self.sb_table = sb_table
        self.date = date

    def run(self,single_sb_sim_cls=SingleSBSimulation):
        self.sb_simulation_summaries = []
        start = datetime.datetime.now()
        for i, row in enumerate(self.sb_table.data.itertuples(), start=1):
            logging.info(f"Running SB {i}/{len(self.sb_table.data)}: {row.sbname} {row.code}")
            single_sb_sim = single_sb_sim_cls(
                                 project_code=row.code, sb_name=row.sbname,
                                 date=self.date,
                                 array_config_12m=self.sb_table.array_config_12m)
            sim_summary = single_sb_sim.simulate()
            sim_summary.analyse()
            self.sb_simulation_summaries.append(sim_summary)
        end = datetime.datetime.now()
        self.run_time = end-start


class CampaignResultFormatter:

    sb_table_keys_for_output = ("code","sbname","sb_uid","p2g_account","sb_state",
                                "sb_state_flag")
    writer = TableWriter()
    p2g_columns = TableWriter.column_order.copy()
    p2g_columns.remove("traceback_of_unexpected_error")

    def __init__(self,campaign):
        self.campaign = campaign
        self.check_consistency()

    def check_consistency(self):
        if not len(self.campaign.sb_table.data) == len(self.campaign.sb_simulation_summaries):
            raise RuntimeError

    def create_master_dataframe(self):
        logging.info("going to create master dataframe")
        rows = []
        for sb_sim_summary, sb_table_row in zip(self.campaign.sb_simulation_summaries,
                                                self.campaign.sb_table.data.itertuples()):
            row = self.build_row(sb_sim_summary=sb_sim_summary, sb_table_row=sb_table_row)
            rows.append(row)
        self.master_table = pd.DataFrame(rows)
        self.master_table.sort_values(by='code',inplace=True)

    def build_row(self, sb_sim_summary, sb_table_row):
        sb = sb_sim_summary.sb
        row = {}
        self.add_general_info(row=row, sb_table_row=sb_table_row)
        self.add_simulation_info(row=row, sb_sim_summary=sb_sim_summary)
        if sb_sim_summary.sb is not None:
            self.add_note_to_aod(row=row, sb=sb)
            self.add_hardcoded_calibrators(row=row, sb=sb)
            self.add_individual_calibrators(row=row,sb=sb)
            self.add_HA_limits(row=row, sb=sb)
        return row

    def add_general_info(self,row,sb_table_row):
        for key in self.sb_table_keys_for_output:
            row[key] = getattr(sb_table_row,key)

    @staticmethod
    def add_note_to_aod(row,sb):
        row["note_to_AoD"] = sb.metadata["note_to_AoD"]

    @staticmethod
    def add_hardcoded_calibrators(row, sb):
        hardcoded_calibrators = [c.cal_type for c in sb.calibrators if
                                 c.is_hardcoded]
        row["hardcoded_calibrators"] = ", ".join(hardcoded_calibrators)

    @staticmethod
    def add_individual_calibrators(row,sb):
        for cal_type in CALIBRATOR_TYPES:
            if cal_type in sb.cal_types:
                row[cal_type] = sb.get_calibrator(cal_type).source_name
            else:
                row[cal_type] = None

    @staticmethod
    def add_HA_limits(row, sb):
        for lim in ("min", "max"):
            row[f"{lim}_HA_OT"] = sb.OT_allowed_HA[lim].hour
        HA_DSA = sb.get_DSA_HA_limits()
        for lim in ("min", "max"):
            row[f"{lim}_HA_DSA"] = HA_DSA[lim].hour
        row["HA_is_restricted"] = (
            sb.OT_allowed_HA["min"] > HA_DSA["min"]
            or sb.OT_allowed_HA["max"] < HA_DSA["max"])

    def add_simulation_info(self, row, sb_sim_summary):
        if sb_sim_summary.unexpected_error is None:
            row["simulated_config"] = sb_sim_summary.config
            row["simulated_date"] = str(sb_sim_summary.date)
            row["simulations"] = self.build_per_HA_summary_string(
                                    simulation_results=sb_sim_summary.simulation_results,
                                    HAs=sb_sim_summary.HAs)
        else:
            unexpected_error = sb_sim_summary.unexpected_error
            row["simulations"] = unexpected_error["error_message"]
            row["traceback_of_unexpected_error"] = unexpected_error["full_traceback"]
        row["inspection_reasons"] = "; ".join(sb_sim_summary.analysis_result.inspection_reasons)
        row["should_be_Waiting"] = sb_sim_summary.analysis_result.should_be_Waiting

    @staticmethod
    def build_per_HA_summary_string(simulation_results, HAs):
        lines = []
        for result, HA in zip(simulation_results, HAs):
            if result.success:
                lines.append(f"{HA.hour:.3g}: success")
            else:
                lines.append(f"{HA.hour:.3g}: {result.fail_reason.error_summary}")
        return "\n".join(lines)

    def write_master_table(self,out_format,output_dir="."):
        filename = f"master_table_{self.campaign.name}.{out_format}"
        logging.info(f"going to write master table to disk (filename: {filename})")
        self.writer.write_table_to_disk(dataframe=self.master_table,filename=filename,
                                        output_dir=output_dir)

    def need_inspection_master_table_selection(self):
        return self.master_table["inspection_reasons"] != ""

    def write_table_for_P2G(self,out_format,output_dir="."):
        need_inspection = self.need_inspection_master_table_selection()
        out = self.master_table[self.p2g_columns][need_inspection]
        filename = f"p2g_table_{self.campaign.name}.{out_format}"
        logging.info(f"going to write P2G table to disk (filename: {filename})")
        self.writer.write_table_to_disk(dataframe=out,filename=filename,
                                        output_dir=output_dir)

    def summarize_number_of_simulated_SBs(self):
        n_simulated_SBs = len(self.campaign.sb_table.data)
        n_simulated_SBs_by_array = self.campaign.sb_table.get_number_of_12m_and_7m_SBs()
        n_12m = n_simulated_SBs_by_array["12m"]
        n_7m = n_simulated_SBs_by_array["7m"]
        return f"nb. of simulated SBs: {n_simulated_SBs} (12m: {n_12m}; 7m: {n_7m})\n"

    def summarize_unexpected_errors(self):
        unexpected_err = [sim_summary.unexpected_error for sim_summary in
                          self.campaign.sb_simulation_summaries]
        unexpected_err = [e for e in unexpected_err if e is not None]
        n_unexpected_err = len(unexpected_err)
        return f"{n_unexpected_err} SBs yielded an unexpected error\n"

    def summarize_number_of_inspections(self):
        need_inspection = self.need_inspection_master_table_selection()
        n_need_inspection = need_inspection.sum()
        return f"{n_need_inspection} SBs need inspection\n"

    @staticmethod
    def get_count_lines(reasons):
        counts = Counter(reasons)
        return "\n".join([f"{reason}: {count}" for reason,count in counts.most_common()]) + "\n"

    def summarize_inspection_reasons(self):
        inspection_reasons = []
        for sim_summary in self.campaign.sb_simulation_summaries:
            inspection_reasons += sim_summary.analysis_result.inspection_reasons
        return "inspection reasons:\n" + self.get_count_lines(inspection_reasons)

    @staticmethod
    def get_fail_reasons_per_SB(sb_simulation_summaries):
        fail_reasons = []
        for sim_summary in sb_simulation_summaries:
            sb_fail_reasons = []
            #iterate over HAs:
            for result in sim_summary.simulation_results:
                if not result.success:
                    if result.fail_reason.category == "missing calibrator":
                        sb_fail_reasons.append(result.fail_reason.error_summary)
                    else:
                        sb_fail_reasons.append(result.fail_reason.category)
            fail_reasons += list(set(sb_fail_reasons))
        return fail_reasons

    def summarize_fail_reasons(self):
        sim_summaries = {"all SBs": self.campaign.sb_simulation_summaries}
        need_inspection_summaries = [summary for summary in
                                     self.campaign.sb_simulation_summaries
                                     if summary.analysis_result.inspection_reasons]
        sim_summaries["SBs needing inspection"] = need_inspection_summaries
        out = ""
        for key,summaries in sim_summaries.items():
            out += f"number of SBs failing with reason X at at least one hour angle ({key}):\n"
            fail_reasons = self.get_fail_reasons_per_SB(summaries)
            out += self.get_count_lines(fail_reasons)
        return out

    def write_campaign_summary(self,output_dir="."):
        lines = [f"summary of campaign '{self.campaign.name}'\n",
                 f"simulated observing date: {self.campaign.date}\n",
                 f"campaign run time: {self.campaign.run_time.total_seconds()/3600} h\n"]
        lines.append(self.summarize_number_of_simulated_SBs())
        lines.append(self.summarize_unexpected_errors())
        lines.append(self.summarize_number_of_inspections())
        lines.append(self.summarize_inspection_reasons())
        lines.append(self.summarize_fail_reasons())
        filename = Path(f"campaign_{self.campaign.name}_summary.txt")
        with open(output_dir/filename, "w") as out:
            out.writelines(lines)


class CampaignRunner:

    def __init__(self,name, sb_table, date):
        self.name = name
        self.sb_table = sb_table
        self.date = date

    def run(self):
        campaign = SimulationCampaign(name=self.name, sb_table=self.sb_table,
                                      date=self.date)
        campaign.run()
        self.result = CampaignResultFormatter(campaign=campaign)
        self.result.create_master_dataframe()