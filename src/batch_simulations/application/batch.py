#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:53:28 2026

@author: gianni
"""

from pathlib import Path
import logging
from batch_simulations.domain.job_planner import JobPlanner
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
from batch_simulations.application.simulation_runner import SimulationRunner
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.infrastructure.table_writer import TableWriter
from batch_simulations.infrastructure.xml_downloader import XMLDownloader
from astropy.coordinates import Angle
from astropy import units as u
import traceback
import pandas as pd


class SingleSBSimulation:

    default_HA_step = Angle(1*u.hour)
    fine_HA_step = Angle(0.25*u.hour)    

    def __init__(self,project_code,sb_name,date,array_config_12m):
        self.project_code = project_code
        self.sb_name = sb_name
        self.date = date
        self.array_config_12m = array_config_12m
        self.xml_filepath = self.build_xml_filepath()

    def build_xml_filepath(self):
        filename = f"{self.project_code}_{self.sb_name}.xml".replace(" ", "_")
        return Path.cwd() / filename

    def get_array_config(self, sb):
        if not sb.is_7m:
            return self.array_config_12m
        return "aca.cm10.pm3.cfg" if sb.requires_TP else "7m"

    def download_xml(self):
        downloader = XMLDownloader(project_code=self.project_code,
                                   sb_name=self.sb_name,
                                   filepath=self.xml_filepath)
        downloader.download_xml()

    def prepare_sb(self):
        self.download_xml()
        sb = BuildSBFromXML.build(self.xml_filepath)
        return sb

    def simulate_HAs(self, planner, array_config, runner_cls=SimulationRunner):
        def run_grid(step):
            jobs = planner.HA_jobs(xml_filepath=self.xml_filepath,
                                   array_config=array_config,date=self.date,
                                   step=step)
            runner = runner_cls()
            return jobs, runner.run_jobs(jobs)
        jobs, results = run_grid(self.default_HA_step)
        if not all(r.success for r in results):
            logging.info("Retrying with finer HA grid")
            jobs, results = run_grid(self.fine_HA_step)
        return [j.HA for j in jobs], results

    def simulate(self):
        try:
            sb = self.prepare_sb()
            array_config = self.get_array_config(sb)
            planner = JobPlanner(sb=sb,default_HA_step=self.default_HA_step,
                                 fine_HA_step=self.fine_HA_step)
            HAs, results = self.simulate_HAs(
                                planner=planner,array_config=array_config)
            return SingleSBSimulationSummary(sb=sb,date=self.date,
                                             HAs=HAs,
                                             config=array_config,
                                             simulation_results=results,
                                             unexpected_error=None)
        except Exception as e:
            logging.exception("Unexpected error while simulating"
                              +f" {self.sb_name} ({self.project_code})\n"
                              +"full traceback:")
            full_trace = traceback.format_exc()
            logging.exception(full_trace)
            unexpected_error = {"error_message":repr(e),
                                "full_traceback":full_trace}
            return SingleSBSimulationSummary(sb=None,date=self.date,
                                             HAs=None,config=None,
                                             simulation_results=None,
                                             unexpected_error=unexpected_error)

    def __del__(self):
        self.xml_filepath.unlink(missing_ok=True)


class AnalysisResult:

    def __init__(self,inspection_reasons=None,should_be_Waiting=False):
        self.inspection_reasons = [] if inspection_reasons is None else inspection_reasons
        self.should_be_Waiting = should_be_Waiting

    def merge(self,other):
        inspection_reasons = list(set(self.inspection_reasons+other.inspection_reasons))
        should_be_Waiting = self.should_be_Waiting or other.should_be_Waiting
        return AnalysisResult(inspection_reasons=inspection_reasons,
                              should_be_Waiting=should_be_Waiting)


class SingleSBSimulationSummary:

    min_runnable_HA_amount = 1

    def __init__(self,sb,config,date,HAs,simulation_results,unexpected_error):
        self.sb = sb
        self.config = config
        self.date = date
        self.HAs = HAs
        self.simulation_results = simulation_results
        self.unexpected_error = unexpected_error

    def HA_widths(self):
        if len(self.HAs) == 1:
            return [0,]
        widths = []
        for i in range(len(self.HAs)):
            if i == 0:
                widths.append((self.HAs[1]-self.HAs[0]).hour / 2)
            elif i == len(self.HAs)-1:
                widths.append((self.HAs[-1]-self.HAs[-2]).hour / 2)
            else:
                widths.append((self.HAs[i]-self.HAs[i-1]).hour/2
                              + (self.HAs[i+1]-self.HAs[i]).hour/2)
        return widths

    def runnable_HA_amount(self):
        HA_widths = self.HA_widths()
        runnable_HA_widths = [HA_width for HA,HA_width,result in
                              zip(self.HAs,HA_widths,self.simulation_results)
                              if (result.success and HA>=self.sb.OT_allowed_HA["min"]
                                  and HA<=self.sb.OT_allowed_HA["max"])]
        return sum(runnable_HA_widths)

    def analyse_runnable_HA_range(self):
        runnable_HA_amount = self.runnable_HA_amount()
        logging.info(f"runnable HA amount: {runnable_HA_amount:.3g} hours")
        if runnable_HA_amount < self.min_runnable_HA_amount:
            return AnalysisResult(inspection_reasons=["runnable HA range is small"],
                                  should_be_Waiting=False)
        return AnalysisResult()

    def analyse_HA_restriction(self):
        if self.sb.any_calibrator_hardcoded():
            #if a simulation is successful, I can only conclude that it is indeed fine
            #if there are no hardcoded calibrators; so for SBs with hardcoded
            #calibrators, I cannot say that HA restrictions should be lifted
            #even if simulation runs fine
            return AnalysisResult()
    
        unnecessarily_restricted_HAs = []
        for HA,result in zip(self.HAs,self.simulation_results):
            HA_is_excluded = HA < self.sb.OT_allowed_HA["min"] or\
                             HA > self.sb.OT_allowed_HA["max"]
            if result.success and HA_is_excluded:
                unnecessarily_restricted_HAs.append(HA)
    
        if not unnecessarily_restricted_HAs:
            logging.info("HA is not unnecessarily restricted")
            return AnalysisResult()
    
        ha_str = ", ".join(str(HA.hour) for HA in unnecessarily_restricted_HAs)
        logging.info(f"HA is unnecessarily restricted at {ha_str}")
        return AnalysisResult(inspection_reasons=[f"unnecessarily restricted HAs: {ha_str}"],
                              should_be_Waiting=False)

    def HA_is_allowed(self,HA):
        return HA >= self.sb.OT_allowed_HA["min"] and HA <= self.sb.OT_allowed_HA["max"]

    def analyse_simulation_failures(self):
        #logic implemented here:
        # - any error at disallowed HA does not trigger anything
        # - "unobservable" never triggers anything
        # - "server error" triggers inspection if occuring at allowed HA
        # - any other error occuring at allowed HA triggers inspection and Waiting
        server_error = False
        failed_and_needs_action = False
        should_be_Waiting = False
    
        for HA, result in zip(self.HAs, self.simulation_results):
            if result.success:
                continue
            
            if not self.HA_is_allowed(HA):
                continue
            
            if result.fail_reason.category == "unobservable":
                continue
    
            if result.fail_reason.category == "server error":
                #server error needs inspection, but does not need to be waiting
                server_error = True
                continue
    
            failed_and_needs_action = True
            should_be_Waiting = True
    
        inspection_reasons = []
        if server_error:
            inspection_reasons.append("server error")
        if failed_and_needs_action:
            inspection_reasons.append("simulation failure")
        return AnalysisResult(inspection_reasons=inspection_reasons,
                              should_be_Waiting=should_be_Waiting)

    def analyse(self):
        if self.unexpected_error is not None:
            self.analysis_result = AnalysisResult(["unexpected error"],
                                                  should_be_Waiting=True)
        else:
            self.analyse_no_unexpected_error()

    def analyse_no_unexpected_error(self):
        self.analysis_result = AnalysisResult()
        for ana in (self.analyse_runnable_HA_range(),
                    self.analyse_HA_restriction(),
                    self.analyse_simulation_failures()):
            self.analysis_result = self.analysis_result.merge(ana)


class SimulationCampaign:

    def __init__(self, name, sb_table, date):
        self.name = name
        self.sb_table = sb_table
        self.date = date

    def run(self,single_sb_sim_cls=SingleSBSimulation):
        self.sb_simulation_summaries = []
        for i, row in enumerate(self.sb_table.data.itertuples(), start=1):
            logging.info(f"Running SB {i}/{len(self.sb_table.data)}: {row.sbname} {row.code}")
            single_sb_sim = single_sb_sim_cls(
                                 project_code=row.code, sb_name=row.sbname,
                                 date=self.date,
                                 array_config_12m=self.sb_table.array_config_12m)
            sim_summary = single_sb_sim.simulate()
            sim_summary.analyse()
            self.sb_simulation_summaries.append(sim_summary)


class CampaignResultFormatter:

    sb_table_keys_for_output = ("code","sbname","sb_uid","p2g_account","sb_state",
                                "sb_state_flag")
    writer = TableWriter()
    p2g_columns = TableWriter.column_order.copy()
    p2g_columns.remove("traceback of unexpected error")

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
            row["traceback of unexpected error"] = unexpected_error["full_traceback"]
        row["inspection reasons"] = "; ".join(sb_sim_summary.analysis_result.inspection_reasons)
        row["should be Waiting"] = sb_sim_summary.analysis_result.should_be_Waiting

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

    def write_table_for_P2G(self,out_format,output_dir="."):
        out = self.master_table[self.p2g_columns]
        need_inspection = out["inspection reasons"] != ""
        filename = f"p2g_table_{self.campaign.name}.{out_format}"
        logging.info(f"going to write P2G table to disk (filename: {filename})")
        self.writer.write_table_to_disk(dataframe=out[need_inspection],filename=filename,
                                        output_dir=output_dir)


class CampaignRunner:

    def __init__(self,name, sb_table, date):
        self.name = name
        self.sb_table = sb_table
        self.date = date

    def run(self):
        self.campaign = SimulationCampaign(name=self.name, sb_table=self.sb_table,
                                           date=self.date)
        self.campaign.run()
        self.campaign_result = CampaignResultFormatter(campaign=self.campaign)
        self.campaign_result.create_master_dataframe()