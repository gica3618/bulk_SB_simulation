#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:53:28 2026

@author: gianni
"""

from pathlib import Path
import logging
from batch_simulations.domain.job_planner import JobPlanner
from batch_simulations.domain.sb import SB
from batch_simulations.domain.calibrator import CALIBRATOR_TYPES
from batch_simulations.application.simulation_runner import SimulationRunner
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.infrastructure.xml_downloader import XMLDownloader
from astropy.coordinates import Angle
from astropy import units as u
from dataclasses import dataclass
import traceback
import datetime
import pandas as pd


@dataclass
class SingleSBSimulationSummary:
    sb: SB | None
    config: str | None
    date: datetime.date
    HAs: list | None
    simulation_results: list | None
    unexpected_error: None | str


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
        XMLDownloader.download_xml(project_code=self.project_code,
                                   sb_name=self.sb_name,
                                   filepath=self.xml_filepath)

    def prepare_sb(self):
        self.download_xml()
        sb = BuildSBFromXML.build(self.xml_filepath)
        return sb

    def simulate_HAs(self, planner, array_config, runner_cls=SimulationRunner):
        def run_grid(step):
            jobs = planner.jobs(xml_filepath=self.xml_filepath,
                                array_config=array_config,date=self.date,step=step)
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
            planner = JobPlanner(sb=sb,default_step=self.default_HA_step,
                                 fine_step=self.fine_HA_step)
            HAs, results = self.simulate_HAs(
                                planner=planner,array_config=array_config)
            return SingleSBSimulationSummary(sb=sb,date=self.date,
                                             HAs=HAs,
                                             config=array_config,
                                             simulation_results=results,
                                             unexpected_error=None)
        except Exception:
            logging.exception("Unexpected error while simulating"
                              +f" {self.sb_name} ({self.project_code})")
            full_trace = traceback.format_exc()
            logging.exception(full_trace)
            return SingleSBSimulationSummary(sb=None,date=self.date,
                                             HAs=None,config=None,
                                             simulation_results=None,
                                             unexpected_error=full_trace)


class SimulationCampaign:

    def __init__(self, name, sb_table, date):
        self.name = name
        self.sb_table = sb_table
        self.date = date

    def run(self,single_sb_sim_cls=SingleSBSimulation):
        self.sb_simulation_summaries = []
        for i, row in enumerate(self.sb_table.data.itertuples(), start=1):
            logging.info("Running SB {i}/{len(self.sb_table)}: {row.sbname} {row.code}")
            single_sb_sim = single_sb_sim_cls(
                                 project_code=row.code, sb_name=row.sbname,
                                 date=self.date,
                                 array_config_12m=self.sb_table.array_config_12m)
            sim_summary = single_sb_sim.simulate()
            self.sb_simulation_summaries.append(sim_summary)


class CampaignResultFormatter:

    sb_table_keys_for_output = ("code","sbname","sb_uid","p2g_account","sb_state",
                                "sb_state_flag")

    def __init__(self,campaign):
        self.campaign = campaign
        self.check_consistency()

    def check_consistency(self):
        if not len(self.campaign.sb_table) == len(self.campaign.sb_simulation_summaries):
            raise RuntimeError

    def create_master_dataframe(self):
        rows = []
        for sb_sim_summary, sb_table_row in zip(self.campaign.sb_simulation_summaries,
                                                self.campaign.sb_table.itertuples()):
            row = self.build_row(sb_sim_summary=sb_sim_summary, sb_table_row=sb_table_row)
            rows.append(row)
        self.master_table = pd.DataFrame(rows)

    def build_row(self, sb_sim_summary, sb_table_row):
        sb = sb_sim_summary.sb
        row = {}
        self.add_general_info(row=row, sb_table_row=sb_table_row)
        self.add_note_to_aod(row=row, sb=sb)
        self.add_hardcoded_calibrators(row=row, sb=sb)
        self.add_individual_calibrators(row=row,sb=sb)
        self.add_HA_limits(row=row, sb=sb)
        self.add_simulation_info(row=row, sb_sim_summary=sb_sim_summary)
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
        row["simulated_config"] = sb_sim_summary.config
        row["simulated_date"] = str(self.campaign.date)
        row["simulations"] = self.build_simulation_result_summary(
                                          simulation_results=sb_sim_summary.simulation_results,
                                          HAs=sb_sim_summary.HAs)

    @staticmethod
    def build_simulation_result_summary(simulation_results, HAs):
        lines = []
        for result, HA in zip(simulation_results, HAs):
            if result.success:
                lines.append(f"{HA.hour:.3g}: success")
            else:
                lines.append(f"{HA.hour:.3g}: {result.fail_reason}")
        return "\n".join(lines)