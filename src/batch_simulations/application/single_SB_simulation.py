#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:33:34 2026

@author: gianni
"""

from astropy.coordinates import Angle
from astropy import units as u
from pathlib import Path
import logging
import traceback
from batch_simulations.infrastructure.xml_downloader import XMLDownloader
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.application.simulation_runner import SimulationRunner
from batch_simulations.domain.job_planner import JobPlanner


def get_unnecessarily_restricted_HAs(has_hardcoded_cals,HAs,results,OT_allowed_HA):
    unnecessarily_restricted_HAs = []
    if has_hardcoded_cals:
        #if a simulation is successful, I can only conclude that it is indeed fine
        #if there are no hardcoded calibrators; so for SBs with hardcoded
        #calibrators, I cannot say that HA restrictions should be lifted
        #even if simulation runs fine
        logging.info("SB has hardcoded calibrators, cannot determine if HAs are "
                     +"unnecessarily restricted")
        return unnecessarily_restricted_HAs
    for HA,result in zip(HAs,results):
        HA_is_excluded = (HA < OT_allowed_HA["min"]) or (HA > OT_allowed_HA["max"])
        if result.success and HA_is_excluded:
            unnecessarily_restricted_HAs.append(HA)
    if not unnecessarily_restricted_HAs:
        logging.info("HA range is not unnecessarily restricted")
    else:
        logging.info("unnecessarily restricted HA(s): "
                     +f"{[HA.hour for HA in unnecessarily_restricted_HAs]}")
    return unnecessarily_restricted_HAs


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

    @staticmethod
    def simulate_again_with_fine_HA_step(results,HAs,has_hardcoded_cals,
                                         OT_allowed_HA):
        if not all(r.success for r in results):
            logging.info("should run again with finer HA step because some "
                         +"simulations failed")
            return True
        unnec_res_HAs = get_unnecessarily_restricted_HAs(
                            has_hardcoded_cals=has_hardcoded_cals,
                            HAs=HAs, results=results, OT_allowed_HA=OT_allowed_HA)
        if unnec_res_HAs:
            logging.info("should run again with finer HA step because of "
                         +"unnecessarily restricted HA(s)")
            return True
        return False

    def simulate_HAs(self, planner, array_config, sb, runner_cls=SimulationRunner):
        def run_grid(job_creator):
            HAs,jobs = job_creator(xml_filepath=self.xml_filepath,
                                   array_config=array_config,date=self.date)
            runner = runner_cls()
            return HAs, runner.run_jobs(jobs)
        HAs, results = run_grid(job_creator=planner.HA_jobs_default_HA_step)
        should_run_again = self.simulate_again_with_fine_HA_step(
                               results=results, HAs=HAs,
                               has_hardcoded_cals=sb.any_calibrator_hardcoded(),
                               OT_allowed_HA=sb.OT_allowed_HA)
        if should_run_again:
            logging.info("Simulating again with finer HA grid")
            HAs, results = run_grid(job_creator=planner.HA_jobs_fine_HA_step)
        return HAs, results

    def simulate(self):
        try:
            sb = self.prepare_sb()
            array_config = self.get_array_config(sb)
            planner = JobPlanner(sb=sb,default_HA_step=self.default_HA_step,
                                 fine_HA_step=self.fine_HA_step)
            HAs, results = self.simulate_HAs(
                                planner=planner,array_config=array_config,sb=sb)
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
        unnecessarily_restricted_HAs = get_unnecessarily_restricted_HAs(
                                          has_hardcoded_cals=self.sb.any_calibrator_hardcoded(),
                                          HAs=self.HAs, results=self.simulation_results,
                                          OT_allowed_HA=self.sb.OT_allowed_HA)
        if not unnecessarily_restricted_HAs:
            logging.info("HA is not unnecessarily restricted")
            return AnalysisResult()
        ha_str = ", ".join(str(HA.hour) for HA in unnecessarily_restricted_HAs)
        logging.info(f"HA is unnecessarily restricted at following HA(s): {ha_str}")
        return AnalysisResult(inspection_reasons=[f"unnecessarily restricted HA(s): {ha_str}"],
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
