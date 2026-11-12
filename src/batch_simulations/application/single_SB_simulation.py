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
import itertools
from batch_simulations.infrastructure.xml_downloader import XMLDownloader
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.application.simulation_runner import SimulationRunner
from batch_simulations.domain.job_planner import JobPlanner
import batch_simulations.utils.angle_comparisons as ac



        # logging.info("should run again with finer HA step because some "
        #              +"simulations failed")


# def get_unnecessarily_restricted_HAs(has_hardcoded_cals,HAs,results,OT_allowed_HA):
#     unnecessarily_restricted_HAs = []
#     if has_hardcoded_cals:
#         #if a simulation is successful, I can only conclude that it is indeed fine
#         #if there are no hardcoded calibrators; so for SBs with hardcoded
#         #calibrators, I cannot say that HA restrictions should be lifted
#         #even if simulation runs fine
#         logging.info("SB has hardcoded calibrators, cannot determine if HAs are "
#                      +"unnecessarily restricted")
#         return unnecessarily_restricted_HAs
#     #tol = 1e-10 #rad; needed because of machine precision rounding errors
#     for HA,result in zip(HAs,results):
#         HA_is_excluded = not HA_is_allowed(HA=HA,OT_allowed_HA=OT_allowed_HA)
#         # HA_is_excluded = (HA.rad + tol < OT_allowed_HA["min"].rad)\
#         #                  or (HA.rad-tol > OT_allowed_HA["max"].rad)
#         if result.success and HA_is_excluded:
#             unnecessarily_restricted_HAs.append(HA)
#     if not unnecessarily_restricted_HAs:
#         logging.info("HA range is not unnecessarily restricted")
#     else:
#         logging.info("unnecessarily restricted HA(s): "
#                      +f"{[HA.hour for HA in unnecessarily_restricted_HAs]}")
#     return unnecessarily_restricted_HAs


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

    def build_sb(self):
        self.download_xml()
        sb = BuildSBFromXML.build(self.xml_filepath)
        return sb

    @staticmethod
    def any_simulation_failed(results):
        return any( [(not r.success) for r in results] )

    def should_run_again_with_fine_step(self,step,results):
        return ac.greater(step,self.fine_HA_step) and self.any_simulation_failed(results)

    #planner_cls and runner_cls are used to make testing easier
    def simulate_HAs(self, array_config, sb, planner_cls=JobPlanner,
                     runner_cls=SimulationRunner):
        planner = planner_cls(DSA_HA_limits=sb.get_DSA_HA_limits())
        #dsa = DSAPolicy()
        def run_grid(step):
            HAs,jobs = planner.HA_jobs(xml_filepath=self.xml_filepath,
                                       array_config=array_config,date=self.date,
                                       step=step)
            runner = runner_cls()
            results = []
            for HA,job in zip(HAs,jobs):
                # reasons_for_not_showing_up_in_DSA = dsa.reasons_for_not_showing_up(
                #                                        sb=sb,start_HA=HA)
                # if reasons_for_not_showing_up_in_DSA:
                #     #reason for doing this: OSS checks calibrator availability before
                #     #the elevation of science target. If calibrator is missing for an HA
                #     #where science target is anyway out of elevation limits, this leads
                #     #to false positives: the SB is marked as needing P2G action,
                #     #but actually no action is needed since SB would never show
                #     #up at that HA since science target is not visible
                #     message = ("simulation skipped, SB would not show up in DSA: "
                #                +"; ".join(reasons_for_not_showing_up_in_DSA))
                #     fail_reason = FailReason(error_message=message,
                #                              error_summary=message,
                #                              category="unobservable")
                #     result = SimulationResult(executed_command=None,
                #                               output_folder=None,
                #                               xml_filename=Path(self.xml_filepath).name,
                #                               success=False,
                #                               fail_reason=fail_reason)
                # else:
                #     result = runner.run(job)
                result = runner.run(job)
                results.append(result)
            return HAs, results
        if sb.HA_is_restricted():
            initial_step = self.fine_HA_step
        else:
            initial_step = self.default_HA_step
        HAs, results = run_grid(step=initial_step)
        if self.should_run_again_with_fine_step(step=initial_step,results=results):
            logging.info("Simulating again with finer HA grid")
            HAs, results = run_grid(step=self.fine_HA_step)
        return HAs, results

    def simulate(self):
        try:
            sb = self.build_sb()
            array_config = self.get_array_config(sb)
            HAs, results = self.simulate_HAs(array_config=array_config,sb=sb)
            return SingleSBSimulationAnalysis(sb=sb,HAs=HAs,
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
            return SingleSBSimulationAnalysis(sb=None,HAs=None,
                                             simulation_results=None,
                                             unexpected_error=unexpected_error)

    def __del__(self):
        self.xml_filepath.unlink(missing_ok=True)


class AnalysisResult:

    state_hierarchy = ("WaitingForP2G","Waiting","Ready",None)

    def __init__(self):
        self.inspection_reasons = set()
        self.recommended_HA_limits = None
        self.recommended_state = self.state_hierarchy[-1]
    
    def suggest_state(self,state):
        if state not in self.state_hierarchy:
            raise ValueError(f"state {state} cannot be suggested")
        if self.state_hierarchy.index(state) < self.state_hierarchy.index(self.recommended_state):
            self.recommended_state = state


class SingleSBSimulationAnalysis:

    concerningly_small_runnable_HA_amount = Angle(1*u.hour)
    max_removable_successful_HA_amount = Angle(1*u.hour)
    min_remaining_HA_amount_when_excluding_HAs = Angle(2*u.hour)

    def __init__(self,sb,HAs,simulation_results,unexpected_error):
        self.sb = sb
        self.HAs = HAs
        self.simulation_results = simulation_results
        self.unexpected_error = unexpected_error

    def analyse(self):
        #I need to consider different categories of errors
        #for example, if at any HA the SB fails because of >2h, it needs to go
        #to WaitingForP2G regardless of what happens at other HAs
        #similarly if it fails with category "other". Actually,
        #in those cases, I should not attempt to restrict the HA, because I do not
        #know what is happening at those HAs
        #about state changes, I will only change Ready/WaitingForCalibrator -> WaitingForP2G,
        #and WaitingForCalibrator -> Ready (in particular, I do not want to change
        #the state of SBs in WaitingForP2G and WaitingForJAO)
        #Also, remember that HA restrictions cannot be relaxed if there are
        #any hardcoded calibrators. they can in principle be tightned, but
        #maybe that is not a good idea, because once tightened, they cannot be
        #relaxed anymore in an automated way.
        
        #cases that I should recommend Ready ((Ready,WaitingForCalibrator -> Ready;
        # WaitingForP2G -> WaitingForP2G):
        # - no simulation failure and no hardcoded calibrator
        # - SBs where HA restrictions can fix the SB, and no hardcoded calibrators
        #cases that I should recommend WaitingForP2G (i.e. need action;
        #(Ready, WaitingForCalibrator -> WaitingForP2G)):
        # - unexpected error, or SB failed at any HA with "other", "server error" or "exceeds 2h limit"
        #cases that I should recommend Waiting; (Ready -> WaitingForP2G,
        #WaitingForCalibrator -> WaitingForCalibrator, WaitingForP2G -> WaitingForP2G)
        # - missing calibrator at any allowed HA and (has hardcoded calibrator or
        # HA restriction cannot fix it) (I do not want to touch HA restrictions for
        #SBs with hardcoded calibrators)

        self.analysis_result = AnalysisResult()
        if self.unexpected_error is not None:
            self.analyse_unexpected_error()
        else:
            self.analyse_regular()

    def analyse_unexpected_error(self):
        self.analysis_result.inspection_reasons.add("unexpected error")
        self.analysis_result.suggest_state("WaitingForP2G")

    def analyse_regular(self):
        self.add_inspection_reasons()

        fail_categories = [result.fail_reason.category for result in 
                           self.simulation_results if not result.success]
        fail_categories = set(fail_categories)
        
        #case without failures, or only "unobservable" type failures
        #if PolCal is hardcoded, I allow suggesting Ready, but if any other calibrator
        #is hardcoded, I do not allow it
        if fail_categories.issubset(set(["unobservable"]))\
                                and not self.sb.any_calibrator_except_PolCal_hardcoded():
            self.analysis_result.recommended_HA_limits = {"min":Angle(-12*u.hour),
                                                          "max":Angle(12*u.hour)}
            self.analysis_result.suggest_state("Ready")

        #cases needing inspection / action from P2G
        for fail_cat in ("exceeds 2h limit","server error","other"):
            if fail_cat in fail_categories:
                self.analysis_result.inspection_reasons.add(fail_cat)
                self.analysis_result.suggest_state("WaitingForP2G")

        #SBs with missing calibrator
        if "missing calibrator" in fail_categories:
            self.analyse_missing_calibrator()

    def add_inspection_reasons(self):
        if self.runnable_HA_range_is_small():
            self.analysis_result.inspection_reasons.add("runnable HA range is small")
        if self.extreme_HA_restriction_imposed_by_P2G():
            self.analysis_result.inspection_reasons.add(
                                    "extreme HA restrictions imposed by P2G")

    def analyse_missing_calibrator(self):
        if self.sb.any_calibrator_hardcoded():
            self.analyse_missing_calibrator_hardcoded()
        else:
            self.analyse_missing_calibrator_no_hardcoded()

    def analyse_missing_calibrator_hardcoded(self):
        #I do not attempt to change HA restrictions of SBs with
        #hardcoded calibrators
        if self.missing_calibrator_at_allowed_HA():
            self.analysis_result.inspection_reasons.add("missing calibrator")
            self.analysis_result.suggest_state("Waiting")

    def analyse_missing_calibrator_no_hardcoded(self):
        recommended_HA_limits = self.get_recommended_HA_limits()
        if recommended_HA_limits is None:
            #cannot fix by changing HA limits
            self.analysis_result.inspection_reasons.add("missing calibrator")
            if self.missing_calibrator_at_allowed_HA():
                self.analysis_result.suggest_state("Waiting")
        else:
            self.analysis_result.recommended_HA_limits = recommended_HA_limits
            self.analysis_result.suggest_state("Ready")

    def HA_is_allowed(self,HA):
        return (ac.geq(HA,self.sb.OT_allowed_HA["min"])
                  and ac.leq(HA,self.sb.OT_allowed_HA["max"]))

    def get_missing_calibrator_HAs(self):
        missing_calibrator_HAs = []
        for HA,result in zip(self.HAs,self.simulation_results):
            if not result.success:
                if result.fail_reason.category == "missing calibrator":
                    missing_calibrator_HAs.append(HA)
        return missing_calibrator_HAs

    def missing_calibrator_at_allowed_HA(self):
        missing_calibrator_HAs = self.get_missing_calibrator_HAs()
        return any([self.HA_is_allowed(HA=HA) for HA in missing_calibrator_HAs])

    def HA_widths(self):
        if len(self.HAs) == 1:
            return [Angle(0*u.hour),]
        widths = []
        for i in range(len(self.HAs)):
            if i == 0:
                widths.append((self.HAs[1]-self.HAs[0]) / 2)
            elif i == len(self.HAs)-1:
                widths.append((self.HAs[-1]-self.HAs[-2]) / 2)
            else:
                widths.append((self.HAs[i]-self.HAs[i-1])/2
                              + (self.HAs[i+1]-self.HAs[i])/2)
        return widths

    def runnable_HA_amount(self):
        HA_widths = self.HA_widths()
        runnable_HA_widths = [HA_width for HA_width,result in
                              zip(HA_widths,self.simulation_results)
                              if result.success]
        if not runnable_HA_widths:
            return Angle(0*u.hour)
        else:
            return sum(runnable_HA_widths)

    def runnable_HA_range_is_small(self):
        runnable_HA_amount = self.runnable_HA_amount()
        logging.info(f"runnable HA amount: {runnable_HA_amount:.3g} hours")
        return ac.less(runnable_HA_amount, self.concerningly_small_runnable_HA_amount)

    def extreme_HA_restriction_imposed_by_P2G(self):
        allowed_width = self.sb.OT_allowed_HA["max"] - self.sb.OT_allowed_HA["min"]
        return ac.less(allowed_width,self.concerningly_small_runnable_HA_amount)

    def get_sim_results(self,HA_limits):
        return [r for r,HA in zip(self.simulation_results,self.HAs)
                if ac.geq(HA,HA_limits["min"]) and ac.leq(HA,HA_limits["max"])]

    def HA_range_contains_missing_calibrator(self,HA_limits):
        results = self.get_sim_results(HA_limits=HA_limits)
        for r in results:
            if r.fail_reason is not None:
                if r.fail_reason.category == "missing calibrator":
                    return True
        return False

    # def HA_range_has_issues(self,HA_limits):
    #     #TODO test needed
    #     results = self.get_sim_results(HA_limits=HA_limits)
    #     if not results:
    #         return
    #     for result in results:
    #         if not result.success:
    #             if result.fail_reason.category != "unobservable":
    #                 return True
    #     return False

    # def determine_HA_limits_necessary_to_fix_SB(self):
    #     best = {}
    #     for HA_min,HA_max in itertools.product(self.HAs,self.HAs):
    #         #TODO machine precision issue?
    #         if HA_min >= HA_max:
    #             continue
    #         HA_limits = {"min": HA_min, "max": HA_max}
    #         has_issues = self.HA_range_has_issues(HA_limits=HA_limits)
    #         if has_issues is None:
    #             raise RuntimeError("for the HA limits given here, there should always "
    #                                +"be some simulated HAs")
    #         if has_issues:
    #             continue
    #         #TODO machine precision issue?
    #         if not best or HA_max - HA_min > best["max"] - best["min"]:
    #             best = HA_limits
    #     return best 

    def width_of_successful_HAs_that_are_excluded(self,HA_limits):
        HA_widths = self.HA_widths()
        excluded_width = 0
        for result,HA,width in zip(self.simulation_results,self.HAs,HA_widths):
            HA_is_excluded = ac.less(HA,HA_limits["min"]) or ac.greater(HA, HA_limits["max"])
            if HA_is_excluded and result.success:
                excluded_width = excluded_width + width
        return excluded_width

    def determine_HA_limits_necessary_for_no_missing_calibrators(self):
        best = {"min":Angle(0*u.hour),"max":Angle(0*u.hour)}
        for HA_min,HA_max in itertools.product(self.HAs,self.HAs):
            if ac.geq(HA_min,HA_max):
                continue
            HA_limits = {"min": HA_min, "max": HA_max}
            if self.HA_range_contains_missing_calibrator(HA_limits=HA_limits):
                continue
            if ac.greater(HA_max-HA_min, best["max"]-best["min"]):
                best = HA_limits
        return best 

    def get_recommended_HA_limits(self):
        HA_limits = self.determine_HA_limits_necessary_for_no_missing_calibrators()

        if ac.geq(HA_limits["min"],HA_limits["max"]):
            logging.info("no HA limits found to exclude all missing calibrator instances")
            return

        excluded_success_width = self.width_of_successful_HAs_that_are_excluded(
                                          HA_limits=HA_limits)
        if ac.greater(excluded_success_width, self.max_removable_successful_HA_amount):
            logging.info("suggested HA limits exclude too large successful HA range")
            return

        width = (HA_limits["max"]-HA_limits["min"])
        if ac.less(width, self.min_remaining_HA_amount_when_excluding_HAs):
            logging.info("suggested HA limits are too tight")
            return

        dsa_limits = self.sb.get_DSA_HA_limits()
        if ac.leq(HA_limits["min"], dsa_limits["min"]):
            HA_limits["min"] = Angle(-12*u.hour)
        if ac.geq(HA_limits["max"], dsa_limits["max"]):
            HA_limits["max"] = Angle(12*u.hour)
        return HA_limits
        


    # def analyse_HA_restriction(self):
    #     unnecessarily_restricted_HAs = get_unnecessarily_restricted_HAs(
    #                                       has_hardcoded_cals=self.sb.any_calibrator_hardcoded(),
    #                                       HAs=self.HAs, results=self.simulation_results,
    #                                       OT_allowed_HA=self.sb.OT_allowed_HA)
    #     if not unnecessarily_restricted_HAs:
    #         logging.info("HA is not unnecessarily restricted")
    #     else:
    #         ha_str = ", ".join(str(HA.hour) for HA in unnecessarily_restricted_HAs)
    #         logging.info(f"HA is unnecessarily restricted at following HA(s): {ha_str}")
    #         self.analysis_result.inspection_reasons.add(f"unnecessarily restricted HA(s): {ha_str}")

    # def analyse_simulation_failures(self):
        #logic implemented here:
        # - any error at disallowed HA does not trigger anything
        # - "unobservable" never triggers anything
        # - "server error" triggers inspection if occuring at allowed HA
        # - any other error occuring at allowed HA triggers inspection and Waiting


        # server_error = False
        # failed_and_needs_action = False
        # should_be_Waiting = False
    
        # for HA, result in zip(self.HAs, self.simulation_results):
        #     if result.success:
        #         continue
            
        #     if not HA_is_allowed(HA=HA,OT_allowed_HA=self.sb.OT_allowed_HA):
        #         continue
            
        #     if result.fail_reason.category == "unobservable":
        #         continue
    
        #     if result.fail_reason.category == "server error":
        #         #server error needs inspection, but does not need to be waiting
        #         self.analysis_result.inspection_reasons.add("server error")
        #         continue
    
        #     self.analysis_result.inspection_reasons.add("simulation failure")
        #     self.analysis_result.should_be_Waiting = True
        #     failed_and_needs_action = True
        #     should_be_Waiting = True
    
        # inspection_reasons = []
        # if server_error:
        #     inspection_reasons.append("server error")
        # if failed_and_needs_action:
        #     inspection_reasons.append("simulation failure")
        # return AnalysisResult(inspection_reasons=inspection_reasons,
        #                       should_be_Waiting=should_be_Waiting)

    # def analyse_no_unexpected_error(self):
        # self.analysis_result = AnalysisResult()
        # for ana in (self.analyse_runnable_HA_range(),
        #             self.analyse_HA_restriction(),
        #             self.analyse_simulation_failures()):
        #     self.analysis_result = self.analysis_result.merge(ana)
