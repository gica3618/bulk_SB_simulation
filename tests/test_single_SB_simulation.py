#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:37:57 2026

@author: gianni
"""

import pytest
import datetime
from pathlib import Path
from types import SimpleNamespace
from astropy.coordinates import Angle,SkyCoord
from astropy import units as u
from batch_simulations.application.single_SB_simulation import SingleSBSimulation,\
        AnalysisResult, SingleSBSimulationAnalysis
from batch_simulations.domain.simulation_result import FailReason
from batch_simulations.domain.dsa_policy import DSAPolicy
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
from batch_simulations.domain.sb import SB
from dataclasses import dataclass


# def test_get_unnecessarily_restricted_HAs():
#     class Result:
#         def __init__(self, success):
#             self.success = success
#     with_fail = [Result(False), Result(True)]
#     all_success = [Result(True), Result(True)]
#     HAs = [Angle(-1*u.hour),Angle(1*u.hour)]
#     OT_allowed_HA_wide = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
#     OT_allowed_HA_narrow = {"min":Angle(0*u.hour),"max":Angle(12*u.hour)}
#     #hardcoded, so cannot decide if unnecessary or not:
#     for OT_allowed_HA,results in itertools.product(
#             (OT_allowed_HA_wide,OT_allowed_HA_narrow),(with_fail,all_success)):
#         res_HAs = get_unnecessarily_restricted_HAs(
#                          has_hardcoded_cals=True,HAs=HAs,results=results,
#                          OT_allowed_HA=OT_allowed_HA)
#         assert res_HAs == []
#     #success, but HA not excluded
#     res_HAs = get_unnecessarily_restricted_HAs(
#                      has_hardcoded_cals=False,HAs=HAs,results=all_success,
#                      OT_allowed_HA=OT_allowed_HA_wide)
#     assert res_HAs == []
#     #HA excluded, but not success
#     res_HAs = get_unnecessarily_restricted_HAs(
#                      has_hardcoded_cals=False,HAs=HAs,results=with_fail,
#                      OT_allowed_HA=OT_allowed_HA_narrow)
#     assert res_HAs == []
#     #HA excluded and success
#     res_HAs = get_unnecessarily_restricted_HAs(
#                      has_hardcoded_cals=False,HAs=HAs,results=all_success,
#                      OT_allowed_HA=OT_allowed_HA_narrow)
#     assert res_HAs == [Angle(-1*u.hour)]


class TestSingleSBSimulation:

    general_sb_simulation = SingleSBSimulation(
                               project_code="2025.1.01422.S",sb_name="LaPequen_a_10_TM1",
                               date=datetime.date(2024,1,2),array_config_12m="c43-4")

    def test_build_xml_filepath(self):
        expected_filename = "2025.1.01422.S_LaPequen_a_10_TM1.xml"
        assert self.general_sb_simulation.build_xml_filepath()\
              == Path.cwd()/expected_filename

    def test_get_array_config(self):
        class FakeSB:
            is_7m = None
            requires_TP = None
        sb = FakeSB()
        sb.is_7m = False
        assert self.general_sb_simulation.get_array_config(sb)\
                               == self.general_sb_simulation.array_config_12m
        sb.is_7m = True
        sb.requires_TP = True
        assert self.general_sb_simulation.get_array_config(sb) == "aca.cm10.pm3.cfg"
        sb.requires_TP = False
        assert self.general_sb_simulation.get_array_config(sb) == "7m"

    # def test_should_simulate_again_with_fine_HA_step(self):
    #     class Result:
    #         def __init__(self, success):
    #             self.success = success
    #     with_fail = [Result(False), Result(True)]
    #     all_success = [Result(True), Result(True)]
    #     HAs = [Angle(-1*u.hour),Angle(1*u.hour)]
    #     OT_allowed_HA = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
    #     sim_again = SingleSBSimulation.should_simulate_again_with_fine_HA_step
    #     #case where HA is not unnecessarily restricted:
    #     for hard in (True,False):
    #         assert not sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=hard,
    #                              OT_allowed_HA=OT_allowed_HA)
    #         assert sim_again(results=with_fail,HAs=HAs,has_hardcoded_cals=hard,
    #                          OT_allowed_HA=OT_allowed_HA)
    #     OT_allowed_HA["min"] = Angle(0*u.hour)
    #     #unnecessarily restricted HA, no hardcoding
    #     assert sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=False,
    #                      OT_allowed_HA=OT_allowed_HA)
    #     ##unnecessarily restricted HA, but hardcoding, so cannot say if it really
    #     #is unnecessary
    #     assert not sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=True,
    #                          OT_allowed_HA=OT_allowed_HA)

    def test_any_simulation_failed(self):
        class Result:
            def __init__(self, success):
                self.success = success
        results = [Result(True),Result(True)]
        assert not SingleSBSimulation.any_simulation_failed(results)
        results.append(Result(False))
        assert SingleSBSimulation.any_simulation_failed(results)

    def test_should_run_again_with_fine_step(self):
        class FakeSingleSBSimulation:
            fine_HA_step = SingleSBSimulation.fine_HA_step
            should_run_again_with_fine_step = SingleSBSimulation.should_run_again_with_fine_step
            @staticmethod
            def any_simulation_failed(results):
                return SingleSBSimulation.any_simulation_failed(results=results)
        @dataclass
        class Result:
            success: bool
        fake = FakeSingleSBSimulation()
        #test that if step is already fine, should not run again, regardless of
        #results:
        for results in [[Result(True),Result(True)],
                        [Result(False),Result(True)]]:
            assert not fake.should_run_again_with_fine_step(step=fake.fine_HA_step,
                                                            results=results)
        #test that it should run again if something failed
        assert fake.should_run_again_with_fine_step(step=fake.fine_HA_step+Angle(0.1*u.hour),
                                                    results=[Result(False),Result(True)])

    def test_simulate_HAs(self,monkeypatch):
        #inspired by ChatGPT
        class FakeSingleSBSimulation:
            xml_filepath = "/tmp/file.xml"
            date = datetime.date(1978,3,3)
            default_HA_step = SingleSBSimulation.default_HA_step
            fine_HA_step = SingleSBSimulation.fine_HA_step
            should_run_again_with_fine_step = SingleSBSimulation.should_run_again_with_fine_step
            simulate_HAs = SingleSBSimulation.simulate_HAs
            @staticmethod
            def any_simulation_failed(results):
                return SingleSBSimulation.any_simulation_failed(results=results)

        @dataclass
        class Result:
            success: bool

        class FakePlanner:
            def __init__(self,DSA_HA_limits):
                pass
            def HA_jobs(self, xml_filepath, array_config, date, step):
                HAs = [Angle(1*u.hour),Angle(1*u.hour)+step]
                jobs = [SimpleNamespace(HA=HA) for HA in HAs]
                return HAs,jobs

        class FakeRunnerFail:
            def run(self,job):
                return Result(False)
    
        class FakeRunnerSuccess:
            def run(self, job):
                return Result(True)
    
        class FakeSB:
            
            rep_coord = SkyCoord('05h47m17.0876901s', '-51d03m59.441135s', frame='icrs')
            HA_is_restricted = SB.HA_is_restricted

            def __init__(self,any_hardcoded,OT_allowed_HA):
                self.any_hardcoded = any_hardcoded
                self.OT_allowed_HA = OT_allowed_HA
            def any_calibrator_hardcoded(self):
                return self.any_hardcoded
            def get_DSA_HA_limits(self):
                return {"min":Angle(-4*u.hour),"max":Angle(3*u.hour)}

        #test without skipped simulation:
        single_sim = FakeSingleSBSimulation()
        kwargs = {"array_config":"c43-5","planner_cls":FakePlanner}
        std_sb = FakeSB(any_hardcoded=False, OT_allowed_HA={"min":Angle(-12*u.hour),
                                                            "max":Angle(12*u.hour)})
        default_step = FakeSingleSBSimulation.default_HA_step
        fine_step = FakeSingleSBSimulation.fine_HA_step

        def check(sb,runner_cls,expected_step,all_HA_should_succeed=False,
                  all_HA_should_fail=False):
            HAs,results = single_sim.simulate_HAs(**kwargs,sb=sb,runner_cls=runner_cls)
            assert HAs[1]-HAs[0] == expected_step
            if all_HA_should_succeed:
                assert all([r.success for r in results])
            if all_HA_should_fail:
                assert not any([r.success for r in results])

        #test success and failure:
        check(sb=std_sb,runner_cls=FakeRunnerSuccess,
              expected_step=default_step,all_HA_should_succeed=True)
        check(sb=std_sb,runner_cls=FakeRunnerFail,
              expected_step=fine_step,all_HA_should_fail=True)

        #test that fine step is used for SB with restricted HA
        restricted_sb = FakeSB(any_hardcoded=False, OT_allowed_HA={"min":Angle(-12*u.hour),
                                                                   "max":Angle(0*u.hour)})
        check(sb=restricted_sb,runner_cls=FakeRunnerSuccess,
              expected_step=fine_step,all_HA_should_succeed=True)

        #test skipping simulations
        # def fake_reasons(self,sb,start_HA):
        #     return ["rassel"]
        # monkeypatch.setattr(DSAPolicy, "reasons_for_not_showing_up",fake_reasons)
        # HAs,results = single_sim.simulate_HAs(**kwargs,sb=std_sb,
        #                                       runner_cls=FakeRunnerSuccess)
        # assert HAs[1]-HAs[0] == SingleSBSimulation.fine_HA_step
        # assert not any([r.success for r in results])
        # for r in results:
        #     assert r.fail_reason.error_message == "simulation skipped, SB would not show up in DSA: rassel"

    def test_simulate(self):
        class FakeJob:
            def __init__(self,HA):
                self.HA = HA
        class FakeSingleSBSimulation:
            default_HA_step = SingleSBSimulation.default_HA_step
            fine_HA_step = SingleSBSimulation.fine_HA_step
            date = datetime.date(2026,1,3)
            HAs = [Angle(1*u.hour),Angle(1.5*u.hour)]
            sb_name = "test_SB"
            project_code = "2025.1.00001.S"
            simulate = SingleSBSimulation.simulate
            def build_sb(self):
                return
            def get_array_config(self,sb):
                return "c43-3"
            def simulate_HAs(self,array_config,sb):
                results = [f"result_{HA}" for HA in self.HAs]
                return self.HAs,results
        fake = FakeSingleSBSimulation()
        sim_analysis = fake.simulate()
        assert sim_analysis.HAs == FakeSingleSBSimulation.HAs
        assert sim_analysis.simulation_results\
                   == [f"result_{HA}" for HA in FakeSingleSBSimulation.HAs]
        assert sim_analysis.unexpected_error is None
        #now the case where an unexpected error is thrown:
        error_message = "something went wrong"
        def build_sb_malfunction():
            raise RuntimeError(error_message)
        fake.build_sb = build_sb_malfunction
        sim_analysis = fake.simulate()
        assert sim_analysis.HAs is None
        assert sim_analysis.simulation_results is None
        for key in ("error_message","full_traceback"):
            assert "RuntimeError" in sim_analysis.unexpected_error[key]
            assert error_message in sim_analysis.unexpected_error[key]


class TestAnalysisResult:

    def test_suggest_state(self):
        analysis_result = AnalysisResult()
        for state in ("Ready","Waiting","WaitingForP2G"):
            analysis_result.suggest_state(state)
            assert analysis_result.recommended_state == state
        for state in ("Ready","Waiting"):
            analysis_result.suggest_state(state)
            assert analysis_result.recommended_state == "WaitingForP2G"
        analysis_result.recommended_state = "Waiting"
        analysis_result.suggest_state("Ready")
        assert analysis_result.recommended_state == "Waiting"
        with pytest.raises(ValueError):
            analysis_result.suggest_state("WaitingForJAO")
        

class TestSingleSBSimulationAnalysis:

    def generate_fake_analysis(self,result_codes,HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)],
                               hardcoded=False,hardcoded_except_PolCal=False,
                               OT_allowed_HA={"min":Angle(-4*u.hour),"max":Angle(3*u.hour)},
                               unexpected_error=None,
                               DSA_HA_limits={"min":Angle(-4*u.hour),"max":Angle(3*u.hour)}):
        if unexpected_error is not None:
            return SingleSBSimulationAnalysis(sb=None,HAs=None,simulation_results=None,
                                              unexpected_error=unexpected_error)
        if len(result_codes) != len(HAs):
            raise ValueError

        class FakeSB:
            def __init__(self,OT_allowed_HA,hardcoded,hardcoded_except_PolCal,
                         DSA_HA_limits):
                self.OT_allowed_HA = OT_allowed_HA
                self.hardcoded = hardcoded
                self.hardcoded_except_PolCal = hardcoded_except_PolCal
                self.DSA_HA_limits = DSA_HA_limits                
            def any_calibrator_hardcoded(self):
                return self.hardcoded
            def any_calibrator_except_PolCal_hardcoded(self):
                return self.hardcoded_except_PolCal
            def get_DSA_HA_limits(self):
                return self.DSA_HA_limits

        sb = FakeSB(OT_allowed_HA=OT_allowed_HA,hardcoded=hardcoded,
                    hardcoded_except_PolCal=hardcoded_except_PolCal,
                    DSA_HA_limits=DSA_HA_limits)

        simulation_results = []
        for code in result_codes:
            if code == "success":
                simulation_results.append(SimpleNamespace(success=True,fail_reason=None))
            else:
                fail_reason = SimpleNamespace(category=code)
                simulation_results.append(SimpleNamespace(success=False,fail_reason=fail_reason))

        return SingleSBSimulationAnalysis(sb=sb,HAs=HAs,simulation_results=simulation_results,
                                          unexpected_error=unexpected_error)

    def test_analyse_with_unexpected_error(self):
        unexpected_error = {"error_message":"kabumm",
                            "full_traceback":"the full traceback"}
        fake_analysis = self.generate_fake_analysis(HAs=None, result_codes=None,
                                                    unexpected_error=unexpected_error)
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.inspection_reasons == set(["unexpected error"])
        assert fake_analysis.analysis_result.recommended_state == "WaitingForP2G"
        assert fake_analysis.analysis_result.recommended_HA_limits is None

    def test_analyse_regular_successful(self):
        for result_codes in [["success",]*3,["success","unobservable","success"]]:
            fake_analysis = self.generate_fake_analysis(result_codes=result_codes)
            fake_analysis.analyse()
            assert fake_analysis.analysis_result.recommended_state == "Ready"
            assert fake_analysis.analysis_result.recommended_HA_limits\
                           == {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
            assert fake_analysis.analysis_result.inspection_reasons == set()
            #case where a calibrator except PolCal is hardcoded:
            fake_analysis = self.generate_fake_analysis(
                                        result_codes=result_codes,
                                        hardcoded_except_PolCal=True)
            fake_analysis.analyse()
            assert fake_analysis.analysis_result.recommended_state is None
            assert fake_analysis.analysis_result.recommended_HA_limits is None
            assert fake_analysis.analysis_result.inspection_reasons == set()

    def test_analyse_regular_WaitingForP2G(self):
        for other_result_codes in [["success"]*2,["other"]*2,
                                   ["success","unobservable"]]:
            for result_code in ("exceeds 2h limit","server error","other"):
                fake_analysis = self.generate_fake_analysis(
                                     result_codes=other_result_codes+[result_code])
                fake_analysis.analyse()
                assert fake_analysis.analysis_result.recommended_state == "WaitingForP2G"
                assert result_code in fake_analysis.analysis_result.inspection_reasons
                assert fake_analysis.analysis_result.recommended_HA_limits is None

    def test_analyse_regular_missing_calibrator_hardcoded_nonallowed_HA(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["success","success","missing calibrator"],
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)],
                             hardcoded=True,
                             OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(0.5*u.hour)})
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.recommended_state is None
        assert fake_analysis.analysis_result.recommended_HA_limits is None
        assert fake_analysis.analysis_result.inspection_reasons == set()

    def test_analyse_regular_missing_calibrator_hardcoded_allowed_HA(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["success","success","missing calibrator"],
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)],
                             hardcoded=True,
                             OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(1*u.hour)})
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.recommended_state == "Waiting"
        assert fake_analysis.analysis_result.recommended_HA_limits is None
        assert fake_analysis.analysis_result.inspection_reasons == set(["missing calibrator"])

    def test_analyse_regular_missing_calibrator_no_hardcoding_fixable(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["success","success","missing calibrator"],
                             HAs=[Angle(-2*u.hour),Angle(0*u.hour),Angle(2*u.hour)])
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.recommended_state == "Ready"
        assert fake_analysis.analysis_result.recommended_HA_limits \
            == {"min":Angle(-2*u.hour),"max":Angle(0*u.hour)}
        assert fake_analysis.analysis_result.inspection_reasons == set()

    def test_analyse_regular_missing_calibrator_no_hardcoding_not_fixable_allowed_HA(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["missing calibrator"]*3)
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.recommended_state == "Waiting"
        assert fake_analysis.analysis_result.recommended_HA_limits is None
        assert "missing calibrator" in fake_analysis.analysis_result.inspection_reasons

    def test_analyse_regular_missing_calibrator_no_hardcoding_not_fixable_disallowed_HA(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["missing calibrator"]*3,
                             OT_allowed_HA={"min":Angle(-5*u.hour),"max":Angle(-4*u.hour)})
        fake_analysis.analyse()
        assert fake_analysis.analysis_result.recommended_state is None
        assert fake_analysis.analysis_result.recommended_HA_limits is None
        assert "missing calibrator" in fake_analysis.analysis_result.inspection_reasons

    def test_HA_is_allowed(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["success"]*3,
                             OT_allowed_HA={"min":Angle(-2*u.hour),"max":Angle(2*u.hour)})
        for HA in (Angle(0*u.hour),Angle(2*u.hour)):
            assert fake_analysis.HA_is_allowed(HA)
        for HA in (Angle(-3*u.hour),Angle(2.1*u.hour)):
            assert not fake_analysis.HA_is_allowed(HA)

    def test_get_missing_calibrator_HAs(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["missing calibrator","success","missing calibrator"],
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)])
        assert fake_analysis.get_missing_calibrator_HAs()\
                                       == [Angle(-1*u.hour),Angle(1*u.hour)]
        fake_analysis = self.generate_fake_analysis(result_codes=["success",]*3)
        assert fake_analysis.get_missing_calibrator_HAs() == []

    def test_missing_calibrator_at_allowed_HA(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["missing calibrator","success","success"],
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)],
                             OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(1*u.hour)})
        assert fake_analysis.missing_calibrator_at_allowed_HA()
        fake_analysis.sb.OT_allowed_HA = {"min":Angle(0*u.hour),"max":Angle(1*u.hour)}
        assert not fake_analysis.missing_calibrator_at_allowed_HA()

    def test_HA_widths(self):
        HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
        result_codes = ["success"]*len(HAs)
        fake_analysis = self.generate_fake_analysis(result_codes=result_codes,HAs=HAs)
        assert fake_analysis.HA_widths() == [Angle(0.5*u.hour),
                                             Angle((0.5+0.75)*u.hour),
                                             Angle((0.75+1.75)*u.hour),
                                             Angle(1.75*u.hour)]
        fake_analysis = self.generate_fake_analysis(result_codes=["succes"],
                                                    HAs=[Angle(1*u.hour)])
        assert fake_analysis.HA_widths() == [Angle(0*u.hour)]

    def test_runnable_HA_amount(self):
        fake_analysis = self.generate_fake_analysis(result_codes=["missing calibrator"]*3)
        assert fake_analysis.runnable_HA_amount() == Angle(0*u.hour)
        fake_analysis = self.generate_fake_analysis(
                            result_codes=["missing calibrator","success","success"],
                            HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(0.5*u.hour)])
        assert fake_analysis.runnable_HA_amount() == Angle((0.5+0.25)*u.hour)+Angle(0.25*u.hour)

    def test_extreme_HA_restriction_imposed_by_P2G(self):
        fake_analysis = self.generate_fake_analysis(
                             result_codes=["success"]*3,
                             OT_allowed_HA={"min":Angle(1*u.hour),"max":Angle(1.5*u.hour)})
        assert fake_analysis.extreme_HA_restriction_imposed_by_P2G()
        fake_analysis.sb.OT_allowed_HA = {"min":Angle(1*u.hour),"max":Angle(3*u.hour)}
        assert not fake_analysis.extreme_HA_restriction_imposed_by_P2G()

    def test_get_sim_results(self):
        result_codes = ["success","other","missing calibrator"]
        fake_analysis = self.generate_fake_analysis(
                             result_codes=result_codes,
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)])
        sim_results = fake_analysis.get_sim_results(
                          HA_limits={"min":Angle(-3*u.hour),"max":Angle(-2*u.hour)})
        assert sim_results == []
        sim_results = fake_analysis.get_sim_results(
                          HA_limits={"min":Angle(-1*u.hour),"max":Angle(3*u.hour)})
        assert len(sim_results) == 3
        assert sim_results[0].success
        assert sim_results[1].fail_reason.category == "other"
        assert sim_results[2].fail_reason.category == "missing calibrator"
        sim_results = fake_analysis.get_sim_results(
                          HA_limits={"min":Angle(0.5*u.hour),"max":Angle(1.5*u.hour)})
        assert len(sim_results) == 1
        assert sim_results[0].fail_reason.category == "missing calibrator"

    def test_HA_range_contains_missing_calibrator(self):
        result_codes = ["success","other","missing calibrator"]
        fake_analysis = self.generate_fake_analysis(
                             result_codes=result_codes,
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)])
        assert fake_analysis.HA_range_contains_missing_calibrator(
                HA_limits={"min":Angle(0.5*u.hour),"max":Angle(1*u.hour)})
        assert not fake_analysis.HA_range_contains_missing_calibrator(
                HA_limits={"min":Angle(-1*u.hour),"max":Angle(0*u.hour)})

    def test_width_of_successful_HAs_that_are_excluded(self):
        result_codes = ["success","success","missing calibrator"]
        fake_analysis = self.generate_fake_analysis(
                             result_codes=result_codes,
                             HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour)])
        HA_limits = {"min":Angle(-1*u.hour),"max":Angle(2*u.hour)}
        assert fake_analysis.width_of_successful_HAs_that_are_excluded(HA_limits=HA_limits)\
                == Angle(0*u.hour)
        HA_limits = {"min":Angle(0*u.hour),"max":Angle(2*u.hour)}
        assert fake_analysis.width_of_successful_HAs_that_are_excluded(HA_limits=HA_limits)\
                == Angle(0.5*u.hour)
        HA_limits = {"min":Angle(1*u.hour),"max":Angle(2*u.hour)}
        assert fake_analysis.width_of_successful_HAs_that_are_excluded(HA_limits=HA_limits)\
                == Angle(0.5*u.hour+1*u.hour)

    def test_determine_HA_limits_necessary_for_no_missing_calibrators(self):
        HAs=[Angle(-1*u.hour),Angle(0*u.hour),Angle(1*u.hour),Angle(2*u.hour)]
        def check(result_codes,expected_limits):
            fake_analysis = self.generate_fake_analysis(result_codes=result_codes,HAs=HAs)
            assert fake_analysis.determine_HA_limits_necessary_for_no_missing_calibrators()\
                    == expected_limits
        result_codes = ["success"]*4
        expected_limits = {"min":Angle(-1*u.hour),"max":Angle(2*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["success","success","missing calibrator","missing calibrator"]
        expected_limits = {"min":Angle(-1*u.hour),"max":Angle(0*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["missing calibrator","success","missing calibrator","missing calibrator"]
        expected_limits = {"min":Angle(0*u.hour),"max":Angle(0*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["missing calibrator"]*4
        expected_limits = {"min":Angle(0*u.hour),"max":Angle(0*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["success","success","missing calibrator","success"]
        expected_limits = {"min":Angle(-1*u.hour),"max":Angle(0*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["missing calibrator","success","missing calibrator","success"]
        expected_limits = {"min":Angle(0*u.hour),"max":Angle(0*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)

    def test_get_recommended_HA_limits(self):
        HAs=[Angle(i*u.hour) for i in range(-3,4)]
        DSA_HA_limits = {"min":min(HAs),"max":max(HAs)}
        def check(result_codes,expected_limits):
            fake_analysis = self.generate_fake_analysis(
                               result_codes=result_codes,HAs=HAs,
                               DSA_HA_limits=DSA_HA_limits)
            assert fake_analysis.get_recommended_HA_limits() == expected_limits
        result_codes = ["success"]*len(HAs)
        expected_limits = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)
        result_codes = ["missing calibrator"]*len(HAs)
        expected_limits = None
        check(result_codes=result_codes,expected_limits=expected_limits)
        #trigger case where too much success is removed:
        result_codes = ["succes","success","success","missing calibrator",
                        "success","success","success"]
        expected_limits = None
        check(result_codes=result_codes,expected_limits=expected_limits)
        #trigger too tight HA limits
        result_codes = ["missing calibrator","missing calibrator","success","success",
                        "missing calibrator","missing calibrator","missing calibrator"]
        expected_limits = None
        check(result_codes=result_codes,expected_limits=expected_limits)
        #"normal" case
        result_codes = ["missing calibrator","successs","success","success",
                        "success","missing calibrator","missing calibrator"]
        expected_limits = {"min":Angle(-2*u.hour),"max":Angle(1*u.hour)}
        check(result_codes=result_codes,expected_limits=expected_limits)


"""
    def test_HA_is_allowed(self):
        OT_allowed_HA = {"min":Angle(-2.1*u.hour),"max":Angle(3*u.hour)}
        class FakeSummary:
            HA_is_allowed = SingleSBSimulationAnalysis.HA_is_allowed
            sb = SimpleNamespace(OT_allowed_HA=OT_allowed_HA)
        fake = FakeSummary()
        for HA in (OT_allowed_HA["min"], OT_allowed_HA["max"],
                   Angle(-2*u.hour),Angle(0*u.hour),Angle(1.1*u.hour)):
            assert fake.HA_is_allowed(HA=HA)
        for HA in (Angle(-2.2*u.hour),Angle(3.1*u.hour)):
            assert not fake.HA_is_allowed(HA=HA)

    def test_HA_widths(self):
        class FakeSummary:
            HA_widths = SingleSBSimulationSummary.HA_widths
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
        fake_summary = FakeSummary()
        widths = fake_summary.HA_widths()
        assert widths == [0.5,0.5+0.75,0.75+1.75,1.75]

    def test_runnable_HA_amount(self):
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
            HA_widths = SingleSBSimulationSummary.HA_widths
            simulation_results = [SimpleNamespace(success=True),SimpleNamespace(success=False),
                                  SimpleNamespace(success=True),SimpleNamespace(success=True)]
            runnable_HA_amount = SingleSBSimulationSummary.runnable_HA_amount
        fake_summary = FakeSummary()
        widths = fake_summary.HA_widths()
        usable_HA_range = fake_summary.runnable_HA_amount()
        expected_range = widths[0]+widths[2]+widths[3]
        assert usable_HA_range == expected_range

    def test_analyse_runnable_HA_range(self):
        class FakeSummary:
            analyse_runnable_HA_range\
                       = SingleSBSimulationSummary.analyse_runnable_HA_range
            def __init__(self):
                self.analysis_result = AnalysisResult()
            def runnable_HA_amount(self):
                return 2
        fake_summary = FakeSummary()
        fake_summary.min_runnable_HA_amount = 1
        fake_summary.analyse_runnable_HA_range()
        assert fake_summary.analysis_result.inspection_reasons == set()
        assert fake_summary.analysis_result.recommended_state is None
        
        fake_summary = FakeSummary()
        fake_summary.min_runnable_HA_amount = 3.3333
        fake_summary.analyse_runnable_HA_range()
        assert fake_summary.analysis_result.inspection_reasons\
                                       == set(["runnable HA range is small"])
        assert fake_summary.analysis_result.recommended_state is None

    def test_analyse_HA_restriction(self):
        class FakeSB:
            OT_allowed_HA = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
            def __init__(self,has_hardcoded):
                self.has_hardcoded = has_hardcoded
            def any_calibrator_hardcoded(self):
                return self.has_hardcoded
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(-1*u.hour),Angle(0.5*u.hour),Angle(4*u.hour)]
            simulation_results = [SimpleNamespace(success=True),SimpleNamespace(success=False),
                                  SimpleNamespace(success=True),SimpleNamespace(success=True)]
            analyse_HA_restriction = SingleSBSimulationSummary.analyse_HA_restriction
            def __init__(self,has_hardcoded):
                self.sb = FakeSB(has_hardcoded=has_hardcoded)
                self.analysis_result = AnalysisResult()

        fake_summary = FakeSummary(has_hardcoded=True)
        fake_summary.analyse_HA_restriction()
        assert fake_summary.analysis_result.inspection_reasons == set()
        assert not fake_summary.analysis_result.should_be_Waiting

        fake_summary = FakeSummary(has_hardcoded=False)
        fake_summary.analyse_HA_restriction()
        assert fake_summary.analysis_result.inspection_reasons == set()
        assert not fake_summary.analysis_result.should_be_Waiting

        fake_summary = FakeSummary(has_hardcoded=False)
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(12*u.hour)}
        fake_summary.analyse_HA_restriction()
        assert fake_summary.analysis_result.inspection_reasons\
                                    == set(["unnecessarily restricted HA(s): -2.0"])
        assert not fake_summary.analysis_result.should_be_Waiting

        fake_summary = FakeSummary(has_hardcoded=False)
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(2*u.hour)}
        fake_summary.analyse_HA_restriction()
        assert fake_summary.analysis_result.inspection_reasons\
                           == set(["unnecessarily restricted HA(s): -2.0, 4.0"])
        assert not fake_summary.analysis_result.should_be_Waiting

    def test_analyse_simulation_failures(self):
        def sim_result(success,category=None):
            fail_reason = None if success else\
                   FailReason(error_message="a",error_summary="b",category=category)
            return SimpleNamespace(success=success,fail_reason=fail_reason)
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(0.5*u.hour)]
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-4*u.hour),"max":Angle(4*u.hour)})
            analyse_simulation_failures = SingleSBSimulationSummary.analyse_simulation_failures
            def __init__(self):
                self.analysis_result = AnalysisResult()

        def assert_correct_analysis(simulation_results,expected_reasons,expected_waiting,
                                    OT_allowed_HA=None):
            fake_summary = FakeSummary()
            if OT_allowed_HA is not None:
                fake_summary.sb.OT_allowed_HA = OT_allowed_HA
            fake_summary.simulation_results = simulation_results
            fake_summary.analyse_simulation_failures()
            assert fake_summary.analysis_result.inspection_reasons\
                                      == expected_reasons
            assert fake_summary.analysis_result.should_be_Waiting == expected_waiting

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=True)],
                                expected_reasons=set(),
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=set(["server error"]),
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=set(["server error"]),
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=True)],
                                expected_reasons=set(),
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=set(["server error", "simulation failure"]),
                                expected_waiting=True)
        
        #making the missing calibrator failure to occur outside of allowed HA:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=set(["server error"]),
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="exceeds 2h limit")],
                                expected_reasons=set(["simulation failure"]),
                                expected_waiting=True)

        #server error outisde of allowed HA, should not trigger inspection:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="server error"),
                                                    sim_result(success=True)],
                                expected_reasons=set(),
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

    def test_analyse_no_unexpected_error(self):
        class FakeSimResult:
            def __init__(self,success,category=None):
                self.success = success
                self.fail_reason = None if category is None\
                                    else FailReason(error_message="bla",
                                                    error_summary="bla",
                                                    category=category)
        #this SB has min allowed HA -12H and max 12H
        sb = BuildSBFromXML.build("tests/xmls/2025.1.01279.S_general_SB.xml")
        summary_kwargs = {"sb":sb,"config":"c43-9","date":datetime.date(1974,2,22),
                          "unexpected_error":None}

        #all success:
        HAs = [Angle(-12*u.hour),Angle(-2*u.hour),Angle(-0*u.hour),Angle(2*u.hour),
               Angle(4*u.hour),]
        simulation_results = [FakeSimResult(True)]*len(HAs)
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse()
        assert summary.analysis_result.inspection_reasons == set()
        assert not summary.analysis_result.should_be_Waiting
        
        #all failures:
        simulation_results = [FakeSimResult(False,category="missing calibrator"),]*5
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse()
        assert summary.analysis_result.inspection_reasons\
                  == set(["simulation failure","runnable HA range is small"])
        assert summary.analysis_result.should_be_Waiting

        #trigger unnecessary restricted HA
        HAs[0] = Angle(-14*u.hour)
        simulation_results = [FakeSimResult(True)]*len(HAs)
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse()
        assert summary.analysis_result.inspection_reasons\
                     == set(['unnecessarily restricted HA(s): -14.0'])
        assert not summary.analysis_result.should_be_Waiting
"""