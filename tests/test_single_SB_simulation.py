#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:37:57 2026

@author: gianni
"""

from unittest.mock import patch
import datetime
from pathlib import Path
from types import SimpleNamespace
from astropy.coordinates import Angle,SkyCoord
from astropy import units as u
from batch_simulations.application.single_SB_simulation import SingleSBSimulation,\
    AnalysisResult, SingleSBSimulationSummary,get_unnecessarily_restricted_HAs
from batch_simulations.domain.simulation_result import FailReason
from batch_simulations.domain.calibrator import Calibrator
from batch_simulations.infrastructure.OT_xml import BuildSBFromXML
import itertools
from scipy import constants


def test_get_unnecessarily_restricted_HAs():
    class Result:
        def __init__(self, success):
            self.success = success
    with_fail = [Result(False), Result(True)]
    all_success = [Result(True), Result(True)]
    HAs = [Angle(-1*u.hour),Angle(1*u.hour)]
    OT_allowed_HA_wide = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
    OT_allowed_HA_narrow = {"min":Angle(0*u.hour),"max":Angle(12*u.hour)}
    #hardcoded, so cannot decide if unnecessary or not:
    for OT_allowed_HA,results in itertools.product(
            (OT_allowed_HA_wide,OT_allowed_HA_narrow),(with_fail,all_success)):
        res_HAs = get_unnecessarily_restricted_HAs(
                         has_hardcoded_cals=True,HAs=HAs,results=results,
                         OT_allowed_HA=OT_allowed_HA)
        assert res_HAs == []
    #success, but HA not excluded
    res_HAs = get_unnecessarily_restricted_HAs(
                     has_hardcoded_cals=False,HAs=HAs,results=all_success,
                     OT_allowed_HA=OT_allowed_HA_wide)
    assert res_HAs == []
    #HA excluded, but not success
    res_HAs = get_unnecessarily_restricted_HAs(
                     has_hardcoded_cals=False,HAs=HAs,results=with_fail,
                     OT_allowed_HA=OT_allowed_HA_narrow)
    assert res_HAs == []
    #HA excluded and success
    res_HAs = get_unnecessarily_restricted_HAs(
                     has_hardcoded_cals=False,HAs=HAs,results=all_success,
                     OT_allowed_HA=OT_allowed_HA_narrow)
    assert res_HAs == [Angle(-1*u.hour)]


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

    def test_should_simulate_again_with_fine_HA_step(self):
        class Result:
            def __init__(self, success):
                self.success = success
        with_fail = [Result(False), Result(True)]
        all_success = [Result(True), Result(True)]
        HAs = [Angle(-1*u.hour),Angle(1*u.hour)]
        OT_allowed_HA = {"min":Angle(-12*u.hour),"max":Angle(12*u.hour)}
        sim_again = SingleSBSimulation.should_simulate_again_with_fine_HA_step
        #case where HA is not unnecessarily restricted:
        for hard in (True,False):
            assert not sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=hard,
                                 OT_allowed_HA=OT_allowed_HA)
            assert sim_again(results=with_fail,HAs=HAs,has_hardcoded_cals=hard,
                             OT_allowed_HA=OT_allowed_HA)
        OT_allowed_HA["min"] = Angle(0*u.hour)
        #unnecessarily restricted HA, no hardcoding
        assert sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=False,
                         OT_allowed_HA=OT_allowed_HA)
        ##unnecessarily restricted HA, but hardcoding, so cannot say if it really
        #is unnecessary
        assert not sim_again(results=all_success,HAs=HAs,has_hardcoded_cals=True,
                             OT_allowed_HA=OT_allowed_HA)

    def test_get_targets_for_elevation_check(self):
        class FakeSB:
            science_targets = [{"name":"miraculix",
                                "coordinates":SkyCoord('05h47m17.08s', '-51d03m59.44s')},
                               {"name":"weissnix",
                                "coordinates":SkyCoord('10h47m17.08s', '-71d03m59.44s')}]
            calibrators = [Calibrator(name="Phase",source_name="J1326-5256",
                                      cal_type="Phase",is_hardcoded=True,
                                      coordinates=SkyCoord('05h00m00s', '-53d00m00s')),
                           Calibrator(name="Bandpass calibrator",source_name="query",
                                      cal_type="Bandpass",is_hardcoded=False,
                                      coordinates=None),
                           Calibrator(name="DGC",source_name="asterix",
                                      cal_type="DGC",is_hardcoded=True,
                                      coordinates=SkyCoord('08h00m00s', '-20d00m00s'))]
        sb = FakeSB()
        targets = SingleSBSimulation.get_targets_for_elevation_check(sb)
        assert len(targets) == 4
        for i in (0,1):
            assert targets[i] == {"name":FakeSB.science_targets[i]["name"],
                                  "DEC":FakeSB.science_targets[i]["coordinates"].dec,"DGC":False}
        assert targets[2] == {"name":"Phase","DEC":FakeSB.calibrators[0].coordinates.dec,"DGC":False}
        assert targets[3] == {"name":"DGC","DEC":FakeSB.calibrators[2].coordinates.dec,"DGC":True}        

    def test_get_targets_beyond_elevation_limits(self,monkeypatch):
        targets = [{"name":"dgc miraculix","DEC":SkyCoord('00h00m00s', '-10d00m00s').dec,
                    "DGC":True},#this should not be limited by elevation
                   {"name":"J12345678","DEC":SkyCoord('00h00m00s', '80d00m00s').dec,
                    "DGC":False},#this one is unobservable by ALMA
                   {"name":"majestix","DEC":SkyCoord("00h00m0.0s","-23d46m56.772s").dec,
                    "DGC":True},#above 85 deg for HA between -0.36h and 0.36h
                   {"name":"majestix","DEC":SkyCoord("00h00m0.0s","-23d46m56.772s").dec,
                    "DGC":False}#above 85 deg for HA between -0.36h and 0.36h
                   ]
        output = SingleSBSimulation.get_targets_beyond_elevation_limits(
                          start_HA=Angle(-0.8*u.hour),execution_time=0.5*constants.hour,
                          targets=targets)
        assert len(output) == 2
        assert output[0] == "J12345678"
        assert output[1] == "majestix"

    def test_simulate_HAs(self):
        #inspired by ChatGPT
        class FakeSingleSBSimulation:
            xml_filepath = "/tmp/file.xml"
            date = datetime.date(1978,3,3)
            simulate_HAs = SingleSBSimulation.simulate_HAs
            @staticmethod
            def should_simulate_again_with_fine_HA_step(results,HAs,has_hardcoded_cals,
                                                 OT_allowed_HA):
                return SingleSBSimulation.should_simulate_again_with_fine_HA_step(
                       results=results,HAs=HAs,has_hardcoded_cals=has_hardcoded_cals,
                       OT_allowed_HA=OT_allowed_HA)
            @staticmethod
            def get_targets_beyond_elevation_limits(start_HA,execution_time,targets):
                return SingleSBSimulation.get_targets_beyond_elevation_limits(
                          start_HA=start_HA,execution_time=execution_time,
                          targets=targets)

        class FakeSingleSBSimulationNoElevationCheck(FakeSingleSBSimulation):
            @staticmethod
            def get_targets_for_elevation_check(sb):
                return []

        class FakeSingleSBSimulationWithElevationCheck(FakeSingleSBSimulation):
            @staticmethod
            def get_targets_for_elevation_check(sb):
                return [{"name":"asterix und underlix",
                         "DEC":SkyCoord("00h00m0.0s","80d46m56.772s").dec,#not observable by ALMA
                         "DGC":False}]

        class Result:
            def __init__(self, success):
                self.success = success

        class FakePlanner:
            default_HA_step = SingleSBSimulation.default_HA_step
            fine_HA_step = SingleSBSimulation.fine_HA_step
            def HA_jobs(self, xml_filepath, array_config, date, step):
                HAs = [Angle(1*u.hour),Angle(1*u.hour)+step]
                jobs = [SimpleNamespace(HA=HA) for HA in HAs]
                return HAs,jobs
            def HA_jobs_default_HA_step(self,**kwargs):
                return self.HA_jobs(**kwargs,step=self.default_HA_step)
            def HA_jobs_fine_HA_step(self,**kwargs):
                return self.HA_jobs(**kwargs,step=self.fine_HA_step)

        class FakeRunnerFail:
            def run(self,job):
                return Result(False)
    
        class FakeRunnerSuccess:
            def run(self, job):
                return Result(True)
    
        class FakeSB:
            
            rep_coord = SkyCoord('05h47m17.0876901s', '-51d03m59.441135s', frame='icrs')
            
            def __init__(self,any_hardcoded,OT_allowed_HA):
                self.any_hardcoded = any_hardcoded
                self.OT_allowed_HA = OT_allowed_HA
            def any_calibrator_hardcoded(self):
                return self.any_hardcoded
            def single_execution_time(self):
                return 0.5*constants.hour

        planner = FakePlanner()

        #test without elevation check:
        single_sim = FakeSingleSBSimulationNoElevationCheck()
        #test success and failure:
        std_sb = FakeSB(any_hardcoded=False, OT_allowed_HA={"min":Angle(-12*u.hour),
                                                        "max":Angle(12*u.hour)})
        kwargs = {"planner":planner,"array_config":"c43-5","sb":std_sb}
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerSuccess)
        assert HAs[1]-HAs[0] == SingleSBSimulation.default_HA_step
        assert all([r.success for r in results])
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerFail)
        assert HAs[1]-HAs[0] == SingleSBSimulation.fine_HA_step
        assert not any([r.success for r in results])
        #test unnecessarily restricted HA (fine HA step):
        restricted_sb = FakeSB(any_hardcoded=False, OT_allowed_HA={"min":Angle(-12*u.hour),
                                                        "max":Angle(0*u.hour)})
        kwargs["sb"] = restricted_sb
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerSuccess)
        assert HAs[1]-HAs[0] == SingleSBSimulation.fine_HA_step
        assert all([r.success for r in results])

        #test elevation check
        single_sim = FakeSingleSBSimulationWithElevationCheck()
        kwargs["sb"] = std_sb
        HAs,results = single_sim.simulate_HAs(**kwargs,runner_cls=FakeRunnerSuccess)
        assert HAs[1]-HAs[0] == SingleSBSimulation.fine_HA_step
        assert not any([r.success for r in results])
        for r in results:
            assert r.fail_reason.error_message == "simulation skipped, elevation outside limits for: asterix und underlix"

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
            def prepare_sb(self):
                return
            def get_array_config(self,sb):
                return "c43-3"
            def simulate_HAs(self,planner,array_config,sb):
                results = [f"result_{HA}" for HA in self.HAs]
                return self.HAs,results
        with patch("batch_simulations.application.single_SB_simulation.JobPlanner") as MockJobPlanner:
            MockJobPlanner.return_value = None
            fake = FakeSingleSBSimulation()
            sim_summary = fake.simulate()
            assert sim_summary.date == FakeSingleSBSimulation.date
            assert sim_summary.HAs == FakeSingleSBSimulation.HAs
            assert sim_summary.simulation_results\
                       == [f"result_{HA}" for HA in FakeSingleSBSimulation.HAs]
            assert sim_summary.unexpected_error is None
            #now the case where an unexpected error is thrown:
            error_message = "something went wrong"
            def prepare_sb_malfunction():
                raise RuntimeError(error_message)
            fake.prepare_sb = prepare_sb_malfunction
            result = fake.simulate()
            assert result.date == FakeSingleSBSimulation.date
            assert result.HAs is None
            assert result.simulation_results is None
            for key in ("error_message","full_traceback"):
                assert "RuntimeError" in result.unexpected_error[key]
                assert error_message in result.unexpected_error[key]


class TestAnalysisResult:

    def test_merge(self):
        empty_analysis = AnalysisResult()
        assert empty_analysis.inspection_reasons == []
        assert not empty_analysis.should_be_Waiting
        analysis1 = AnalysisResult(inspection_reasons=["problem"],should_be_Waiting=True)
        analysis2 = AnalysisResult(inspection_reasons=["problem","catastrophy"],
                                   should_be_Waiting=True)
        analysis3 = AnalysisResult(inspection_reasons=["minor"],should_be_Waiting=False)
        for ana in (analysis1,analysis2,analysis3):
            merged = ana.merge(empty_analysis)
            assert len(merged.inspection_reasons) > 0
            assert sorted(merged.inspection_reasons) == sorted(ana.inspection_reasons)
            assert merged.should_be_Waiting == ana.should_be_Waiting
        merged = analysis1.merge(analysis2)
        assert sorted(merged.inspection_reasons) == sorted(["problem","catastrophy"])
        assert merged.should_be_Waiting 
        merged = analysis1.merge(analysis3)
        assert sorted(merged.inspection_reasons) == sorted(["problem","minor"])
        assert merged.should_be_Waiting
        merged = analysis2.merge(analysis3)
        assert sorted(merged.inspection_reasons) == sorted(["problem","catastrophy","minor"])
        assert merged.should_be_Waiting


class TestSingleSBSimulationSummary:

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
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-2*u.hour),"max":Angle(2*u.hour)})
            runnable_HA_amount = SingleSBSimulationSummary.runnable_HA_amount
            HA_is_allowed = SingleSBSimulationSummary.HA_is_allowed
        fake_summary = FakeSummary()
        widths = fake_summary.HA_widths()
        usable_HA_range = fake_summary.runnable_HA_amount()
        expected_range = widths[0]+widths[2]
        assert usable_HA_range == expected_range

    def test_analyse_runnable_HA_range(self):
        class FakeSummary:
            def runnable_HA_amount(self):
                return 2
            analyse_runnable_HA_range = SingleSBSimulationSummary.analyse_runnable_HA_range
        fake_summary = FakeSummary()
        fake_summary.min_runnable_HA_amount = 1
        ana = fake_summary.analyse_runnable_HA_range()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary.min_runnable_HA_amount = 3.3333
        ana = fake_summary.analyse_runnable_HA_range()
        assert ana.inspection_reasons == ["runnable HA range is small"]
        assert not ana.should_be_Waiting

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
        fake_summary = FakeSummary(has_hardcoded=True)
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary = FakeSummary(has_hardcoded=False)
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == []
        assert not ana.should_be_Waiting
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(12*u.hour)}
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == ["unnecessarily restricted HA(s): -2.0"]
        assert not ana.should_be_Waiting
        fake_summary.sb.OT_allowed_HA = {"min":Angle(-1*u.hour),"max":Angle(2*u.hour)}
        ana = fake_summary.analyse_HA_restriction()
        assert ana.inspection_reasons == ["unnecessarily restricted HA(s): -2.0, 4.0"]
        assert not ana.should_be_Waiting

    def test_HA_is_allowed(self):
        class FakeSummary:
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-2.1*u.hour),"max":Angle(3*u.hour)})
            HA_is_allowed = SingleSBSimulationSummary.HA_is_allowed
        fake_summary = FakeSummary()
        for HA in (fake_summary.sb.OT_allowed_HA["min"], fake_summary.sb.OT_allowed_HA["max"],
                   Angle(-2*u.hour),Angle(0*u.hour),Angle(1.1*u.hour)):
            assert fake_summary.HA_is_allowed(HA)
        for HA in (Angle(-2.2*u.hour),Angle(3.1*u.hour)):
            assert not fake_summary.HA_is_allowed(HA)

    def test_analyse_simulation_failures(self):
        def sim_result(success,category=None):
            fail_reason = None if success else\
                   FailReason(error_message="a",error_summary="b",category=category)
            return SimpleNamespace(success=success,fail_reason=fail_reason)
        class FakeSummary:
            HAs = [Angle(-2*u.hour),Angle(0.5*u.hour)]
            sb = SimpleNamespace(OT_allowed_HA={"min":Angle(-4*u.hour),"max":Angle(4*u.hour)})
            analyse_simulation_failures = SingleSBSimulationSummary.analyse_simulation_failures
            HA_is_allowed = SingleSBSimulationSummary.HA_is_allowed

        def assert_correct_analysis(simulation_results,expected_reasons,expected_waiting,
                                    OT_allowed_HA=None):
            fake_summary = FakeSummary()
            if OT_allowed_HA is not None:
                fake_summary.sb.OT_allowed_HA = OT_allowed_HA
            fake_summary.simulation_results = simulation_results
            ana = fake_summary.analyse_simulation_failures()
            assert sorted(ana.inspection_reasons) == sorted(expected_reasons)
            assert ana.should_be_Waiting == expected_waiting

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=True),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error"],
                                expected_waiting=False)

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error"],
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="unobservable"),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False)
        
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error", "simulation failure"],
                                expected_waiting=True)
        
        #making the missing calibrator failure to occur outside of allowed HA:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="server error")],
                                expected_reasons=["server error",],
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="missing calibrator"),
                                                    sim_result(success=False,
                                                               category="exceeds 2h limit")],
                                expected_reasons=["simulation failure"],
                                expected_waiting=True)

        #server error outisde of allowed HA, should not trigger inspection:
        assert_correct_analysis(simulation_results=[sim_result(success=False,
                                                               category="server error"),
                                                    sim_result(success=True)],
                                expected_reasons=[],
                                expected_waiting=False,
                                OT_allowed_HA={"min":Angle(-1*u.hour),"max":Angle(4*u.hour)})

    def test_analyse(self):
        class FakeSummary:
            analyse = SingleSBSimulationSummary.analyse
            unexpected_error = {"error_message":"kabumm",
                                "full_traceback":"the full traceback"}
        fake_summary = FakeSummary()
        fake_summary.analyse()
        assert fake_summary.analysis_result.inspection_reasons == ["unexpected error"]
        assert fake_summary.analysis_result.should_be_Waiting

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
        summary.analyse_no_unexpected_error()
        assert summary.analysis_result.inspection_reasons == []
        assert not summary.analysis_result.should_be_Waiting
        
        #all failures:
        simulation_results = [FakeSimResult(False,category="missing calibrator"),]*5
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse_no_unexpected_error()
        assert sorted(summary.analysis_result.inspection_reasons)\
                  == sorted(["simulation failure","runnable HA range is small"])
        assert summary.analysis_result.should_be_Waiting

        #trigger unnecessary restrikted HA
        HAs[0] = Angle(-14*u.hour)
        simulation_results = [FakeSimResult(True)]*len(HAs)
        summary = SingleSBSimulationSummary(HAs=HAs,simulation_results=simulation_results,
                                            **summary_kwargs)
        summary.analyse_no_unexpected_error()
        assert summary.analysis_result.inspection_reasons == ['unnecessarily restricted HA(s): -14.0']
        assert not summary.analysis_result.should_be_Waiting