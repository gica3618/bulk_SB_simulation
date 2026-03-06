#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 12:06:38 2026

@author: gianni
"""

import os
import sys
sys.path.append("..")
import batch_SB_simulations
import shutil
from pathlib import Path
import pytest
import numpy as np
from astropy import units as u
from astropy.coordinates import Angle
import pandas as pd
import pickle
import datetime

def get_xml_filepath(filename):
    xml_folder = 'tests/test_xmls'
    return os.path.join(xml_folder,filename)


class TestSimulationResult:

    output_folder = "tests/simulateSB_outputs"

    def get_simulation_output(self,ID):
        filepath = os.path.join(self.output_folder,ID,f"simulateSB_output_{ID}.pkl")
        with open(filepath, "rb") as f:
            return pickle.load(f)

    def test_wrong_xml(self):
        simulation_output = self.get_simulation_output(ID="success")
        with pytest.raises(AssertionError):
            batch_SB_simulations.SimulationResult(
               executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
               simulation_output=simulation_output,output_folder=os.path.join(self.output_folder,"success"),
                         xml_filename="abc.xml")

    def test_fail_from_summary(self):
        filepath = os.path.join(self.output_folder,"success_but_summaryfile_reports_failed",
                                "simulateSB_output_success.pkl")
        with open(filepath, "rb") as f:
            simulation_output = pickle.load(f)
        result = batch_SB_simulations.SimulationResult(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
                       simulation_output=simulation_output,
                       output_folder=os.path.join(self.output_folder,"success_but_summaryfile_reports_failed"),
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert not result.success
        assert result.fail_reason == 'summary file does not report success'

    def test_success(self):
        simulation_output = self.get_simulation_output(ID="success")
        result = batch_SB_simulations.SimulationResult(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT,2026-02-20 -C c43-3",
                       simulation_output=simulation_output,
                       output_folder=os.path.join(self.output_folder,"success"),
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert result.success
        assert not result.server_error
        assert result.fail_reason is None

    def test_failed(self):
        simulation_output = self.get_simulation_output(ID="failed")
        result = batch_SB_simulations.SimulationResult(
                       executed_command="simulateSB.py 2025.1.00378.S_SchedBlock0.xml TRANSIT-4H,2026-02-20 -C c43-3",
                       simulation_output=simulation_output,
                       output_folder=os.path.join(self.output_folder,"failed"),
                       xml_filename="2025.1.00378.S_SchedBlock0.xml")
        assert not result.success
        assert not result.server_error
        assert result.fail_reason == "no check"


class TestSBSimulator:

    filename = "example_polarisation_2023.1.00013.S.xml"
    xml_filepath = get_xml_filepath(filename)
    simulator = batch_SB_simulations.SBSimulator(xml_filepath=xml_filepath)
    epoch = "TRANSIT-1h"

    def test_init(self):
        assert self.simulator.xml_filename == self.filename
        assert self.simulator.xml_name == "example_polarisation_2023.1.00013.S"

    def test_work_folder_path(self):
        self.simulator.epoch = self.epoch
        array_config = "test.cfg"
        self.simulator.array_config = array_config
        expected_path = os.path.join(os.getcwd(), f'simulation_output_{self.simulator.xml_name}_{self.epoch}_{array_config}')
        assert expected_path == self.simulator.get_work_folder_path()

    def test_create_folder(self):
        self.simulator.epoch = self.epoch
        array_config = "test.cfg"
        self.simulator.array_config = array_config
        folderpath = self.simulator.get_work_folder_path()
        self.simulator.create_work_folder()
        assert os.path.isdir(folderpath)
        with pytest.raises(FileExistsError):
            self.simulator.create_work_folder()
        os.rmdir(folderpath)

    def test_prepare_simulations(self):
        def get_foldername(epoch,array_config):
            return f'simulation_output_{self.simulator.xml_name}_{epoch}_{array_config}'
        cwd = os.getcwd()
        array_config = "c43-6"
        self.simulator.prepare_simulation(epoch=self.epoch,array_config=array_config)
        expected_work_folder = os.path.join(cwd,get_foldername(
                                epoch=self.epoch, array_config=array_config))
        print(expected_work_folder)
        assert os.path.isdir(expected_work_folder)
        os.rmdir(expected_work_folder)
        #case with config file:
        array_config = "welrkj.cfg"
        array_config_filepath = os.path.join(cwd,array_config)
        Path(array_config_filepath).touch()
        expected_work_folder = os.path.join(cwd,get_foldername(
                                  epoch=self.epoch, array_config=array_config))
        self.simulator.prepare_simulation(epoch=self.epoch,array_config=array_config)
        assert os.path.isdir(expected_work_folder)
        assert os.path.isfile(os.path.join(expected_work_folder,array_config))
        shutil.rmtree(expected_work_folder)
        os.remove(array_config_filepath)

    def test_get_command(self):
        self.simulator.epoch = self.epoch
        array_config = "test.cfg"
        self.simulator.array_config = array_config
        expected_command = f'simulateSB.py {self.xml_filepath} {self.epoch}'\
                         +f' -C {array_config}'
        assert expected_command == self.simulator.get_command()

    def test_copy_cfg_file(self):
        self.simulator.epoch = self.epoch
        array_config = "test.cfg"
        array_config_filepath = os.path.join(os.getcwd(),array_config)
        Path(array_config_filepath).touch()
        self.simulator.array_config = array_config
        folderpath = self.simulator.get_work_folder_path()
        self.simulator.create_work_folder()
        self.simulator.copy_cfg_file()
        assert os.path.isfile(os.path.join(folderpath,array_config))
        shutil.rmtree(folderpath)
        os.remove(array_config_filepath)


class TestSB:

    @staticmethod
    def test_init():
        
        def assert_is_not_XX(sb,exluded):
            for is_XX in ("is_7m","is_VLBI","is_Polarisation","is_B2B","is_Solar"):
                if is_XX in exluded:
                    continue
                assert not getattr(sb, is_XX)
        
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_polarisation_2023.1.00013.S.xml"))
        assert sb.is_Polarisation
        assert_is_not_XX(sb=sb, exluded=["is_Polarisation",])
        
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_VLBI_2022.1.01268.V.xml"))
        assert sb.is_VLBI
        assert_is_not_XX(sb=sb, exluded=["is_VLBI",])

        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_solar_2022.1.01544.S.xml"))
        assert sb.is_Solar
        assert_is_not_XX(sb=sb, exluded=["is_Solar",])

        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01389.S_B2B.xml"))
        assert sb.is_B2B
        assert_is_not_XX(sb=sb, exluded=["is_B2B",])

        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_NoteToAoD_2023.1.00578.S.xml"))
        assert sb.is_7m
        assert_is_not_XX(sb=sb, exluded=["is_7m",])

        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        assert_is_not_XX(sb=sb,exluded=[])
        
    @staticmethod
    def test_has_no_PolCal():
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_polarisation_2023.1.00013.S.xml"))
        assert not sb.has_no_PolCal()
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        assert sb.has_no_PolCal()

    @staticmethod
    def test_has_at_least_one_pol_cal():
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_polarisation_2023.1.00013.S.xml"))
        assert sb.has_at_least_one_PolCal()
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        assert not sb.has_at_least_one_PolCal()

    @staticmethod
    def test_get_calibrator():
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_polarisation_2023.1.00013.S.xml"))
        pol_cal = sb.get_calibrator("Polarization")
        assert pol_cal.name == "Polarization calibrator"
        assert pol_cal.source_name == "J0522-3627"
        bandpass = sb.get_calibrator("Bandpass")
        assert bandpass.name == "Bandpass"
        assert bandpass.source_name == "J0538-4405"
        phase = sb.get_calibrator("Phase")
        assert phase.name == "Phase"
        assert phase.source_name == "query"
        #attempting to get PolCal from non-polarization SB, should fail:
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        with pytest.raises(AssertionError):
            sb.get_calibrator("Polarization")

    @staticmethod
    def test_consistency_checks():
        batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_several_polcal_only_one_in_obsgroups.xml"))
        with pytest.raises(ValueError):
            batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_several_polcal_in_obsgroups.xml"))
        with pytest.raises(ValueError):
            batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_without_PolCal.xml"))
        with pytest.raises(ValueError):
            batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_PolCal_not_hardcoded.xml"))
        with pytest.raises(ValueError):
            batch_SB_simulations.SB(xml_filepath=get_xml_filepath("NonPolarisation_with_PolCal.xml"))

    @staticmethod
    def test_determine_minmax_HA_considered_by_DSA():
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == 42.19541666666667
        assert sb.min_HA_DSA.hour == -3
        assert sb.max_HA_DSA.hour == 2
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("general_DEC-15.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == -15.39436944
        assert sb.min_HA_DSA.hour == -4
        assert sb.max_HA_DSA.hour == 3
        #example with leading PolCal:
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("example_polarisation_2023.1.00013.S.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == -5.376798611111111
        assert sb.xml.get_representative_coordinates().ra.deg == 83.80885625
        assert sb.get_PolCal().coordinates.ra.deg == 80.7416026833
        #this is a case where the min needs to be adjusted, since without PolCal,
        #it would be -4
        assert sb.min_HA_DSA.hour == -3 - (83.80885625-80.7416026833)/360*24
        assert sb.max_HA_DSA.hour == 3
        #example with leading PolCal around RA=0
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_trailing_PolCal_around_RA0deg.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == -59.52381111
        assert sb.xml.get_representative_coordinates().ra.deg == 1
        assert sb.get_PolCal().coordinates.ra.deg == 359.471941933
        delta_ra = (360+1)-359.471941933
        assert np.isclose(sb.min_HA_DSA.hour, -3 - delta_ra/360*24, atol=0, rtol=1e-8)
        assert sb.max_HA_DSA.hour == 3
        #example with trailing PolCal
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_leading_PolCal.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == -59.52381111
        assert sb.xml.get_representative_coordinates().ra.deg == 165.0
        assert sb.get_PolCal().coordinates.ra.deg == 165.814104429
        assert sb.min_HA_DSA.hour == -3 + (165.814104429-165)/360*24
        assert sb.max_HA_DSA.hour == 3
        #example with trailing PolCal around RA=0
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisatoin_trailing_PolCal_around_RA0.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == -59.52381111
        assert sb.xml.get_representative_coordinates().ra.deg == 355.0
        assert sb.get_PolCal().coordinates.ra.deg == 1.14856453333
        delta_ra = (360+1.14856453333) - 355
        assert sb.min_HA_DSA.hour == -3 + delta_ra/360*24
        assert sb.max_HA_DSA.hour == 3
        #leading PolCal, no adjustement needed (anyway min is -3 because of high DEC)
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_leading_PolCal_highDEC.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        assert sb.xml.get_representative_coordinates().dec.deg == 0.5238108333333333
        assert sb.xml.get_representative_coordinates().ra.deg == 150
        assert sb.get_PolCal().coordinates.ra.deg == 139.933497596
        assert sb.min_HA_DSA.hour == -3
        assert sb.max_HA_DSA.hour == 2
        #PolCal leading too much
        with pytest.raises(ValueError):
            sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_PolCal_leads_too_much.xml"))
        #PolCal trailing too much
        with pytest.raises(ValueError):
            sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("Polarisation_PolCal_trails_too_much.xml"))

    @staticmethod
    def test_determine_HA_to_simulate():
        sb = batch_SB_simulations.SB(xml_filepath=get_xml_filepath("2025.1.01279.S_general_SB.xml"))
        sb.determine_minmax_HA_considered_by_DSA()
        sim_HA = sb.determine_HAs_to_simulate(HA_step=1)
        assert np.all([HA.hour for HA in sim_HA] == np.arange(-3,2.1))
        sim_HA = sb.determine_HAs_to_simulate(HA_step=0.25)
        assert np.all([HA.hour for HA in sim_HA] == np.arange(-3,2.1,0.25))
        with pytest.raises(AssertionError):
            sb.determine_HAs_to_simulate(HA_step=-1)
        sb.max_HA_DSA = Angle(2.1*u.hour)
        with pytest.raises(AssertionError):
            sb.determine_HAs_to_simulate(HA_step=1)
        sb.min_HA_DSA = Angle(2*u.hour)
        sb.max_HA_DSA = Angle(1*u.hour)
        with pytest.raises(AssertionError):
            sb.determine_HAs_to_simulate(HA_step=1)
        sb.min_HA_DSA = Angle(-3.2*u.hour)
        sb.max_HA_DSA = Angle(3*u.hour)
        sim_HA = sb.determine_HAs_to_simulate(HA_step=0.5)
        assert np.all([HA.hour for HA in sim_HA] == np.arange(-3,3.1,0.5))

    @staticmethod
    def test_get_epoch():
        HA_0 = Angle(0*u.deg)
        epoch = batch_SB_simulations.SB.get_epoch(HA=HA_0, date=None)
        assert epoch == "TRANSIT"
        epoch = batch_SB_simulations.SB.get_epoch(HA=HA_0, date=datetime.date(2024,4,1))
        assert epoch == "TRANSIT,2024-04-01"
        HA_1 = Angle(1*u.hour)
        epoch = batch_SB_simulations.SB.get_epoch(HA=HA_1, date=None)
        assert epoch == "TRANSIT+1h"
        HA_minus25 = Angle(-2.5*u.hour)
        epoch = batch_SB_simulations.SB.get_epoch(HA=HA_minus25, date=None)
        assert epoch == "TRANSIT-2.5h"
        epoch = batch_SB_simulations.SB.get_epoch(HA=HA_minus25, date=datetime.date(2020,11,2))
        assert epoch == "TRANSIT-2.5h,2020-11-02"


class TestBatchSimulations:
    
    input_folder = "tests/test_input_tables"
    input_filepaths = {"12m":os.path.join(input_folder,"configuration_lookup_table_cycle_12.csv"),
                       "7m":os.path.join(input_folder,"2026-02-06_schedBlockList_7M.csv")}
    obs_dates = [datetime.date(2025,1,3),datetime.date(2025,2,14)]
    
    def test_read_12m_array_config_number(self):
        kwargs = {"input_filepaths":self.input_filepaths,"obs_dates":self.obs_dates}
        for array_config_12m,number in zip(("c43-3","c43-10"),(3,10)):
            batch_sim = batch_SB_simulations.BatchSimulations(
                                    **kwargs,array_config_12m=array_config_12m)
            assert batch_sim.array_config_number_12m == number
        for array_config_12m in ("c-1","c43-0","c43-11","c42-5"):
            with pytest.raises(AssertionError):
                batch_sim = batch_SB_simulations.BatchSimulations(
                               **kwargs,array_config_12m=array_config_12m)

    def test_read_12m_input(self):
        input_filepaths = {"12m":self.input_filepaths["12m"],"7m":None}
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m="c43-5")
        data = pd.read_csv(input_filepaths["12m"])
        data_selection = ~data["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))\
                         & data["selected_c5"]
        assert len(batch_sim.input_data) == len(data[data_selection])
        assert sorted(batch_sim.input_data["sb_uid"])\
               == sorted(data[data_selection]["sb_uid"])
        assert sorted(batch_sim.input_data.columns)\
                  == sorted(batch_SB_simulations.BatchSimulations.input_columns)
        #explicit check using the manually (with Excel Filter) extracted uids:
        uids = pd.read_csv(os.path.join(self.input_folder,"c5_uids.csv"))
        assert sorted(batch_sim.input_data["sb_uid"]) == sorted(uids["sb_uid"])
        

    def test_read_7m_input(self):
        input_filepaths = {"12m":None,"7m":self.input_filepaths["7m"]}
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m=None)
        data = pd.read_csv(input_filepaths["7m"])
        assert len(batch_sim.input_data) == len(data)
        assert sorted(batch_sim.input_data.columns)\
                  == sorted(batch_SB_simulations.BatchSimulations.input_columns)

    def test_read_input(self):
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=self.input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m="c43-4")
        print("columns: ",batch_sim.input_data.columns)
        data12m = pd.read_csv(self.input_filepaths["12m"])
        data_selection_12m = data12m["selected_c4"]\
                & ~data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        data7m = pd.read_csv(self.input_filepaths["7m"])
        assert len(batch_sim.input_data) == len(data7m) + len(data12m[data_selection_12m])
        assert sorted(batch_sim.input_data["sb_uid"])\
               == sorted(list(data12m[data_selection_12m]["sb_uid"]) + list(data7m["SB UID"]))
        assert sorted(batch_sim.input_data.columns)\
                  == sorted(batch_SB_simulations.BatchSimulations.input_columns)

    def test_filter(self):
        p2g = "gianni"
        def p2g_filter(row):
            return row["p2g_account"] == p2g
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=self.input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m="c43-4",SB_filter=p2g_filter)
        data12m = pd.read_csv(self.input_filepaths["12m"])
        data_selection_12m = data12m["selected_c4"] & (data12m["p2g_account"]==p2g)\
                            &  ~data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        data12m_filtered = data12m[data_selection_12m]
        assert len(data12m_filtered) > 0
        data7m = pd.read_csv(self.input_filepaths["7m"])
        data7m_filtered = data7m[data7m["P2G"]==p2g]
        assert len(data7m_filtered) > 0
        assert len(batch_sim.input_data) == len(data12m_filtered) + len(data7m_filtered)
        assert sorted(batch_sim.input_data["sb_uid"])\
               == sorted(list(data12m_filtered["sb_uid"]) + list(data7m_filtered["SB UID"]))
        code = "2025.1.00240.S"
        def project_filter(row):
            return row["code"] == code
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=self.input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m="c43-4",SB_filter=project_filter)
        data_selection_12m = data12m["selected_c4"] &  (data12m["code"]==code)\
                            &  ~data12m["sb_state"].isin( ('FullyObserved','ObservingTimedOut'))
        data12m_filtered = data12m[data_selection_12m]
        assert len(data12m_filtered) == 0
        data7m_filtered = data7m[data7m["Project Code"]==code]
        assert len(data7m_filtered) > 0
        assert len(batch_sim.input_data) == len(data12m_filtered) + len(data7m_filtered)
        assert sorted(batch_sim.input_data["sb_uid"])\
               == sorted(list(data12m_filtered["sb_uid"]) + list(data7m_filtered["SB UID"]))

    def test_get_array_config(self):
        config_12m = "c43-4"
        xml_folder = "tests/test_xmls"
        batch_sim = batch_SB_simulations.BatchSimulations(
                       input_filepaths=self.input_filepaths,obs_dates=self.obs_dates,
                       array_config_12m=config_12m)
        sb = batch_SB_simulations.SB(
               xml_filepath=os.path.join(xml_folder,"example_NoteToAoD_2023.1.00578.S.xml"))
        assert batch_sim.get_array_config(sb=sb) == "aca.cm10.pm3.cfg"
        sb = batch_SB_simulations.SB(xml_filepath=os.path.join(xml_folder,"2025.1.01279.S_general_SB.xml"))
        assert batch_sim.get_array_config(sb=sb) == config_12m
        sb = batch_SB_simulations.SB(xml_filepath=os.path.join(xml_folder,"example_general_7M_Band3.xml"))
        assert batch_sim.get_array_config(sb=sb) == "7m"

    # def test_add_sb_info_to_output(self):
    #     output = {"a":1,"b":2}
    #     batch_sim = batch_SB_simulations.BatchSimulations(
    #                    input_filepaths=self.input_filepaths,array_config_12m="c43-4")
    #     sb = batch_SB_simulations.SB(xml_filepath="tests/test_xmls/example_NoteToAoD_2023.1.00578.S.xml")
    #     output = batch_sim.add_sb_info_to_output(sb=sb, output=output)
    #     expected_output = {"a":1,"b":2,"note_to_AoD":"test Note to AoD",
    #                        "hardcoded_calibrator(s)":False,'Polarization':None,
    #                        'Bandpass':"query",'Phase':"query",'Check':"query",
    #                        'Amplitude':None,'DGC':None,"min_HA_OT":-3,"max_HA_OT":2.25,
    #                        "min_HA_DSA":-4,"max_HA_DSA":3,
    #                        "simulated_config":"aca.cm10.pm3.cfg"}
    #     assert output == expected_output