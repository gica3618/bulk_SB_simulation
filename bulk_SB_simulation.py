#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 09:05:57 2023

@author: gianni
"""

import os
import shutil
import xml.etree.ElementTree as ET
import csv
import itertools
import numpy as np
import subprocess
import time
import logging
from astropy import units as u
from astropy.coordinates import Angle,SkyCoord
import sys


class Calibrator():

    calibrator_types = ['Polarization','Bandpass','Phase','Check','Amplitude']
    calibrator_keywords = {cal:[cal,cal.lower()] for cal in calibrator_types}
    
    def __init__(self,name,cal_type,xml_data):
        self.name = name
        assert cal_type in self.calibrator_types
        self.cal_type = cal_type
        self.xml_data = xml_data
        
    def is_query(self):
        is_query = self.xml_data.findtext('sbl:isQuery',namespaces=OT_XML_File.namespaces)
        if is_query == 'false':
            return False
        elif is_query == 'true':
            return True
        else:
            raise RuntimeError(f'unknown xml value for isQuery: {is_query}')


class OT_XML_File():
    namespaces = {'sbl':'Alma/ObsPrep/SchedBlock',
                  'prj':"Alma/ObsPrep/ObsProject",
                  'val':"Alma/ValueTypes"}
    long_lat_keys = {'longitude':'ra','latitude':'dec'}

    def __init__(self,filepath,xml_str=None):
        if filepath is not None:
            assert xml_str is None
            tree = ET.parse(filepath)
            self.root = tree.getroot()
        else:
            self.root = ET.fromstring(xml_str)
        self.determine_calibrators()

    #I took the following method directly from simulateSB.py and simplified it
    @classmethod
    def from_download(cls,project_code,SB):
        scriptGetSB = "/groups/science/scripts/P2G/getsb/getsb.py"
        if not os.path.isfile(scriptGetSB):
            scriptGetSB = "/users/ahirota/AIV/science/scripts/P2G/getsb/getsb.py"
        try:
            import cx_Oracle
            serverName = None
        except:
            logging.info("cx_Oracle is not available, and thus will run getsb.py"+
                         " on red-osf")
            serverName = "red-osf.osf.alma.cl"
        if serverName:
            userName = os.getenv("USER")
            cmd = ["ssh"]
            cmd.append(f"{userName}@{serverName}")
            scriptName = "PYTHONPATH=/users/ahirota/local/lib64/python2.6/"\
                          +"site-packages/cx_Oracle-5.2.1-py2.6-linux-x86_64.egg"\
                          +f":$PYTHONPATH {scriptGetSB}"
            cmd.append(f"{scriptName} -p '{project_code}' -s '{SB}'")
            cmd.append("-S ora.sco.alma.cl:1521/ONLINE.SCO.CL")
        else:
            cmd = [scriptGetSB]
            cmd.extend(["-p", project_code, "-s", SB])
            cmd.extend(["-S", "ora.sco.alma.cl:1521/ONLINE.SCO.CL"])
        logging.info("# Retrieving SB xml with the following command [%s]"\
                      % (" ".join(cmd)))
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,env=None)
        xml_str, _ = p.communicate()
        if hasattr(sys.stdout,"detach"):
            # Py3
            xml_str = xml_str.decode("utf-8")
        return cls(filepath=None,xml_str=xml_str)

    def find_unique_element(self,tag):
        elements = self.root.findall(tag,namespaces=self.namespaces)
        n_elements = len(elements)
        assert n_elements == 1, f'found {n_elements} matching elements for {tag}'
        return elements[0]

    def read_allowed_HA(self):
        allowed_HA = {}
        for key in ('minAllowedHA','maxAllowedHA'):
            tag = f'sbl:Preconditions/prj:{key}'
            element = self.find_unique_element(tag)
            HA = float(element.text)
            xml_unit = element.attrib['unit']
            if xml_unit == 'deg':
                unit = u.deg
            elif xml_unit == 'h':
                unit = u.hour
            else:
                raise RuntimeError(f'unknown unit {unit} for {key}')
            HA = Angle(HA,unit=unit)
            output_key = key[:3]
            allowed_HA[output_key] = HA
        assert allowed_HA['min'] < allowed_HA['max']
        return allowed_HA

    def read_RequiresTPAntenna(self):
        text = self.root.findtext('sbl:SchedulingConstraints/sbl:sbRequiresTPAntennas',
                                  namespaces=self.namespaces)
        if text is None:
            return None
        if text == 'true':
            return True
        elif text == 'false':
            return False
        else:
            raise RuntimeError(f'unknown xml content: {text}')

    def read_coordinates(self,coord_data):
        coord = {}
        for xml_key,output_key in self.long_lat_keys.items():
            element = coord_data.find(f'val:{xml_key}',namespaces=self.namespaces)
            assert element.attrib['unit'] == 'deg'
            coord[output_key] = float(element.text)
        return SkyCoord(ra=coord['ra']*u.deg,dec=coord['dec']*u.deg)

    def get_representative_coordinates(self):
        tag = 'sbl:SchedulingConstraints/sbl:representativeCoordinates'
        coord_data = self.find_unique_element(tag)
        return self.read_coordinates(coord_data=coord_data)

    def read_modeName(self):
        return self.root.findtext('sbl:modeName',namespaces=self.namespaces)

    def is_Polarisation(self):
        return self.read_modeName() == 'Polarization Interferometry'

    def is_VLBI(self):
        return 'VLBI' in self.read_modeName()

    def is_Solar(self):
        return 'Solar' in self.read_modeName()

    def get_nominal_config(self):
        return self.root.findtext('sbl:SchedulingConstraints/sbl:nominalConfiguration',
                                  namespaces=self.namespaces)

    def is_7m(self):
        return self.get_nominal_config() == '7M'

    def get_FieldSource_names(self):
        field_sources = self.root.findall("sbl:FieldSource",
                                          namespaces=self.namespaces)
        names = [f.find('sbl:name',namespaces=self.namespaces).text for f in
                 field_sources]
        return names

    def determine_calibrators(self):
        calibrators = []
        field_source_names = self.get_FieldSource_names()
        for source in field_source_names:
            candidate_cal_types = []
            for cal,keywords in Calibrator.calibrator_keywords.items():
                if any([keyword in source for keyword in keywords]):
                    candidate_cal_types.append(cal)
            if len(candidate_cal_types) == 0:
                logging.info(f'Source {source} is not a calibrator')
                continue
            elif len(candidate_cal_types) == 1:
                cal_type = candidate_cal_types[0]
                tag = f"sbl:FieldSource[sbl:name='{source}']"
                cal_elements = self.root.findall(tag,namespaces=self.namespaces)
                for xml_data in cal_elements:
                    calibrators.append(Calibrator(name=source,cal_type=cal_type,xml_data=xml_data))
            else:
                raise RuntimeError(f'could not uniquely identify cal type of source {source}')
        self.calibrators = calibrators

    def get_PolCal(self):
        pol_calibrators = [cal for cal in self.calibrators if
                           cal.cal_type=='Polarization']
        assert len(pol_calibrators) == 1
        return pol_calibrators[0]

    def PolCal_is_fixed(self):
        pol_cal = self.get_PolCal()
        return (not pol_cal.is_query())

    def get_PolCal_coordinates(self):
        pol_cal = self.get_PolCal()
        PolCal_data = pol_cal.xml_data
        coord_data = PolCal_data.find('sbl:sourceCoordinates',
                                      namespaces=self.namespaces)
        return self.read_coordinates(coord_data=coord_data)

    def get_NotetoAoD(self):
        return self.root.findtext('prj:note',namespaces=self.namespaces)


class SingleSimulation():

    def __init__(self,project_code,SB,epoch,array_config):
        self.project_code = project_code
        self.SB = SB
        self.epoch = epoch
        self.array_config = array_config
        self.parent_directory = os.getcwd()
        self.output_folder = os.path.join(self.parent_directory,
                                          f'simulation_output_{project_code}_{SB}')
        if os.path.isdir(self.output_folder):
            self.delete_output_folder()
        self.create_output_folder()

    def run(self):
        self.clean_output_folder()
        #in principle I could test here if the xml is already available from a previous
        #simulation and use it, but the time to get the xml is only ~2s, so for simplicity
        #let's just download the xml every time...
        command = f'simulateSB.py {self.project_code} {self.SB} {self.epoch}'
        if self.array_config is not None:
            command += f' -C {self.array_config}'
            if self.array_config[-4:] == '.cfg':
                logging.info(f'will use cfg file {self.array_config}, copying it'
                             +' to the work folder')
                shutil.copy(src=self.array_config,dst=self.output_folder)
        logging.info(f'going to execute the following command:\n{command}')
        os.chdir(self.output_folder)
        output = subprocess.run(command,shell=True,universal_newlines=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        pipe = {'stdout':output.stdout,'stderr':output.stderr}
        loggers = {'stdout':logging.info,'stderr':logging.error}
        for key,logger in loggers.items():
            logger(f'{key} of simulateSB.py:')
            logger(pipe[key])
        if output.returncode != 0:
            success = False
            pipe_messages = {key:p.split('\n') for key,p in pipe.items()}
            pipe_messages = {key:[m for m in messages if m!=''] for key,messages
                             in pipe_messages.items()}
            found_error_message = False
            #check messages, starting from the latest
            msg_iterator = itertools.zip_longest(pipe_messages['stdout'][::-1],
                                                 pipe_messages['stderr'][::-1],
                                                 fillvalue='')
            max_messages_to_go_back = 5
            for i,(std_msg,error_msg) in enumerate(msg_iterator):
                #give preference to error_msg (i.e. check it first):
                for msg in (error_msg,std_msg):
                    casefolded_msg = msg.casefold()
                    if 'error' in casefolded_msg\
                                         or 'exception' in casefolded_msg:
                        logging.info('identified output error message from failed '
                                     +f'simulation: {msg}')
                        found_error_message = True
                        fail_reason = msg
                        break
                if found_error_message or i >= max_messages_to_go_back-1:
                    break
            if not found_error_message:
                logging.info('did not find error message, will take last output'
                             +' of stdout instead')
                fail_reason = pipe['stdout'][-1]
        else:
            if self.summary_file_reports_success():
                success = True
                fail_reason = None
            else:
                success = False
                fail_reason = 'summary file does not report success'
        os.chdir(self.parent_directory)
        return command,success,fail_reason

    def summary_file_reports_success(self):
        summary_filename = f'log_{self.project_code}_{self.SB}.xml_OSS_summary.txt'
        filepath = os.path.join(self.output_folder,summary_filename)
        with open(filepath) as f:
            contents = f.readlines()
        return 'SUCCESS' in contents[1]

    def create_output_folder(self):
        logging.info(f'going to create {self.output_folder}')
        os.mkdir(self.output_folder)

    def clean_output_folder(self):
        for filename in os.listdir(self.output_folder):
            filepath = os.path.join(self.output_folder,filename)
            logging.info(f'deleting {filepath}')
            os.remove(filepath)

    def delete_output_folder(self):
        logging.info(f'going to delete {self.output_folder}')
        shutil.rmtree(self.output_folder)


class BulkSimulation():

    SB_data_to_extract = ('code','sbname','p2g_account','sb_state',
                          'dc_letter_grade')
    ALMA_site_latitude = np.radians(-23.029)
    min_elevation = np.radians(20)
    max_elevation = np.radians(88)

    def __init__(self,SB_list,custom_SB_filter=None,HAs=None):
        self.SB_list = SB_list
        self.custom_SB_filter = custom_SB_filter
        self.HAs = HAs
        self.read_SB_list()
        self.extract_SB_data_to_simulate()

    def read_SB_list(self):
        with open(self.SB_list,mode='r') as csv_file:
            csv_reader = csv.DictReader(csv_file)
            fieldnames = csv_reader.fieldnames
            data = {fieldname:[] for fieldname in fieldnames}
            for row in csv_reader:
                for fieldname,value in row.items():
                    data[fieldname].append(value)
        self.raw_SB_data = data

    def get_custom_selection(self):
        keys = list(self.raw_SB_data.keys())
        n_raw_SBs = len(self.raw_SB_data[keys[0]])
        if self.custom_SB_filter is None:
            logging.info('no custom selection of SBs')
            custom_selection = [True,]*n_raw_SBs
        else:
            logging.info('doing custom selection of SBs')
            custom_selection = []
            for i in range(n_raw_SBs):
                data = {key:d[i] for key,d in self.raw_SB_data.items()}
                custom_selection.append(self.custom_SB_filter(data))
        return custom_selection

    def extract_SB_data_to_simulate(self):
        raise NotImplementedError

    def determine_min_max_HA_to_simulate(self,OT_xml):
        rep_coord = OT_xml.get_representative_coordinates()
        #DSA will consider the following HA limits:
        if rep_coord.dec.deg >= -5:
            min_HA = Angle(-3*u.hour)
            max_HA = Angle(2*u.hour)
        else:
            min_HA = Angle(-4*u.hour)
            max_HA = Angle(3*u.hour)
        logging.info(f'preliminary HA range: {min_HA.hour} h to {max_HA.hour} h')
        #however, P2G can change the allowed HA in the SB, so we need to adjust
        #the HA if necessary:
        allowed_HA = OT_xml.read_allowed_HA()
        #generally, we expect this to be between -12h and 12h
        #but sometimes it is between 0 and 24h, which needs to be converted
        #to avoid rounding errors, do this work in degrees (Angle(180*u.deg).HA gives
        #12.000000000000002)
        if allowed_HA['max'].deg > 180:
            logging.info('converting HA range to be between -12h and 12h')
            allowed_HA['min'] -= Angle(180*u.deg)
            allowed_HA['max'] -= Angle(180*u.deg)
        for key,HA in allowed_HA.items():
            error_message = f'need allowed_HA ({key}) to wrap at 12h to be'\
                +f' comparable to min_HA and max_HA, but found HA = {HA.deg} deg'
            assert -180 <= HA.deg <= 180,error_message
        assert allowed_HA['min'] < allowed_HA['max']
        min_HA = Angle(max(min_HA.hour,allowed_HA['min'].hour) * u.hour)
        max_HA = Angle(min(max_HA.hour,allowed_HA['max'].hour) * u.hour)
        logging.info(f'HA range after reading xml: {min_HA.hour} h to {max_HA.hour} h')
        #for Pol observations, there is an additional condition on the HA:
        #for the first execution of the SB, Pol Cal needs to have HA between -4 and -0.5
        #but we only need to adjust the lower bound (min_HA). The reason is that if
        #min_HA is increased due to these pol cal restrictions,
        #in the second execution, we have already passed the original min_HA. So
        #there is no possibility that we will have an execution at an HA that is lower
        #than what is allowed by pol cal. On the other hand, we can do the second
        #execution at an HA that is larger than what would be allowed in the
        #first execution. So we should not decrease max_HA
        if OT_xml.is_Polarisation():
            logging.info('polarization SB, going to check if min HA needs'
                         +' to be adjudsted')
            min_pol_cal_HA = Angle(-4*u.hour)
            pol_cal_coord = OT_xml.get_PolCal_coordinates()
            logging.info(f'Pol Cal RA: {pol_cal_coord.ra.deg} deg')
            logging.info(f'representative coord RA: {rep_coord.ra.deg} deg')
            #note that in special cases, delta_ra can be very large although
            #angular separation is small, e.g. if rep_coord.ra=1deg and
            #pol_cal_coord.ra = 359 deg
            delta_ra = rep_coord.ra - pol_cal_coord.ra
            #if delta_ra is positive, then rep coord has larger ra, thus smaller HA
            #because HA = LST - ra; therefore target_HA=pol_cal_HA-delta_ra
            target_HA_for_min_pol_cal_HA = min_pol_cal_HA - delta_ra
            #we need angles between -12 and 12 h in order to compare to min_HA
            target_HA_for_min_pol_cal_HA.wrap_at(12*u.hour,inplace=True)
            min_HA = Angle(max(min_HA.hour,target_HA_for_min_pol_cal_HA.hour)
                           *u.hour)
            logging.info(f'min HA after checking pol cal: {min_HA.hour} h')
        return min_HA,max_HA

    def get_HAs_to_simulate(self,OT_xml):
        if self.HAs is not None:
            return [Angle(HA*u.hour) for HA in self.HAs]
        min_HA,max_HA = self.determine_min_max_HA_to_simulate(OT_xml=OT_xml)
        min_HA,max_HA = min_HA.hour,max_HA.hour
        assert min_HA < max_HA
        HAs = [min_HA,max_HA]
        HA = np.ceil(min_HA).astype(int)
        while HA < max_HA:
            HAs.append(HA)
            HA += 1
        HAs = np.unique(HAs)
        return [Angle(HA*u.hour) for HA in HAs]

    @staticmethod
    def convert_HA_to_str(HA):
        if HA.hour == 0:
            return 'TRANSIT'
        else:
            return f'TRANSIT{HA.hour:+}h'

    def get_array_configs(self,OT_xml):
        raise NotImplementedError

    def get_n_selected_SBs(self):
        return len(list(self.selected_SBs.values())[0])

    def get_elevation_at_ALMA_site(self,DEC,HA):
        sin_elevation = np.sin(self.ALMA_site_latitude)*np.sin(DEC.rad)\
              + np.cos(self.ALMA_site_latitude)*np.cos(DEC.rad)*np.cos(HA.rad)
        return np.arcsin(sin_elevation)

    def run_simulations(self,obs_dates,skip_fixed_cal_SBs=False):
        SB_project_codes = self.selected_SBs['code']
        n_SBs = len(SB_project_codes)
        self.executed_simulations = []
        execution_times = []
        for i,code in enumerate(SB_project_codes):
            start = time.time()
            print(f'simulating SB {i+1}/{n_SBs}')
            SB = self.selected_SBs['sbname'][i]
            print(f'SB name: {SB} ({code})')
            OT_xml = OT_XML_File.from_download(project_code=code,SB=SB)
            if OT_xml.is_Solar():
                logging.info('project is solar, will not simulate')
                continue
            if OT_xml.is_VLBI():
                logging.info('project is VLBI, will not simulate')
                continue
            executed_simulation = {'project_code':code,'SB':SB,
                                   'p2g':self.selected_SBs['p2g_account'][i],
                                   'sb_state':self.selected_SBs['sb_state'][i],
                                   'note_to_aod':OT_xml.get_NotetoAoD(),
                                   'executed_commands_list':[],'epochs':[],
                                   'successful_simulations_list':[],
                                   'failed_simulations_list':[],'fail_reasons_list':[],
                                   'skip_reason':''}
            if OT_xml.is_7m():
                requires_tp = OT_xml.read_RequiresTPAntenna()
                assert requires_tp is not None
                executed_simulation['xml_requires_TP'] = requires_tp
            if OT_xml.is_Polarisation():
                if not OT_xml.PolCal_is_fixed():
                    logging.info('Polarization Calibrator not fixed. Not going to simulate')
                    executed_simulation['skip_reason'] = 'Pol Cal not fixed'
                    self.executed_simulations.append(executed_simulation)
                    continue
                else:
                    logging.info('Pol Cal fixed, as expected')
            if skip_fixed_cal_SBs:
                if any([not cal.is_query() for cal in OT_xml.calibrators]):
                    fixed_calibrators = [cal.cal_type for cal in OT_xml.calibrators
                                         if not cal.is_query()]
                    fixed_calibrators = ','.join(fixed_calibrators)
                    logging.info(f'Calibrators ({fixed_calibrators}) fixed, not going to simulate')
                    executed_simulation['skip_reason'] = f'calibrators ({fixed_calibrators}) fixed'
                    self.executed_simulations.append(executed_simulation)
                    continue
            HAs = self.get_HAs_to_simulate(OT_xml=OT_xml)
            logging.info('HAs to simulate: '+str([f'{HA.hour}h' for HA in HAs]))
            array_configs = self.get_array_configs(OT_xml=OT_xml)
            representative_coord = OT_xml.get_representative_coordinates()
            at_least_1_simulation_ran = False
            for HA,obs_date,array_config in\
                              itertools.product(HAs,obs_dates,array_configs):
                logging.info(f'considering HA={HA.hour}h, config={array_config}')
                #TODO the following is only checking at one point in time (the considered
                #HA); but the observations have a finite duration, so sometimes
                #the elevation is ok at this HA, but not at later HA run during
                #the simulation. Then the simualtion fails; but DSA should consider
                #this, so actually the SB is fine (DSA will not schedule (?))
                #can I detect these cases and skip the corresponding simulations
                #(e.g. read the simulation log file?)
                representative_target_elevation = self.get_elevation_at_ALMA_site(
                                                   DEC=representative_coord.dec,HA=HA)
                if representative_target_elevation < self.min_elevation:
                    logging.info('representative target elevation too low, '
                                 'will not show up in DSA, skipping this HA')
                    continue
                if representative_target_elevation > self.max_elevation:
                    logging.info('representative target elevation too high, '
                                 'will not show up in DSA, skipping this HA')
                    continue
                HA_str = self.convert_HA_to_str(HA=HA)
                epoch = f'{HA_str},{obs_date.isoformat()}'
                logging.info(f'simulating with epoch={epoch}')
                simulation = SingleSimulation(project_code=code,SB=SB,epoch=epoch,
                                              array_config=array_config)
                command,success,fail_reason = simulation.run()
                executed_simulation['executed_commands_list'].append(command)
                executed_simulation['epochs'].append(epoch)
                at_least_1_simulation_ran = True
                if success:
                    executed_simulation['successful_simulations_list'].append(command)
                else:
                    executed_simulation['failed_simulations_list'].append(command)
                    executed_simulation['fail_reasons_list'].append(fail_reason)
                    print(f'simulation failed ({command})')
            if at_least_1_simulation_ran:
                simulation.delete_output_folder()
            else:
                logging.warn('All simulations for {code} were skipped')
                executed_simulation['failed_simulations_list'].append('No simulations were run')
                executed_simulation['fail_reasons_list'].append(
                                  'All simulations skipped, please investigate')
            self.executed_simulations.append(executed_simulation)
            end = time.time()
            execution_times.append(end-start)
            mean_execution_time = np.mean(execution_times)
            n_remaining_projects = n_SBs-i-1
            guess_remaining_time = n_remaining_projects*mean_execution_time
            print(f'remaining time for the {n_remaining_projects} remaining projects:'+
                  f' {int(guess_remaining_time/60)} min')
        #need to do the following outside of the loop (e.g. for case that do not do
        #the whole loop, such as no fixed PolCal)
        for i in range(len(self.executed_simulations)):
            updated_sim = self.executed_simulations[i].copy()
            for key in ('executed_commands_list','successful_simulations_list',
                        'failed_simulations_list','fail_reasons_list'):
                updated_sim[key[:-5]] = '\n'.join(updated_sim[key])
            self.executed_simulations[i] = updated_sim
        self.executed_simulations = sorted(self.executed_simulations,
                                           key=lambda sim: sim['project_code'])

    def write_failed_and_skipped_simulations_to_csvfile(self,filepath):
        fieldnames = ['project_code','SB','p2g','sb_state','failed_simulations',
                      'fail_reasons','skip_reason','note_to_aod']
        with open(filepath,'w',newline='') as csvfile:
            writer = csv.DictWriter(csvfile,fieldnames=fieldnames,
                                    extrasaction='ignore')
            writer.writeheader()
            for sim in self.executed_simulations:
                if len(sim['fail_reasons_list']) == 0 and sim['skip_reason'] == '':
                    continue
                writer.writerow(sim)
        print(f'wrote failed and skipped simulations to {filepath}')

    def get_skip_reasons(self):
        reasons = [sim['skip_reason'] for sim in self.executed_simulations]
        return [r for r in reasons if r!='']

    def get_fail_reasons(self):
        reasons = []
        for sim in self.executed_simulations:
            reasons += sim['fail_reasons_list']
        return reasons

    def get_reason_counts(self,reasons):
        unique_reasons = set(reasons)
        reason_counts = {reason:reasons.count(reason) for reason in unique_reasons}
        #order the dict by the number of counts:
        reason_counts = {reason:count for reason,count in
                         sorted(reason_counts.items(),key=lambda item:-item[1])}
        return reason_counts

    def write_statistics(self,filepath):
        simulations_with_failure = [sim for sim in self.executed_simulations
                                    if len(sim['failed_simulations_list'])>0]
        failed_project_codes = [sim['project_code'] for sim in
                                simulations_with_failure]
        failed_project_codes = set(failed_project_codes)
        all_project_codes = set(self.selected_SBs['code'])
        n_simulations = sum([len(sim['executed_commands_list']) for sim in
                             self.executed_simulations])
        n_failed_simulations = sum([len(sim['failed_simulations_list']) for sim in
                                    simulations_with_failure])
        n_skipped_SBs = len([sim for sim in self.executed_simulations
                             if sim['skip_reason'] != ''])
        skip_reasons = self.get_skip_reasons()
        fail_reasons = self.get_fail_reasons()
        skip_reason_counts = self.get_reason_counts(reasons=skip_reasons)
        fail_reason_counts = self.get_reason_counts(reasons=fail_reasons)
        with open(filepath, "w") as file:
            file.write(f'total number of projects: {len(all_project_codes)}\n')
            file.write('number of projects with at least one failed simulation: '+
                       f'{len(failed_project_codes)}\n')
            file.write(f'total number of SBs: {self.get_n_selected_SBs()}\n')
            file.write('total number of SBs considered: '
                       +f'{len(self.executed_simulations)}\n')
            file.write(f'total number of skipped SBs: {n_skipped_SBs}\n')
            file.write(f'skip reason counts: {skip_reason_counts}\n')
            file.write('number of SBs with at least one failed simulation: '+
                       f'{len(simulations_with_failure)}\n')
            file.write(f'total number of simulations: {n_simulations}\n')
            file.write(f'total number of failed simulations: {n_failed_simulations}\n')
            file.write('fail reason counts:\n')
            for fail_reason,count in fail_reason_counts.items():
                file.write(f'- "{fail_reason}": {count}\n')
        logging.info(f'wrote statistics to {filepath}')


class BulkSimulation12m(BulkSimulation):

    excluded_states = ('FullyObserved','ObservingTimedOut')

    def __init__(self,SB_list,array_config,support_arcs=None,custom_SB_filter=None,
                 HAs=None):
        self.support_arcs = support_arcs
        self.array_config = array_config
        BulkSimulation.__init__(self,SB_list=SB_list,
                                custom_SB_filter=custom_SB_filter,HAs=HAs)

    def get_array_configs(self,OT_xml):
        return [self.array_config,]

    def get_array_config_number(self):
        antennas,config_number = self.array_config.split('-')
        assert antennas == 'c43'
        config_number = int(config_number)
        assert config_number in np.arange(1,11,dtype=int)
        return config_number

    def extract_SB_data_to_simulate(self):
        if self.support_arcs is not None:
            logging.info('going to select following support ARCs: '
                         +str(self.support_arcs))
            support_arc_selection = [(sup_arc in self.support_arcs) for sup_arc in
                                     self.raw_SB_data['support_arc']]
        else:
            logging.info('no support ARC filtering')
            support_arc_selection = [True,]*len(self.raw_SB_data['support_arc'])
        state_selection = [(state not in self.excluded_states) for state in
                           self.raw_SB_data['sb_state']]
        config_number = self.get_array_config_number()
        config_key = f'selected_c{config_number}'
        array_config_selection = [entry=='TRUE' for entry in
                                  self.raw_SB_data[config_key]]
        custom_selection = self.get_custom_selection()
        assert len(support_arc_selection) == len(state_selection)\
                 == len(array_config_selection) == len(custom_selection)
        all_selections = [support_arc_selection,state_selection,
                          array_config_selection,custom_selection]
        selection = list(map(all, zip(*all_selections)))
        self.selected_SBs = {}
        for fieldname in self.SB_data_to_extract:
            self.selected_SBs[fieldname] = list(itertools.compress(
                                           self.raw_SB_data[fieldname],selection))
        logging.info(f'of {len(selection)} SBs in the list, {sum(selection)}'
                     +' were selected')


class BulkSimulation7m(BulkSimulation):

    project_tracker_fieldnames = {'code':'Project Code','sbname':'SB Name',
                                  'p2g_account':'P2G','sb_state':'State',
                                  'dc_letter_grade':'Grade'}

    def __init__(self,SB_list,array_configs=None,custom_SB_filter=None,HAs=None):
        self.array_configs = array_configs
        BulkSimulation.__init__(self,SB_list=SB_list,
                                custom_SB_filter=custom_SB_filter,HAs=HAs)

    def get_array_configs(self,OT_xml):
        if self.array_configs is not None:
            return self.array_configs
        TP_config = ['aca.cm10.pm3.cfg',]
        all_configs = ['7m','aca.cm11.cfg'] + TP_config
        requires_TP = OT_xml.read_RequiresTPAntenna()
        if requires_TP is None:
            #this happens for SBs before cycle 10
            logging.info('did not find RequiresTPAntenna in the xml,'+
                         ' will test all 7m array configs')
            return all_configs
        else:
            if requires_TP:
                logging.info(f'Requires TP, will only test {TP_config}')
                return TP_config
            else:
                logging.info('does not require TP, will test all 7m array configs')
                return all_configs

    def extract_SB_data_to_simulate(self):
        #the list extracted from PT should already be filtered for support ARC,
        #SB state etc., so we only apply a custom filter, if any
        self.selected_SBs = {}
        custom_selection = self.get_custom_selection()
        for key,PT_key in self.project_tracker_fieldnames.items():
            self.selected_SBs[key] = list(itertools.compress(
                                         self.raw_SB_data[PT_key],custom_selection))
        #PT includes the substates, so let's add it to the state:
        sub_states = list(itertools.compress(
                                     self.raw_SB_data['SubState'],custom_selection))
        self.selected_SBs['sb_state'] = [state+substate for state,substate in
                                         zip(self.selected_SBs['sb_state'],sub_states)]
        logging.info(f'of {len(custom_selection)} SBs in the list, '
                     +f'{self.get_n_selected_SBs()} were selected')


class CheckTP_for_7m_SBs():

    '''Check which SBs need TP'''
    config7m_ID = '-C 7m'
    configTP_ID = '-C aca.cm10.pm3.cfg'

    def __init__(self,SB_list,obs_date,custom_SB_filter=None):
        self.SB_list = SB_list
        self.obs_date = obs_date
        self.custom_SB_filter = custom_SB_filter
        self.array_configs = ['aca.cm10.pm3.cfg','7m']

    def check_TP(self,check_results_filepath,statistics_filepath):
        self.general_sim = BulkSimulation7m(
                              SB_list=self.SB_list,array_configs=self.array_configs,
                              HAs=[-1,],custom_SB_filter=self.custom_SB_filter)
        #if any calibrator is fixed, OSS will not determine the SNR of the calibrator
        #using the antenna configuration
        #in that case, we cannot determine whether TP is needed or not, so we skip those SBs
        self.general_sim.run_simulations(obs_dates=[self.obs_date,],skip_fixed_cal_SBs=True)
        self.general_sim.write_statistics(statistics_filepath)
        self.general_sim.write_failed_and_skipped_simulations_to_csvfile(
                                filepath='failed_and_skipped_simulations_7m_TPcheck.csv')
        self.run_additional_simulations()
        self.merge_simulations()
        self.determine_if_TP_is_needed()
        self.sort_executed_simulations()
        self.write_needsTP_to_csvfile(filepath=check_results_filepath)

    def get_failed_simulations(self):
        failed_simulations = []
        for sim in self.general_sim.executed_simulations:
            #If both with TP and without TP fail, then we do additional simulations
            if len(sim['failed_simulations_list']) == 2:
                failed_simulations.append(sim)
        return failed_simulations

    @staticmethod
    def get_sim_ID(sim):
        return '_'.join([sim['project_code'],sim['SB']])

    def run_additional_simulations(self):
        '''For those SBs that fail, we run all HAs'''
        failed_simulations = self.get_failed_simulations()
        self.additional_sim_IDs = [self.get_sim_ID(sim) for sim in
                                   failed_simulations]
        logging.info('going to run additional simulations for '+
                     f'{len(failed_simulations)} SBs')
        def SB_filter(data):
            if self.custom_SB_filter is None:
                custom = True
            else:
                custom = self.custom_SB_filter(data)
            ID = '_'.join([data['Project Code'],data['SB Name']])
            return (ID in self.additional_sim_IDs) and custom
        self.additional_sim = BulkSimulation7m(SB_list=self.SB_list,
                                               array_configs=self.array_configs,
                                               HAs=None,custom_SB_filter=SB_filter)
        self.additional_sim.run_simulations(obs_dates=[self.obs_date,])

    def merge_simulations(self):
        self.executed_simulations = []
        for sim in self.general_sim.executed_simulations:
            general_sim_ID = self.get_sim_ID(sim)
            if general_sim_ID in self.additional_sim_IDs:
                candidates = [additional_sim for additional_sim in
                              self.additional_sim.executed_simulations
                              if self.get_sim_ID(additional_sim)==general_sim_ID]
                assert len(candidates) == 1
                self.executed_simulations.append(candidates[0])
            else:
                self.executed_simulations.append(sim)

    def check_TP_general_sim(self,general_sim):
        logging.info(f'checking TP for {general_sim["SB"]} ({general_sim["project_code"]})')
        if len(general_sim['executed_commands_list']) == 0:
            return 'unable_to_determine'
        else:
            assert len(general_sim['executed_commands_list']) == 2
        if len(general_sim['successful_simulations_list']) == 2:
            logging.info('both simulations successful, no need for TP')
            return False
        elif len(general_sim['successful_simulations_list']) == 0:
            logging.info('both simulations failed, not possible to determine'+
                         ' if TP is needed')
            return 'unable_to_determine'
        else:
            assert len(general_sim['successful_simulations_list'])\
                          == len(general_sim['failed_simulations_list']) == 1
            successful_sim = general_sim['successful_simulations_list'][0]
            failed_sim = general_sim['failed_simulations_list'][0]
            if self.config7m_ID in successful_sim:
                assert self.configTP_ID in failed_sim
                logging.info('7m success, 7m+TP failed. Unexpected, not possible to determine'
                             +'if TP is needed')
                return 'unable_to_determine'
            elif self.configTP_ID in successful_sim:
                assert self.config7m_ID in failed_sim
                logging.info('7m+TP success, 7m failed. This needs TP.')
                return True

    def check_TP_additional_sim(self,additional_sim):
        considered_epochs = set(additional_sim['epochs'])
        needs_TP_per_epoch = []
        for epoch in considered_epochs:
            successful_commands = [command for command in
                                   additional_sim['successful_simulations_list']
                                   if epoch in command]
            assert len(successful_commands) <= 2
            if len(successful_commands) == 2:
                for config_ID in (self.config7m_ID,self.configTP_ID):
                    assert config_ID in ''.join(successful_commands)
                needs_TP_per_epoch.append(False)
            elif len(successful_commands) == 1:
                failed_commands = [command for command in
                                   additional_sim['failed_simulations_list']
                                   if epoch in command]
                assert len(failed_commands) == 1
                successful_command = successful_commands[0]
                failed_command = failed_commands[0]
                if self.config7m_ID in failed_command:
                    assert self.configTP_ID in successful_command
                    needs_TP_per_epoch.append(True)
                else:
                    needs_TP_per_epoch.append('unable_to_determine')
            else:
                needs_TP_per_epoch.append('unable_to_determine')
        if True in needs_TP_per_epoch and False not in needs_TP_per_epoch:
            #at least for one epoch, it is True, and others are unable_to_determine
            return True
        elif False in needs_TP_per_epoch and True not in needs_TP_per_epoch:
            #at least for one epoch if False, and others are unable_to_determine
            return False
        else:
            #all entries unable_to_determine, or both True and False are present
            assert (True in needs_TP_per_epoch and False in needs_TP_per_epoch)\
                  or (True not in needs_TP_per_epoch and False not in needs_TP_per_epoch)
            return 'unable_to_determine'

    def determine_if_TP_is_needed(self):
        updated_simulations = []
        logging.info('determining if TP is needed')
        for sim in self.executed_simulations:
            #need to take a copy such that the original sim does not get modified
            output_sim = sim.copy()
            output_sim['TP requested'] = sim['xml_requires_TP']
            sim_ID = self.get_sim_ID(sim)
            if sim_ID in self.additional_sim_IDs:
                needs_TP = self.check_TP_additional_sim(additional_sim=sim)
            else:
                needs_TP = self.check_TP_general_sim(general_sim=sim)
            output_sim['TP needed'] = needs_TP
            #be careful here; can't simple say "if sim['TP needed']", because sim['TP needed']
            #can be 'unable_to_determine'
            if output_sim['TP needed'] == True and not output_sim['TP requested']:
                output_sim['P2G action'] = 'activate TP'
            if output_sim['TP needed'] == False and output_sim['TP requested']:
                output_sim['P2G action'] = 'deactivate TP'
            else:
                output_sim['P2G action'] = ''
            updated_simulations.append(output_sim)
        self.executed_simulations = updated_simulations

    def sort_executed_simulations(self):
        p2g_action = []
        no_action = []
        for sim in self.executed_simulations:
            no_action.append(sim) if sim['P2G action']=='' else p2g_action.append(sim)
        sort_key = lambda sim: sim['project_code']
        p2g_action.sort(key=sort_key)
        no_action.sort(key=sort_key)
        self.executed_simulations = p2g_action + no_action

    def write_needsTP_to_csvfile(self,filepath):
        fieldnames = ['project_code','SB','p2g','sb_state','failed_simulations',
                      'fail_reasons','skip_reason','successful_simulations',
                      'TP requested','TP needed','P2G action']
        with open(filepath,'w',newline='') as csvfile:
            writer = csv.DictWriter(csvfile,fieldnames=fieldnames,
                                    extrasaction='ignore')
            writer.writeheader()
            for sim in self.executed_simulations:
                writer.writerow(sim)
        print(f'wrote TP assessment to {filepath}')

    def write_needsTP_statistics(self,filepath):
        all_needs_TP = [sim['TP needed'] for sim in self.executed_simulations]
        with open(filepath,"w") as file:
            file.write(f'total number of SBs: {len(all_needs_TP)}\n')
            file.write('needs TP:\n')
            for value in (True,False,'unable_to_determine'):
                n_entries = all_needs_TP.count(value)
                file.write(f'{value}: {n_entries}\n')
        logging.info(f'wrote statistics for needs_TP to {filepath}')


if __name__ == '__main__':
    from datetime import date
    logging.basicConfig(format='%(levelname)s: %(message)s',level=logging.INFO,
                        stream=sys.stdout)
    
    # test_xml = OT_XML_File('example_xml_files/example_polarisation_2022.1.01477.S_polRA_5deg_RepRA_355deg.xml')
    # #test_xml = OT_XML_File('example_xml_files/example_polarisation_2022.1.01477.S_polRA_356deg_RepRA_5deg.xml')
    # test_sim = BulkSimulation12m(SB_list='lookup_table_for_test.csv',
    #                               support_arcs=['EA',],array_config='c43-9')
    # HA_range = test_sim.determine_min_max_HA_to_simulate(OT_xml=test_xml)
    # print(HA_range)

    # def custom_filter(dat):
    #     return dat['dc_letter_grade'] == 'B'
    # test_sim = BulkSimulation12m(SB_list='../lookup_table_12mSBs_20230828.csv',
    #                               support_arcs=['EA',],array_config='c43-9',
    #                               custom_SB_filter=custom_filter)

    # test_sim = BulkSimulation7m(SB_list='../7m_SBs_2023-10-16.csv',
    #                             custom_SB_filter=None)
    
    # test_sim = BulkSimulation7m(SB_list='SBs_7m_test_elevation.csv',
    #                             custom_SB_filter=None)
    # test_sim.run_simulations(obs_dates=[date(year=2023,month=10,day=1),])
    # test_sim.write_failed_simulations_to_csvfile(filepath='test_output.csv')
    # test_sim.print_statistics()
    
    #test_xml = OT_XML_File('../example_xml_files/example_polarisation_2023.1.00013.S.xml')
    #test_xml = OT_XML_File('../example_xml_files/example_polarisation_2022.1.01477.S.xml')
    #test_xml = OT_XML_File('../example_xml_files/example_cycle10_7m_2023.1.01099.S.xml')
    #test_xml = OT_XML_File('../example_xml_files/example_solar_2022.1.01544.S.xml')
    #test_xml = OT_XML_File('../example_xml_files/example_VLBI_2022.1.01268.V.xml')
    #test_xml = OT_XML_File('../example_xml_files/2023.1.00908.S_IRC+1021_a_07_7M.xml')
    test_xml = OT_XML_File('../example_xml_files/example_NoteToAoD_2023.1.00578.S.xml')