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


class OT_XML_File():
    sbl = 'Alma/ObsPrep/SchedBlock'
    prj = "Alma/ObsPrep/ObsProject"
    val = "Alma/ValueTypes"
    long_lat_keys = {'longitude':'ra','latitude':'dec'}

    def __init__(self,filepath,xml_str=None):
        if filepath is not None:
            assert xml_str is None
            tree = ET.parse(filepath)
            self.root = tree.getroot()
        else:
            self.root = ET.fromstring(xml_str)

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
        elements = self.root.findall(tag)
        n_elements = len(elements)
        assert n_elements == 1, f'found {n_elements} matching elements for {tag}'
        return elements[0]

    def read_allowed_HA(self):
        allowed_HA = {}
        for key in ('minAllowedHA','maxAllowedHA'):
            tag = f'{{{self.sbl}}}Preconditions/{{{self.prj}}}{key}'
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
        text = self.root.findtext(f'{{{self.sbl}}}SchedulingConstraints/'+
                                  f'{{{self.sbl}}}sbRequiresTPAntennas')
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
            element = coord_data.find(f'{{{self.val}}}{xml_key}')
            assert element.attrib['unit'] == 'deg'
            coord[output_key] = float(element.text)
        return SkyCoord(ra=coord['ra']*u.deg,dec=coord['dec']*u.deg)

    def get_representative_coordinates(self):
        tag = f'{{{self.sbl}}}SchedulingConstraints'\
              +f'/{{{self.sbl}}}representativeCoordinates'
        coord_data = self.find_unique_element(tag)
        return self.read_coordinates(coord_data=coord_data)

    def read_modeName(self):
        return self.root.findtext(f'{{{self.sbl}}}modeName')

    def is_Polarisation(self):
        return self.read_modeName() == 'Polarization Interferometry'

    def is_VLBI(self):
        return 'VLBI' in self.read_modeName()

    def is_Solar(self):
        return 'Solar' in self.read_modeName()

    def get_PolCal_data(self):
        tag = f"{{{self.sbl}}}FieldSource[{{{self.sbl}}}name='Polarization']"
        return self.find_unique_element(tag)
        
    def PolCal_is_fixed(self):
        PolCal_data = self.get_PolCal_data()
        is_query = PolCal_data.findtext(f'{{{self.sbl}}}isQuery')
        if is_query == 'false':
            return True
        elif is_query == 'true':
            return False
        else:
            raise RuntimeError(f'unknown xml value for isQuery: {is_query}')

    def get_PolCal_coordinates(self):
        PolCal_data = self.get_PolCal_data()
        coord_data = PolCal_data.find(f'{{{self.sbl}}}sourceCoordinates')
        return self.read_coordinates(coord_data=coord_data)


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
        logging.info('going to execute the following command:')
        logging.info(command)
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

    def __init__(self,SB_list,custom_SB_filter=None):
        self.SB_list = SB_list
        self.custom_SB_filter = custom_SB_filter
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

    def determine_HA_range_to_simulate(self,ot_xml):
        rep_coord = ot_xml.get_representative_coordinates()
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
        allowed_HA = ot_xml.read_allowed_HA()
        for key,HA in allowed_HA.items():
            error_message = f'need allowed_HA ({key}) to wrap at 12h to be'\
                             +' comparable to min_HA and max_HA'
            #because of rounding errors, I need to make the following comparison
            #in degrees, not in hour... (Angle(180*u.deg).HA gives
            #12.000000000000002)
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
        if ot_xml.is_Polarisation():
            logging.info('polarization SB, going to check if min HA needs'
                         +' to be adjudsted')
            min_pol_cal_HA = Angle(-4*u.hour)
            pol_cal_coord = ot_xml.get_PolCal_coordinates()
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

    def get_transits(self,min_HA,max_HA):
        min_HA,max_HA = min_HA.hour,max_HA.hour
        assert min_HA < max_HA
        transit_offsets = [min_HA,max_HA]
        trans = np.ceil(min_HA).astype(int)
        while trans < max_HA:
            transit_offsets.append(trans)
            trans += 1
        transit_offsets = np.unique(transit_offsets)
        transits = []
        for transit_offset in transit_offsets:
            if transit_offset == 0:
                transits.append('TRANSIT')
            else:
                sign = np.sign(transit_offset)
                if sign == 1:
                    sign_str = '+'
                else:
                    assert sign == -1
                    sign_str = '-'
                transits.append(f'TRANSIT{sign_str}{np.abs(transit_offset)}h')
        return transits

    def get_array_configs(self,ot_xml):
        raise NotImplementedError

    def run_simulations(self,obs_dates):
        project_codes = self.selected_SBs['code']
        n_projects = len(project_codes)
        self.executed_simulations = []
        execution_times = []
        for i,code in enumerate(project_codes):
            start = time.time()
            print(f'simulating SB {i+1}/{n_projects}')
            SB = self.selected_SBs['sbname'][i]
            print(f'SB name: {SB} ({code})')
            ot_xml = OT_XML_File.from_download(project_code=code,SB=SB)
            if ot_xml.is_Solar():
                logging.info('project is solar, will not simulate')
                continue
            if ot_xml.is_VLBI():
                logging.info('project is VLBI, will not simulate')
                continue
            executed_simulation = {'project_code':code,'SB':SB,
                                   'p2g':self.selected_SBs['p2g_account'][i],
                                   'sb_state':self.selected_SBs['sb_state'][i],
                                   'failed_commands':[],'fail_reasons':[]}
            if ot_xml.is_Polarisation():
                if not ot_xml.PolCal_is_fixed():
                    logging.info('Pol Cal is not fixed. Not going to simulate')
                    executed_simulation['failed_commands'] = 'N/A (no simulation)'
                    executed_simulation['fail_reasons'] = 'Pol Cal is not fixed'
                    self.executed_simulations.append(executed_simulation)
                    continue
                else:
                    logging.info('Pol Cal is fixed, as expected')
            min_HA,max_HA = self.determine_HA_range_to_simulate(ot_xml=ot_xml)
            transits = self.get_transits(min_HA=min_HA,max_HA=max_HA)
            logging.info('transits to simulate: '+str(transits))
            array_configs = self.get_array_configs(ot_xml=ot_xml)
            for transit,obs_date,array_config in\
                              itertools.product(transits,obs_dates,array_configs):
                epoch = f'{transit},{obs_date.isoformat()}'
                logging.info(f'simulating with epoch={epoch}')
                simulation = SingleSimulation(project_code=code,SB=SB,epoch=epoch,
                                              array_config=array_config)
                command,success,fail_reason = simulation.run()
                if not success:
                    executed_simulation['failed_commands'].append(command)
                    executed_simulation['fail_reasons'].append(fail_reason)
                    print(f'simulation failed ({command})')
            simulation.delete_output_folder()
            self.executed_simulations.append(executed_simulation)
            end = time.time()
            execution_times.append(end-start)
            mean_execution_time = np.mean(execution_times)
            n_remaining_projects = n_projects-i-1
            guess_remaining_time = n_remaining_projects*mean_execution_time
            print(f'remaining time for the {n_remaining_projects} remaining projects:'+
                  f' {int(guess_remaining_time/60)} min')

    def write_failed_simulations_to_csvfile(self,filepath):
        fieldnames = list(self.executed_simulations[0].keys())
        with open(filepath,'w',newline='') as csvfile:
            writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
            writer.writeheader()
            for sim in self.executed_simulations:
                if len(sim['failed_commands']) == 0:
                    continue
                for key in ('fail_reasons','failed_commands'):
                    sim[key] = '\n'.join(sim[key])
                writer.writerow(sim)
        print(f'wrote failed simulations to {filepath}')


class BulkSimulation12m(BulkSimulation):

    excluded_states = ('FullyObserved','ObservingTimedOut')

    def __init__(self,SB_list,array_config,support_arcs,custom_SB_filter=None):
        self.support_arcs = support_arcs
        self.array_config = array_config
        BulkSimulation.__init__(self,SB_list=SB_list,
                                custom_SB_filter=custom_SB_filter)

    def get_array_configs(self,ot_xml):
        return [self.array_config,]

    def get_array_config_number(self):
        antennas,config_number = self.array_config.split('-')
        assert antennas == 'c43'
        config_number = int(config_number)
        assert config_number in np.arange(1,11,dtype=int)
        return config_number

    def extract_SB_data_to_simulate(self):
        support_arc_selection = [(sup_arc in self.support_arcs) for sup_arc in
                                 self.raw_SB_data['support_arc']]
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
        logging.info(f'of {len(selection)} SBs in the list, {selection.sum()}'
                     +' were selected')


class BulkSimulation7m(BulkSimulation):

    project_tracker_fieldnames = {'code':'Project Code','sbname':'SB Name',
                                  'p2g_account':'P2G','sb_state':'State',
                                  'dc_letter_grade':'Grade'}

    def get_array_configs(self,ot_xml):
        TP_config = ['aca.cm10.pm3.cfg',]
        all_configs = ['7m','aca.cm11.cfg'] + TP_config
        requires_TP = ot_xml.read_RequiresTPAntenna()
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
        #SB state,... so we only apply a custom filter, if any
        self.selected_SBs = {}
        custom_selection = self.get_custom_selection()
        for key,PT_key in self.project_tracker_fieldnames.items():
            self.selected_SBs[key] = list(itertools.compress(
                                         self.raw_SB_data[PT_key],custom_selection))
        n_selected_SBs = len(list(self.selected_SBs.values())[0])
        logging.info(f'of {len(custom_selection)} SBs in the list, '
                     +f'{n_selected_SBs} were selected')


if __name__ == '__main__':
    from datetime import date
    logging.basicConfig(format='%(levelname)s: %(message)s',level=logging.INFO,
                        stream=sys.stdout)
    
    # test_xml = OT_XML_File('example_xml_files/example_polarisation_2022.1.01477.S_polRA_5deg_RepRA_355deg.xml')
    # #test_xml = OT_XML_File('example_xml_files/example_polarisation_2022.1.01477.S_polRA_356deg_RepRA_5deg.xml')
    # test_sim = BulkSimulation12m(SB_list='lookup_table_for_test.csv',
    #                              support_arcs=['EA',],array_config='c43-9')
    # HA_range = test_sim.determine_HA_range_to_simulate(ot_xml=test_xml)
    # print(HA_range)

    # def custom_filter(dat):
    #     return dat['dc_letter_grade'] == 'B'
    # test_sim = BulkSimulation12m(SB_list='../lookup_table_12mSBs_20230828.csv',
    #                               support_arcs=['EA',],array_config='c43-9',
    #                               custom_SB_filter=custom_filter)

    def custom_filter(dat):
        return dat['Project Code'][:4] != '2023'
    test_sim = BulkSimulation7m(SB_list='../SBs_7m_2023-08-28.csv',
                                custom_SB_filter=None)

    #test_xml = OT_XML_File('example_VLBI_2022.1.01268.V.xml')