#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 09:00:06 2024

@author: gianni
"""

import subprocess
import logging
import datetime


class Calibrator():
    calibrator_types = ['Polarization','Bandpass','Phase','Check','Amplitude',
                        'DGC']
    calibrator_keywords = {cal:[cal,cal.lower()] for cal in calibrator_types}
    
    def __init__(self,name,source_name,cal_type,is_query,coordinates):
        self.name = name
        self.source_name = source_name
        assert cal_type in self.calibrator_types
        self.cal_type = cal_type
        self.is_query = is_query
        self.coordinates = coordinates


class GetCalibratorCandidates():
    
    def __init__(self,xml_filepath):
        self.xml_filepath = xml_filepath

    def run(self,integration_time,array_config,search_radius=None,epoch=None,
            calibrator_type=None,spectral_spec=None,src=None,no_spwavg=False):
        command = f'getCalibratorCandidates.py {self.xml_filepath} -t {integration_time}'\
                 +f' -C {array_config}'
        if search_radius is not None:
            command += f' -r {search_radius}'
        if epoch is not None:
            command += f' -e {epoch}'
        if calibrator_type is not None:
            command += f' -c {calibrator_type}'
        if spectral_spec is not None:
            command += f' --spectralSpec={spectral_spec}'
        if src is not None:
            command += f' --src={src}'
        if no_spwavg:
            #getCalibratorCandidates.py -h tells me that this option is only for
            #phase query, but in the BP, it is used for bandpass.
            #TODO find out who is correct...
            assert calibrator_type == 'phase','no_spwavg is allowed for phase cal query only'
            command += ' --no_spwavg'
        output = subprocess.run(command,shell=True,universal_newlines=True,
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        return GetCalibratorCandidatesResult(output=output)


class GetCalibratorCandidatesResult():

    def __init__(self,output):
        self.output = output
        if self.output.returncode != 0:
            logging.error('getCalibratorCandidates.py crashed')
            logging.error(self.output.stderr)
            self.success = False
        else:
            self.success = True
        self.read_calibrator_candidates()

    def read_calibrator_candidates(self):
        #TODO I tested the case where not calibrator candidate is found, works fine
        #but what if just one is found, does it work?
        stdout = self.output.stdout.splitlines()
        self.calibrator_candidates = []
        for i,line in enumerate(stdout):
            if '-> Listing ranked candidate list...' in line:
                first_calibrator_candidate = i+4
            if '6th col:' in line:
                last_calibrator_candidate = i-2
            if 'sources passed the selection criteria' in line:
                #example: "3 sources passed the selection criteria: [J0529-0519, J0532-0307, J0541-0541]"
                 n_valid_sources = int(line.split('sources')[0])
                 source_names = line.split('[')[1].replace(']','')
                 source_names = source_names.split(',')
                 source_names = [sn.strip() for sn in source_names]
        for line in stdout[first_calibrator_candidate:last_calibrator_candidate+1]:
            calibrator_data = line.split('|')
            assert len(calibrator_data) == 22
            assert calibrator_data[0] == calibrator_data[-1] == ''
            #SourceName (and type IDs)   |    Az|    El|  eRa| eDec|  Spec Index|  RCQ|Nobs|LastDate|Estimated Flux|  SNR (200 sec)|  Sep|isObs|   fShadow| fRes|dDays|   UVmax|  UVmin|Score|   Reason
            #|[J0529-0519] 1   |  38.2|  68.1| 0.14| 0.23|-0.70+- 0.15|-1.00|  40|20240530| 0.107+- 0.008|   49.1+-   3.6|  1.4| True|0.00(0.00)|  0.0|   54|-16969.3|    nan| 5.00|         |'
            print(calibrator_data[1])
            source_name,type_ID = calibrator_data[1].split(']')
            source_name = source_name.replace('[','').strip()
            type_ID = str(type_ID).strip()
            spec_index_data = calibrator_data[6].split('+-')
            datestr =  calibrator_data[9].strip()
            assert len(datestr) == 8
            flux_data = calibrator_data[10].split('+-')
            snr_data = calibrator_data[11].split('+-')
            isObs_data = calibrator_data[13].strip()
            if isObs_data == 'True':
                isObs = True
            else:
                assert isObs_data == 'False',isObs_data
                isObs = False
            calibrator = {'SourceName':source_name,'type_ID':type_ID,'Az':float(calibrator_data[2]),
                          'El':float(calibrator_data[3]),'eRa':float(calibrator_data[4]),
                          'eDec':float(calibrator_data[5]),
                          'specIndex':float(spec_index_data[0]),
                          'specIndex_error':float(spec_index_data[1]),
                          'reduced_chi2':float(calibrator_data[7]),
                          'Nobs':int(calibrator_data[8]),
                          'LastDate':datetime.date(year=int(datestr[:4]),
                                                   month=int(datestr[4:6]),
                                                   day=int(datestr[6:])),
                          'EstimatedFlux':float(flux_data[0]),
                          'EstimatedFluxError':float(flux_data[1]),
                          'SNR':float(snr_data[0]),
                          'SNRError':float(snr_data[1]),'Sep':float(calibrator_data[12]),
                          'isObservable':isObs,
                          'fShadow':calibrator_data[14].strip(),#don't know how to read this
                          'fRes':float(calibrator_data[15]),
                          'dDays':int(calibrator_data[16]),
                          'UVmax':float(calibrator_data[17]),
                          'UVmin':float(calibrator_data[18]),
                          'Score':float(calibrator_data[19]),
                          'Reason':calibrator_data[20].strip()
                          }
            self.calibrator_candidates.append(calibrator)
        


if __name__ == '__main__':
    test_get_cal_cand = GetCalibratorCandidates(
                         xml_filepath='/users/gcataldi/new_bulk_simulations/2023.1.01430.S_MMS_1_a_08_TM2.xml')
    test_result = test_get_cal_cand.run(integration_time=200,array_config='c43-3',
                                        search_radius=5)