#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:11:14 2026

@author: gianni
"""
import datetime
from dataclasses import dataclass


@dataclass
class CalibratorCandidate:
    source_name: str
    type_IDs: str
    Az:float
    El:float
    eRa:float
    eDec:float
    specIndex:float
    specIndex_error:float
    reduced_chi2:float
    Nobs:int
    LastDate:datetime.date
    EstimatedFlux:float
    EstimatedFluxError:float
    SNR:float
    SNRError:float
    Sep:float
    isObservable:bool
    fShadow:float
    fShadow_noncritical:float
    fRes:float
    dDays:int
    UVmax:float
    UVmin:float
    Score:float
    Reason:str


def parse_calibrator_candidate_data(line):
    calibrator_data = line.split('|')
    if len(calibrator_data) != 22:
        raise ValueError("Unexpected line format")
    if not (calibrator_data[0] == calibrator_data[-1] == ''):
        raise ValueError("Unexpected line format")
    #SourceName (and type IDs)   |    Az|    El|  eRa| eDec|  Spec Index|  RCQ|Nobs|LastDate|Estimated Flux|  SNR (200 sec)|  Sep|isObs|   fShadow| fRes|dDays|   UVmax|  UVmin|Score|   Reason
    #|[J0529-0519] 1   |  38.2|  68.1| 0.14| 0.23|-0.70+- 0.15|-1.00|  40|20240530| 0.107+- 0.008|   49.1+-   3.6|  1.4| True|0.00(0.00)|  0.0|   54|-16969.3|    nan| 5.00|         |'
    source_name,type_IDs = calibrator_data[1].split(']')
    source_name = source_name.replace('[','').strip()
    type_IDs = str(type_IDs).strip()
    spec_index_data = calibrator_data[6].split('+-')
    datestr =  calibrator_data[9].strip()
    if len(datestr) != 8:
        raise ValueError("unexpected format of date string")
    flux_data = calibrator_data[10].split('+-')
    snr_data = calibrator_data[11].split('+-')
    isObs_data = calibrator_data[13].strip()
    if isObs_data not in ("True", "False"):
        raise ValueError(f"unexpected value of isObs: {isObs_data}")
    isObs = isObs_data == "True"
    fshadow, fshadow_noncritical = calibrator_data[14].split("(")
    fshadow_noncritical = fshadow_noncritical.rstrip(")")
    candidate = CalibratorCandidate(
                   source_name=source_name,
                   type_IDs=type_IDs,
                   Az=float(calibrator_data[2]),
                   El=float(calibrator_data[3]),
                   eRa=float(calibrator_data[4]),
                   eDec=float(calibrator_data[5]),
                   specIndex=float(spec_index_data[0]),
                   specIndex_error=float(spec_index_data[1]),
                   reduced_chi2=float(calibrator_data[7]),
                   Nobs=int(calibrator_data[8]),
                   LastDate=datetime.datetime.strptime(datestr, "%Y%m%d").date(),
                   EstimatedFlux=float(flux_data[0]),
                   EstimatedFluxError=float(flux_data[1]),
                   SNR=float(snr_data[0]),
                   SNRError=float(snr_data[1]),
                   Sep=float(calibrator_data[12]),
                   isObservable=isObs,
                   fShadow=float(fshadow),
                   fShadow_noncritical=float(fshadow_noncritical),
                   fRes=float(calibrator_data[15]),
                   dDays=int(calibrator_data[16]),
                   UVmax=float(calibrator_data[17]),
                   UVmin=float(calibrator_data[18]),
                   Score=float(calibrator_data[19]),
                   Reason=calibrator_data[20].strip())
    return candidate