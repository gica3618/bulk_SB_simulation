#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:18:15 2026

@author: gianni
"""

import datetime
import math
from batch_simulations.infrastructure.calibrator_query import parse_calibrator_candidate_data


def test_parse_calibrator_candidate_data():
    line = ("|[J0529-0519] 1   |  38.2|  68.1| 0.14| 0.23|-0.70+- 0.15|-1.00|  "
            +"40|20240530| 0.107+- 0.008|   49.1+-   3.6|  1.4| True|0.00(0.00)"
            +"|  0.0|   54|-16969.3|    nan| 5.00|         |")
    # line = ("|[J0429+2724] 1              |   2.4|  39.5| 0.11| 0.10|-0.70+- 0.15|-1.00| "
    #         +" 26|20260311| 0.052+- 0.009|    6.3+-   1.1|  3.2| True|0.12(0.12)|  0.0| "
    #         +"  47|-16978.5|    nan|-1.00|       TW|")
    candidate = parse_calibrator_candidate_data(line)
    assert candidate.source_name == "J0529-0519"
    assert candidate.type_IDs == "1"
    assert candidate.Az == 38.2
    assert candidate.El == 68.1
    assert candidate.eRa == 0.14
    assert candidate.eDec == 0.23
    assert candidate.specIndex == -0.7
    assert candidate.specIndex_error == 0.15
    assert candidate.reduced_chi2 == -1
    assert candidate.Nobs == 40
    assert candidate.LastDate == datetime.date(2024,5,30)
    assert candidate.EstimatedFlux == 0.107
    assert candidate.EstimatedFluxError == 0.008
    assert candidate.SNR == 49.1
    assert candidate.SNRError == 3.6
    assert candidate.Sep == 1.4
    assert candidate.isObservable
    assert candidate.fShadow == 0
    assert candidate.fShadow_noncritical == 0
    assert candidate.fRes == 0
    assert candidate.dDays == 54
    assert candidate.UVmax == -16969.3
    assert math.isnan(candidate.UVmin)
    assert candidate.Score == 5
    assert candidate.Reason == ""