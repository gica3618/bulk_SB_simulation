#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:52:29 2026

@author: gianni
"""

import logging
import time


def retry(operation, max_trials, sleep_time, retry_condition):
    for trial in range(1, max_trials + 1):
        logging.info(f"Trial {trial}/{max_trials}")
        result = operation()
        if not retry_condition(result):
            return result
        logging.info(f"Retry condition triggered, sleeping {sleep_time}s")
        time.sleep(sleep_time)
    logging.info("Reached maximum number of trials")
    return result