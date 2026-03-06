#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:40:09 2026

@author: gianni
"""

from dataclasses import dataclass
import os
import itertools
import logging


@dataclass
class SimulationResult:
    executed_command: str
    output_folder: str
    xml_filename: str
    success: bool
    fail_reason: str | None
    server_error: bool

    @classmethod
    def from_completed_process(cls, executed_command, process, output_folder,
                               xml_filename):
        if xml_filename not in executed_command:
            raise ValueError("executed command expected to contain xml filename")
        if process.returncode != 0:
            success = False
            pipe = {'stdout':process.stdout,'stderr':process.stderr}
            raw_fail_reason = cls.get_raw_fail_reason(pipe=pipe)
        else:
            summary_filename = f'log_{xml_filename}_OSS_summary.txt'
            filepath = os.path.join(output_folder,summary_filename)
            if cls.summary_file_reports_success(filepath=filepath):
                success = True
                raw_fail_reason = None
            else:
                success = False
                raw_fail_reason = 'summary file does not report success'
        server_error = cls.server_error(raw_fail_reason)
        if raw_fail_reason is not None:
            fail_reason = cls.shorten_fail_reason(raw_fail_reason)
        return cls(executed_command=executed_command,output_folder=output_folder,
                   xml_filename=xml_filename,success=success,
                   fail_reason=fail_reason,server_error=server_error)

    @staticmethod
    def server_error(raw_fail_reason):
        if raw_fail_reason is None:
            return False
        else:
            return ("unexpected response from the source catalogue" in raw_fail_reason)\
                       or ("socket.gaierror" in raw_fail_reason)

    @staticmethod
    def shorten_fail_reason(raw_fail_reason):
        missing_cal_prefix = "Exception: Although 1 source(s) requested, got only 0 for "
        if missing_cal_prefix in raw_fail_reason:
            suffix = " query."
            missing_cal = raw_fail_reason.removeprefix(missing_cal_prefix).removesuffix(suffix)
            return f"no {missing_cal}"
        exceeds_2hlimit_prefix = "Exception: Refusing the SB execution as it will "\
                                 +"exceed the limit (2.00 hours) by "
        if exceeds_2hlimit_prefix in raw_fail_reason:
            excess = raw_fail_reason.removeprefix(exceeds_2hlimit_prefix)
            return f"SB exceeds 2h limit by {excess}"
        return raw_fail_reason

    @staticmethod
    def get_raw_fail_reason(pipe):
        pipe_messages = {key:p.split('\n') for key,p in pipe.items()}
        pipe_messages = {key:[m for m in messages if m!=''] for key,messages
                         in pipe_messages.items()}
        #check messages, starting from the latest
        msg_iterator = itertools.zip_longest(pipe_messages['stdout'][::-1],
                                             pipe_messages['stderr'][::-1],
                                             fillvalue='')
        max_messages_to_go_back = 5
        for i,(std_msg,error_msg) in enumerate(msg_iterator):
            #give preference to error_msg (i.e. check it first):
            for msg in (error_msg,std_msg):
                casefolded_msg = msg.casefold()
                if 'error' in casefolded_msg or 'exception' in casefolded_msg:
                    logging.info(f'identified error message: {msg}')
                    return msg
            if i+1 >= max_messages_to_go_back:
                break
        logging.info('did not find error message, will take last output'
                     +' of stdout instead')
        return pipe_messages['stdout'][-1]

    @staticmethod
    def summary_file_reports_success(filepath):
        with open(filepath) as f:
            contents = f.readlines()
        return 'SUCCESS' in contents[1]

    def read_calibrator_queries(self):
        raise NotImplementedError