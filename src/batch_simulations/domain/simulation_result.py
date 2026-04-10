#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:40:09 2026

@author: gianni
"""

from dataclasses import dataclass
from pathlib import Path
import itertools
import logging


@dataclass
class SimulationResult:
    executed_command: str
    output_folder: Path
    xml_filename: str
    success: bool
    fail_reason: str | None
    server_error: bool

    max_messages_to_go_back = 5

    @classmethod
    def from_completed_process(cls, executed_command, process, output_folder,
                               xml_filename):
        if xml_filename not in executed_command:
            raise ValueError("executed command expected to contain xml filename")
        if process.returncode != 0:
            success = False
            pipe = {'stdout':process.stdout,'stderr':process.stderr}
            fail_reason = cls.get_fail_reason(pipe=pipe)
        else:
            summary_filename = Path(f'log_{xml_filename}_OSS_summary.txt')
            filepath = output_folder / summary_filename
            if cls.summary_file_reports_success(filepath=filepath):
                success = True
                fail_reason = None
            else:
                success = False
                fail_reason = 'summary file does not report success'
        server_error = cls.server_error(fail_reason)
        if fail_reason is not None:
            fail_reason = cls.shorten_fail_reason(fail_reason)
        return cls(executed_command=executed_command,output_folder=output_folder,
                   xml_filename=xml_filename,success=success,
                   fail_reason=fail_reason,server_error=server_error)

    @staticmethod
    def server_error(fail_reason):
        if fail_reason is None:
            return False
        else:
            return ("unexpected response from the source catalogue" in fail_reason)\
                       or ("socket.gaierror" in fail_reason)

    @staticmethod
    def shorten_fail_reason(fail_reason):
        missing_cal_prefix = "Exception: Although 1 source(s) requested, got only 0 for "
        if missing_cal_prefix in fail_reason:
            suffix = " query."
            missing_cal = fail_reason.removeprefix(missing_cal_prefix).removesuffix(suffix)
            return f"no {missing_cal}"
        exceeds_2hlimit_prefix = "Exception: Refusing the SB execution as it will "\
                                 +"exceed the limit (2.00 hours) by "
        if exceeds_2hlimit_prefix in fail_reason:
            excess = fail_reason.removeprefix(exceeds_2hlimit_prefix)
            return f"SB exceeds 2h limit by {excess}"
        return fail_reason

    @classmethod
    def get_fail_reason(cls,pipe):
        pipe_messages = {key:p.split('\n') for key,p in pipe.items()}
        pipe_messages = {key:[m for m in messages if m!=''] for key,messages
                         in pipe_messages.items()}
        #check messages, starting from the latest
        msg_iterator = itertools.zip_longest(pipe_messages['stdout'][::-1],
                                             pipe_messages['stderr'][::-1],
                                             fillvalue='')
        for i,(std_msg,error_msg) in enumerate(msg_iterator):
            #give preference to error_msg (i.e. check it first):
            for msg in (error_msg,std_msg):
                casefolded_msg = msg.casefold()
                if 'error' in casefolded_msg or 'exception' in casefolded_msg:
                    logging.info(f'identified error message: {msg}')
                    return msg
            if i+1 >= cls.max_messages_to_go_back:
                break
        logging.info('did not find error message, will take last output'
                     +' of stdout instead')
        return pipe_messages['stdout'][-1]

    @staticmethod
    def summary_file_reports_success(filepath):
        with open(filepath) as f:
            contents = f.readlines()
        return contents[1].split()[-1] == 'SUCCESS'

    def read_calibrator_queries(self):
        raise NotImplementedError