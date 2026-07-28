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


class FailReason:

    CATEGORIES = ("missing calibrator",
                  "exceeds 2h limit",
                  "unobservable",
                  "server error",
                  "other")

    def __init__(self,error_message,error_summary,category):
        self.error_message = error_message
        self.error_summary = error_summary
        if not category in self.CATEGORIES:
            raise ValueError("invalid category '{category}'")
        self.category = category

    @classmethod
    def from_error_message(cls,error_message):
        prefix = "Exception: Although 1 source(s) requested, got only 0 for"
        if error_message.startswith(prefix):
            missing_cal = error_message.removeprefix(prefix).removesuffix("query.").strip()
            error_summary = f"no {missing_cal}"
            return cls(error_message=error_message,error_summary=error_summary,
                       category="missing calibrator")

        prefix = "Exception: Refusing the SB execution as it will exceed the limit (2.00 hours) by"
        if error_message.startswith(prefix):
            excess = error_message.removeprefix(prefix).strip()
            error_summary = f"SB exceeds 2h limit by {excess}"
            return cls(error_message=error_message,error_summary=error_summary,
                       category="exceeds 2h limit")

        prefix = "Observation.SBExecutionMode.SBExecutionError:"
        if error_message.startswith(prefix):
            message = error_message.removeprefix(prefix).strip()
    
            if message.startswith("No visible science target"):
                error_summary = message
    
            elif message.startswith("All science targets in") and message.endswith("are unobservable"):
                error_summary = message
    
            elif (message.startswith("Although execution of")
                                     and "it is not observable" in message):
                parts = message.split()
                calibrator_type = parts[3].replace("'", "")
                calibrator_name = parts[5]
                error_summary = f"{calibrator_type} {calibrator_name} not observable"
            else:
                error_summary = error_message
            return cls(error_message=error_message,error_summary=error_summary,
                       category="unobservable")

        if error_message.startswith("Exception: Specified elevation") and\
                                   error_message.endswith("is out of range."):
            return cls(error_message=error_message,error_summary="elevation out of range",
                       category="unobservable")

        if ("unexpected response from the source catalogue" in error_message)\
                   or ("socket.gaierror" in error_message):
            return cls(error_message=error_message,error_summary=error_message,
                       category="server error")

        return cls(error_message=error_message,error_summary=error_message,
                   category="other")


@dataclass
class SimulationResult:
    executed_command: str | None
    output_folder: Path | None
    xml_filename: str
    success: bool
    fail_reason: FailReason | None

    max_messages_to_go_back = 5

    @classmethod
    def from_completed_process(cls, executed_command, process, output_folder,
                               xml_filename):
        if xml_filename not in executed_command:
            raise ValueError("executed command expected to contain xml filename")
        if process.returncode != 0:
            success = False
            pipe = {'stdout':process.stdout,'stderr':process.stderr}
            error_message = cls.get_error_message(pipe=pipe)
            fail_reason = FailReason.from_error_message(error_message=error_message)
        else:
            summary_filename = Path(f'log_{xml_filename}_OSS_summary.txt')
            filepath = output_folder / summary_filename
            if cls.summary_file_reports_success(filepath=filepath):
                success = True
                fail_reason = None
            else:
                success = False
                fail_reason = FailReason(
                                error_message=None,
                                error_summary='summary file does not report success',
                                category="other")
        return cls(executed_command=executed_command,output_folder=output_folder,
                   xml_filename=xml_filename,success=success,
                   fail_reason=fail_reason)

    @classmethod
    def get_error_message(cls,pipe):
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