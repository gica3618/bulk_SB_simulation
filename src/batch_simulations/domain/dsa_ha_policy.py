#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 15:39:13 2026

@author: gianni
"""

from astropy.coordinates import Angle
from astropy import units as u
import logging


class DSAHourAnglePolicy:

    @staticmethod
    def compute(sb):
        #Note: this method does not take into account that DSA does not show SBs
        #if science target or hardcoded calibrators are unobservable (below/above)
        #elevation limit
        #DSA will consider the following HA limits:
        if sb.rep_coord.dec.deg >= -5:
            min_HA_DSA = Angle(-3*u.hour)
            max_HA_DSA = Angle(2*u.hour)
        else:
            min_HA_DSA = Angle(-4*u.hour)
            max_HA_DSA = Angle(3*u.hour)
        logging.info(f'preliminary DSA HA range: {min_HA_DSA.hour} h to {max_HA_DSA.hour} h')
        #Polarization SBs are to be executed sequentially (usually for 2 or 3 executions)
        #to cover enough parallactic angles for calibration.
        #So, for Pol observations, there is an additional condition on the HA:
        #For the first execution of the SB, Pol Cal needs to have HA between
        #-3h and -0.5h.
        #So we might need to adjust the lower bound (min_HA).
        #On the other hand, we can do the second execution at an HA that is larger
        #than what would be allowed in the first execution. So we should not decrease max_HA,
        #because DSA might consider up to max_HA at least for second/third executions
        if sb.is_Polarisation:
            logging.info('polarization SB, going to check if min HA considered'
                         +' by DSA needs to be adjudsted')
            min_pol_cal_HA = Angle(-3*u.hour) #see Best Practices
            pol_cal_coord = sb.get_PolCal().coordinates
            logging.info(f'Pol Cal RA: {pol_cal_coord.ra.deg} deg')
            logging.info(f'representative coord RA: {sb.rep_coord.ra.deg} deg')
            #note that in special cases, delta_ra can be very large although
            #angular separation is small, e.g. if rep_coord.ra=1deg and
            #pol_cal_coord.ra = 359 deg
            delta_ra = sb.rep_coord.ra - pol_cal_coord.ra
            #if delta_ra is positive, then rep coord has larger ra, thus smaller HA
            #because HA = LST - ra; therefore target_HA=pol_cal_HA-delta_ra
            target_HA_at_min_pol_cal_HA = min_pol_cal_HA - delta_ra
            #we need angles between -12 and 12 h in order to compare to min_HA_DSA
            target_HA_at_min_pol_cal_HA.wrap_at(12*u.hour,inplace=True)
            if target_HA_at_min_pol_cal_HA > max_HA_DSA:
                raise ValueError("once Pol cal is at min HA, target is already unobservable")
            #next is just a consistency check, i.e. the Pol Cal is not leading too much
            max_pol_cal_HA = Angle(-0.5*u.hour)
            target_HA_at_max_pol_cal_HA = max_pol_cal_HA - delta_ra
            target_HA_at_max_pol_cal_HA.wrap_at(12*u.hour,inplace=True)
            if target_HA_at_max_pol_cal_HA < min_HA_DSA:
                raise ValueError("once target becomes observable, Pol Cal is "
                                 +"already outside of -3h to -0.5h range")
            #finally, update the min HA considered by DSA
            min_HA_DSA = max(min_HA_DSA,target_HA_at_min_pol_cal_HA)
            logging.info(f'min HA considered by DSA after checking pol cal: {min_HA_DSA.hour} h')
        if min_HA_DSA >= max_HA_DSA:
            raise RuntimeError("invalid DSA HA limits computed "
                               +f"(min: {min_HA_DSA.hour}, max={max_HA_DSA.hour})")
        return {"min":min_HA_DSA,"max":max_HA_DSA}