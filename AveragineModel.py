#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 09:32:10 2017

@author: sgiese
"""
import numpy as np
import pyopenms as oms
from scipy import stats


def averagine_model(mass, n_peaks=5, only_intensity=False):
    """
    Give a mass(!) of a peptide returns the expected relative intensity
    distribution by averagine modelling.
    
    Parameters
    ----------
    mass:
        Mass of a peptide [Da]
        
    n_peaks:
        Number of isotope peaks to compute.
    
    only_intensity:
        bool, if True, only returns the intensities (no mass column)

    Returns
    -------
    np.array
        Array of floats representing the abundance of the isotope peaks and
        their mass values as computed by OpenMS.


    Example:
    -------
        >>> pepmass = 1500
        >>> averagine_model(pepmass)

    """
    isotope_dist = oms.IsotopeDistribution()
    isotope_dist.setMaxIsotope(n_peaks)
    isotope_dist.estimateFromPeptideWeight(mass)
    peaks = np.array(isotope_dist.getContainer())
    
    if only_intensity:
        return(peaks[:,1])
    else:
        return(peaks)


def score_cluster(mz_temp, int_temp, charge_temp):
    """
    Scores an isotope cluster using simple pearson correlation.
    The lightest peak from the mz array  will be used to compute the
    averagine intensities. All peaks are used for the scoring.

    :param mz_temp:
    :param int_temp:
    :param charge_temp:
    :return:
    """
    mass = mz_temp[0]*charge_temp
    peaks = averagine_model(mass, len(mz_temp))
    peaks[:,0] = mz_temp
#    print(peaks[:,1])
#    print(int_temp/int_temp.sum())
#    print(stats.pearsonr(peaks[:,1], (int_temp/int_temp.sum())))
    return(stats.pearsonr(peaks[:,1], (int_temp/int_temp.sum()))[0], peaks)


