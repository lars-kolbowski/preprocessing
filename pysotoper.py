#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:48:37 2017

@author: sgiese
"""

import Deisotoper as de
import sys
import argparse
import textwrap


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
         prog='pysotoper',
         formatter_class=argparse.RawDescriptionHelpFormatter,
         description=textwrap.dedent('''\
         Command line tool for deisotoping of high resolution mass spectra.
         
         --------------------------------
         Pysotoper can be used to annotate fragment peaks in MGF spectra
         with their charge state. High-resolution data is necessary for this.
         
         Ambiguous isotope cluster are supported and overlapping isotope 
         distributions can be resolved.'''))
    
    
    parser.add_argument('infile', metavar='infile', type=str, nargs='+',
                       help='Location of the input file.')
        
    parser.add_argument('--ppm', dest='ppm_tolerance', default=20.0,
                        help='tolerance in ppm for mass errors in the isotope peaks  (default: 20.0)')
    
    # =========================================================================
    #     Param reading.
    # =========================================================================
    args = parser.parse_args()
    
    infile = args.infile[0]
    ppm_tolerance = args.ppm_tolerance    
    
    # =========================================================================
    #    Run script.
    # =========================================================================
    print(args)
    
    outfile = infile.replace(".mgf", "_deisotoped.mgf")
    dt = de.Deisotoper()
    print("Deisotoping spectra.")
    dt_spectra = dt.deisotope_spectra(infile, n_jobs=-1)
    
    print("Writing to file.")
    with open(outfile, 'w') as fobj:
        fobj.write("\r\n".join([spectrum.to_mgf() for spectrum in dt_spectra]))
