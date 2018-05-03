#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:04:31 2017

@author: sgiese
"""

import Deisotoper as de
import glob
import os




def get_files(basepath):
    """
    Returns the list of files that was found in the basepath directory.
    """
    files = []
    for filename in glob.iglob(basepath, recursive=True):
        files.append(filename)
    #ensure reproduceable order
    files = sorted(files)
    return(files)


#import sys
#sys.exit()

def process_files(basepath, outpath, n_jobs=-1):
    """
    Worker function to process all files.
    """
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    files = get_files(basepath)
    
    #init logger
    logger = open(basepath.split("/")[-2] + "_log2.txt", "w")
    
    print ("Files to go: {}".format(len(files)))
    #count successful processing files
    counter = 0
    for infile in files:
        basename = os.path.basename(infile)
        outfile = outpath + "/" + basename + "_deisotoped.mgf"
    
        print ("Analyzing file {} / {}".format(counter, len(files)))
        print(basename)
        
        #file does not exist?
        if not os.path.exists(outfile):
            dt = de.Deisotoper()
            
            print("Deisotoping spectra.")
            dt_spectra = dt.deisotope_spectra(infile, n_jobs=n_jobs)
                        
            print("Writing to file.")
            with open(outfile, 'w') as fobj:
                fobj.write("\r\n".join([spectrum.to_mgf() for spectrum in dt_spectra]))
            
            
            logger.write("{} Done.".format(outfile))
        
        else:
            pass
        counter += 1
        
    logger.close()    
    print ("Done.")
    print ("Files: {}".format(len(files)))
    print ("Processed: {}".format(counter))

# =============================================================================
# start stuff
# =============================================================================
basepath = "/home/sgiese/data/SwantjeTest/test/chaetomium_new/process/*/*/*"
basepath = "/home/sgiese/data/SwantjeTest/test/process/HSA/*.mgf"
basepath = "/home/sgiese/data/SwantjeTest/test/process/Chaetomium/*.mgf"

#set the input directory and the output directory
basepath = "/home/sgiese/data/Myco_deiso/*.mgf"
outpath = "/home/sgiese/data/Myco_deiso/processed/"


if __name__ == "__main__":
    process_files(basepath, outpath, n_jobs=-1)


# =============================================================================
# debug
# =============================================================================
#import ProteoFileReader as PFR
#n_jobs = 1
#infile = "Error_MGF_7" 
#dt = de.Deisotoper()
#
#MGF_file = PFR.MGF_Reader()
#MGF_file.load(infile)
#for ii, spectrum in enumerate(MGF_file):
#    break
#
#
#mz = spectrum.getMZ()
#intensity = spectrum.getIntensities()
#G = dt.spec2graph(mz, intensity, 20, 1, 7)
#cluster_ar = dt.extract_isotope_cluster(G, True)
#cluster_res = dt.resolve_ambiguous(cluster_ar, mz, intensity, 0.6, 0.25, 0.3, True)
#cluster_df = dt.assemble_spectrum(cluster_res, mz, intensity, spectrum.getTitle())
#
#
#dt_spectra = dt.deisotope_spectra(infile, n_jobs=1)
