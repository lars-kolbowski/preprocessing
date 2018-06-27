#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 16:04:31 2017

@author: sgiese
"""

import Deisotoper as de
import glob
import os
import numpy as np
import pandas as pd


def get_files(basepath, filetype):
    """
    Returns the list of files that was found in the basepath directory.
    """
    files = []
    for filename in glob.iglob(basepath+"/*{}".format(filetype), recursive=True):
        files.append(filename)
    #ensure reproduceable order
    files = sorted(files)
    return(files)

def to_mgf(spectrum):
        mgf_str="""
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (spectrum.title, spectrum.RT, spectrum.pepmz, spectrum.pepint, 
        spectrum.charge, "\r\n".join(["%s %s" % (i, j) for i,j in zip(np.ravel(spectrum.peaks[:,0]),
                                      np.ravel(spectrum.peaks[:,1]))]))
        return(mgf_str)
        
#import sys
#sys.exit()
        
def deconvolute(spectrum):
    """
    Deconvolutes a spectrum, requires charge annotations!
    """
    pass
    
def process_files(basepath, outpath, filetype="mzML", return_type="df",
                  n_jobs=-1):
    """
    Worker function to process all files.
    
    n_jobs=-1
    filetype="mzML"
    return_type="df"
    """
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    files = get_files(basepath, filetype=filetype)
    
    #init logger
    logger = open(basepath.split("/")[-2] + "_log2.txt", "w")
    
    print ("Files to go: {}".format(len(files)))
    #count successful processing files
    counter = 0
    if filetype == "mzML":
        meta_df = mzml_to_meta(basepath, filetype)
        meta_df.to_csv(outpath +  "scan_summary_from_mzml.csv")
    
    for infile in files:
        basename = os.path.basename(infile)
        outfile = outpath + "/" + basename + "_deisotoped.mgf"
    
        print ("Analyzing file {} / {}".format(counter, len(files)))
        print(basename)
        
        #file does not exist?
        if not os.path.exists(outfile):
            dt = de.Deisotoper()
            
            print("Deisotoping spectra.")
            
            
            if return_type == "df":
                df_spectra = dt.deisotope_spectra(infile, n_jobs=n_jobs, 
                                                  in_type=filetype,
                                                  return_type=return_type)
                df_spectra["infile"] = basename
                df_spectra.to_csv(outfile+"_df.csv")
                
            else:
                dt_spectra = dt.deisotope_spectra(infile, n_jobs=n_jobs, 
                                                  in_type=filetype)
                
                print("Writing to file.")
                with open(outfile, 'w') as fobj:                      
                    #fobj.write("\r\n".join([spectrum.to_mgf() for spectrum in dt_spectra]))
                    fobj.write("\r\n".join([to_mgf(spectrum) for spectrum in dt_spectra]))

        else:
            pass
        counter += 1
        
    logger.close()    
    print ("Done.")
    print ("Files: {}".format(len(files)))
    print ("Processed: {}".format(counter))



def mzml_to_meta(basepath, filetype):
    """
    """
    import re
    import pyopenms as oms
    def get_scan(scanheader):
        try:
            return(re.search("scan=(\d+)", scanheader).groups()[0])
        except:
            return(re.search("index=(\d+)", scanheader).groups()[0])
        
    files = get_files(basepath, filetype=filetype)
    meta_dic = {"file":[],
                "scan":[],
                "header":[],
                "pmz":[],
                "pRT":[],
                "pint":[],
                "pz":[]}
    
    for file in files:
        mzml_file = oms.MzMLFile()
        exp = oms.MSExperiment()
        mzml_file.load(file, exp)
        for spectrum in exp:
            if spectrum.getMSLevel() == 2:
                meta_dic["file"].append(file)
                # print (spectrum.getNativeID())
                meta_dic["header"].append(spectrum.getNativeID())
                meta_dic["pmz"].append(spectrum.getPrecursors()[0].getMZ())
                meta_dic["pint"].append(spectrum.getPrecursors()[0].getIntensity())
                meta_dic["pz"].append(spectrum.getPrecursors()[0].getCharge())
                meta_dic["pRT"].append(spectrum.getRT())
                meta_dic["scan"].append(get_scan(str(spectrum.getNativeID())))
    mzml_meta_df = pd.DataFrame(meta_dic)
    return(mzml_meta_df)
# =============================================================================
# start stuff
# =============================================================================
basepath = "/home/sgiese/data/SwantjeTest/test/chaetomium_new/process/*/*/*"
basepath = "/home/sgiese/data/SwantjeTest/test/process/HSA/*.mgf"
basepath = "/home/sgiese/data/SwantjeTest/test/process/Chaetomium/*.mgf"
basepath = "/data/sgiese/Projects/deisotoping"
basepath = "/home/sgiese/data/4swantje/"
outpath = "/home/sgiese/data/4swantje/"

#set the input directory and the output directory
#basepath = "/home/sgiese/data/Myco_deiso/*.mgf"
#outpath = "/home/sgiese/data/Myco_deiso/processed/"


if __name__ == "__main__":
    process_files(basepath, outpath, filetype="mzML", n_jobs=-1, return_type="df")


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
