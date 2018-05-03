#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 10:43:33 2017

@author: sgiese
"""
import Deisotoper as de
import ProteoFileReader as PFR
import pandas as pd
import AveragineModel as AM
import numpy as np
import time


if __name__ == "__main__":
    infile = "data/mscon_PF_20_100_0_B160803_02.mgf"
    outfile = infile.replace(".mgf", "_deisotoped.mgf")
    dt = de.Deisotoper()
    print("Deisotoping spectra.")
    dt_spectra = dt.deisotope_spectra(infile, show_progress=True,
                                              n_jobs=-1)
    
    print("Writing to file.")
    with open(outfile, 'w') as fobj:
        fobj.write("\r\n".join([spectrum.to_mgf() for spectrum in dt_spectra]))
        
def test():
    #%%
    infile = "data/mscon_PF_20_100_0_B160803_02.mzML"
    exp = PFR.mzMLReader(infile)
    print (exp)
    
    infile = "data/test.mgf"
    infile = "data/mscon_PF_20_100_0_B160803_02.mgf"

    deisotoper = de.Deisotoper()
    deisotoped_spectra = deisotoper.deisotope_spectra(infile, show_progress=True,
                                                      n_jobs=-1)

    #%%
    spectrum = pd.read_csv("data/test.dta2d", sep="\t")
    mz = spectrum["MZ"].values
    intensity = spectrum["INT"].values
    deisotoper = de.Deisotoper()
    G = deisotoper.spec2graph(mz, intensity)
    cluster_ar = deisotoper.extract_isotope_cluster(G, verbose=False)
    cluster_dic_resolved = deisotoper.resolve_ambiguous(cluster_ar, mz,
                                                        intensity)
    cluster_df = de.deisotoper.assemble_spectrum(cluster_dic_resolved,
                                              mz, intensity, 1337)
    n_accepted = sum([i[0] for i in cluster_dic_resolved.values()])
    n_rejected = len(cluster_dic_resolved) - n_accepted
    

    #%%
    #plt.bar(mz, intensity)
    verbose = False
    deisotoper = de.Deisotoper()

    # =========================================================================
    #
    # =========================================================================
    #example 1 - pseudo overlapping - charge 2 but 1 also possible
    mz = np.array([0,100,250, 700, 700.5, 701, 701.5, 800])
    z_test = 2
    exp1 = AM.averagine_model(700*z_test, n_peaks=4,  only_intensity=True)*100
    intensity = np.array([10, 10, 10,
                          exp1[0],
                          exp1[1],
                          exp1[2],
                          exp1[3],
                          10])

    GT = {}
    GT["ratio"] = "1:0"
    GT["Case"] = "test 1 - only one true isotope cluster"
    GT["TestID"] = "1"
    G = deisotoper.spec2graph(mz, intensity)
    de.plot_graph(G)
    cluster_dic = deisotoper.extract_isotope_cluster(G, verbose)
    cluster_dic_resolved = deisotoper.resolve_ambiguous(cluster_dic, mz, intensity)

    #plt.bar(mz, intensity, width=0.2)
    #plt.xlim(699, 702)
    #%%
    # =========================================================================
    # example 2 - overlapping charge 1 and 2
    # =========================================================================
    verbose = True
    deisotoper = de.Deisotoper()
    
    a = time.time()
    mz = np.array([0,100,250, 300,500, 501, 501.5, 502, 503])
    exp1 = AM.averagine_model(500, n_peaks=4, only_intensity=True)*100  * 3 #+ #np.random.normal(0, 10, 4)
    exp2 = AM.averagine_model(501.5*2, n_peaks=3, only_intensity=True)*100 * 6 #+ np.random.normal(0, 10, 3)
    intensity = np.array([10,10,10,10,
                          exp1[0],
                          exp1[1]+exp2[0],
                                  exp2[1],
                          exp1[2]+exp2[2],
                          exp1[3]])

    GT = {}
    GT["ratio"] = "1:2"
    GT["Case"] = "test 2 - Overlapping different charge"
    GT["TestID"] = "2"
    G = deisotoper.spec2graph(mz, intensity)
    de.plot_graph(G)
    cluster_ar = deisotoper.extract_isotope_cluster(G, verbose)
    cluster_dic_resolved = deisotoper.resolve_ambiguous(cluster_ar, mz, intensity, verbose=verbose, GT=GT)
    print (cluster_ar)
    b = time.time()
    took = (b-a) / 60.
    print (took*30000)
