import os
import re
import numpy as np
import subprocess
from multiprocessing import Pool
import sys

sys.path.append('D:\\software\\wxpython\\wx-3.0-msw')
from gooey import Gooey, GooeyParser

# import pyopenms as oms

import re
from pyteomics import mzml


def split_mzml(mzml_file):
    """
    function to split a mzML file into dict of MS2_Spectra objects (can be written to mgf format)
    by fragmentation method

    Parameters:
    -----------------------------------------
    mzml_file: str,
            path to mzML file

    Return: dict {fragMethod: list(MS2_spectrum)

    """

    mzml_reader = mzml.read(mzml_file)
    ordered_ms2_spectra = {
        "CID": [],
        "HCD": [],
        "ETD": [],
        "ETciD": [],
        "EThcD": []
    }

    for spectrum in mzml_reader:
        if spectrum['ms level'] == 2:
            try:
                groups = re.search("@([A-z]+)([0-9.]+)@?([A-z]+)?([0-9.]+)?", spectrum['scanList']['scan'][0]['filter string']).groups()
            except:
                print spectrum['scanList']['scan'][0]['filter string']
            title = os.path.split(mzml_file)[1].split('mzML')[0] + spectrum['id']
            rt = spectrum['scanList']['scan'][0]['scan start time'] * 60
            precursor = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            pre_mz = precursor['selected ion m/z']
            try:
                pre_int = precursor['peak intensity']
            except KeyError:
                pre_int = 0
            pre_z = precursor['charge state']
            peaks = zip(spectrum['m/z array'], spectrum['intensity array'])

            ms2class_spectrum = MS2_spectrum(title, rt, pre_mz, pre_int, pre_z, peaks)

            if "etd" in groups:
                if "cid" in groups:
                    ordered_ms2_spectra['ETciD'].append(ms2class_spectrum)
                elif "hcd" in groups:
                    ordered_ms2_spectra['EThcD'].append(ms2class_spectrum)
                else:
                    ordered_ms2_spectra['ETD'].append(ms2class_spectrum)
            elif "cid" in groups:
                ordered_ms2_spectra['CID'].append(ms2class_spectrum)
            elif "hcd" in groups:
                ordered_ms2_spectra['HCD'].append(ms2class_spectrum)

    return {k: v for k, v in ordered_ms2_spectra.items() if len(v) > 0}


class MS2_spectrum():
    """
    Class container for MS2 spectra.
    We need the following input parameters:
    title, RT, pepmass, pepint, charge, peaks, peakcharge=[]

    Parameters:
    -----------------------------------------
    title: str,
            title of the spectrum
    RT: float,
        retention time of the precursor
    pepmass: float,
              mass of the precursor
    charge: int,
             charge of the precursor
    peaks: [(float, float)],
           mass intensity
    peakcharge: arr,
                charge array for the peaks

    """
    def __init__(self, title, RT, pepmass, pepint, charge, peaks, peakcharge=[]):
        self.title = title
        self.RT = RT
        self.pepmass = pepmass
        self.pepint = pepint
        self.charge = charge
        self.peaks = peaks
        self.peakcharge = peakcharge

    def getPrecursorMass(self):
        """
        Returns the precursor mass
        """
        return(self.pepmass)

    def getPrecursorIntensity(self):
        """
        Returns the precursor mass
        """
        return(self.pepint)

    def getRT(self):
        """
        Returns the precursor mass
        """
        return(self.RT)

    def getTitle(self):
        """
        Returns the precursor mass
        """
        return(self.title)

    def getPeaks(self):
        """
        Returns the precursor mass
        """
        return(self.peaks)

    def getMasses(self):
        """
        TODO:
        Returns the precursor mass
        """
        return(self.peaks[:,0])

    def getIntensities(self):
        """
        TODO:
        Returns the precursor mass
        """
        return(self.peaks[:,1])

    def getUnchargedMass(self):
        """
        Computs the uncharged mass of a fragment:
        uncharged_mass = (mz * z ) - z
        """
        return( (self.pepmass * self.charge) -  self.charge)

    def printf(self):
        print ("Title, RT, PEPMASS, PEPINT, CHARGE")
        print (self.title, self.RT, self.pepmass, self.pepint, self.charge)

    def to_mgf(self):
        # need dummy values in case no peak charges are in the data
        if len(self.peakcharge) == 0:
            self.peakcharge = [""]*self.peaks.shape[0]
        mgf_str="""
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (self.title, self.RT, self.pepmass, self.pepint, self.charge, "\r\n".join(["%s %s %s" % (i[0], i[1], j, ) for i,j in zip(self.peaks, self.peakcharge)]))
        return(mgf_str)

#==============================================================================
# File Reader
#==============================================================================
class MGF_Reader():
    """A MGF_Reader is associated with a FASTA file or an open connection
    to a file-like object with content in FASTA format.
    It can generate an iterator over the sequences.

    Usage:
    --------------------------
    >>reader = MGF_Reader() \r\n
    >>reader.load(infile) \r\n
    >>#do something \r\n
    >>reader.store(outfile, outspectra) \r\n
    """
    def load(self, infile, getpeakcharge=False):
        """
        Function to set the input file for the MGF file.

        Parameters:
        -----------------------------
        infile: str,
                file location



        """
        self.infile = infile
        self.peakcharge = getpeakcharge

    def __iter__(self):
        mgf_file = open(self.infile, "r")
        found_ions = False
        for line in mgf_file:
            if len(line.strip()) == 0:
                continue
            if line.startswith("BEGIN IONS"):
               # init arrays for peaks
                found_ions = True
                mass = []
                intensity = []
                peakcharge = []
            elif line.startswith("TITLE"):
                title = re.search("TITLE=(.*)", line).groups()[0]

            elif line.startswith("RTINSECONDS"):
                RT = float(re.search("RTINSECONDS=(.*)", line).groups()[0])

            elif line.startswith("PEPMASS"):
                precursor = re.search("PEPMASS=(.*)", line).groups()[0].split()
                pep_mass = float(precursor[0])
                try:
                    pep_int = float(precursor[1])
                except:
                    pep_int = -1.0

            elif line.startswith("CHARGE"):
                charge = float(re.search("CHARGE=(\d)", line).groups()[0])

            elif "=" in line:
                print ("unhandled paramter: %s" % (line))

            elif line.startswith("END IONS"):
                ms = MS2_spectrum(title, RT, pep_mass, pep_int, charge, np.array(zip(mass, intensity)), peakcharge)
                yield ms
            else:
                if found_ions is True:
                    peak = line.split()
                    mass.append(float(peak[0]))
                    intensity.append(float(peak[1]))
                    if self.peakcharge:
                        if len(peak) >2:
                            peakcharge.append(peak[2])


def add_relaxation_mgf(args):
    # mgf_file, outdir, differences = argls[0], argls[1], argls[2]
    mass_diff = 1.00335483
    if '.mgf' in args['file']:
        filename = os.path.split(args['file'])[1].replace('mzML', 'mgf')
        # read mgf
        spectra = MGF_Reader()
        spectra.load(args['mgffile'])
        out_writer = open(os.path.join(args['outdir'], filename), "w")
        for spectrum in spectra:
            # calculate mass (neglect proton bec. later on difference used)
            mass = spectrum.getPrecursorMass() * spectrum.charge
            spectra_add_mip = [str((mass + x * mass_diff) / spectrum.charge) for x in args['differences']]
            stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={}
CHARGE={}+
RTINSECONDS={}
ADDITIONALMZ={}
{}
END IONS     """.format(spectrum.getTitle(),
                                spectrum.getPrecursorMass(),
                                int(spectrum.charge), spectrum.getRT(),
                                ';'.join(spectra_add_mip),
                                "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks]))
            out_writer.write(stavrox_mgf)

    else:
        # read mzml, sort spectra into acquisition modes
        # output as MS2spectrum object
        splitted_spectra = split_mzml(args['file'])
        orig_name = os.path.split(args['file'])[1].replace('mzML', 'mgf')

        for acq in splitted_spectra:
            filename = acq + '_' + orig_name
            # add line to info
            # iterate over spectra and write MGF to new file
            out_writer = open(os.path.join(args['outdir'], filename), "w")
            for spectrum in splitted_spectra[acq]:
                # calculate mass (neglect proton bec. later on difference used)
                mass = spectrum.getPrecursorMass() * spectrum.charge
                spectra_add_mip = [str((mass + x * mass_diff) / spectrum.charge) for x in args['differences']]
                stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={}
CHARGE={}+
RTINSECONDS={}
ADDITIONALMZ={}
{}
END IONS     """.format(spectrum.getTitle(),
                                        spectrum.getPrecursorMass(),
                                        int(spectrum.charge), spectrum.getRT(),
                                        ';'.join(spectra_add_mip),
                                        "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks]))
                out_writer.write(stavrox_mgf)


def mscon_cmd(rawdir, outdir, settings, mgf):
    # create list with all files in directory
    files_list = os.listdir(rawdir)

    # subselect all raw files and join to full path
    raw_files_list = [x for x in files_list if '.' in x]
    raw_files_list = [os.path.join(rawdir, x) for x in raw_files_list
                      if not x in os.listdir(outdir)] #.replace('raw', 'mzML')

    filter_formatted = []
    for i in range(len(settings)):
        filter_formatted.append('--filter')
        filter_formatted.append(settings[i])

    if mgf:
        cmd_list = [[rawfile, '--mgf', '-o', outdir] + filter_formatted for rawfile in raw_files_list]
    else:
        cmd_list = [[rawfile, '-o', outdir] + filter_formatted for rawfile in raw_files_list]

    return cmd_list


def mscon_wrap(cmds):
    mscon = ['C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.9740\\msconvert.exe']
    convert = subprocess.Popen(mscon + cmds)
    convert.communicate()


def main():
    raw_dir = 'D:/user/Swantje/Lars' #'//130.149.167.198/rappsilbergroup/users/lswantje/Lars'
    split_acq = True
    nthr = 2

    # MS2 peak filtering - can be left like this
    mscon_settings = ["MS2Denoise 20 100 false"] #

    # monoisotopic shifts to add - can also be left like this (exact mass difference is set in function)
    # up to -3 Da led to 5% more PSMs for ribo data set - might be worth it, since it doesn't increase search time much
    da_range = [-1, -2, -3]

    filt_dir = os.path.join(raw_dir, 'preprocessed', 'mgf_filtered')
    rel_dir = os.path.join(raw_dir, 'preprocessed', 'relaxed')

    if not os.path.exists(filt_dir):
        os.makedirs(filt_dir)

    if split_acq:
        conv_cmds = mscon_cmd(rawdir=raw_dir, outdir=filt_dir, settings=mscon_settings, mgf=False)
    else:
        conv_cmds = mscon_cmd(rawdir=raw_dir, outdir=filt_dir, settings=mscon_settings, mgf=True)

    pool = Pool(processes=nthr)
    pool.map(mscon_wrap, conv_cmds)
    pool.close()
    pool.join()
    print 'msconvert finished'

    if not os.path.exists(rel_dir):
        os.makedirs(rel_dir)

    relaxation_arg = [{'file': os.path.join(filt_dir, f), 'outdir': rel_dir, 'differences': da_range} for f in os.listdir(filt_dir)]
    pool = Pool(processes=nthr)
    pool.map(add_relaxation_mgf, relaxation_arg)
    pool.close()
    pool.join()
    print 'relaxation added'


if __name__ == '__main__':
    main()
