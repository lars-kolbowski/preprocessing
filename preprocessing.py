import os
import numpy as np
import subprocess
from multiprocessing import Pool
import sys
import re
import getopt
from pyteomics import mzml
from functools import partial


def read_cmdline():
    try:
        opts, args = getopt.getopt(sys.argv[1:], '', ['input=', 'config=', 'outpath='])
    except getopt.GetoptError:
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file>')
        sys.exit()

    for opt, arg in opts:
        if opt == '--input':
            input_arg = arg
        elif opt == '--outpath':
            outdir = arg
        elif opt == '--config':
            config = arg

    if 'input_arg' not in locals() or 'config' not in locals():
        print('preprocessing.py --input <folder or single file to process> '
              '--outpath <directory for output, default is separate folder in input directory> '
              '--config <path to config file>')
        sys.exit()
    # if no outdir defined use location of input
    if 'outdir' not in locals() and os.path.isdir(input_arg):
        outdir = os.path.join(input_arg, 'processed')
    elif 'outdir' not in locals() and not os.path.isdir(input_arg):
        outdir = os.path.join(os.path.split(input_arg)[0], 'processed')

    return input_arg, outdir, config


def split_mzml(mzml_file, detector="all"):
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
        "EThcD": [],
        "unknown": []
    }

    n = 0
    for spectrum in mzml_reader:
        if spectrum['ms level'] == 2:
            n += 1
            filter_str = spectrum['scanList']['scan'][0]['filter string']
            try:
                detector_str = re.search("^(FT|IT)", filter_str).groups()[0]
                frag_groups = re.findall("@([A-z]+)([0-9.]+)", filter_str)
            except AttributeError:
                raise StandardError("filter string parse error: %s" % filter_str)

            if not detector == "all":
                if not detector == detector_str:
                    continue

            title = os.path.split(mzml_file)[1].split('.mzML')[0] + " " + spectrum['id']
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

            frag_methods = [f[0] for f in frag_groups]

            if "etd" in frag_methods:
                if "cid" in frag_methods:
                    ordered_ms2_spectra['ETciD'].append(ms2class_spectrum)
                elif "hcd" in frag_methods:
                    ordered_ms2_spectra['EThcD'].append(ms2class_spectrum)
                else:
                    ordered_ms2_spectra['ETD'].append(ms2class_spectrum)
            elif "cid" in frag_methods:
                ordered_ms2_spectra['CID'].append(ms2class_spectrum)
            elif "hcd" in frag_methods:
                ordered_ms2_spectra['HCD'].append(ms2class_spectrum)
            else:
                ordered_ms2_spectra['unknown'].append(ms2class_spectrum)
    if len(ordered_ms2_spectra['unknown']) > 0:
        raise Warning("The fragmentation method of %i spectra could not be identified" % len(ordered_ms2_spectra['unkown']))

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
        return (self.pepmass)

    def getPrecursorIntensity(self):
        """
        Returns the precursor mass
        """
        return (self.pepint)

    def getRT(self):
        """
        Returns the precursor mass
        """
        return (self.RT)

    def getTitle(self):
        """
        Returns the precursor mass
        """
        return (self.title)

    def getPeaks(self):
        """
        Returns the precursor mass
        """
        return (self.peaks)

    def getMasses(self):
        """
        TODO:
        Returns the precursor mass
        """
        return (self.peaks[:, 0])

    def getIntensities(self):
        """
        TODO:
        Returns the precursor mass
        """
        return (self.peaks[:, 1])

    def getUnchargedMass(self):
        """
        Computs the uncharged mass of a fragment:
        uncharged_mass = (mz * z ) - z
        """
        return ((self.pepmass * self.charge) - self.charge)

    def printf(self):
        print ("Title, RT, PEPMASS, PEPINT, CHARGE")
        print (self.title, self.RT, self.pepmass, self.pepint, self.charge)

    def to_mgf(self):
        # need dummy values in case no peak charges are in the data
        if len(self.peakcharge) == 0:
            self.peakcharge = [""] * self.peaks.shape[0]
        mgf_str = """
BEGIN IONS
TITLE=%s
RTINSECONDS=%s
PEPMASS=%s %s
CHARGE=%s
%s
END IONS
        """ % (self.title, self.RT, self.pepmass, self.pepint, self.charge,
               "\r\n".join(["%s %s %s" % (i[0], i[1], j,) for i, j in zip(self.peaks, self.peakcharge)]))
        return (mgf_str)


# ==============================================================================
# File Reader
# ==============================================================================
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
                        if len(peak) > 2:
                            peakcharge.append(peak[2])


def mscon_cmd(filepath, outdir, settings, mgf):
    filename = os.path.split(filepath)[1]

    if (filename[:filename.rfind('.')] in [x[:x.rfind('.')] for x in os.listdir(outdir)]) or (os.path.isdir(filepath)):
        return []

    filter_formatted = []
    for i in range(len(settings)):
        filter_formatted.append('--filter')
        filter_formatted.append(settings[i])

    if mgf:
        cmd_list = [filepath, '--mgf', '-o', outdir] + filter_formatted
    else:
        cmd_list = [filepath, '-o', outdir] + filter_formatted
    # TODO modify title via msconvert
    return cmd_list


def write_mgf(spectra, outfile):
    out_writer = open(os.path.join(outfile), "w")
    for spectrum in spectra:
        stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={}
CHARGE={}+
RTINSECONDS={}
{}
END IONS     """.format(spectrum.getTitle(),
                            spectrum.getPrecursorMass(),
                            int(spectrum.charge), spectrum.getRT(),
                            "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks if i[1] > 0]))
        out_writer.write(stavrox_mgf)


def process_file(filepath, outdir, mscon_settings, split_acq, detector_filter, mscon_exe):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    conv_cmds = mscon_cmd(filepath=filepath, outdir=outdir, settings=mscon_settings, mgf=not split_acq)

    if len(conv_cmds) > 0:
        msconvert = subprocess.Popen([mscon_exe] + conv_cmds)
        msconvert.communicate()

    if split_acq:
        filename = os.path.split(filepath)[1]
        mzml_file = os.path.join(outdir, filename[:filename.rfind('.')]+'.mzML')
        splitted_spectra = split_mzml(mzml_file, detector_filter)

        for acq in splitted_spectra:
            write_mgf(spectra=splitted_spectra[acq],
                      outfile=os.path.join(outdir, acq + '_' + os.path.split(mzml_file)[1]))


if __name__ == '__main__':
    # read cmdline arguments / get deafult values
    input_arg, outdir, config_path = read_cmdline()
    execfile(config_path)

    # get files in directory
    if os.path.isdir(input_arg):
        full_paths = [os.path.join(input_arg, rawfile) for rawfile in os.listdir(input_arg)]
        full_paths = [x for x in full_paths if not os.path.isdir(x)]
    # if single file given reformat to list
    # TODO allow txt file with filepath
    else:
        full_paths = [input_arg]

    pool = Pool(processes=nthr)
    pool.map(partial(process_file, outdir=outdir, mscon_settings=mscon_settings, split_acq=split_acq,
                     detector_filter=detector_filter, mscon_exe=msconvert_exe), full_paths)
    pool.close()
    pool.join()
