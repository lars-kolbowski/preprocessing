import os
import numpy as np
import subprocess
from multiprocessing import Pool
import sys
import re
import getopt
from pyteomics import mzml
from functools import partial
import ProteoFileReader
import mass_recal
import zipfile


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

            ms2class_spectrum = ProteoFileReader.MS2_spectrum(title, rt, pre_mz, pre_int, pre_z, peaks)

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


def mscon_cmd(filepath, outdir, settings, mgf):
    filename = os.path.split(filepath)[1]

    if (filename[:filename.rfind('.')] in [x[:x.rfind('.')] for x in os.listdir(outdir)]) or (os.path.isdir(filepath)):
        print('File ' + filename + ' already existing in output directory, not done again.')
        return []

    filter_formatted = []
    for i in range(len(settings)):
        filter_formatted.append('--filter')
        filter_formatted.append(settings[i])

    if mgf:
        cmd_list = [filepath, '--mgf', '-o', outdir] + filter_formatted
    else:
        cmd_list = [filepath, '-o', outdir] + filter_formatted
    return cmd_list


def write_mgf(spectra, outfile):
    out_writer = open(os.path.join(outfile), "w")
    for spectrum in spectra:
        title = re.match('(B|E)[0-9]{6}_[0-9]{2}.+?( )', spectrum.getTitle()).group(0)[:-1]
        scan = re.search('scan=[0-9]*', spectrum.getTitle()).group(0)[5:]
        stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
{}
END IONS     """.format('.'.join([title, scan, scan, str(int(spectrum.charge))]),
                        spectrum.getPrecursorMass(),
                        spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                        int(spectrum.charge), spectrum.getRT(),
                        "\n".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks if i[1] > 0]))
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
                      outfile=os.path.join(outdir, acq + '_' + filename[:filename.rfind('.')]+'.mgf'))


if __name__ == '__main__':
    # read cmdline arguments / get deafult values
    input_arg, outdir, config_path = read_cmdline()
    try:
        execfile(config_path)
    except NameError:
        exec(open(config_path).read())

    # get files in directory
    if os.path.isdir(input_arg):
        full_paths = [os.path.join(input_arg, rawfile) for rawfile in os.listdir(input_arg)]
        full_paths = [x for x in full_paths if not os.path.isdir(x)]
    # if single file given reformat to list
    # TODO allow txt file with filepath
    else:
        full_paths = [input_arg]

    print("""file input:
{}
""".format('\n'.join(full_paths)))

    # start msconvert for conversion and peak filtering
    pool = Pool(processes=nthr)
    pool.map(partial(process_file, outdir=outdir, mscon_settings=mscon_settings, split_acq=split_acq,
                     detector_filter=detector_filter, mscon_exe=msconvert_exe), full_paths)
    pool.close()
    pool.join()

    recal_in = [os.path.join(outdir, x) for x in os.listdir(outdir) if '.mgf' in x]
    if mass_recalibration:
        # pool = Pool(processes=nthr)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        output = zipfile.ZipFile(outdir + '/recalibrated_files.zip', 'w', zipfile.ZIP_DEFLATED)
        # TODO change to parallel with manual input of error
        for inputfile in recal_in:
            mass_recal.main(fasta=database, xi_cnf=xi_recal_config, outpath=outdir,
                            xi_dir=xi_offline, mgf=inputfile, threads=str(nthr))
            output.write(os.path.join(outdir, 'recal_' + os.path.split(inputfile)[1]),
                         arcname='recal_' + os.path.split(inputfile)[1])
            # pool.map(partial(mass_recal.main, fasta=database, xi_cnf=xi_recal_config, outpath=outdir + '/recal',
        #                  xi_jar=xi_offline), recal_in)
        # pool.close()
        # pool.join()
