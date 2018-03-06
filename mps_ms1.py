import os
import pyopenms as oms
import numpy as np
from ProteoFileReader import MGF_Reader


def add_relaxation_mgf(mgf, mps, outfile, create_comparison=False):
    # mgf_file, outdir, differences = argls[0], argls[1], argls[2]
    mass_diff = 1.00335483
    # if '.mgf' in args['file']:
    filename = os.path.split(mgf)[1]
    # read mgf
    spectra = MGF_Reader()
    spectra.load(mgf)
    out_writer = open(outfile, "w")
    for spectrum in spectra:
        # calculate mass (neglect proton bec. later on difference used)
        scan = (spectrum.getTitle().split('.')[-2])
        # try:
        if not create_comparison:
            differences = mps[mps[:, 0] == scan, 1:]
            differences = [-i for i in range(len(differences)) if differences[i] == 1]
        else:
            differences = [0, -1, -2, -3, -4]
        mass = spectrum.getPrecursorMass() * spectrum.charge
        spectra_add_mip = [str((mass + x * mass_diff) / spectrum.charge) for x in differences if x != 0]
        if 0 in differences:
            prec_mz = spectrum.getPrecursorMass()
        else:
            prec_mz = spectra_add_mip[0]
            spectra_add_mip = spectra_add_mip[1:]
        # except KeyError:
        #     differences = [-2, -1, 0]

        stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={}
CHARGE={}+
RTINSECONDS={}
ADDITIONALMZ={}
{}
END IONS     """.format(spectrum.getTitle(), prec_mz, int(spectrum.charge), spectrum.getRT(),
                        ';'.join(spectra_add_mip),
                "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks]))
        out_writer.write(stavrox_mgf)


def read_mzml(infile):
    # init variables
    mzml_file = oms.MzMLFile()
    exp = oms.MSExperiment()

    # load spectra into exp
    mzml_file.load(infile, exp)
    return (exp)


def get_error(mz1, mz2, charge=None, ppm=True):
    if ppm:
        return abs(mz1 - mz2) / mz2 * 1e6
    else:
        if charge is not None:
            return (mz1 - mz2) * charge


def ms1_peaks(exp, tolerance=6, mps_range=[-1, -2, -3, -4]):
    # loop through spectra and make count if precursor in MS2 spectrum
    peaks_found = []
    n_ms2 = 0
    for i, spectrum in enumerate(exp):

        if spectrum.getMSLevel() == 1:
            continue
        n_ms2 += 1
        # if iidone % 5000 == 0:
        #     print "{}/{} done..".format(iidone, nspectra)

        # if iidone not in matched_scans:
        #     continue
        precursor = spectrum.getPrecursors()[0]
        prec_mz, prec_charge, prec_int = precursor.getMZ(), precursor.getCharge(), precursor.getIntensity()
        ms1_prev = 0
        for j in range(i, 0, -1):
            tmp_ms1 = exp[j]
            if tmp_ms1.getMSLevel() == 1:
                ms1_prev += 1
                res = tmp_ms1[tmp_ms1.findNearest(prec_mz)]
                if abs(res.getMZ() - prec_mz) / prec_mz <= tolerance:
                    theo_mip = np.array([prec_mz + (mip_i * 1.00335483) / prec_charge for mip_i in mps_range])
                    mip_nearest = [tmp_ms1.findNearest(x) for x in theo_mip]
                    error = np.array([get_error(tmp_ms1[mea].getMZ(), expi, ppm=True) for expi, mea in
                                      zip(theo_mip, mip_nearest)])
                    range_found = [True if x <= tolerance else False for x in error]
                    if sum(range_found) == 0:
                        # TODO: try if sensible to not mps search these
                        peaks_found.append(
                            [i + 1, True] + [True] * len(mps_range)
                            # [i + 1, True] + [False] * len(mps_range)
                        )
                        break
                    else:
                        # check for continous peaks except -1 peak
                        # TODO: allow gap?
                        found = False
                        lightest_peak = 1
                        for i_mip in range(len(range_found), 1, -1):
                            if range_found[i_mip - 1] & (i_mip > lightest_peak):
                                lightest_peak = i_mip
                            if sum(range_found[:i_mip]) == len(range_found[:i_mip]):
                                # takes lightest 2 continous + existing lighter peaks, excludes heaviar
                                sel = [False] * (i_mip - 1) + range_found[i_mip - 1:] # 2
                                peaks_found.append([i + 1, False] + sel)
                                found = True
                                break
                        # if no continous found take all
                        if not found:
                            peaks_found.append(
                                # [i + 1, True] + [True] * len(mps_range)
                                [i + 1, True] + [True] * lightest_peak + [False] * (len(range_found) - lightest_peak)
                            )
                            break
                        else:
                            found = False
                            break
                else:
                    if ms1_prev == 3:
                        break
                    continue

    return np.array(peaks_found)


if __name__ == '__main__':
    isotope_diff = 1.00335483
    # mzml_dir = 'D:/user/Swantje/data/PC/mzML/'
    mzml_dir = 'D:/user/Swantje/data/Chaetomium/frac3_6/mzML/'
    # mgf_filtered_dir = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/lars_PC_4frag_BS3_Lumos/All_prepro_peakfiles/mscon_PF_20_100_0/'
    mgf_filtered_dir = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/fr3_6/All_prepro_peakfiles/mscon_PF_20_100_0/'
    setting_name = 'ms1_sel4_4'
    # mgf_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/lars_PC_4frag_BS3_Lumos/All_prepro_peakfiles/' + setting_name
    # info_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/lars_PC_4frag_BS3_Lumos/relaxation_tbls/' + setting_name
    mgf_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/fr3_6/All_prepro_peakfiles/' + setting_name
    info_out = 'D:/user/Swantje/projects/pipeline_prepro_xi_fdr/chaetomium/fr3_6/relaxation_tbls/' + setting_name

    if not os.path.exists(mgf_out):
        os.makedirs(mgf_out)
    if not os.path.exists(info_out):
        os.makedirs(info_out)

    for mzml_file in [x for x in os.listdir(mzml_dir) if '.mzML' in x]:
        exp = read_mzml(mzml_dir + mzml_file)
        exp_id = mzml_file[:10]

        mps_df = ms1_peaks(exp)
        np.savetxt(info_out + '/%s_%s.csv' % (setting_name, mzml_file), mps_df, delimiter=',')
        corresponding_mgf = [x for x in os.listdir(mgf_filtered_dir) if exp_id in x][0]
        add_relaxation_mgf(mgf=mgf_filtered_dir + corresponding_mgf, mps=mps_df,
                           outfile=mgf_out + '/%s_' % setting_name + corresponding_mgf)

    # for mgf_file in os.listdir(mgf_filtered_dir):
    #     add_relaxation_mgf(mgf=mgf_filtered_dir + mgf_file, mps=[], create_comparison=True,
    #                        outfile='D:/user/Swantje/projects/pipeline_prepro_xi_fdr/lars_PC_4frag_BS3_Lumos/All_prepro_peakfiles/' + 'decoy_4_' + mgf_file)
