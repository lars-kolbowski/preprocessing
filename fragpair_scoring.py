import numpy as np
from joblib import Parallel, delayed
import pandas as pd
import ProteoFileReader
import itertools
import pickle
import re

proton = 1.00727647


def get_mass(mz, z, z_max, i):
    """
    Returns mass from m/z and charge. If unknown charge (=0) then all possible masses up to precursor z are returned.
    """
    if float(z) == 0:
        z = range(1, int(z_max) + 1)
    else:
        z = [z]
    return [((float(z_) * float(mz)) - float(z_)*proton, i) for z_ in z]


def uncharged_masses(ms2_spec, charge_max, df=False):
    """
    Return peak masses for a spectrum.
    """
    if not df:
        try:
            ms2_mz = ms2_spec.peaks[:,0]
            ms2_int = ms2_spec.peaks[:,1]
        except IndexError:
            return
        ms2_charge = ms2_spec.peakcharge
        if len(ms2_charge) == 0:
            ms2_charge = [0] * len(ms2_mz)
    else:
        ms2_spec_mod = ms2_spec.loc[ms2_spec.groupby('clusterIds_proc')['mz'].idxmin()].copy()
        ms2_mz = ms2_spec_mod.mz
        ms2_int = ms2_spec_mod.intensity
        ms2_charge = ms2_spec_mod.charge_pre_match
    # currently not used, was used to only combine higher intense peaks
    try:
        int_cutoff = sorted(ms2_int)[0]
    except IndexError:
        int_cutoff = 0
    ms2_int = ms2_int/max(ms2_int)
    mass_list = [get_mass(mzi, zi, charge_max, ii) for mzi, zi, ii in zip(ms2_mz, ms2_charge, ms2_int)] #ii>=int_cutoff
    mass_list_low = [] # (get_mass(mzi, zi, charge_max), ii)for mzi, zi, ii in zip(ms2_mz, ms2_charge, ms2_int) if ii<int_cutoff

    return mass_list, mass_list_low


def sum_frags(peak_masses, peak_masses_low):
    """
    Sums up combinations of masses and forms mean of intensities.
    """
    if peak_masses is None:
        return -1
    comb, comb_int = [], []

    for el in list(itertools.combinations(peak_masses, 2)):
        comb += list(itertools.product(*el))

    return np.array([(frag_comb[0][0] + frag_comb[1][0], np.mean([frag_comb[0][1], frag_comb[1][1]])) for frag_comb in comb])


def score_masses(possible_masses, paired_masses, ms2err):
    """
    For each possible (monoisotopic) mass a score is returned.
    """
    count = []
    mass, inten = zip(*paired_masses)
    mass = np.array(mass)
    inten = np.array(inten)
    for isotope_peak in possible_masses:
        resdiff = np.abs(mass - isotope_peak)
        res_err = resdiff / isotope_peak *10**6

        score = [inten[i]**2 for i in range(len(res_err)) if res_err[i] <= ms2err]
        count += [sum(score)]
    return count


def precursor_matches(spec, ms2err, df=False):
    """

    Returns scores for each possible (monoisotopic) precursor of a scan.

    """
    if not df:
        scan = spec.getTitle().split('.')[-2]
        prec_charge = spec.charge
        prec_mz = spec.pepmass
        prec_mass = (prec_charge * prec_mz) - prec_charge*proton
        possible_masses = [prec_mass + (mip_i * 1.00335483) for mip_i in [-4, -3, -2, -1, 0]]
    else:
        scan = spec[1].scan.unique()[0]
        id = re.match('(B|E)[0-9]{6}_[0-9]{2}', spec[1].run[spec[1]['PSMID'] == spec[0]].unique()[0]).group(0)
        scan = str(id) + '_' + str(scan)
        prec_charge = spec[1]['PrecursorZ'].unique()
        prec_mz = spec[1].npp_mz.unique() #PrecursorMZ.unique()
        prec_mass = (prec_charge * prec_mz) - prec_charge*proton
        possible_masses = [prec_mass + (mip_i * 1.00335483) for mip_i in [-4, -3, -2, -1, 0]]

    if df:
        peak_masses, peak_masses_low = uncharged_masses(spec[1], prec_charge, df=True)
    else:
        try:
            peak_masses, peak_masses_low = uncharged_masses(spec, prec_charge, df=False)
        except TypeError:
            return
    paired_masses = sum_frags(peak_masses, peak_masses_low)
    count_matched = score_masses(possible_masses, paired_masses, ms2err)

    return [scan] + count_matched


def test_monoisotopic_mass(results_store, ms2err, df=False):
    """

    Adds all the peaks up and checks if the monoisotopic mass is contained.

    """

    # corrected = [precursor_matches(spectrum, ms2err, df=False) for spectrum in results_store]
    if df:
        corrected = Parallel(n_jobs=20)(delayed(precursor_matches)(x, ms2err, df=True) for x in results_store.groupby('PSMID'))
    else:
        corrected = Parallel(n_jobs=20)(delayed(precursor_matches)(spectrum, ms2err) for spectrum in results_store)

    return corrected


def max_frag_count(x):
    a = x[x == x.max()]
    if len(a) > 1:
        return 1
    else:
        return int(a.index.tolist()[0])


if __name__ == '__main__':
    base_dir = 'D:/user/Swantje/projects/fragmass_sum/'
    mgf_deiso = base_dir + '/deisotoped/hits_mscon_PF_20_100_0_B160803_02._deisotoped.mgf'
    file_prefix = 'B160803_02_'

    prec_infos = pd.DataFrame.from_csv(base_dir + './prec_infos.csv')
    hits_mip = prec_infos['Xi_mip']
    ms2_tol = 10

    exp = ProteoFileReader.MGF_Reader()
    exp.load(mgf_deiso, getpeakcharge=True)

    b = test_monoisotopic_mass(exp, ms2_tol, df=False)
    # pickle.dump(b, open(base_dir + '/frag_counts_2Da_annotator_new.p', 'wb'))

    frag_count = pd.DataFrame(b)
    frag_count.columns = ['scan'] + ['-' + str(x) for x in range(4, -1, -1)]
    frag_count['scan'] = file_prefix + frag_count['scan']
    frag_count.index = frag_count.scan
    del frag_count['scan']

    frag_count['max_frag'] = frag_count.apply(max_frag_count, axis=1)
    frag_count = pd.concat([frag_count, hits_mip], axis=1, join='inner')

    frag_count.to_csv(base_dir + '/fragment_scores.csv')

