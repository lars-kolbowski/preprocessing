import os
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import ProteoFileReader


def xi_wrapper(arguments):
    xi = subprocess.Popen(arguments)
    xi.communicate()
    return


def run_xi_lin(peakfile, fasta, cnf, outpath, xipath, threads='5'):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    xi_cmds = ['java', '-cp', xipath, 'rappsilber.applications.Xi',
               '--fasta=' + fasta,
               '--xiconf=UseCPUs:' + threads,
               '--peaks=' + peakfile,
               '--config=' + cnf,
               '--output=' + outpath + '/xi_' + os.path.split(peakfile)[1].replace('.mgf', '.csv')]

    xi = subprocess.Popen(xi_cmds)
    xi.communicate()


def get_ppm_error(xi_df, outfile):
    xi_df = xi_df[(xi_df.decoy == 0) & (xi_df['match score'] > 7)]
    median_err = np.median(xi_df['Precoursor Error'])

    fig, ax = plt.subplots()
    sns.distplot(xi_df['Precoursor Error'], norm_hist=False, kde=False)
    ax.axvline(median_err)
    plt.savefig(outfile)
    plt.close()

    return median_err


def adjust_prec_mz(mgf_file, error, outpath):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    exp = ProteoFileReader.MGF_Reader()
    exp.load(mgf_file)
    outfile = os.path.join(outpath, 'adj_' + os.path.split(mgf_file)[1])

    out_writer = open(os.path.join(outfile), "w")
    for spectrum in exp:
        prec_mz_new = spectrum.getPrecursorMass()/(1-error/10.**6)
        stavrox_mgf = """
MASS=Monoisotopic
BEGIN IONS
TITLE={}
PEPMASS={} {}
CHARGE={}+
RTINSECONDS={}
{}
END IONS     """.format(spectrum.getTitle(),
                        prec_mz_new, spectrum.getPrecursorIntensity() if spectrum.getPrecursorIntensity() > 0 else 0,
                            int(spectrum.charge), spectrum.getRT(),
                            "\r".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks if i[1] > 0]))
        out_writer.write(stavrox_mgf)


def mass_recal(mgf, fasta, xi_cnf, outpath, xi_jar):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    filename = os.path.split(mgf)[1]
    # linear small search in Xi
    run_xi_lin(peakfile=mgf, fasta=fasta, cnf=xi_cnf, outpath=os.path.join(outpath), xipath=xi_jar)

    # evaluate results, get mean / median
    ms1_err = get_ppm_error(xi_df=pd.DataFrame.from_csv(os.path.join(outpath, 'xi_' + filename.replace('.mgf', '.csv'))),
                            outfile=os.path.join(outpath, 'MS1_err_' + filename + '.png'))

    # shift all old m/z by value
    adjust_prec_mz(mgf_file=mgf, error=ms1_err, outpath=os.path.join(outpath, 'adj_mscon_PF_20'))


if __name__ == '__main__':
    base_dir = '//130.149.167.198/rappsilbergroup/users/lswantje/prepro_mito/trypnew/processed'
    mgf_dir = base_dir
    fasta = '//130.149.167.198/rappsilbergroup/users/MitoProject/For Swantje/20170109_uniprot_mitoIDrun_FASTA.fasta'
    outpath = base_dir + '/error_shift'
    for in_file in os.listdir(mgf_dir):
        if '.mgf' in in_file:
            mass_recal(mgf=os.path.join(mgf_dir, in_file),
                       fasta=fasta,
                       xi_cnf='D:/user/Swantje/projects/pipeline_prepro_xi_fdr/resources/xi_linear.conf',
                       outpath=outpath,
                       xi_jar='D:/user/Swantje/XiSearch_1.6.731/*')
