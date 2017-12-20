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


def run_xi_lin(peakfile, fasta, cnf, outpath, xipath, threads='1'):
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    elif os.path.isfile( outpath + '/xi_' + os.path.split(peakfile)[1].replace('.mgf', '.csv')):
        return

    xi_cmds = ['java', '-cp', xipath + '/fastutil-8.1.0.jar;' + xipath + '/XiSearch.jar', 'rappsilber.applications.Xi',
               '--fasta=' + fasta,
               '--xiconf=UseCPUs:' + threads,
               '--peaks=' + peakfile,
               '--config=' + cnf,
               '--output=' + outpath + '/xi_' + os.path.split(peakfile)[1].replace('.mgf', '.csv')]

    xi = subprocess.Popen(xi_cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    xi.communicate()


def get_ppm_error(xi_df, outfile):
    xi_df = xi_df[(xi_df.decoy == 0) & (xi_df['match score'] > 7)]
    if len(xi_df) < 75:
        print os.path.split(outfile)[1] + ': not enough data to shift'
        err = raw_input('Enter error to correct by (0 for no correction):\n')
        try:
            err = float(err)
            if (err != 0):
                return err
        except ValueError:
            return None
    median_err = np.median(xi_df['Precoursor Error'])

    fig, ax = plt.subplots()
    sns.distplot(xi_df['Precoursor Error'], norm_hist=False, kde=False)
    ax.axvline(median_err)
    plt.savefig(outfile)
    plt.close()

    return median_err


def adjust_prec_mz(mgf_file, error, outpath):
    outfile = os.path.join(outpath, 'recal_' + os.path.split(mgf_file)[1])
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    elif os.path.isfile(outfile):
        return
    exp = ProteoFileReader.MGF_Reader()
    exp.load(mgf_file)

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
                            "\n".join(["%s %s" % (i[0], i[1]) for i in spectrum.peaks if i[1] > 0]))
        out_writer.write(stavrox_mgf)


def main(mgf, fasta, xi_cnf, outpath, xi_dir, threads):
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    filename = os.path.split(mgf)[1]
    # linear small search in Xi
    run_xi_lin(peakfile=mgf, fasta=fasta, cnf=xi_cnf, outpath=os.path.join(outpath), xipath=xi_dir, threads=threads)

    # evaluate results, get median ms1 error
    ms1_err = get_ppm_error(xi_df=pd.DataFrame.from_csv(os.path.join(outpath, 'xi_' + filename.replace('.mgf', '.csv'))),
                            outfile=os.path.join(outpath, 'MS1_err_' + filename + '.png'))

    error_file = open(outpath + '/ms1_err.csv', 'a')
    error_file.write(filename + ',' + str(ms1_err) + '\n')
    error_file.close()

    if ms1_err is not None: # shift all old m/z by value
        adjust_prec_mz(mgf_file=mgf, error=ms1_err, outpath=os.path.join(outpath))


if __name__ == '__main__':
    base_dir = '//130.149.167.198/rappsilbergroup/users/lswantje/prepro_mito/tryp/new/mgf'
    mgf_dir = base_dir
    fasta = '//130.149.167.198/rappsilbergroup/users/MitoProject/For Swantje/20170109_uniprot_mitoIDrun_FASTA.fasta'
    outpath = base_dir + '/error_shift'
    for in_file in os.listdir(mgf_dir):
        if '.mgf' in in_file:
            main(mgf=os.path.join(mgf_dir, in_file),
                 fasta=fasta,
                 xi_cnf='D:/user/Swantje/projects/pipeline_prepro_xi_fdr/resources/xi_linear_by_tryp.conf',
                 outpath=outpath,
                 xi_dir='D:/user/Swantje/XiSearch_1.6.731')
