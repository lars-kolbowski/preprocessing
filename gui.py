import os
import preprocessing
from multiprocessing import Pool
from functools import partial
import sys

sys.path.append('D:\\software\\wxpython\\wx-3.0-msw')
from gooey import Gooey, GooeyParser


@Gooey()
# if __name__ == '__main__':
def main():
    parser = GooeyParser()
    parser.add_argument('msconvert_exe', widget="FileChooser",
                        default='C:/Program Files/ProteoWizard/ProteoWizard 3.0.9740/msconvert.exe')
    parser.add_argument('mscon_settings',
                        default=['MS2Denoise 20 100 false', 'titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState>'])
    parser.add_argument('split_acq', default='False')
    parser.add_argument('nthr', default=2)
    parser.add_argument('input', widget="DirChooser")
    parser.add_argument('output', widget="DirChooser")

    args = parser.parse_args()

    full_paths = [os.path.join(args.input, rawfile) for rawfile in os.listdir(args.input)]
    full_paths = [x for x in full_paths if not os.path.isdir(x)]

    if args.split_acq == 'True':
        split_acq = True
    elif args.split_acq == 'False':
        split_acq = False

    pool = Pool(processes=int(args.nthr))
    pool.map(partial(preprocessing.process_file, outdir=args.output, mscon_settings=args.mscon_settings, split_acq=split_acq,
                     mscon_exe=args.msconvert_exe), full_paths)
    pool.close()
    pool.join()

if __name__ == '__main__':
    main()
