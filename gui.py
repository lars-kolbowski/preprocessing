import os
import preprocessing
from multiprocessing import Pool
from functools import partial
import sys

sys.path.append('D:/software/wxpython/wx-3.0-msw')
from gooey import Gooey, GooeyParser


@Gooey(
    program_description="Converts raw files to mgfs whith optional MS2 denoising.",
    default_size=(610, 400))
# if __name__ == '__main__':
def main():
    parser = GooeyParser()
    parser.add_argument('input', widget="DirChooser")
    parser.add_argument('output', widget="DirChooser")
    parser.add_argument('config', widget="FileChooser")

    args = parser.parse_args()

    full_paths = [os.path.join(args.input, rawfile) for rawfile in os.listdir(args.input)]
    full_paths = [x for x in full_paths if not os.path.isdir(x)]

    if args.input == args.output:
        args.output = os.path.join(args.output, 'processed')

    execfile(args.config, globals())

    pool = Pool(processes=int(nthr))
    pool.map(partial(preprocessing.process_file, outdir=args.output, mscon_settings=mscon_settings, split_acq=split_acq,
                     detector_filter=detector_filter, mscon_exe=msconvert_exe), full_paths)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
