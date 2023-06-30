'''
Simulate an input file for tumor fraction prediction
'''
import os
import sys
import argparse

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cftf.simulation import draw_allele

def main(
    output_file: str, 
    num_measurements: int,
    tumor_fraction: float, 
    vaf_min: float, 
    vaf_max: float,
    vaf_precision: float,
    ):

    vaf_range = np.arange(vaf_min, vaf_max, vaf_precision)
    vafs = np.random.choice(vaf_range, size = num_measurements, replace = True)

    true_prob = vafs * tumor_fraction

    alleles = draw_allele(true_prob)

    fp = len(str(vaf_precision).split('.')[1])

    with open(output_file, 'w') as handle:

        handle.write('#num_measurements={}\n'.format(num_measurements))
        handle.write('#tumor_fraction={}\n'.format(tumor_fraction))
        handle.write('#vaf_min={}\n'.format(vaf_min))
        handle.write('#vaf_max={}\n'.format(vaf_max))
        handle.write('#vaf_precision={}\n'.format(vaf_precision))

        for v, a in zip(vafs, alleles):
            handle.write(f"{float(v):.{fp}f}\t{int(a)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Description of your program")

    parser.add_argument("-o", "--output-file", type=str, help="Output file path")
    parser.add_argument("-n", "--num-measurements", type=int, help="Number of measurements")
    parser.add_argument("-t", "--tumor-fraction", type=float, help="Tumor fraction")
    parser.add_argument("--vaf-min", type=float, default = 0.01, help="Minimum VAF value")
    parser.add_argument("--vaf-max", type=float, default = 0.5, help="Maximum VAF value")
    parser.add_argument("--vaf-precision", type=float, default = 0.01, help="VAF precision")

    args = parser.parse_args()

    main(
        args.output_file,
        args.num_measurements,
        args.tumor_fraction,
        args.vaf_min,
        args.vaf_max,
        args.vaf_precision,
    )

