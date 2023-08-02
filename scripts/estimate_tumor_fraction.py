'''Script to estimate tumor fraction based on a monte-carlo approach
'''

import os
import sys
import argparse

import numpy as np
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from cftf.simulation import draw_allele

def main(input_file, output_file, num_simulations, tf_min, tf_max, tf_precision, overwrite):

    if output_file is None:
        output_file = input_file + '.sim'

    alternative_alleles = 0
    with open(input_file, 'r') as handle:
        c = 0
        vafs = []
        for line in handle:
            if line.startswith('#'):
                continue
            line = line.strip('\n').split('\t')
            vaf, ra = float(line[0]), int(line[1])
            vafs.append( vaf )
            c += 1
            alternative_alleles += ra
    vafs = np.array(vafs)
    # assert np.array_equal(vafs, vafs2), "not equal"
    print('Number of measurements: {}'.format(len(vafs)))
    print('Number of alternative alleles: {}'.format(alternative_alleles))

    tfs = np.arange(tf_min, tf_max, tf_precision)

    header = np.arange(0, 10000, 10)


    if os.path.isfile(output_file) and overwrite:
        os.remove(output_file)

    if not os.path.isfile(output_file):
        with open(output_file, 'w') as handle:
            # f'#Number of measurements: {len(vafs)}\n#Number of alternative alleles: {alternative_alleles}\n'+
            handle.write(
                '\t'.join(header.astype(int).astype(str).tolist())+'\n'
            )

        with tqdm(total=num_simulations * len(tfs)) as pbar:
            for tf in tfs:
                sim_alternative_alleles = np.zeros((num_simulations, ), dtype=int)
                for sim_num in range(num_simulations):
                    sim_alternative_alleles[sim_num] = draw_allele(vafs*tf).sum()
                    pbar.update(1)

                sim_alternative_hist, header = np.histogram(
                    sim_alternative_alleles, 
                    bins = 1000,
                    range = (0, 10000)
                )
                
                with open(output_file, 'a') as handle:
                    handle.write(
                        '\t'.join(sim_alternative_hist.astype(str).tolist())+'\n'
                    )

    with open(output_file, 'r') as handle:
        for line_num, line in enumerate(handle):
            if line.startswith("#"):
                continue
            elif line_num == 0:
                header = np.array(line.strip('\n').split('\t')).astype(int)
                
                col_num = np.searchsorted(header, alternative_alleles) - 1
                continue

            tf = tfs[line_num-1]

            counts = np.array(line.strip('\n').split('\t')).astype(int)
            counts = counts/counts.sum()
            #print(counts)

            if counts[col_num] == 0:
                continue

            else:
                print('Tumor fraction {}: {}%'.format(round(tf,4), counts[col_num]*100))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Description of your program")

    parser.add_argument("--input-file", type=str, help="Input file path")
    parser.add_argument("--output-file", type=str, default = None)
    parser.add_argument("--num-simulations", type=int, default = 10000, help="Number of simulations")
    parser.add_argument("--tf-min", type=float, default = 0.00, help="Minimum tumor fraction")
    parser.add_argument("--tf-max", type=float, default = 0.20, help="Maximum tumor fraction")
    parser.add_argument("--tf-precision", type=float, default = 0.01, help="Tumor fraction precision")
    parser.add_argument("--overwrite", action='store_true')

    args = parser.parse_args()

    main(
        args.input_file,
        args.output_file,
        args.num_simulations,
        args.tf_min,
        args.tf_max,
        args.tf_precision,
        args.overwrite,
    )