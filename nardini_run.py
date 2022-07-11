#!/usr/bin/python3

import sys
import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
from localcider.sequenceParameters import SequenceParameters

parser = argparse.ArgumentParser(description= "")
parser.add_argument("-seq", "--seq", help="Sequence to analyze. Can be single sequence or a fasta file")
parser.add_argument("-name", "--name", help="Sequence name")
parser.add_argument("-write", "--write", action='store_true', help="Whether or not to write output to file")
parser.add_argument("-fasta", "--fasta", action='store_true', help="Whether or not input is a fasta file")
parser.add_argument("-plot", "--plot", action='store_true', help="Whether or not to generate a plot of the z-scores")
parser.add_argument("-iter", "--iter", default=100000, help="The number of random iterations to generate")
parser.add_argument("-seed", "--seed", default=None, help="The random seed for randomization")
args = parser.parse_args()

print("Running NARDINI...")
print("\n")

if args.fasta is False:
    seqObj = SequenceParameters(str(args.seq))
    zscores = seqObj.calculate_zscore(num_scrambles=int(args.iter), random_seed=args.seed)
    keys = list(zscores.keys())
    for i in keys:
        arr = zscores[i][3]
        zmat = np.matrix(arr)
    if args.write is False:
        print(zmat)
    else:
        with open(str(args.name) + "_nardini_output.txt", 'w') as f:
            for line in zmat:
                np.savetxt(f, line, fmt='%.5f')
else:
    with open(str(args.seq)) as f:
        tmp = []
        for record in SeqIO.parse(f, 'fasta'):
            sequence = str(record.seq)
            seq_id = str(record.id)
            seqObj = SequenceParameters(sequence)
            zscores = seqObj.calculate_zscore(num_scrambles=int(args.iter), random_seed=args.seed)
            keys = list(zscores.keys())
            for i in keys:
                arr = zscores[i][3]
            flatarr = arr.flatten()
            listarr = flatarr.tolist()
            indexes = [8,16,17,24,25,26,32,33,34,35,40,41,42,43,44,48,49,50,51,52,53,56,57,58,59,60,61,62]
            for index in sorted(indexes, reverse=True):
                del listarr[index]
            df = pd.DataFrame(listarr, columns = [seq_id])
            tmp.append(df)
        zmat = pd.concat(tmp, axis=1)

    if args.write is False:
        print(zmat)
    else:
        if args.fasta is False:
            with open(str(args.name) + "_nardini_output.txt", 'w') as f:
                print(zmat, file=f)
            f.close()
        else:
            zmat.to_csv(str(args.name) + "_nardini_output.txt", index=False, sep="\t")

if args.plot is True:
    if args.write is False:
        sys.exit("To plot the output, the results must first be written out. Use the -write argument to write the NARDINI output to a file.")
    else:
        import subprocess
        print("\n")
        print("Plotting output...")
        subprocess.call(['Rscript', 'nardini_plots.R', args.name])

print("Done!")
