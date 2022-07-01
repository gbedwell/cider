#!/usr/bin/python3

import argparse
from tabulate import tabulate
from Bio import SeqIO
from localcider.sequenceParameters import SequenceParameters

parser = argparse.ArgumentParser(description= "")
parser.add_argument("-seq", "--seq", help="Sequence to analyze. Can be single sequence or a fasta file")
parser.add_argument("-name", "--name", help="Sequence name")
parser.add_argument("-write", "--write", action='store_true', help="Whether or not to write output to file")
parser.add_argument("-fasta", "--fasta", action='store_true', help="Whether or not input is a fasta file")
args = parser.parse_args()

if args.fasta is False:
    seqObj = SequenceParameters(str(args.seq))
    table = [ ["Kappa", seqObj.get_kappa()],
              ["Net Charge per Residue", seqObj.get_NCPR(pH=None)],
              ["Fraction of Charged Residues", seqObj.get_FCR(pH=None)],
              ["Fraction of Negative Residues", seqObj.get_fraction_negative()],
              ["Fraction of Positive Residues", seqObj.get_fraction_positive()] ]
    output = tabulate(table, headers=["Parameter", "Value"])
else:
    with open(str(args.seq)) as f:
        tmp = []
        for record in SeqIO.parse(f, 'fasta'):
            sequence = str(record.seq)
            seq_id = str(record.id)

            seqObj = SequenceParameters(sequence)
            table = [ ["Kappa", seqObj.get_kappa()],
                      ["Net Charge per Residue", seqObj.get_NCPR(pH=None)],
                      ["Fraction of Charged Residues", seqObj.get_FCR(pH=None)],
                      ["Fraction of Negative Residues", seqObj.get_fraction_negative()],
                      ["Fraction of Positive Residues", seqObj.get_fraction_positive()] ]
            tmp.append(seq_id + " results" + "\n" + "\n" + tabulate(table, headers=["Parameter", "Value"]) + "\n" + "\n" + "\n")
            output = "".join(tmp)

if args.write is False:
    print(output)
else:
    with open(str(args.name) + "_cider_output.txt", 'w') as f:
        print(output, file=f)
