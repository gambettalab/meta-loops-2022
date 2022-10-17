#!/usr/bin/env python

import argparse
from Bio import Seq, SeqIO, SeqRecord

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Soft-mask a FASTA file, given its hard-masked version,')
    parser.add_argument('input_unmasked', help='input unmasked FASTA file')
    parser.add_argument('input_hardmasked', help='input hard-masked FASTA file')
    parser.add_argument('output_softmasked', help='output soft-masked FASTA file')
    return parser.parse_args(args)


def main(args=None):
    args = parse_arguments(args)

    input_unmasked = SeqIO.parse(open(args.input_unmasked), 'fasta')
    input_hardmasked = SeqIO.parse(open(args.input_hardmasked), 'fasta')

    with open(args.output_softmasked, 'w') as output_softmasked:
        for seq, seq_masked in zip(input_unmasked, input_hardmasked):
            assert(seq.id == seq_masked.id)
            assert(len(seq.seq) == len(seq_masked.seq))
            seq_softmasked = ''.join(lu if lu == lh else lu.lower()
                for lu, lh in zip(seq.seq, seq_masked.seq))
            seq = SeqRecord.SeqRecord(Seq.Seq(seq_softmasked), id=seq.id, description='')
            SeqIO.write(seq, output_softmasked, 'fasta')


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
