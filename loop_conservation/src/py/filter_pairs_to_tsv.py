#!/usr/bin/env python
"""Extract multi-way interactions overapping a given genomic region from .pairs.gz file.
"""

import argparse
import gzip
import re


def parse_arguments(args=None):
    """Parse command-line arguments.

    Args:
        args (None, optional): Command-line arguments

    Returns:
        Namespace: Populated namespace
    """
    parser = argparse.ArgumentParser(
        description='Extract multi-way interactions overapping a given genomic region from .pairs.gz file.')
    parser.add_argument('input_pairs', help='input pairs.gz file')
    parser.add_argument(
        'output_tsv', help='output .tsv.gz file (read_id, chrom, start 0-based inclusive, end 0-based non-inclusive')
    parser.add_argument('output_stats', help='output .tsv.gz file')
    parser.add_argument('chrom', help='chromosome')
    parser.add_argument('start', type=int, help='start position, 0-based inclusive')
    parser.add_argument('end', type=int, help='end position, 0-based non-inclusive')
    return parser.parse_args(args)


def main(args=None):
    """Main function.

    Args:
        args (None, optional): Command-line arguments
    """
    args = parse_arguments(args)

    input_pairs = gzip.open(args.input_pairs, "rt")
    output_tsv = gzip.open(args.output_tsv, "wt")

    header_dict = {}
    alignment_buffer = []
    read_id = None

    global all_interactions, overlapping_interactions
    all_interactions = 0
    overlapping_interactions = 0

    def add_to_buffer(alignment_buffer, chrom, start, end):
        """Add a tuple (chromosome, start, end) to `alignment_buffer`. Stored coordinates are in BED format (zero-based half-open intervals).

        Args:
            alignment_buffer (list): List of tuples
            chrom (str): Chromosome
            start (int): Start coordinate, one-based
            end (int): End coordinate, one-based, can be smaller than the start coordinate
        """
        if start < end:
            alignment_tuple = (chrom, start - 1, end)
        else:
            alignment_tuple = (chrom, end - 1, start)

        if not (alignment_buffer and alignment_buffer[-1] == alignment_tuple):
            if alignment_tuple[0] != '!':
                alignment_buffer.append(alignment_tuple)

    def flush_buffer():
        """Flush the alignments from `alignment_buffer` to the `output_tsv` file.
        """
        global all_interactions, overlapping_interactions

        all_interactions += 1

        this_overlapping = False
        for alignment in alignment_buffer:
            if alignment[0] == args.chrom and alignment[1] < args.end and args.start < alignment[2]:
                this_overlapping = True

        if this_overlapping:
            overlapping_interactions += 1
            for alignment in alignment_buffer:
                output_tsv.write(str(all_interactions) + '\t' + alignment[0] + '\t' +
                                 str(alignment[1]) + '\t' + str(alignment[2]) + '\n')

        alignment_buffer.clear()

    for line in input_pairs:
        if line[0] == '#':
            if line.startswith('#columns: '):
                header = re.sub('^#columns: ', '', line.strip('\n'))
                header_dict = {v: index for index,
                               v in enumerate(header.split())}
                output_tsv.write('read_id\tchrom\tstart\tend\n')
            continue

        line_split = line.strip('\n').split('\t')
        if alignment_buffer and line_split[0] != read_id:
            flush_buffer()

        # assert(line_split[header_dict['pos1']] == line_split[header_dict['pos51']]
        #        or line_split[header_dict['pos5']] == line_split[header_dict['pos31']])
        # assert(line_split[header_dict['pos2']] == line_split[header_dict['pos52']]
        #        or line_split[header_dict['pos5']] == line_split[header_dict['pos32']])

        add_to_buffer(alignment_buffer, line_split[header_dict['chrom1']],
                      int(line_split[header_dict['pos51']]), int(line_split[header_dict['pos31']]))
        add_to_buffer(alignment_buffer, line_split[header_dict['chrom2']],
                      int(line_split[header_dict['pos52']]), int(line_split[header_dict['pos32']]))
        read_id = line_split[0]

    flush_buffer()

    input_pairs.close()
    output_tsv.close()

    output_stats = open(args.output_stats, "wt")
    output_stats.write('all_interactions\toverlapping_interactions\n')
    output_stats.write(f'{all_interactions}\t{overlapping_interactions}\n')
    output_stats.close()


if __name__ == '__main__':
    import sys
    main(sys.argv[1:])
