#!/usr/bin/env python
import argparse


def read_bed(file_name):
    bed_list = []
    with open(file_name, 'r') as file:
        for line in file:
            bed_list.append((line.split()[0], int(line.split()[1]), int(line.split()[2])))
    return bed_list


def write_bed(bed, file_name):
    with open(file_name, 'w') as file:
        file.writelines(f'{name}\t{start}\t{end}\n' for name, start, end in bed)


def read_fasta(file_name):
    fasta = {}
    with open(file_name, 'r') as file:
        for line in file:
            if line[0] == '>':
                chrom = line[1:].strip()
                fasta[chrom] = ''
            else:
                fasta[chrom] += line.strip()
    return fasta


def sort(bed):
    return sorted(bed)


def subtract_interval(bed1, subtr):
    subtracted = []
    for interval in sorted(bed1):
        if subtr[0] != interval[0] or subtr[2] <= interval[1] or interval[2] <= subtr[1]:
            subtracted.append(interval)
            continue
        if subtr[1] <= interval[1] and interval[2] <= subtr[2]:
            continue
        if subtr[1] <= interval[1] and subtr[2] < interval[2]:
            subtracted.append((interval[0], subtr[2], interval[2]))
            continue
        if interval[1] < subtr[1] and interval[2] <= subtr[2]:
            subtracted.append((interval[0], interval[1], subtr[1]))
            continue
        subtracted.append((interval[0], interval[1], subtr[1]))
        subtracted.append((interval[0], subtr[2], interval[2]))
    return subtracted


def subtract(bed1, bed2):
    for interval in bed2:
        bed1 = subtract_interval(bed1, interval)
    return bed1


def merge(bed, gap):
    bed = sorted(bed)
    merged = []
    i = 0
    while i < len(bed):
        mergering = bed[i]
        while True:
            if i + 1 == len(bed) or bed[i][0] != bed[i + 1][0] or bed[i + 1][1] - bed[i][2] > gap:
                i += 1
                break
            else:
                mergering = (mergering[0], mergering[1], bed[i + 1][2])
                i += 1
        merged.append(mergering)
    return merged


def intersect(bed1, bed2):
    return subtract(bed1, subtract(bed1, bed2))


def getfasta(fasta, bed, output_name):
    with open(output_name, 'w') as file:
        for interval in bed:
            print(fasta[interval[0]][interval[1]:interval[2] + 1])
            try:
                file.write(f'>{interval[0]}:{interval[1]}{interval[2]}\n')
                file.write(f'{fasta[interval[0]][interval[1]:interval[2] + 1]}\n')
            except (IndexError, KeyError):
                pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='The simple script for processing BED files (available in Windows).'
                                                 'Type "python3 goodtools.py [name of function] -h" for further help.'
                                                 'For example "python3 goodtool.py sort -h"')
    subparsers = parser.add_subparsers(dest='method')

    sort_parser = subparsers.add_parser('sort')
    sort_parser.add_argument('--input', type=str, required=True, help='Your file in BED format')
    sort_parser.add_argument('--output', type=str, required=False, help='The name of output file')

    merge_parser = subparsers.add_parser('merge')
    merge_parser.add_argument('--input', type=str, required=True, help='Your file in BED format')
    merge_parser.add_argument('--gap', type=int, default=0, required=False,
                              help='The maximum distance between two intervals, allowing to merge them')
    merge_parser.add_argument('--output', type=str, required=False, help='The name of output file')

    subtract_parser = subparsers.add_parser('subtract')
    subtract_parser.add_argument('--input', type=str, nargs='+', required=True,
                                 help='Requires two files in BED format. The intervals from the second file'
                                      ' will be subtracted from the first file')
    subtract_parser.add_argument('--output', type=str, required=False, help='The name of output file')

    getfasta_parser = subparsers.add_parser('getfasta')
    getfasta_parser.add_argument('--input', type=str, nargs='+', required=True,
                                 help='Requires two files: one in BED format, the second in fasta format.'
                                      'First specify the BED file, than the fasta file')
    getfasta_parser.add_argument('--output', type=str, required=False, help='The name of output file in fasta format')

    intersect_parser = subparsers.add_parser('intersect')
    intersect_parser.add_argument('--input', type=str, nargs='+', required=True,
                                  help='Requires two files in BED format.'
                                       ' Output will contain the intervals presented in both BED files')
    intersect_parser.add_argument('--output', type=str, required=False, help='The name of output file')

    args = parser.parse_args()

    if not args.output:
        args.output = 'output'

    if args.method == 'sort':
        write_bed(sort(read_bed(args.input)), args.output)
    if args.method == 'subtract':
        write_bed(subtract(read_bed(args.input[0]), read_bed(args.input[1])), args.output)
    if args.method == 'merge':
        write_bed(merge(read_bed(args.input), args.gap), args.output)
    if args.method == 'intersect':
        write_bed(intersect(read_bed(args.input[0]), read_bed(args.input[1])), args.output)
    if args.method == 'getfasta':
        getfasta(read_bed(args.input[0]), read_fasta(args.input[1]), args.output)
