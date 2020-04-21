#!/usr/bin/env python
import argparse


def symbols2values(symbols):
    return [ord(el) - ord('!') for el in symbols]


def calculate_gc_content(read):
    if len(read) == 0:
        return 0
    return len([i for i in read if i in "GC"]) * 100 / len(read)


def filter(read, min_length, gc_min, gc_max):
    if len(read) < min_length:
        return 'wrong length'
    gc_content = calculate_gc_content(read)
    if gc_min <= gc_content <= gc_max:
        return 'good'
    else:
        return 'wrong gc_content'


def reverse_tuple(a, b):
    return a[::-1], b[::-1]


def sliding_window(read, quality, threshold, window):
    quality_val = symbols2values(quality)
    for i in range(len(quality_val) - window + 1):
        if sum(quality_val[i:i + window]) < threshold * window:
            return read[:i], quality[:i]
    return read, quality


def leading(read, quality, threshold):
    quality_val = symbols2values(quality)
    thresh_index = -1
    for i in range(len(quality_val)):
        if quality_val[i] < threshold:
            thresh_index = i
        else:
            break
    return read[thresh_index + 1:], quality[thresh_index + 1:]


def trailing(read, quality, threshold):
    return reverse_tuple(*leading(*reverse_tuple(read, quality), threshold))


def crop(read, quality, length):
    if len(read) <= length:
        raise ValueError
    else:
        return reverse_tuple(*headcrop(*reverse_tuple(read, quality), length))


def headcrop(read, quality, length):
    if len(read) <= length:
        raise ValueError
    else:
        return read[length:], quality[length:]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='This is a Python3 script working in two modes: filtering and trimming'
                                                 ' for .fastq reads')
    subparsers = parser.add_subparsers(dest='mode')

    filter_parser = subparsers.add_parser('filter',
                                          help="Type 'fastq_barbershop.py filter -h' for further help")
    filter_parser.add_argument('--min_length', type=int, required=True, default=0,
                               help='Minimal read length, default = 0')
    filter_parser.add_argument('--keep_filtered', action="store_true",
                               help='This flag keep filtered reads in distinct output file')
    filter_parser.add_argument('--gc_bounds', type=int, nargs='+', required=True,
                               help='The threshold for GC rate in read')
    filter_parser.add_argument('--output_base_name', type=str, required=False,
                               help='The prefix name for output files')
    filter_parser.add_argument('fastq', type=str,
                               help='Your input file in .fastq format')

    trimmer_parser = subparsers.add_parser('trimmer', help="Type 'fastq_barbershop.py trimmer -h' for further help")
    trimmer_subparser = trimmer_parser.add_subparsers(dest='method')

    sliding_window_parser = trimmer_subparser.add_parser('SLIDINGWINDOW',
                                                                  help='Sliding window scanning the read from '
                                                                  'the start to the end and removes the nucleotides '
                                                                  'from the end of the read'
                                                                  'when the mean quality of nucleotides inside '
                                                                  'the window of specified size drops below specified '
                                                                  'threshold. '
                                                                  'Requires two arguments: '
                                                                  '--threshold - the threshold of quality '
                                                                  '--window - the size of sliding window')
    sliding_window_parser.add_argument('--threshold', type=int, required=True)
    sliding_window_parser.add_argument('--window', type=int, required=True,)

    leading_parser = trimmer_subparser.add_parser('LEADING',
                                                  help='Cutting nucleotides with quality below threshold '
                                                       'from the start of read. '
                                                       'Requires one argument: --threshold - the threshold of quality')
    leading_parser.add_argument('--threshold', type=int, required=True)

    trailing_parser = trimmer_subparser.add_parser('TRAILING',
                                                   help='Cutting nucleotides with quality below '
                                                        ' threshold from the end of read. '
                                                        'Requires one argument: --threshold '
                                                        ' - the threshold of quality')
    trailing_parser.add_argument('--threshold', type=int, required=True)

    crop_parser = trimmer_subparser.add_parser('CROP',
                                               help='Cutting the specified number of nucleotides '
                                                    'from the end of the read. '
                                                    'Requires one argument: --length - the number'
                                                    ' of nucleotides to be cutted')
    crop_parser.add_argument('--length', type=int, required=True)

    headcrop_parser = trimmer_subparser.add_parser('HEADCROP',
                                                   help='Cutting the specified number of nucleotides '
                                                        'from the start of the read. '
                                                        'Requires one argument: --length - the number'
                                                        ' of nucleotides to be cutted')
    headcrop_parser.add_argument('--length', type=int, required=True)

    trimmer_parser.add_argument('fastq', type=str, help='Your input file in .fastq format')

    args = parser.parse_args()

    with open(args.fastq, 'r') as file:
        fastq = [line.strip() for line in file.readlines()]

    if args.mode == 'trimmer':

        if args.method == 'SLIDINGWINDOW':
            with open('output', 'w') as file:
                for i in range(0, len(fastq), 4):
                    file.write(f'{fastq[i]}\n')
                    file.write(f'{sliding_window(fastq[i + 1], fastq[i + 3], args.threshold, args.window)[0]}\n')
                    file.write(f'{fastq[i + 2]}\n')
                    file.write(f'{sliding_window(fastq[i + 1], fastq[i + 3], args.threshold, args.window)[1]}\n')

        if args.method == 'LEADING':
            with open('output', 'w') as file:
                for i in range(0, len(fastq), 4):
                    file.write(f'{fastq[i]}\n')
                    file.write(f'{leading(fastq[i + 1], fastq[i + 3], args.threshold)[0]}\n')
                    file.write(f'{fastq[i + 2]}\n')
                    file.write(f'{leading(fastq[i + 1], fastq[i + 3], args.threshold)[1]}\n')

        if args.method == 'TRAILING':
            with open('output', 'w') as file:
                for i in range(0, len(fastq), 4):
                    file.write(f'{fastq[i]}\n')
                    file.write(f'{trailing(fastq[i + 1], fastq[i + 3], args.threshold)[0]}\n')
                    file.write(f'{fastq[i + 2]}\n')
                    file.write(f'{trailing(fastq[i + 1], fastq[i + 3], args.threshold)[1]}\n')

        if args.method == 'CROP':
            with open('output', 'w') as file:
                for i in range(0, len(fastq), 4):
                    file.write(f'{fastq[i]}\n')
                    file.write(f'{crop(fastq[i + 1], fastq[i + 3], args.length)[0]}\n')
                    file.write(f'{fastq[i + 2]}\n')
                    file.write(f'{crop(fastq[i + 1], fastq[i + 3], args.length)[1]}\n')

        if args.method == 'HEADCROP':
            with open('output', 'w') as file:
                for i in range(0, len(fastq), 4):
                    file.write(f'{fastq[i]}\n')
                    file.write(f'{headcrop(fastq[i + 1], fastq[i + 3], args.length)[0]}\n')
                    file.write(f'{fastq[i + 2]}\n')
                    file.write(f'{headcrop(fastq[i + 1], fastq[i + 3], args.length)[1]}\n')

    if args.mode == 'filter':
        assert args.min_length > 0, "min_length < 0."
        assert 1 <= len(args.gc_bounds) <= 2, "The wrong number of boundaries"

        gc_min = args.gc_bounds[0]
        if len(args.gc_bounds) == 2:
            gc_max = args.gc_bounds[1]
        else:
            gc_max = 100

        if args.output_base_name is None:
            output_path = args.fastq[:args.fastq.rindex('.')]
        else:
            output_path = args.output_base_name

        pass_path = output_path + '__passed.fastq'
        fail_path = output_path + '__failed.fastq'

        if args.keep_filtered:
            failed = open(fail_path, 'w')
            print(failed)

        with open(pass_path, 'w') as passed:
            dropped_gc = 0
            dropped_len = 0
            for i in range(0, len(fastq), 4):
                filter_result = filter(fastq[i + 1], args.min_length, gc_min, gc_max)
                if filter_result == 'good':
                    passed.write(f'{fastq[i]}\n{fastq[i + 1]}\n{fastq[i + 2]}\n{fastq[i + 3]}\n')
                elif args.keep_filtered:
                    failed.write(f'{fastq[i]}\n{fastq[i + 1]}\n{fastq[i + 2]}\n{fastq[i + 3]}\n')
                    if filter_result == 'wrong length':
                        dropped_len += 1
                    else:
                        dropped_gc += 1
                elif filter_result == 'wrong length':
                        dropped_len += 1
                else:
                    dropped_gc += 1

        if args.keep_filtered:
            print(failed)
            failed.close()

        print(f'Number of dropped reads with wrong length: {dropped_len}\n'
              f'Number of dropped reads with wrong gc: {dropped_gc}')