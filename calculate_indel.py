# Supplementary Data 4: calculate_indel.py

import difflib
import numpy as np
from scipy.stats import binom
from Bio import SeqIO
import argparse
import itertools


READ_TRUNCATION = 20
HASH_READ_TRUNCATION = 0
MIN_READ_LENGTH = 56
MAX_AMBIGUOUS_BASES = 5
MAX_INDEL_MISMATCH = 6
ERROR_TOLERANCE_THRESHOLD = 0.15

INITIAL_SEARCH_WINDOW = 20  # 20 works well
SEARCH_INCREMENT = 3
MAX_SEARCH_WINDOW = 50

KMER_SIZE = 15

SINGLE_FILE_STRUCTURE = '{}_out.csv'


def find_loc(guide, target):
    loc = target.find(guide)
    return (loc, loc + len(guide))


def rc(seq):
    base_pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(base_pairs[i] for i in seq[::-1])


def generate_hash(seq):
    kmer_index = {}
    for i in range(len(seq) - KMER_SIZE):
        kmer = seq[i:i + KMER_SIZE]
        if kmer in kmer_index:
            kmer_index[kmer] = None
        else:
            kmer_index[kmer] = i
    return kmer_index


def calc_mle(total_reads, indel_counts, background):
    indel_range = np.array(range(indel_counts))
    distrib = binom.pmf(indel_counts - indel_range, total_reads -
                        indel_range, background)
    if len(distrib) == 0:
        mle_freq = 0
    else:
        mle_freq = distrib.argmax() / float(total_reads)

    z = 1.96
    upper_bound = (total_reads * mle_freq + z**2 / 2 + z * np.sqrt(total_reads *
                                                                   mle_freq * (1 - mle_freq) + z**4 / 4)) / (total_reads + z**2)
    lower_bound = (total_reads * mle_freq + z**2 / 2 - z * np.sqrt(total_reads *
                                                                   mle_freq * (1 - mle_freq) + z**4 / 4)) / (total_reads + z**2)

    return mle_freq, lower_bound, upper_bound


def write_mle(sample_sheet, output_file, verbose, quiet):
    if not quiet:
        print 'Applying MLE correction'
    with open(output_file) as start_output_file:
        file_read = start_output_file.read().split('\n')
        output_header = file_read[0]
        output_text = file_read[1:-1]
    controls = []
    with open(sample_sheet) as in_handle:
        for i, l in enumerate(in_handle):
            if len(l.strip().split(',')) < 5:
                print 'Sample and Control flags not detected'
                break
            elif l.strip().split(',')[4][0].upper().strip() == 'C':
                controls.append(i)
    background_list = [float(output_text[i].split(',')[
                             7]) / 100 for i in controls]
    background = sum(background_list) / len(background_list)

    with open(output_file, 'w') as out_handle:
        out_handle.write(output_header+'\n')
        for i, l in enumerate(output_text):
            if i in controls:
                out_handle.write('{},{},{},{}\n'.format(l, 'NA', 'NA', 'NA'))
            else:
                samp_data = [int(l.split(',')[i]) for i in (1, 2, 5, 6)]
                mle_percentage, lower_bound, upper_bound = calc_mle(
                    sum(samp_data), samp_data[1], background)

                out_handle.write('{},{},{},{}\n'.format(
                    l, mle_percentage, lower_bound, upper_bound))


def op_ver(opcodes):
    '''
    Designed to parse the opcodes from difflib.SequenceMatcher to generate edits. Detects if there are an odd number of edits
    and if there are edits with intervening equal regions.
    '''
    ops = [x[0][0] for x in opcodes]
    if len(ops) % 2:
        # assumes read is longer than target
        if not (ops[0] == 'd' and ops[-1] == 'd' and set(ops[1::2]) == set(['e'])):
            return False
        else:
            proc_ops = [(x[0][0], x[3], x[4] - x[3], x[1], x[1] - x[2])
                        for x in opcodes[2:-2:2]]
            return proc_ops
    else:
        return False


def indel_calc_window_hash(seq_handle, target):
    '''
     Iterates through a SeqRecord iterator and calculates statistics about each read for a given window with hash algorithm
    '''
    perf_total, indel_total, err_total, rejected_total, miscall_total, replace_total = (
        0,) * 6
    target_index = generate_hash(target)
    for readout in seq_handle:
        read = str(readout.seq)[HASH_READ_TRUNCATION:]
        if len(read) < MIN_READ_LENGTH or read.count('N') > MAX_AMBIGUOUS_BASES:  # filtering for junk
            rejected_total += 1
        elif target in read:
            perf_total += 1
        else:
            read_index = generate_hash(read)
            mapping = {}
            for kmer in read_index:
                if read_index[kmer] is not None and kmer in target_index and target_index[kmer] is not None:
                    mapping[read_index[kmer]] = target_index[kmer]
            if len(mapping) == 0:
                err_total += 1
            else:
                index_diff = (
                    mapping[i] - i if i in mapping else None for i in range(len(read) + KMER_SIZE + 1))
                collapsed_dif = [[k, len(list(g))]
                                 for k, g in itertools.groupby(index_diff)]

                start = True
                indels = 0
                sing_mismatch = 0
                mult_mismatch = 0
                offset = 0
                if collapsed_dif[-1][0] is not None:
                    err_total += 1
                else:
                    for el in collapsed_dif[:-1]:

                        if start:
                            # advance to first non nan location (trim back from
                            # start of read to first alignment)
                            if el[0] is not None:
                                offset = el[0]
                                start = False
                        if el[0] is not None:
                            doff = el[0] - offset
                            # append indel start loc to iloc and length of indel to ilen
                            # insertion deletion combinations are summarized as follows for computaitonal simplicity
                            # insertion deletion with len(ins)>len(del) = insertion
                            # insertion deletion with len(ins)<len(del) = deletion
                            # insertion deletion with len(ins)==len(del) =
                            # mismatches (currently not considered indel)
                            if doff != 0:
                                indels += 1
                        else:
                            if el[1] < (KMER_SIZE + 1):
                                sing_mismatch += 1
                            elif el[1] > (KMER_SIZE):
                                mult_mismatch += 1
                    if indels > 0:
                        indel_total += 1
                    elif mult_mismatch > 0:
                        replace_total += 1
                        # print collapsed_dif
                    elif sing_mismatch > 0:
                        miscall_total += 1
                    else:
                        err_total += 1

    return (perf_total, indel_total, err_total, rejected_total, miscall_total, replace_total)


def indel_calc_window(seq_handle, target):
    '''
    Iterates through a SeqRecord iterator and calculates statistics about each read for a given window
    '''
    perf_total, indel_total, err_total, rejected_total, miscall_total, replace_total = (
        0,) * 6
    for readout in seq_handle:
        read = str(readout.seq)[READ_TRUNCATION:]
        if len(read) < MIN_READ_LENGTH or read.count('N') > MAX_AMBIGUOUS_BASES:  # filtering for junk
            rejected_total += 1
        elif target in read:
            perf_total += 1
        else:
            opcodes = difflib.SequenceMatcher(
                None, read, target, autojunk=False).get_opcodes()
            # filter out any reads with more than allowed indels + mismatches
            if len(opcodes) > 3 + MAX_INDEL_MISMATCH * 2:
                err_total += 1
            else:
                # if there are not an odd number of edits, try to shift
                # sequence and reattempt
                if not len(opcodes) % 2:
                    opcodes = difflib.SequenceMatcher(
                        None, read, target[1:-1], autojunk=False).get_opcodes()
                indel_list = op_ver(opcodes)
                if not indel_list:
                    err_total += 1
                else:
                    # check if only single mismatched bases, interpreted as
                    # miscalled bases
                    miscall = set.union(set(x[2] for x in indel_list), set(
                        x[4] for x in indel_list), set(x[0] for x in indel_list)) == set(['r', 1, -1])
                    # check for larger replacement regions (not
                    # insertions/deletions)
                    mismatch = set(x[0] for x in indel_list) == set('r')

                    if miscall:
                        miscall_total += 1
                    elif mismatch:
                        replace_total += 1
                    else:
                        indel_total += 1
    return (perf_total, indel_total, err_total, rejected_total, miscall_total, replace_total)


def file_calc(f_name, guide_loc, target, file_type, hash_flag):
    '''
    Attempts different windows to pass error threshold
    '''
    error_flag = True
    window_size = INITIAL_SEARCH_WINDOW
    min_error = 100
    min_total = []
    note = ''

    if hash_flag:
        algorithm = indel_calc_window_hash
    else:
        algorithm = indel_calc_window

    while error_flag:  # attempt windows while above threshold
        target_window = target[guide_loc[0] -
                               window_size:guide_loc[1] + window_size]
        with open(f_name, 'rU') as f_handle:
            total_list = algorithm(
                SeqIO.parse(f_handle, file_type), target_window)

        err_total = total_list[2]
        rejected_total = total_list[3]
        error_percentage = float(err_total) / \
            (sum(total_list) - rejected_total) * 100

        if error_percentage < min_error:  # check if better than previously achieved
            min_error = error_percentage
            min_total = total_list

        error_flag = (error_percentage > ERROR_TOLERANCE_THRESHOLD) and (
            window_size > MAX_SEARCH_WINDOW)
        window_size += SEARCH_INCREMENT

    if error_percentage > ERROR_TOLERANCE_THRESHOLD:
        note = 'Error threshold not met returning best attempt'
    return min_total, note


def prep_entry(f_name, guide, target, file_type, hash_flag):
    '''
    Finds guide location
    '''

    if guide in target:
        total_list, note = file_calc(
            f_name, find_loc(guide, target), target, file_type, hash_flag)
    elif rc(guide) in target:
        total_list, note = file_calc(f_name, find_loc(
            rc(guide), target), target, file_type, hash_flag)
    else:
        total_list = (0,) * 6
        note = 'Guide not found in target sequence'
    return total_list, note


def whole_file_read(sample_sheet, file_type, output_file, hash_flag, mle, verbose, quiet):
    '''
    Reads through a complete file and constructs corresponding output file
    '''
    if not quiet:
        print 'Reading input sheet from {}'.format(sample_sheet)
    if mle:
        mle_string = ''
    else:
        mle_string = ', MLE corrected rate, lower bound, upper bound'
    with open(sample_sheet) as in_handle, open(output_file, 'w') as out_handle:
        out_handle.write(
            'sample,perfect matches,indels,misaligned reads,reads below threshold, reads with miscalled bases, reads with replacements,indel percentage, notes{}\n'.format(mle_string))
        for l in in_handle:
            sample_name, file_name, guide, target = l.strip().split(',')[:4]
            if verbose:
                print 'Analyzing sample {} from {}'.format(sample_name, file_name)
            guide = guide.upper().strip()
            target = target.upper().strip()
            total_list, note = prep_entry(
                file_name, guide, target, file_type, hash_flag)
            indel_total = total_list[1]
            rejected_total = total_list[2] + total_list[3]
            indel_rate = float(indel_total) / \
                (sum(total_list) - rejected_total)
            total_list_string = ','.join(str(s) for s in total_list)
            out_handle.write('{},{},{},{}\n'.format(
                sample_name, total_list_string, indel_rate, note))
    if not mle:
        write_mle(sample_sheet, output_file, verbose, quiet)


def single_entry_read(sample_sheet, file_type, input_name, hash_flag, verbose, quiet):
    '''
    Reads through a single sample
    '''

    with open(sample_sheet) as in_handle:
        for l in in_handle:
            sample_name, file_name, guide, target = l.strip().split(',')[:4]
            if sample_name.strip() == input_name.strip():
                with open(SINGLE_FILE_STRUCTURE.format(input_name.strip()), 'w') as out_handle:
                    guide = guide.upper().strip()
                    target = target.upper().strip()
                    total_list, note = prep_entry(
                        file_name, guide, target, file_type, hash_flag)

                    indel_total = total_list[1]
                    rejected_total = total_list[2] + total_list[3]
                    indel_rate = float(indel_total) / \
                        (sum(total_list) - rejected_total)

                    total_list_string = ','.join(str(s) for s in total_list)

                    out_handle.write('{},{},{},{}\n'.format(
                        sample_name, total_list_string, indel_rate, note))


def combine_files(sample_sheet, file_type, output_file, mle, verbose, quiet):
    '''
    Combines separately processed files
    '''
    if mle:
        mle_string = ''
    else:
        mle_string = ', MLE corrected rate, lower bound, upper bound'

    with open(sample_sheet) as in_handle, open(output_file, 'w') as out_handle:
        out_handle.write(
            'sample,perfect matches,indels,misaligned reads,reads below threshold, reads with miscalled bases, reads with replacements,indel percentage, notes{}\n'.format(mle_string))
        for l in in_handle:
            sample_name, file_name, guide, target = l.strip().split(',')[:4]
            with open(SINGLE_FILE_STRUCTURE.format(sample_name.strip()), 'w') as samp_handle:
                out_handle.write(samp_handle.readline())
    if not mle:
        write_mle(sample_sheet, output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Analyze sequencing data for the presence of indels')
    combine = parser.add_mutually_exclusive_group()
    verbosity = parser.add_mutually_exclusive_group()
    combine.add_argument(
        '-c', '--combine', help='combines files generated by individual samples', action='store_true')
    parser.add_argument(
        '-f', '--fasta', help='reads fasta files (default is fastq)', action='store_true')
    parser.add_argument(
        '-no-m', '--no-mle', dest='nomle', help='does not calculate MLE', action='store_true')
    parser.add_argument('-o', '--output', dest='output_file',
                        help='output file name', default='calc_indel_out.csv')
    parser.add_argument(
        '-a', '--hash', help='uses alternative hashing algorithm', action='store_true')
    parser.add_argument('-i', '--input', dest='sample_sheet',
                        help='input file name', default='sample_sheet.csv')
    combine.add_argument('-s', '--sample', dest='input_name',
                         help='sample name for running in single sample mode')
    verbosity.add_argument(
        '-v', '--verbose', help='outputs verbose', action='store_true')
    verbosity.add_argument(
        '-q', '--quiet', help='supresses output', action='store_true')

    args = parser.parse_args()

    file_type = 'fasta' if args.fasta else 'fastq'
    if args.combine:
        combine_files(args.sample_sheet, file_type,
                      args.output_file, args.nomle, args.verbose, args.quiet)
    elif args.input_name:
        single_entry_read(args.sample_sheet, file_type,
                          args.input_name, args.hash, args.verbose, args.quiet)
    else:
        whole_file_read(args.sample_sheet, file_type,
                        args.output_file, args.hash, args.nomle, args.verbose, args.quiet)
