#!/usr/bin/env python3

# map orthologs with halLiftover and postprocess with HALPER 
# modified from https://github.com/pfenninglab/halLiftover-postprocessing/blob/master/halper_map_peak_orthologs.sh

from argparse import ArgumentParser
import sys,os
from subprocess import PIPE, Popen, run

def process_peak(input):
    # reformat input peaks, if only 3 columns, add 4th name column, otherwise just take the first 4 columns
    out_peak_file = os.path.join(output_dir, os.path.basename(input) + ".reformat")
    
    if not os.path.isfile(out_peak_file):
        with open(input) as f:
            line = f.readline().split("\t")
            n_cols = len(line)
        if n_cols == 3:
            awk_cmd = f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2, $3, $1 ":" $2 "-" $3}}\' {input} > {out_peak_file}'
            run(awk_cmd, shell = True)
        elif n_cols == 4:
            out_peak_file = input
        else:
            cut_cmd = f'cut -f1-4 {input} > {out_peak_file}'
            run(cut_cmd, shell = True)
    return out_peak_file


def gen_summit_file(input):
    summit_file = os.path.join(output_dir, os.path.basename(input).split(".")[0] + "_summit.bed")

    if not os.path.isfile(summit_file):
        # if summit file does not already exist
        print("Generating summit file", file = sys.stderr)
        # read one line and count number of columns
        with open(input) as f:
            line = f.readline().split("\t")
            n_cols = len(line)
        if n_cols == 10:
            # assumes its narrowPeak format, 10th column is summit position
            print("Assuming narrowPeak file, using the 10th column as summit position", file = sys.stderr)
            awk_cmd = f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2+$10, $2+$10+1, $4, $5, $6}}\' {input} > {summit_file}'
            
        else:
            # otherwise, take the middle point as summit
            print("No summits found, taking the middle point as the summits.", fiile = sys.stderr)
            awk_cmd =  f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $4, $5, $6}}\' {input} > {summit_file}'

        run(awk_cmd, shell = True) 
    return summit_file

def run_halLiftover(input, to_run = False):
    output = os.path.join(output_dir, os.path.basename(input).split(".")[0] + "_" + OUTPUT_LABEL + "_halLiftout.bed")
    cmd = f'halLiftover {args.hal_file} {args.src_genome} {input} {args.dest_genome} {output}'
    if not to_run:
        return output,cmd
    else:
        # submit job, hard coded 
        qsub_cmd = f'echo "{cmd}"|queue_commands_slurm.py -i - -p hal -qsub_memory 20g -t 5:00:00'
        print(qsub_cmd)
        run(qsub_cmd, shell = True)
        return output



def run_HALPER(qFile, tFile, sFile, oFile, max_len, min_len, protect, threads, to_run = False):
    output = os.path.join(output_dir, os.path.basename(oFile))
    cmd = f'HALPER.R -q {qFile} -t {tFile} -s {sFile} -o {output} --max_len {max_len} --min_len {min_len} --protect_dist {protect} -T {threads}'
    if to_run:
        run(cmd, shell = True)
    else:
        return cmd




if __name__ == '__main__':
    parser = ArgumentParser(description='Performs liftover of genomic regions between two genomes in a HAL file and postprocess with HALPER.')
    parser.add_argument('hal_file', help="Input HAL file.")
    parser.add_argument('src_genome', help="Source genome name.")
    parser.add_argument('dest_genome', help="Destination genome name.")
    parser.add_argument('qFile', help="input bed file, 1st 4 columns must be in standard bed format, peak names must be unique")
    parser.add_argument('oFile', help="HALPER output file name")

    parser.add_argument("--output_dir", help = "output directory for halLiftover results and HALPER output")
    parser.add_argument('--sFile', help = "input summit file. 4th column expected to be unique peak names matched with qFile")

    parser.add_argument('-max_len', default = 1500, help='maximum number of base pairs of the ortholog. if < 10, treated as max fraction. Default 1500')
    parser.add_argument('-min_len', default = 100, help='minimum number of base pairs of the ortholog; if < 1, treated as min fraction. Default 100')
    parser.add_argument('-protect_dist',help='summit protection distance', default=25)
    parser.add_argument('-T', '--threads', help = "number of threads used for HALPER")
    parser.add_argument('--submit', action = 'store_true', help = "if set, submit halLiftover commands as jobs")

    args = parser.parse_args()
    
    input_file = args.qFile
    hal_file = args.hal_file
    output_dir = args.output_dir

    OUTPUT_LABEL = args.src_genome + '-' + args.dest_genome

    # reformat peak
    input_peak = process_peak(input_file)
    
    if args.sFile is not None:
        # if a summit file is provided
        summit_file = args.sFile
    else:
        # if summit file is not provided, generate a summit file
        summit_file = gen_summit_file(input_file)

    if args.submit:
        peak_liftover_file = run_halLiftover(input_peak, to_run = True)
        summit_liftover_file = run_halLiftover(summit_file, to_run = True)
    else:
        peak_liftover_file, peak_liftover_cmd = run_halLiftover(input_peak)
        summit_liftover_file, summit_liftover_cmd = run_halLiftover(summit_file)

    halper_cmd = run_HALPER(qFile = input_peak, tFile = peak_liftover_file, sFile = summit_liftover_file, oFile = args.oFile, max_len = args.max_len, min_len = args.min_len, protect = args.protect_dist, threads = args.threads)

    #print(halper_cmd)

    if args.submit:
        print(halper_cmd + "\n")
    else:
        print(peak_liftover_cmd + "\n", summit_liftover_cmd + "\n", halper_cmd + "\n")

    