#!/usr/bin/env python3


'''
Wrapper script for epigenetics comparaitve analysis using halLiftover, Ensembl compara, or UCSC liftover. 
It will build the directory structure and appropriately name files for subsequent analysis.
It will generate text files with commands that can subsequently be submitted to the cluster
'''

from argparse import ArgumentParser
import sys,os
from subprocess import PIPE, Popen, run
import pandas as pd
import re, fnmatch
from itertools import permutations

def reformat_peaks(line):
    # line: one line in input table, default to species, factor, peak_file
    species, factor, peak_file = line
    peaks = pd.read_table(peak_file, header=None)
    peaks_new = peaks.iloc[:,0:3]
    peaks_new['name'] = peaks[0]+":"+peaks[1].astype(str) + "-" + peaks[2].astype(str)

    new_peak_file = os.path.join(OUTDIR, factor, "peaks", f"{factor}_{species}.bed")
    peaks_new.to_csv(new_peak_file, sep = "\t", header = False, index = False)
    return new_peak_file

def find_chain(s_pair):
    # given two species names,find the corresponding chain file
    s_pattern = f'{s_pair[0]}To{s_pair[1]}*'
    rule = re.compile(fnmatch.translate(s_pattern), re.IGNORECASE)
    chain = [name for name in os.listdir(CHAIN_FILE_DIR) if rule.match(name)]

    if len(chain) == 0:
        sys.exit(f"No chain file found for {s_pair[0]} and {s_pair[1]} in {CHAIN_FILE_DIR}")
    elif len(chain) > 1:
        sys.exit(f"More than one chain file found for {s_pair[0]} and {s_pair[1]} in {CHAIN_FILE_DIR}")
    else:
        return(os.path.join(CHAIN_FILE_DIR, chain[0]))

def gen_liftover_out_file(s_pair):
    return os.path.join(OUTDIR, f'{FACTOR}-{s_pair[0]}-to-{s_pair[1]}-liftover.bed')

def gen_liftover_cmd(s_pair):
    # given two species names, generate liftover commands
    s1, s2 = s_pair
    liftover_out_file = gen_liftover_out_file(s_pair)
    s1_peak = species_peak_dict[s1]
    chain = find_chain(s_pair)
    liftover_cmd = f'liftOver {s1_peak} {chain} {liftover_out_file} {liftover_out_file}.unmapped -minMatch=0.1 -multiple'
    return liftover_cmd


def main():
    parser = ArgumentParser(description='Wrapper for comparative epigenomics analysis. Build directory structure, reformat peaks, and generate subsequent commands.')
    parser.add_argument('input_file', help='Tab delimited input file specifying factor, species, and peak files.')
    parser.add_argument('-m','--method', default = "ucsc", required = True, help = 'Alignment method to use. Choose from hal, compara, and ucsc')
    parser.add_argument('-o','--outdir', default = ".", help="Output directory. Default to current directory")
    parser.add_argument('--hal', help = 'Path to a hal file if method is "hal"')
    parser.add_argument('--chains', default='/hpf/largeprojects/mdwilson/lib/chains', help='Path to a directory with necessary UCSC chain files if method is "ucsc". Default to mdwilson/lib/chains')
    parser.add_argument('--compara_v', default='102', help='Version of Ensembl compara data base to use if method is "compara", default to 102')
    
    args = parser.parse_args()

    input = args.input_file
    input_data=pd.read_table(input, header=None, names = ["species","factor","peak"])

    global FACTOR 
    FACTOR = input_data["factor"][0] # TODO: make it work with multiple factors
    
    global CHAIN_FILE_DIR 
    CHAIN_FILE_DIR = args.chains

    # make output directory 
    global OUTDIR 
    OUTDIR = args.outdir
    run(f"mkdir -p {OUTDIR}/{FACTOR} {OUTDIR}/{FACTOR}/peaks", shell=True)

    # generate output file names
    commands_file = f'{FACTOR}_{args.method}_qsub_commands'
    R_input = f'{FACTOR}_{args.method}_Rscript_input.txt'

    # reformat peaks
    species_list = input_data["species"]
    new_peak_list = input_data.apply(reformat_peaks, axis = 1)
    global species_peak_dict 
    species_peak_dict = dict(zip(species_list, new_peak_list))

    # generate commands for each pair of species
    s_pairs = list(permutations(species_list,2))

    if args.method == "ucsc":
        with open(commands_file, 'w') as file:
            for x in s_pairs:
                file.write(gen_liftover_cmd(x) + "\n")


    # generate Rscript input
    r_input = pd.DataFrame(s_pairs)
    r_input["file"] = [gen_liftover_out_file(x) for x in s_pairs]
    r_input["peaks"] = [species_peak_dict[x] for x in r_input[0]]
    r_input.to_csv(R_input, sep = "\t", header = False, index = False)

    R_cmd = f'./process_comparative_output.R -i {R_input} -f {FACTOR} --mode S -t 4 --type {args.method}'
    #print(R_cmd)
    R_cmd_file = f'{FACTOR}_{args.method}_Rscript_command'
    with open(R_cmd_file,'w') as file:
        file.writelines(R_cmd)

if __name__ == "__main__":
    main()
