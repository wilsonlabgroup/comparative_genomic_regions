#!/usr/bin/env python3


'''
Wrapper script for epigenetics comparaitve analysis using halLiftover, Ensembl compara, or UCSC liftover. 
It will build the directory structure and appropriately name files for subsequent analysis.
It will generate text files with commands that can subsequently be submitted to the cluster

TODO: 
 - flexible species naming (bit tricky to do as it involves too many factors (ensembl, ucsc chains, hal...))
 - Split input peak files to multiple temp files and run hal liftover as an array job. 

'''

from argparse import ArgumentParser
import sys,os
from subprocess import PIPE, Popen, run
import pandas as pd
import re, fnmatch
from itertools import permutations


def gen_summit_file(input, output):
    print("Generating summit file", file = sys.stderr)
    # read one line and count number of columns
    with open(input) as f:
        line = f.readline().split("\t")
        n_cols = len(line)
    if n_cols == 10:
        # assumes its narrowPeak format, 10th column is summit position
        print("Assuming narrowPeak file, using the 10th column as summit position", file = sys.stderr)
        awk_cmd = f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, $2+$10, $2+$10+1, $4, $5, $6}}\' {input} > {output}'
        
    elif n_cols == 3:
        # if 3 columns, take the middle point as summit and merge peak coord as peak name
        print("No summits found, taking the middle point as the summits.", file = sys.stderr)
        awk_cmd =  f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $1":"$2"-"$3}}\' {input} > {output}'
    elif n_cols == 4:
        # if 4 columns, keep the forth column as peak name
        print("No summits found, taking the middle point as the summits.", file = sys.stderr)
        awk_cmd =  f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $4}}\' {input} > {output}'
    else:
        # otherwise, take the middle point as summit
        print("No summits found, taking the middle point as the summits.", file = sys.stderr)
        awk_cmd =  f'awk \'BEGIN {{FS="\t"; OFS="\t"}} {{print $1, int(($2+$3)/2), int(($2+$3)/2)+1, $4, $5, $6}}\' {input} > {output}'

    run(awk_cmd, shell = True) 
    return None


def reformat_peaks(line, summit = True):
    # line: one line in input table, default to species, factor, peak_file
    if line.size == 3:
        species, factor, peak_file = line
    else:
        species, factor, peak_file, summit_file = line
    
    # process peak files
    peaks = pd.read_table(peak_file, header=None)
    peaks_new = peaks.iloc[:,0:3]
    
    # if 4th column exist, add it to the new peak name 
    if peaks.shape[1] == 3:
        # if only three columns, create a 4th column as peak name
        peaks_new['name'] = peaks[0]+":"+peaks[1].astype(str) + "-" + peaks[2].astype(str)
    else:
        # otherwise use the 4th column as peakname
        peaks_new['name'] = peaks[3].astype(str)
    # write out new peak file
    new_peak_file = os.path.join(OUTDIR, factor, "peaks", f"{factor}_{species}.bed")
    peaks_new.to_csv(new_peak_file, sep = "\t", header = False, index = False)

    if summit:
        new_summit_file = os.path.join(OUTDIR, factor, "peaks", f"{factor}_{species}_summits.bed")
    # process summit file
        if os.path.isfile(os.path.expanduser(summit_file)):
            #print("use summit file")
            # if summit file exists
            summits = pd.read_table(summit_file, header=None)
            if summits.shape[1] == 3:
                print(f"summit file {summit_file} does not have peak names")
                sys.exit(1) 
            else:
                # use 4th column as peakname
                # check if peaknames are consistent with names in peak file
                peaknames = peaks_new.iloc[:,3].sort_values()
                summit_peaknames = summits.iloc[:,3].sort_values()
                res = peaknames == summit_peaknames
                if res.all():
                    summits_new = summits.iloc[:,0:4]
                    summits_new.to_csv(new_summit_file, sep = "\t", header = False, index = False)
                else:
                    print(f"peak names don't match between summit {summit_file} and peak file {peak_file}")
                    sys.exit()
                
        else: # summit file doesnt exist, make a summit file from peak file
            gen_summit_file(peak_file, new_summit_file)
        return new_peak_file, new_summit_file
    else:
        # ignore summit file
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

def gen_ortho_out_file(s_pair, ortho_method = "liftover"):
    return os.path.join(OUTDIR, f'{FACTOR}-{s_pair[0]}-to-{s_pair[1]}-{ortho_method}.bed')

def gen_liftover_cmd(s_pair, label):
    # given two species names, generate liftover commands
    s1, s2 = s_pair
    liftover_out_file = gen_ortho_out_file(s_pair, label)
    s1_peak = species_peak_dict[s1]
    chain = find_chain(s_pair)
    liftover_cmd = f'liftOver {s1_peak} {chain} {liftover_out_file} {liftover_out_file}.unmapped -minMatch=0.1 -multiple'
    return liftover_cmd

def gen_hal_cmd(s_pair):
    # given two species names and hal file, generate hal liftover commands
    s1, s2 = s_pair
    hal_out_file = gen_ortho_out_file(s_pair, ortho_method = "hal_raw")
    s1_peak = species_peak_dict[s1]
    hal_cmd = f'halLiftover {args.hal} {s1} {s1_peak} {s2} {hal_out_file}'
    if args.summits:
        s1_summit = species_summit_dict[s1]
        hal_summit_out_file = gen_ortho_out_file(s_pair, ortho_method="hal_summits")
        hal_cmd2 = f'halLiftover {args.hal} {s1} {s1_summit} {s2} {hal_summit_out_file}'
        return(hal_cmd + "\n" + hal_cmd2)
    else:
        return(hal_cmd)
    
def gen_fragMerger_cmd(s_pair, label):
    # given two species names and hal file, generate fragment merger commands
    s1 = s_pair[0]
    hal_out_file = gen_ortho_out_file(s_pair, ortho_method = "hal_raw")
    merger_out_file = gen_ortho_out_file(s_pair, ortho_method = label)
    s1_peak = species_peak_dict[s1]
    if args.summits:
        s1_lifted_summit = gen_ortho_out_file(s_pair, ortho_method="hal_summits")
        merger_cmd = f'halLiftover_fragments_merger.R -q {s1_peak} -t {hal_out_file} -o {merger_out_file} -s {s1_lifted_summit} --min_frac {args.min_frac} --max_frac {args.max_frac}'
    else:
        merger_cmd = f'halLiftover_fragments_merger.R -q {s1_peak} -t {hal_out_file} -o {merger_out_file} --min_frac {args.min_frac} --max_frac {args.max_frac}'
    return(merger_cmd)

def main():
    parser = ArgumentParser(description='Wrapper for comparative epigenomics analysis. Build directory structure, reformat peaks, and generate subsequent commands.')
    parser.add_argument('input_file', help='Tab delimited input file specifying factor, species, peak files, and summit file (optional).')
    parser.add_argument('-m','--method', default = "ucsc", required = True, help = 'Alignment method to use. Choose from hal, compara, and ucsc')
    parser.add_argument('-o','--outdir', default = ".", help="Output directory. Default to current directory")
    parser.add_argument('--hal', help = 'Path to a hal file if method is "hal"')
    parser.add_argument('--chains', default='/hpf/largeprojects/mdwilson/lib/chains', help='Path to a directory with necessary UCSC chain files if method is "ucsc". Default to mdwilson/lib/chains')
   # parser.add_argument('--compara_v', default='102', help='Version of Ensembl compara data base to use if method is "compara", default to 102')
    parser.add_argument('--summits', action='store_true', help = "If set, expects the forth column in the input file to be summits; otherwise make summits from peak file by either taking the real summits (if narrowPeak format) or take the mid point of peaks as summits")
    parser.add_argument('-max_frac', default = 1.25, help='maximum number of base pairs of the ortholog. if < 10, treated as max fraction of original peak. Default 1.25')
    parser.add_argument('-min_frac', default = 0.75, help='minimum number of base pairs of the ortholog; if < 1, treated as min fraction of original peak. Default 0.75')


    global args
    args = parser.parse_args()

    input = args.input_file
    input_data=pd.read_table(input, header=None)
    
    if input_data.shape[1] == 3: 
        input_data.columns = ["species","factor","peak"]
    elif input_data.shape[1] == 4:
        print("summit files exist")
        input_data.columns = ["species","factor","peak","summit"]

    # set a few global paramters 
    global FACTOR 
    FACTOR = input_data["factor"][0] # TODO: make it work with multiple factors
    
    global CHAIN_FILE_DIR 
    CHAIN_FILE_DIR = args.chains

    # make output directory 
    global OUTDIR 
    OUTDIR = args.outdir
    run(f"mkdir -p {OUTDIR}/{FACTOR} {OUTDIR}/{FACTOR}/peaks", shell=True)

    # generate output file names
    commands_file = os.path.join(OUTDIR, f'{FACTOR}_{args.method}_qsub_commands')
    R_input = os.path.join(OUTDIR, f'{FACTOR}_{args.method}_Rscript_input.txt')
    R_cmd_file = os.path.join(OUTDIR,f'{FACTOR}_{args.method}_Rscript_command')

    # reformat peaks
    species_list = input_data["species"]
    if args.summits:
        res = input_data.apply(reformat_peaks, axis = 1)
        new_peak_list, new_summit_list = zip(*res)
    else:
        new_peak_list = input_data.apply(reformat_peaks, axis = 1, summit = False)
    global species_peak_dict 
    species_peak_dict = dict(zip(species_list, new_peak_list))

    if args.summits:
        global species_summit_dict
        species_summit_dict = dict(zip(species_list, new_summit_list))

    # generate commands for each pair of species
    s_pairs = list(permutations(species_list,2))

    if args.method == "ucsc":
        with open(commands_file, 'w') as file:
            for x in s_pairs:
                file.write(gen_liftover_cmd(x, "ucsc") + "\n")

    if args.method == "hal":
        with open(commands_file, 'w') as file:
            for x in s_pairs:
                file.write(gen_hal_cmd(x,) + "\n")
        # generate fragment merger commands
        with open(R_cmd_file, 'w') as file:
            for x in s_pairs:
                file.write(gen_fragMerger_cmd(x, "hal") + "\n")
        

    # generate Rscript input
    r_input = pd.DataFrame(s_pairs)
    r_input["file"] = [gen_ortho_out_file(x, ortho_method = args.method) for x in s_pairs]
    r_input["peaks"] = [species_peak_dict[x] for x in r_input[0]]
    r_input.to_csv(R_input, sep = "\t", header = False, index = False)

    R_cmd = f'process_comparative_output.R -i {R_input} -f {FACTOR} --mode S -t 4 --type {args.method}'
    #print(R_cmd)
    
    with open(R_cmd_file,'a') as file:
        file.write(R_cmd + "\n")

if __name__ == "__main__":
    main()
