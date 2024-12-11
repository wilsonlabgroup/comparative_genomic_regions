This script is used to overlap orthologous regions and peaks in each pair of species.

Usage example:

`process_comparative_output.R -i /hpf/largeprojects/mdwilson/huayun/scripts/hal_liftover/R_script_test/rela_halper_frac_rscript_input.txt -f rela --mode S -t 4 --type ucsc --label 10-20-halper-frac --ol_perc 10 --ol_bp 20`
 

Example command should take ~ 3 mins


Options:

`-i INPUT-FILE`: a tab delimited file with four columns:
- species 1
- species 2
- path to liftover results from  species 1 peaks to species 2
- path to species 1 peaks (assumes first 4 columns are chr, start, end, peakname)

`-f FACTOR`: the factor being analyzed i.e. rela

`-m MAXMIN`: Maxgap:Mininum overlap between orthologous regions and peaks, in bp. Used in findOverlap() command. 

`--ol_perc`: minimum percentage of peaks overlapping orthologs regions, default to 10, which means the overlap betwween peaks and orthologous regions is >= 10% of the peak width

`--ol_bp`:  minimum number of bps of peaks overlapping orthologs regions, default to 20, which means the overlap betwween peaks and orthologous regions is >= 20bp

`--mode`: one of L|P|S, recommend using "S"
- Legacy: liftover from query species to anchor species to overlap with anchor peaks
- Strict: liftover and overlap from both perspectives
- Permissive: liftover and overlap from either perspective

`-t THREADS`: number of threads.The actual number of threads used will be -t x number of species. For example, if there are 4 species, and the available threads is 16, you should set -t to 4. 

`--type`: Which alignment was used to generate the orthologs regions. Choose from 'hal','compara', and 'ucsc'. Set to 'ucsc' for the example command.

`--label`: Additional text labels to add to output files


**output**

Output files are generated in the current directory. Files are named based on the `-i input` filename and `--label`

One output table for each species in this format (from human perspective):

Filename: `rela_max0min1_peakConservation_homo_sapiens_ucsc_10-20-halper-frac_S.txt`

|seqnames|start|end|conservation|bos_taurus|mus_musculus
| ---      | ---       | ---      | ---       | ---      | ---       |
|1|778564|778947|10|bos_taurus.25:27772292-27772614|NA
|1|779009|779287|00|NA|NA
|1|827060|827624|1X|bos_taurus.1:157504740-157504990;bos_taurus.27:6107719-6108107|NA


First 3 columns corerspond to human peaks


Outputs `.RData` file containing the R environment of the run for debugging and accessing detailed information


**NOTES on input files**

`liftover results from  species 1 peaks to species 2`:  is the output format of UCSC liftover or halLiftover. It should be a four column file with:

1-3th.: chr, start, and end of orthologous region

4th.: name of the orginal peak


