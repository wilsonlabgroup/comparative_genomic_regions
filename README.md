# Pipeline for comparing genomic regions across species

Given a list of functional regions (e.g. ChIPseq peaks) from one species, extract orthologous regions from a multiple genome alignment for another species, and optionally overlap them with functional regions in this species.

## Wrapper script to generate halLiftover/uscs liftOver commands and subsequent R commands.

**Usage example:**

```shell
run_comparative_epi.py example_data/chipseq_peak_compare_input.txt -o test -m hal --hal EC_6species_Ensembl102_cactus.hal --summit
```

**Options:**

`input_file`: a 3 or 4 column file specifying species names (need to be consistent with species names in hal file or ucsc chain filenames), factor names, peak files, and optionally peak summit files. See `example_data/chipseq_peak_compare_input.txt`. **NOTE** if summit files with peaks names (4th column) is supplied, the peak names must match the ones in peak files

`-m --method`: method used for getting orthologous regions. `hal` for halLiftover, `ucsc` for UCSC liftover

`-o --outdir`: output directory. All commands and output files will be written there.

See help messages for other options.

**Outputs:**

- A directory named by factor name will contain reformated peaks and summits

- `{factor}_{method}_qsub_commands`: halLiftover/ucsc liftover commands. Run these commands first.
- `TF_hal_Rscript_command`: run this file after the liftover is finished

The specific steps are explained below:

## Getting orthologous regions

### Multiple genome alignment resources

- [Ensembl Compara](https://useast.ensembl.org/info/genome/compara/multiple_genome_alignments.html)

Set up Ensembl mySQL database locally and use the ensembl API to extract orthologous regions. 

- [progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) 


Use [HalLiftover](https://github.com/ComparativeGenomicsToolkit/Hal) to extract orthologous regions. 


```shell
halLiftover hal_file src_genome input_bed_file dest_genome output
```

## Merge fragmented orthologous regions returned from progressive cactus .hal file

halLiftover returns highly fragmented regions. An simple approach is used to merge these fragments: 

1. Fragments within a certain distance are merged. The default distance is 50bp
2. Require the merged orthologous region length to be within a certain range, defined by either fraction of original peaks, or a fixed length. 
3. If lifted over summit file is provided, require the merged orthologous region to overlap lifted-over summits.

**Usage example:** 

```shell
halLiftover_fragments_merger.R -q test/TF/peaks/TF_bos_taurus.bed -t test/TF-bos_taurus-to-homo_sapiens-hal_raw.bed -o test/TF-bos_taurus-to-homo_sapiens-hal.bed -s test/TF-bos_taurus-to-homo_sapiens-hal_summits.bed --min_frac 0.75 --max_frac 1.25
```

change to `--min_frac 100 --max_frac 1000` to restrict orthologous regions to be between 100-1000bp

**Outputs:**

Two outputs are produced:
- the file defined by `-o`: this will be the **merged** orthologous regions. This file is used by default as input to `process_comparative_output.R`.

- output`_fragments.bed`: this will be the **raw fragments** from halLiftover that are merged. Use this file as input to `process_comparative_output.R` to get precise amount of orthologous region overlaps.


## Compare orthologous regions with peaks
`process_comparative_output.R` is used to overlap orthologous regions and peaks in each pair of species and report each peak's conservation in all other species. 

This script works with orthologous regions from different multiple genome alignment sources: Ensembl compara, cactus, and pairwise alignment chain files (UCSC).  

Usage example:

`process_comparative_output.R -i R_script_test/rela_halper_frac_rscript_input.txt -f rela --mode S -t 4 --type ucsc --ol_perc 10 --ol_bp 20`
 

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


