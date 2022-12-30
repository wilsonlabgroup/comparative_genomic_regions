#/bin/bash

halfile=$1
species1=$2
species2=$3
submit=$4

# get bed file
halStats --bedSequences $species1 $halfile > $species1".bed"

# split bed file per chromosome and another file containing all contigs. 7000000 is an arbitary cutoff
awk '{if($3 < 7000000)print}' $species1".bed" > $species1"-contigs.bed"
awk -v s=$species1 '{if($3 > 7000000)print $0> s"-"$1".bed"}' $species1".bed"

rm -f $species1".bed"

# generate 
for f in $species1*bed;do echo "halLiftover --outPSL $halfile $species1 $f $species2 stdout|pslPosTarget stdin `basename $f .bed`"-to-"$species2".psl"";done > $species1_"halliftover_cmds"

# submit
job_prefix=`echo ${species1:0:3}`

if [ $# -eq 4 ]
then
    submit="T"
else
    submit="F"
fi

echo $submit

if [ $submit == "T" ]
then
    
    queue_commands_slurm.py -i $species1_"halliftover_cmds" -p $job_prefix"hal" -qsub_memory 30g -t 20:00:00 
fi