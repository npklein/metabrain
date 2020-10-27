#!/bin/bash

JOB_PREFIX=$1
START=$2
END=$3
QOS_INPUT=$4

QOS_OPTIONS=("leftover" "regular" "priority" "panic" "ds", "dev")
QOS="regular"
if [[ " ${QOS_OPTIONS[@]} " =~ " $QOS_INPUT " ]]; then
    if [[ $QOS_INPUT = "panic" ]]; then
      QOS_INPUT="panic mode"
    fi
    QOS=$QOS_INPUT
fi

JOBS_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-10-12-decon-optimizer/jobs/$JOB_PREFIX
LOGFILE=$JOBS_DIR/start.txt

for I in $(seq $START $END); do
	JOBNAME="${JOB_PREFIX}${I}"
	FILENAME="${LOG_DIR}${JOBNAME}.out"
	if [ -f $FILENAME ]; then
		FINISHED=$(tail $FILENAME | grep 'Shutting down' | wc -l)
		if [ $FINISHED = 0 ]; then
      # echo sbatch --qos=$QOS $FILENAME
      echo -e "$JOBNAME:\t$(sbatch --qos=$QOS $FILENAME)" >> $LOGFILE
		fi
	fi
done