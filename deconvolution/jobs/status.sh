#!/bin/bash

LOG_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/
JOB_PREFIX=CIA
START=0
END=9

echo -e "JOBNAME\tJOBID\tNODE\tELAPSED\tERRORS\tPROGRESS\tPANIC\tFINISHED"
for I in $(seq $START $END); do
	JOBNAME="${JOB_PREFIX}${I}"
	FILENAME="${LOG_DIR}${JOBNAME}.out"
	if [ -f $FILENAME ]; then
		JOBID=$(tail $FILENAME | grep 'Resources consumed' | cut -f5 -d' ')
		NODE=$(tail $FILENAME | grep 'Resources consumed' | cut -f13 -d' ' | sed 's/.$//')
    ELAPSED=$(tail $FILENAME | grep "${JOBID}.batch" | awk -F' ' '{print $2}')
		ERRORS=$(cat $FILENAME | grep -i 'error' | wc -l)
    PROGRESS=$(cat $FILENAME | grep 'Processing' | tail -n 1 | grep -oP '[0-9]{2}.[0-9]{2}%')
    PANIC=$(cat $FILENAME | grep 'Panic!!!' | wc -l)
		FINISHED=$(tail $FILENAME | grep 'Shutting down' | wc -l)
		echo -e "$JOBNAME\t$JOBID\t$NODE\t$ELAPSED\t$ERRORS\t$PROGRESS\t$PANIC\t$FINISHED"
	fi
done