#!/bin/bash

LOG_DIR=/groups/umcg-biogen/tmp03/output/2019-11-06-FreezeTwoDotOne/2020-03-12-deconvolution/jobs/output/
JOB_PREFIX=CIA
START=1
END=2

echo -e "JOBNAME\tJOBID\tNODE\tERRORS\tPROGRESS\tFINISHED"
for I in $(seq $START $END); do
        JOBNAME="${JOB_PREFIX}${I}"
        FILENAME="${LOG_DIR}${JOBNAME}.out"
        if [ -f $FILENAME ]; then
                JOBID=$(cat $FILENAME | grep 'Resources consumed' | cut -f5 -d' ')
                NODE=$(cat $FILENAME | grep 'Resources consumed' | cut -f13 -d' ' | sed 's/.$//')
                ERRORS=$(cat $FILENAME | grep -i 'error' | wc -l)
                PROGRESS=$(cat $FILENAME | grep 'Processing' | tail -n 1)
                FINISHED=$(tail $FILENAME | grep 'Shutting down' | wc -l)
                echo -e "$JOBNAME\t$JOBID\t$NODE\t$ERRORS\t$PROGRESS\t$FINISHED"
        fi
done
