dir=$1
yourfilenames=`ls $dir/*.sh`

for entry in $yourfilenames
do
	echo $entry
	qsub $entry
done