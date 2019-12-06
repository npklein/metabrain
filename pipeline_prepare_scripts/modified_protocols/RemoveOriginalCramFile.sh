#MOLGENIS nodes=1 ppn=1 mem=1gb walltime=0:10:00

### variables to help adding to database (have to use weave)
#string project
#string uniqueID
###
#string alignedBamOrCram
#string cramFileDir

echo "## "$(date)" Start $0"

newCram=${cramFileDir}${uniqueID}.cram


original_cram_size=$(stat -c %s $alignedBamOrCram)
new_cram_size=$(stat -c %s $newCram)
original_cram_size_5_percent=$(($original_cram_size / 20))


if [[ $new_cram_size -le $((original_cram_size + original_cram_size_5_percent)) && $new_cram_size -ge $((original_cram_size - original_cram_size_5_percent)) ]];
then
    echo "$newCram (size: $new_cram_size) within 5% of $alignedBamOrCram (size: $original_cram_size +/- $original_cram_size_5_percent)"
    echo "Removing $alignedBamOrCram"
    rm $alignedBamOrCram
else
    echo "ERROR: $newCram (size: $new_cram_size) not within 5% of $alignedBamOrCram (size: $original_cram_size +/- $original_cram_size_5_percent)"
    echo "$alignedBamOrCram not removed"
    exit 1;
fi

echo
echo
echo "success!"

echo "## "$(date)" ##  $0 Done "
