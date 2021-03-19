# Merge STAR gene count files (counts per sample) to matrix with rows genes, columns samples
# NOTE: Needs Python version 3.5+
import glob
import argparse

parser = argparse.ArgumentParser(description='Merge multiple STAR count files into a matrix.')
parser.add_argument('star_base_path', help='base path from where to search for star *ReadsPerGene.out.tab files')
parser.add_argument('outfile', help='output file name')

args = parser.parse_args()

# save the expression per sample in a dictionary to loop over later
expression_per_sample = {}

# keep list of feature and samples so that we can keep the same order later
list_of_features = []
list_of_samples = []
# also keep a set to check if it's already in (faster than using the list)
set_of_features = set([])
set_of_samples = set([])

# loop over all the STAR expression files
for f in glob.iglob(args.star_base_path+'/**/star/*ReadsPerGene.out.tab', recursive=True):
    print(f)
    # sample_name of the sample is everything before .ReadsPerGene
    sample_name = f.split('/')[-1].split('.ReadsPerGene')[0]
    if sample_name not in set_of_samples:
        list_of_samples.append(sample_name)
    if sample_name in expression_per_sample:
        raise RuntimeError('sample_name should be unique')
    expression_per_sample[sample_name] = {}
    with open(f) as input_file:
        for line in input_file:
            line = line.strip().split('\t')
            feature =  line[0]
            if feature not in set_of_features:
                list_of_features.append(feature)
                set_of_features.add(feature)
            count = line[1]
            if feature in expression_per_sample[sample_name]:
                raise RuntimeError('feature should be unique per sample')
            expression_per_sample[sample_name][feature] = count


# Write the matrix
with open(args.outfile,'w') as out:
    for sample in list_of_samples:
        out.write('\t'+sample)
    out.write('\n')

    for feature in list_of_features:
        out.write(feature)
        for sample in list_of_samples:
            out.write('\t'+expression_per_sample[sample][feature])
        out.write('\n')


print('Output written to:',args.outfile)
