import argparse
import statistics


parser = argparse.ArgumentParser(description='Subset covariate table on samples, then remove rows that have no variance')
parser.add_argument('covar_table',help='table with covariates')
parser.add_argument('sample_list',help='file with list of samples to keep')
parser.add_argument('out_covar',help='outfile with covariates') 

args = parser.parse_args()


with open(args.sample_list) as input_file:
    samples = input_file.read().split('\n')


with open(args.covar_table) as input_file, open(args.out_covar,'w') as out:
    header = input_file.readline().strip().split('\t')
    index_to_print = []
    for index, element in enumerate(header):
        if element in samples:
            index_to_print.append(index)
            out.write('\t'+element)
    out.write('\n')
    for line in input_file:
        line = line.strip().split('\t')
        covar_values = [float(line[x]) for x in index_to_print]
        sd = statistics.stdev(covar_values)
        if sd == 0:
            print('Standard deviation of '+line[0]+' is 0 for these samples, excluding it from covar table')
        else:
            out.write(line[0])
            for index in index_to_print:
                out.write('\t'+line[index])

            out.write('\n')
