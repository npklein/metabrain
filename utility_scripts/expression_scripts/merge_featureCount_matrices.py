import gzip
import glob
import argparse

parser = argparse.ArgumentParser(description='Merge multiple FeatureCount matrices')
parser.add_argument('input_dir')

args = parser.parse_args()


for type in ["transcript.countAll", "transcript.countFraction", "exon.countAll", "exon.countFraction", "metaExon.countAll", "metaExon.countFraction"]:
    print(type)
    out_name_list = []
    for f in glob.glob(args.input_dir+'/*'+type+'*gz'):
        out_name_list.append(f.split('/')[-1].split('.')[0])
    
     
    out_name = '_'.join(out_name_list)
    prev_header = None
    with gzip.open(out_name+'.'+type+'.txt','wt') as out:
        for index, f in enumerate(glob.glob(args.input_dir+'/*'+type+'*gz')):
            print(f)
            with gzip.open(f) as input_file:
                header = input_file.readline().decode('utf-8')
                if index == 0:
                    out.write(header)
                elif prev_header != header:
                    print(prev_header)
                    print(header)
                    raise RuntimeError('headers not the same')
                prev_header = header
                for line in input_file:
                    out.write(line.decode('utf-8'))
