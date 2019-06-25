import datetime
import gzip
import glob
import argparse
import datetime
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Merge multiple FeatureCount matrices')
parser.add_argument('input_dir')
parser.add_argument('output_dir')

args = parser.parse_args()


def merge_files(type):
    print(type)
    out_name_list = []
    for f in glob.glob(args.input_dir+'/*'+type+'*gz')+glob.glob(args.input_dir+'/*'+type+'*txt'):
        # sometimes has date a first in name, other times not. Change split value based on that.
        if f.split('/')[-1].startswith('20'):
            s = 1
        else:
            s = 0
        out_name_list.append(f.split('/')[-1].split('.')[s])
    
    now = datetime.datetime.now()
    out_name = now.strftime("%Y-%m-%d")+'.'+'-'.join(out_name_list)
    print('name:',out_name)
    prev_header = None

    with gzip.open(args.output_dir+'/'+out_name+'.'+type+'.txt.gz','wt') as out:
        for index, f in enumerate(glob.glob(args.input_dir+'/*'+type+'*.gz')+glob.glob(args.input_dir+'/*'+type+'*.txt')):
            if f.endswith('gz'):
                input_file = gzip.open(f)
            else:
                input_file = open(f)
            print(f)

            header = input_file.readline()
            if f.endswith('gz'):
                header = header.decode('utf-8')
            if index == 0:
                out.write(header)
            elif prev_header != header:
                print(prev_header[1:10])
                print(header[1:10])
                raise RuntimeError('headers not the same')
            prev_header = header
            for line in input_file:
                if f.endswith('gz'):
                    line = line.decode('utf-8')
                out.write(line)
            input_file.close()

#list_of_types = ["metaExon.countAll", "metaExon.countFraction", "transcript.countAll", "transcript.countFraction", "exon.countAll", "exon.countFraction", "metaExon.countAll", "metaExon.countFraction"]
list_of_types = ["transcript.countAll", "transcript.countFraction"]
p = Pool(len(list_of_types))
p.map(merge_files, list_of_types)
