import datetime
import gzip
import glob
import argparse

parser = argparse.ArgumentParser(description='Merge multiple FeatureCount matrices')
parser.add_argument('input_dir')
parser.add_argument('output_dir')

args = parser.parse_args()


for type in ["metaExon.countAll", "metaExon.countFraction", "transcript.countAll", "transcript.countFraction"]:#, "exon.countAll", "exon.countFraction", "metaExon.countAll", "metaExon.countFraction"]:
    print(type)
    out_name_list = []
    for f in glob.glob(args.input_dir+'/*'+type+'*gz')+glob.glob(args.input_dir+'/*'+type+'*txt'):
        out_name_list.append(f.split('/')[-1].split('.')[1])
    
    now = datetime.datetime.now()
    out_name = now.strftime("%Y-%m-%d")+'.'+'-'.join(out_name_list)
    prev_header = None

    with open(args.output_dir+'/'+out_name+'.'+type+'.txt','w') as out:
        print(glob.glob(args.input_dir+'/*'+type+'*.txt'))
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
                print(prev_header)
                print(header)
                raise RuntimeError('headers not the same')
            prev_header = header
            for line in input_file:
                if f.endswith('gz'):
                    line = line.decode('utf-8')
                out.write(line)
            input_file.close()
