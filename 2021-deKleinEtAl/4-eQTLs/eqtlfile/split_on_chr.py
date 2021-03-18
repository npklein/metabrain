import gzip

line_per_chr = {}
print('start reading')
with gzip.open('eQTL-Cis-EUR-1Mb-biogenformat.txt.gz', 'rt') as input_file:
    header = input_file.readline()
    for index, line in enumerate(input_file):
        if index % 10000 == 0:
            print(index, 'lines processed')
        chr = line.split('\t')[2]
        if chr not in line_per_chr:
            line_per_chr[chr] = []
        line_per_chr[chr].append(line)
        if index == 10000:
            break

print('Done reading, start writing')
for chr in line_per_chr:
    print('write chr',chr)
    with gzip.open('chr'+chr+'.eQTL-Cis-EUR-1Mb-biogenformat.txt.gz','wt') as out:
        out.write(header)
        for line in line_per_chr[chr]:
            out.write(line)
