import glob

for f in glob.glob('*/parameter_files/parameters.csv'):
    with open(f) as input_file:
        header = inut_file.readline()
        for line in input_file:
            print(line)
            exit()
