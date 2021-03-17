import gzip
import sys

input=sys.argv[1]
output=sys.argv[2]

fh=gzip.open(input,'rt')
fho=gzip.open(output,'wt')


outln=[None]*21

	

ctr = 0
for line in fh:
	line = line.strip()
	elems = line.split("\t")
	if ctr == 0:
		# header line
		fho.write("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tRefAllele\tAltAllele\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta\tMeta-SE\tIncludedDatasetsBeta(SE)\tFDR\n")
	else:
		# copy required stuff
		for i in range(0,8):
			outln[i] = elems[i]
		# split alleles
		alleles = elems[8]
		a = alleles.split("/")
		outln[8]=a[0]		# ref
		outln[9]=a[1]		# alt
		
		outln[10]=elems[9]	# assessed
		outln[11]=elems[10]	# z-score
		outln[12]=elems[11]	# datasets
		outln[13]=elems[12]	# zscores
		outln[14]=elems[13]	# samples
		outln[15]=elems[16]	# hgnc
		outln[16]=elems[17]	# correlations
		try:
			# split beta
			betacol=elems[18]
			b = betacol.split(" ")
			outln[17] = b[0].strip()
			# remove braces from SE
			se = b[1].strip("(").strip(")").strip()
			outln[18] = se
		except:
			print("Error splitting: "+betacol+"\t"+";".join(b))
			exit()	
		outln[19] = elems[19]
		outln[20] = elems[21]
		fho.write("\t".join(outln)+"\n")
	ctr = ctr + 1
	if ctr % 10000 == 0:
		print("\r"+str(ctr)+" lines processed")
fh.close()
fho.close()




