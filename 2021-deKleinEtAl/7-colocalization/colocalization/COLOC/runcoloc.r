
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("At least four arguments must be supplied: inputdir eqtlsamples cases controls", call.=FALSE)
} 

library(coloc)


inputdir=args[1]
eqtlsamples=args[2]
cases=args[3]
controls=args[4]


#n1<-3780 

# MS
#s2<-14802
#n2<-26703+s2
#s2 = s2/n2

# ALS
#s2<-2579
#n2<-2767+s2

# ALZ
#s2<-24087
#n2<-47793+s2

# PKD
#s2<-20184
#n2<-397324+s2

#s2 = s2/n2

# inputdir="D:\\Sync\\SyncThing\\Postdoc2\\2019-BioGen\\data\\2019-04-Freeze2\\2019-07-29-GWAS\\b38\\eqtlmerge\\ALZ\\"

# get a list of files in the inputdir
files <- list.files(path=inputdir, pattern="*.txt.gz", full.names=TRUE, recursive=FALSE)

# iterate over genes for disease in folder inputfolder
lapply(files,function(filename) {

  if(!endsWith(filename, "coloc.gz") && !endsWith(filename,"colocsummary.gz")){
	fileout=paste(filename,"-coloc.gz",sep="")
	fileoutsummary=paste(filename,"-colocsummary.gz",sep="")
	if(file.exists(fileoutsummary)){
	  file.remove(fileoutsummary)  
	}
	print(paste("Reading",filename))
	data<-read.table(gzfile(filename),header=TRUE,sep="\t")
	
	beta1<-data$b1
	varbeta1<-data$se1*data$se1
	beta2<-data$b2
	varbeta2<-data$se2*data$se2
	
	maf<-data$maf
	
	myres <- coloc.abf(dataset1=list(beta=beta1, varbeta=varbeta1, N=eqtlsamples,type="quant", snp=data$SNP),
					   dataset2=list(beta=beta2, varbeta=varbeta2, N=cases+controls,s=controls, type="cc", snp=data$SNP),
					   MAF=maf)
	write.table(myres$results,file=gzfile(fileout),quote=FALSE,sep="\t")
	write.table(myres$summary,file=gzfile(fileoutsummary),quote=FALSE,sep="\t")  
  }
  
})
