#################################################
###
###  
###
###  March 2015
###  Hector Tejero, Translational Bioinformatics Unit, CNIO
#############################################################


#import sys
#export PYTHONPATH=$PYTHONPATH:~/NGS_docs/PyVCF-master/
import sys 
import re
import os

if len(sys.argv)==2:
	inputFile = sys.argv[1]
	print "Filtering file " + inputFile
else: 	
	print "Error: File to filter not specified"
	sys.exit()

import vcf
vcf_reader = vcf.VCFReader(open(inputFile, 'r'))

fileName, fileExtension = os.path.splitext(inputFile)

outputFile = "/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/By_Gene/matrix/" + fileName + "_matrix.tsv"

f = open(outputFile, 'w')
#f = open("/home/htejero/EndocrCNIO/NGS_docs/By_", 'w')

f.write("Snp " + ' '.join(vcf_reader.samples) + "\n")

for record in vcf_reader:
	f.write("_".join([record.CHROM, str(record.POS), str(record.ID)]) + " ")
#	print str(len(record.samples))	
	for sample in record.samples:
		if sample['GT']==".":
			f.write("0\t")
		elif sample['GT']=="0/1":
			f.write("1\t")
		elif sample['GT']=="1/1":
			f.write("2\t")
		else:
			f.write("2\t")
	f.write("\n")


f.close()
