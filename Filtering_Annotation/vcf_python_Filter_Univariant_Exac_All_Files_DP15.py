#################################################
###
###  This script filter single sample vcf files according to the some thresholds
###
###  January 2015
###  Hector Tejero, Translational Bioinformatics Unit, CNIO
#############################################################


#import sys
#export PYTHONPATH=$PYTHONPATH:~/NGS_docs/PyVCF-master/


import pandas as pd
import sys  # para meter el fichero por la linea de comando 
import os  # para gestionar los nombres de fichero 
import re
import vcf
print "### Loading ExAc..."

ExAc = pd.read_csv("/local/htejero/sequencingAP/databases/ExAC.csv", sep='\t')

#ExAc = pd.read_csv("/local/htejero/sequencingAP/databases/ExAC.csv", sep='\t')

#inputFile = '/local/htejero/EndocrCNIO/MySeq/vcf/12S0762_S52.vcf.gz' # '/local/htejero/EndocrCNIO/MySeq/vcf/12S0760_S80.vcf.gz'

if len(sys.argv)>=2:
	Files = sys.argv[1:]
	
else: 	
	print "Error: File to filter not specified"
	sys.exit()

print "### Reading Vcfs..."

for inputFile in Files:
	print "Filtering file " + inputFile
	
	vcf_reader = vcf.VCFReader(open(inputFile, 'r'))

	fileName, fileExtension = os.path.splitext(inputFile)

	if fileExtension=='.gz': 
		realName, ext2 = os.path.splitext(fileName)
		outputFile = realName + "_pyFiltered_ExAc" + ext2 
		outputDiscarded = realName + "_pyDiscarded_ExAc" + ext2 
	elif fileExtension==".vcf":
		outputFile = fileName + "_pyFiltered_ExAc" + fileExtension
		outputDiscarded = fileName + "_pyDiscarded_ExAc" + fileExtension

	vcf_writer = vcf.Writer(open(outputFile, 'w'), vcf_reader)

	vcf_writer_discarded = vcf.Writer(open(outputDiscarded, 'w'), vcf_reader)

	for record in vcf_reader: # for each record 	
		#alt = ",".join(record.ALT)	
		#record.genotype('NA00001')['GT']	
		filters = []  #A list for all the conditions of the filters
		sample = record.samples[0]
		GT = sample['GT']  #Genotype 
		DP = sample['DP']  #Deep read 
		#VF = sample['VF']  #Variant frequency 
		AD = sample['AD']  #Allele depth
		VF = float(sample['AD'][1])/sample['DP']
		
	
		#FILTERS 	
		if record.ID:
			filters.append(True)  #If the variation has an rs ID  pass the filter
		else:
			filters.append(record.QUAL>30)	
			filters.append(DP > 15 )
			if record.num_het==1:
				filters.append(VF>0.3 and VF<0.6)
			elif record.num_hom_alt==1:
				filters.append(VF>0.6)
			else:
				pass
				

		
		if all(filters):   #if all the filters are correct
			vcf_writer.write_record(record)
			#print record.CHROM, record.POS, record.ID, record.REF, record.ALT , record.QUAL, DP, VF,  ','.join(map(str, AD ))


		else:
			#Si no esta miramos si esta en ExAc
			pos = ExAc[(ExAc.CHROMOSOME==int(record.CHROM[3:])) & (ExAc.LOCATION==record.POS)]
			if pos.shape[0]>0:
				if any([x in record.ALT for x in pos.ALT]):
					#print record.CHROM, record.POS, record.ID, record.REF, record.ALT 
					vcf_writer.write_record(record)
				
			else:				
				vcf_writer_discarded.write_record(record)

