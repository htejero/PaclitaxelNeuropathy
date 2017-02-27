#################################################
###
###  This script reads a d vcf file and a case-control tsv file
###  and summarizes the genotype information in a table
###  for the moment, it only gets the snps in which there is just one
###  ALT form with respect to the REF
###
###  January 2015
###  Hector Tejero, Translational Bioinformatics Unit, CNIO
#############################################################


#import sys
#export PYTHONPATH=$PYTHONPATH:~/NGS_docs/PyVCF-master/


MINFREQ = 1 -0.30

import pandas as pd

print "### Reading Case-Control File...."

#df1 = pd.read_csv('/local/htejero/EndocrCNIO/MySeq/Case_control_vcf.tsv', sep='\t')

df1 = pd.read_csv('/local/htejero/EndocrCNIO/MySeq/Case_control_MySeq_Exome_207.csv', sep=',')

#df1 = df1.drop(173) #remove a bad file 14S1685

genes = pd.read_csv('/home/htejero/NGS_docs/genes_in_Gene_Snps_Complete_list.txt')

#print "### Loading ExAc..."


ExAc = pd.read_csv("/local/htejero/sequencingAP/databases/ExAC.csv", sep='\t')



print "### Reading Vcf..."

import vcf
#vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/vcf/MySeq_EndocrCNIO_all_samples.vcf', 'r'))
#vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/All_Samples_pyFiltered_ExAc_SnpEff_Merged_CNIOSampleName.vcf', 'r'))

vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/All_Samples_pyFiltered_ExAc_SnpEff_Merged_CNIOSampleName.vcf', 'r'))

#vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/NonIncluded/All_Samples_pyFiltered_ExAc_SnpEff_Merged_NonIncluded.vcf', 'r'))

samples = vcf_reader.samples # [i.sample for i in record.sampl
idx = [i in samples for i in df1.Codigo]
df1 = df1.loc[idx]  

#dado el 

notox =  df1[df1.Muestra=="NO_TOX_Pacl"].Codigo
tox = df1[df1.Muestra=="TOX_PNP_Pacl"].Codigo
#tox = df1[df1.Muestra!="NO_TOX_Pacl"].Codigo
#tox = df1[df1.Muestra=="TOX_PNP_Bortezomib"].Codigo


print "CHROM", "POS", "ID", "REF", "ALT", "QUAL",  "G00_NO_TOX_PCL", "G01_NO_TOX_PCL", "G11_NO_TOX_PCL", "G00_TOX" , "G01_TOX" , "G11_TOX", "GENE", "CONSEQUENCE", "MUTATION", "NFE_ALLELE_FREQUENCY_ExAc"

for record in vcf_reader: # for each record 	
	#alt = ",".joLST3in(record.ALT)	
	#record.genotype('NA00001')['GT']
	freq = 'NA'	
	for sample in record.samples:
		name = sample.sample
		GT = sample['GT']
		DP = sample['DP']
		#VF = sample['VF']
		
		genotypes_tox = [record.genotype(i)['GT'] for i in tox]
		genotypes_notox = [record.genotype(i)['GT'] for i in notox]
		
	
		nGTox = [genotypes_tox.count('.'),  genotypes_tox.count('0/1'),genotypes_tox.count('1/1')]

		
		nGNoTox = [genotypes_notox.count('.'), genotypes_notox.count('0/1'), genotypes_notox.count('1/1')]
		
	if sum(nGTox)==len(tox) and sum(nGNoTox)==len(notox):  #Filter those with several ALT phenotypes
		
		#FILTERS 		
		

		#FIND EXAC FREQUENCY 

		#A = ExAc.loc[(ExAc.CHROMOSOME==record.CHROM) & (ExAc.LOCATION==record.POS) & (ExAc.ALT=="G")]
		
		# Some mutations fall in diferent genes. We take just that we are interested in.
		
		if 'ANN' in record.INFO.keys():
			if len(record.INFO['ANN'])==1:
					gene = record.INFO['ANN'][0].split('|')[3]
					consequence = record.INFO['ANN'][0].split('|')[1]
					mut = record.INFO['ANN'][0].split('|')[10]
			else:
				for i in range(len(record.INFO['ANN'])):
					gene = record.INFO['ANN'][i].split('|')[3]
					consequence = record.INFO['ANN'][0].split('|')[1]
					mut = record.INFO['ANN'][0].split('|')[10]
					if sum(gene==genes.Gene)>0:
						break
		else:
			gene = 'NA'
			consequence = 'NA'
			mut = 'NA'
					
		#Find the Exac Frequency
		pos = ExAc[(ExAc.CHROMOSOME==int(record.CHROM[3:])) & (ExAc.LOCATION==record.POS)]

		if len(pos)>1:  #is the variant
			pos = pos[pos.ALT==str(record.ALT[0])]
		

		if len(pos)==1:

			if  (record.REF==pos['REF'].iloc[0] and (str(record.ALT[0])==pos['ALT'].iloc[0])): #(record.REF in ExAc['REF'][pos] &  str(record.ALT[0])==ExAc['ALT'][0]):
				freq = pos['NFE ALLELE FREQUENCY'].iloc[0]
			else:
				freq = 'NA'
		else: 
			freq = 'NA'
			


		print record.CHROM, record.POS, record.ID, record.REF, ','.join(map(str, record.ALT )) , record.QUAL,  ' '.join(map(str, nGNoTox)), ' '.join(map(str, nGTox)),  gene, consequence, mut, freq

#	else:
#		gene = record.INFO['ANN'][0].split('|')[3]		
#		print record.CHROM, record.POS, record.ID, record.REF, record.ALT, gene

	 
