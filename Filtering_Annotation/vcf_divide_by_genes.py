#################################################
###
###
###  March 2015
###  Hector Tejero, Translational Bioinformatics Unit, CNIO
#############################################################


#import sys
#export PYTHONPATH=$PYTHONPATH:~/NGS_docs/PyVCF-master/

import vcf
import re
import pandas as pd

#vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/vcf/Filtered/MySeq_pyFiltered_NonSyn.vcf', 'r'))

vcf_reader = vcf.VCFReader(open('/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/All_Samples_pyFiltered_ExAc_SnpEff_Merged_CNIOSampleName.vcf', 'r'))

vcf_writer = vcf.Writer(open('/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/By_Gene/None.vcf', 'w'), vcf_reader)

genes = pd.read_csv('/home/htejero/NGS_docs/genes_in_Gene_Snps_Complete_list.txt')

gene = ""
for record in vcf_reader: # for each record 	
#	if 'GI' in record.INFO.keys():
#		next_gene = list(set(record.INFO['GI']))[0]
	#print record.INFO['ANN']
	if len(record.INFO['ANN'])==1:
		next_gene = record.INFO['ANN'][0].split('|')[3]
		
	else:
		for i in range(len(record.INFO['ANN'])):
			next_gene = record.INFO['ANN'][i].split('|')[3]
			if sum(next_gene==genes.Gene)>0: 
				break

	print gene, next_gene
	if next_gene==gene:
		vcf_writer.write_record(record)
	else:
		vcf_writer.close()
		gene = next_gene
		filename = "/local/htejero/EndocrCNIO/MySeq/allVcf_MySeq_Exomes/Snp_eff_Vcfs/Selected_Muts/By_Gene/" + "All_Samples_pyFiltered_ExAc_SnpEff_Merged_CNIOSampleName_By_Gene_" + gene + ".vcf"
		vcf_writer = vcf.Writer(open(filename, 'w'), vcf_reader)
		vcf_writer.write_record(record)


