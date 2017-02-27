## MYSEQ 

module add NGS/vcftools/0.1.12a

module add NGS/samtools 

module add NGS/tabix

for file in *pyFiltered_ExAc_SnpEff_SelectedMuts.vcf; 
do
filename=$(basename "$file")
filename=${filename%.*}
outfile="$filename""_sorted.vcf"
#echo $outfile

vcf-sort $file > $outfile; 

done

for file in *pyFiltered_ExAc_SnpEff_SelectedMuts_sorted.vcf; do bgzip "$file"; done

for file in  *pyFiltered_ExAc_SnpEff_SelectedMuts_sorted.vcf.gz; do  tabix -fp vcf $file; done

vcf-merge  *pyFiltered_ExAc_SnpEff_SelectedMuts_sorted.vcf.gz > All_Samples_pyFiltered_ExAc_SnpEff_Merged_SelectedMuts.vcf
