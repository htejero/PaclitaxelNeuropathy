metaDir="outputEndocrCNIO1_7-11"
bed_file="/home/htejero/NGS_docs/Paclitaxel_MySeq_Genes_Filtered_Sorted_Collapsed.bed"
module add NGS/bedtools/2.22.0 
for i in $(ls -d */ ); 
do
 dir=${i%%/}
 infile="./$dir/calling/$metaDir""_$dir""_BOTH.annotated.vcf"
 outfile="./$dir/calling/Gene_Panel_$metaDir""_$dir"".vcf"
 
 echo $infile 
 bedtools intersect -a $infile -b $bed_file  -header | awk '$5!="." {print}' > $outfile;
 grep -f /local/htejero/EndocrCNIO/MySeq/Ids_of_SNPS_Gene_Panel.txt $infile  >> $outfile  #Hay una lista de SNPs (con Ids) que 
 cp $outfile /local/htejero/EndocrCNIO/MySeq/fromExome/

done
