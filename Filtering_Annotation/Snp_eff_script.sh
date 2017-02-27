for file in *.vcf;
do

filename=$(basename "$file")
filename=${filename%.*}
outfile="Snp_eff_Vcfs/""$filename""_SnpEff.vcf"
echo $outfile
java -Xmx4g -jar  ~/snpEff/snpEff.jar -canon -no-downstream -no-upstream -no-intergenic -no-utr GRCh37.75 $file > $outfile

# echo $outfile 
# bedtools intersect -a $infile -b $bed_file  -header | awk '$5!="." {print}' > $outfile;
# grep -f /local/htejero/EndocrCNIO/MySeq/Ids_of_SNPS_Gene_Panel.txt $infile  >> $outfile  #Hay una lista de SNPs (con Ids) que 
# cp $outfile /local/htejero/EndocrCNIO/MySeq/fromExome/


 


done


# java -Xmx4g -jar  ~/snpEff/snpEff.jar -canon -no-downstream -no-upstream -no-intergenic -no-utr GRCh37.75 Gene_Panel_outputEndocrCNIO1_1-11_C1VMEACXX_8_1ss_pyFiltered_ExAc.vcf > $outfile
