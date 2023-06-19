SAMPLE_ID=$1
OUT=/run/${SAMPLE_ID}.avinput
IN_VCF=/run/${SAMPLE_ID}.vcf
/annovar/convert2annovar.pl -format vcf $IN_VCF -outfile $OUT -include -withzyg -includeinfo
/annovar/table_annovar.pl $OUT /ref/hg19_annovar_db -buildver hg19 -out /run/${SAMPLE_ID} -remove -otherinfo -protocol refGene,cytoBand,snp138,1000g2015aug_all,cosmic68,clinvar_20221231 -operation g,r,f,f,f,f -nastring NA -csvout