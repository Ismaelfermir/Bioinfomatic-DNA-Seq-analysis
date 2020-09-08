#!/bin/bash
#	Author: Ismael Fernandez
#	Modified year: 2020

echo $1
cd $1

db='/annovar/humandb'
Annovar='/annovar/table_annovar.pl'

for file in *snp.vcf
do
	file_name=$( echo $file | sed 's:.varscan.snp.vcf::g' )
	echo $file_name
	echo "Anotando: " $file_name
	time perl $Annovar --buildver hg19 --remove --outfile $file_name"_snp" --protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,avsnp150,cosmic70,clinvar_20190305 --operation g,r,f,f,f,f,f,f --nastring "." --vcfinput $file $db
	echo "Done!"
done

for file in *indel.vcf
do
	file_name=$( echo $file | sed 's:.varscan.indel.vcf::g' )
	echo $file_name
	echo "Anotando: " $file_name
	time perl $Annovar --buildver hg19 --remove --outfile $file_name"_indel" --protocol refGene,cytoBand,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,avsnp150,cosmic70,clinvar_20190305 --operation g,r,f,f,f,f,f,f --nastring "." --vcfinput $file $db
	echo "Done!"
done
