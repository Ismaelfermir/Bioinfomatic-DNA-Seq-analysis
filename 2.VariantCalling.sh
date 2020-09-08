#!/bin/bash
#	Author: Ismael Fernandez
#	Modified year: 2020

# Samtools ordenar el alineamiento

echo $1
cd $1


for file in *.bam
do
	file_name=$( echo $file | sed 's:_L001_R1_001_BaseRecalibrator.bam::g' )
	echo "$file_name"
	echo "Ordenando alineamiento: " $file_name
	time samtools sort -@ 6 $file > $file_name.sorted.bam
	echo "Done!"
done

for file in *sorted.bam
do
	file_name=$( echo $file | sed 's:sorted.bam::g' )
	echo "Indexando bam para " $file_name
	time samtools index $file > $file.bai
	echo "Done!"
done

# mpileup2017 para archivo de VarScan

for file in *.sorted.bam
do
	file_name=$( echo $file | sed 's:.sorted.bam::g' )
	echo "Creando archivo mpileup para VarScan: " $file_name
	time samtools mpileup -f /reference.fa $file > $file_name.mpileup
	echo "Done!"
done



#                      ***************************************


# VarScan SNPs e INDELs

cd $1

GENOME="/reference.fa"
Varscan_cns='java -jar /herramientas/VarScan-2.3.9/VarScan.v2.3.9.jar mpileup2cns'
Varscan_indel='java -jar /herramientas/VarScan-2.3.9/VarScan.v2.3.9.jar mpileup2indel'
Varscan_snp='java -jar /herramientas/VarScan-2.3.9/VarScan.v2.3.9.jar mpileup2snp'

for file in *.mpileup
do
	file_name=$( echo $file | sed 's:.mpileup::g' )
	echo "Calling SNPs para: "$file_name
	time $Varscan_snp $file --strand-filter 1 \
	--min-var-freq 0 \
	--output-vcf 1 > ./$file_name".varscan.snp.vcf"
	echo "VCF SNPs file for: "$file_name" has been generated"
	echo '----'
	echo "Calling INDELs para: "$file_name
	time $Varscan_indel $file --strand-filter 1 \
	--min-var-freq 0 \
	--output-vcf 1 > ./$file_name".varscan.indel.vcf"
	echo "VCF INDELs file for: "$file_name" has been generated"
	echo "Done!"
	echo '---------------------------------------------------------------------------'
done

