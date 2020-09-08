#!/bin/bash

cp /herramientas/bwa-0.7.17/bwa /usr/bin/

echo 'El directorio de los fastq.gz no debe tener espacios'

cd $1
echo 'Redireccionando al directorio'
echo $1
sleep 2

echo 'Selecionando la referencia'
reference='/reference.fa'



echo 'BWA tool'
for  file1 in *_R1_*
do
	file2=$(echo $file1 | sed 's:_R1_:_R2_:g')
	samplename=$(echo $file1 | sed 's:_R1_001.fastq.gz::g')
	echo 'Comenzando alineamiento de la muestra:' $samplename
	time bwa mem -t 6 -M /reference.fa $file1 $file2 > $samplename'.sam'
done



echo 'Samtools tool'
for  file in *.sam
do
	samplename=$(echo $file | sed 's:.sam::g')
	echo 'Transformando .sam en .bam --' $samplename
	time samtools view -S -b $file > $samplename'.bam'
done



echo 'Picard tool'
for  file in *.bam
do
	samplename=$(echo $file | sed 's:.bam::g')
	echo 'Añadiendo nombre a los grupos de lecturas del .bam --' $samplename
	time java -jar /herramientas/picard.jar AddOrReplaceReadGroups \
	I=$file \
	O=$samplename'_ReadGroups.bam' \
	RGID=4 \
	RGLB=DTS \
	RGPL=ILLUMINA \
	RGPU=unit1 \
	RGSM=$samplename
done



echo 'Samtools tool'
for  file in *_ReadGroups.bam
do
	samplename=$(echo $file | sed 's:_ReadGroups.bam::g')
	echo 'Ordenando .bam --' $samplename
	time samtools sort -@ 6 $file -o $samplename'_sorted.bam'
done



echo 'Picard tool'
for  file in *sorted.bam
do
	samplename=$(echo $file | sed 's:_sorted.bam::g')
	echo 'Marcando duplicados --' $samplename
	time java -jar /herramientas/picard.jar MarkDuplicates I=$file  O=$samplename'_MarkDuplicates.bam' M=$samplename'_MarkDuplicate.txt'
done



echo 'Picard tool'
for  file in *_MarkDuplicates.bam
do
	samplename=$(echo $file | sed 's:_MarkDuplicates.bam::g')
	echo 'Ordenando archivo de alineamiento --' $samplename
	time java -jar /herramientas/picard.jar SortSam I=$file  O=$samplename'_sortv2.bam' SORT_ORDER=coordinate
done



echo 'GATK recalibrator tool'
for  file in *_sortv2.bam
do
	samplename=$(echo $file | sed 's:_sortv2.bam::g')
	echo 'Obteniendo tabla para recalibrar --' $samplename
	time gatk BaseRecalibrator \
	-I $file \
	-R /reference.fa \
	--known-sites SNP.data.base.vcf \
	-O $samplename'.table'
done



echo 'GATK recalibrator tool'
for  file in *_sortv2.bam
do
	samplename=$(echo $file | sed 's:_sortv2.bam::g')
	sampletable=$(echo $file | sed 's:_sortv2.bam:.table:g')
	echo 'Recalibrando el alineamiento a partir de la tabla de recalibrado --' $samplename
	time  gatk ApplyBQSR \
	-R /reference.fa \
	-I $file \
	--bqsr-recal-file $sampletable \
	-O $samplename'_BaseRecalibrator.bam'
done


echo '--------------------------------------------------------------------------------------'
echo 'Archivos alineados'
