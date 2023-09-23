#!/bin/bash
doc1="

#####################################
Pipeline de Metabarcoding para la
caracterización de macroinvertebrados

Versión: 0.1
Autor: Rosa Bermejo
Tutor: Rubén Cañas
Trabajo: TFM bioinformática VIU
Fecha: Septiembre 2023
####################################
"
doc2="
Opciones:

-h : Muestra la ayuda y sale
-i : Directorio de entrada de los fastq de las muestras
-r : Clasificador de referencia // En caso de estar usando BLAST, nombre común para las secuencias y taxonomía. Las secuencias de referencia deberán ser \${REF}_seqs.qza y la taxonomía \${REF}_tax.qza.
-f : Funcionalidad. Escoge entre utilizar un clasificador (1) o utilizar BLAST (2) para la clasificación de tus secuencias.
La salida por defecto del programa se realiza en la carpeta de entrada
"
# Variables globales / archivos / etc
#REFCLASS='~/References/NCBI/metazoan_COI/ncbi_COI_pt_classifier.qza'
SHDIR=${PWD}
ADAPTER='/home/BioinfoShare/References/SARS-CoV-2/referenceGenome/IlluminaAdaptors.fasta'
T=12

sample_table='samples_fastq_table.tsv'
sample_manifest='qiime_manifest.csv'
metadata='metadata.tsv'

# Entrada de Argumentos #
while [ -n "$1" ]; do
	case "$1" in
		-i|--input)
			WD="$2"
			shift
			;;
		-r|--reference)
			REFCLASS="$2"
			shift
			;;
		-f|--funcionality)
			FUNC="$2"
			shift
			;;
		-h|--help)
			echo "$doc1$doc2"
			exit
			;;
		*)
			echo "Opción $1 no reconocida máquina"
			exit
			;;
		esac
		shift
done
echo $doc1
if [[ $FUNC != 1 && $FUNC != 2 ]]; then
	echo "Se ha introducido una funcionalidad desconocida"
	echo "Saliendo del programa..."
	exit
fi
# El primer paso es buscar los archivos en el directorio de trabajo
echo '-----------------------'
echo 'Recogiendo las muestras'
echo '-----------------------'
echo

bash get_samples.sh ${WD}
mkdir -p ${WD}/fastp
mkdir -p ${WD}/QC
mkdir -p ${WD}/qiime
mkdir -p ${WD}/biom
mkdir -p ${WD}/results
mkdir -p ${WD}/raw_data
cd ${SHDIR}

echo '---------------------------------------'
echo 'Realizando preprocesado de las lecturas'
echo '---------------------------------------'
echo 

cat ${WD}/${sample_table} | while read sample fastq1 fastq2;
do
	fastp \
		--in1 ${WD}/${fastq1} \
		--in2 ${WD}/${fastq2} \
		--out1 ${WD}/fastp/${sample}_R1_qc.fastq.gz \
		--out2 ${WD}/fastp/${sample}_R2_qc.fastq.gz \
		--detect_adapter_for_pe \
		--adapter_fasta ${ADAPTER} \
		--length_required 50 \
		--cut_tail \
		--cut_window_size 10 \
		--cut_mean_quality 30 \
		--thread ${T} \
		--max_len1 250 \
		--max_len2 250 \
		--json ${WD}/QC/${sample}_fastp.json \
		--html ${WD}/QC/${sample}_fastp.html
done

multiqc ${WD}/QC -o ${WD}/QC

echo
echo "Resumen de QC almacenado en ${WD}/multiqc_report"
echo

source ~/anaconda3/etc/profile.d/conda.sh
conda activate qiime

echo '---------------------------'
echo 'Realizando análisis de ASVs'
echo '---------------------------'
echo 


# Importación de las secuencias a un artefacto de qiime
qiime tools import \
	--type 'SampleData[PairedEndSequencesWithQuality]' \
	--input-path ${WD}/qiime_manifest.csv \
	--input-format PairedEndFastqManifestPhred33 \
	--output-path ${WD}/qiime/sample_reads.qza

# Resumen de demultiplexado (en realidad para nosotros será simplemente un resumen sobre cómo están nuestras secuencias. Sería equivalente al report del multiqc.
qiime demux summarize \
	--i-data ${WD}/qiime/sample_reads.qza \
	--o-visualization ${WD}/qiime/sample_reads.qzv
qiime tools export \
	--input-path ${WD}/qiime/sample_reads.qzv \
	--output-path ${WD}/QC/summ_demux/

# Eliminación de malas secuencias y quimeras / Generación de ASVs
qiime dada2 denoise-paired \
	--i-demultiplexed-seqs ${WD}/qiime/sample_reads.qza \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-trunc-len-f 0 \
	--p-trunc-len-r 0 \
	--p-n-threads ${T} \
	--o-table ${WD}/qiime/table.qza \
	--o-representative-sequences ${WD}/qiime/rep-seqs.qza \
	--o-denoising-stats ${WD}/qiime/denoising-stats.qza 

# Revisión de los resultados del denoising
qiime feature-table summarize \
	--i-table ${WD}/qiime/table.qza \
	--o-visualization ${WD}/qiime/table.qzv \
	--m-sample-metadata-file ${WD}/${metadata}
qiime feature-table tabulate-seqs \
	--i-data ${WD}/qiime/rep-seqs.qza \
	--o-visualization ${WD}/qiime/rep-seqs.qzv
qiime metadata tabulate \
	--m-input-file ${WD}/qiime/denoising-stats.qza \
	--o-visualization ${WD}/qiime/denoising-stats.qzv

qiime tools export \
	--input-path ${WD}/qiime/table.qzv \
	--output-path ${WD}/QC/tabla/
qiime tools export \
	--input-path ${WD}/qiime/rep-seqs.qzv \
	--output-path ${WD}/QC/rep-seqs-qc
qiime tools export \
	--input-path ${WD}/qiime/denoising-stats.qzv \
	--output-path ${WD}/QC/denoising

# Árbol filogenético de nuestras muestras aunque finalmente no lo utilizamos (mafft-MSA and fasttree)
qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences ${WD}/qiime/rep-seqs.qza \
	--o-alignment ${WD}/qiime/aligned-rep-seqs.qza \
	--o-masked-alignment ${WD}/qiime/masked-aligned-rep-seqs.qza \
	--o-tree ${WD}/qiime/unrooted-tree.qza \
	--o-rooted-tree ${WD}/qiime/rooted-tree.qza

# Clasificación taxonómica

echo '-----------------------------------'
echo 'Realizando clasificación taxonómica'
echo '-----------------------------------'
echo

if [ $FUNC == 1 ]; then

	qiime feature-classifier classify-sklearn \
		--i-classifier ${REFCLASS} \
		--i-reads ${WD}/qiime/rep-seqs.qza \
		--p-n-jobs ${T} \
		--o-classification ${WD}/results/seqs_taxonomy.qza
elif [ $FUNC == 2 ]; then
	qiime feature-classifier classify-consensus-blast \
		--i-query ${WD}/qiime/rep-seqs.qza \
		--p-perc-identity 0.9 \
		--p-evalue 1e-50 \
		--i-reference-reads ${REFCLASS}_seqs.qza \
		--i-reference-taxonomy ${REFCLASS}_tax.qza \
		--o-classification ${WD}/results/seqs-blast_taxonomy.qza
fi

# Realizamos un barplot sobre los resultados
qiime taxa barplot \
	--i-table ${WD}/qiime/table.qza \
	--i-taxonomy ${WD}/results/seqs*_taxonomy.qza \
	--m-metadata-file ${WD}/${metadata} \
	--o-visualization ${WD}/results/taxa-bar-plots.qzv
## Interesante revisarlo en qiime view

echo '---------------'
echo 'Fin del proceso'
echo '---------------'
echo
## A partir de aquí proceder con phyloseq

. ${SHDIR}/qiime2biom.sh \
	-i ${WD}/qiime/table.qza \
	-t ${WD}/results/seqs*_taxonomy.qza \
	-o ${WD}/biom/table-with-tax-json.biom



conda deactivate

mv ${WD}/*fastq.gz ${WD}/raw_data/
