#!/bin/bash

hstring="
Programa para la transformación de la asignación taxonómica y el recuento de OTU/ASV en formato biom a partir de artefactos de qiime (qza).
El programa permitirá obtener el archivo que se pasa al paquete phyloseq de R.

Argumentos:
-i : Input de la tabla de recuento de OTU/ASVs
-t : Input de la clasificación taxonómica
-o : Nombre o ruta para el fichero (.biom) de salida. Si no se introduce nombre se queda como predeterminado 'table-with-tax-json.biom'
-h : Imprime este mensaje y sale del programa
"
OUTNAME=''
while [ -n "$1" ]; do
        case "$1" in
                -i)
                        TABLE="$2"
                        shift
                        ;;
                -t)
                        TAX="$2"
                        shift
                        ;;
                -o)
                        OUTNAME="$2"
                        shift
                        ;;
                -h)
                        echo "$hstring"
                        exit
                        ;;
                *)
                        echo "Opción $1 no reconocida"
                        exit
                        ;;
                esac
                shift
done

if [[ $OUTNAME == '' ]]; then
	OUTNAME='table-with-tax-json.biom'
fi

source ~/anaconda3/etc/profile.d/conda.sh
conda activate qiime

mkdir -p .temp
qiime tools export \
	--input-path ${TABLE} \
	--output-path ./.temp/ \
	> /dev/null 2>&1

qiime tools export \
	--input-path ${TAX} \
	--output-path ./.temp/ \
	> /dev/null 2>&1

echo -e "#OTUID\ttaxonomy\tconfidence" > ./.temp/tax_table.tsv;\
	tail -n +2 ./.temp/taxonomy.tsv >> ./.temp/tax_table.tsv

biom add-metadata \
	-i ./.temp/feature-table.biom \
	-o ./.temp/table-with-tax.biom \
	--observation-metadata-fp ./.temp/tax_table.tsv \
	--sc-separated taxonomy

biom convert \
	-i ./.temp/table-with-tax.biom \
	-o ${OUTNAME} \
	--table-type="OTU table" \
	--to-json

bash ~/Pipelines/Metabarcoding/fix_tax.sh ${OUTNAME}
rm -r ./.temp/
conda deactivate
