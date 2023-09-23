source ~/anaconda3/etc/profile.d/conda.sh
conda activate qiime

qiime rescript get-ncbi-data   --p-query "txid33208[ORGN] AND (cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title] NOT environmental[Title] NOT environmental sample[Title] NOT uncultured[Title] NOT environmental samples[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title])" --p-ranks kingdom subkingdom superphylum phylum subphylum superclass class subclass superorder order suborder superfamily family subfamily genus species --p-rank-propagation --o-sequences metazoan_COI_seqs.qza --o-taxonomy metazoan_COI_tax.qza

conda deactivate
