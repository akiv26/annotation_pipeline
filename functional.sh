#!/bin/bash
gff3=$1
genome=$2
base=$3

gffread ${gff3} -g ${genome} -x transcripts.fasta
python3 /home/Roma/maker_test/EVidenceModeler-1.1.1/simple_example/get_longest_ORF.py transcripts.fasta > proteins.fasta

makeblastdb -in ${base} -dbtype prot
blastp -query proteins.fasta -db ${base} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 > output.blastp
/home/Roma/maker_test/interproscan-5.56-89.0/interproscan.sh -appl pfam -dp -f TSV -goterms -iprlookup -pa -tp -i proteins.fasta -o output.iprscan
maker_map_ids --prefix ${4} --justify 8 ${gff3} > hsap_contig.map

filename=$(basename ${gff3})
columns=$(echo $filename | grep -o '\.' | wc -l)
gffn=$(echo $filename | cut -f1-${columns} -d '.')

cp ${gff3} ${gffn}.renamed.gff3
map_gff_ids hsap_contig.map ${gffn}.renamed.gff3

filename=$(basename proteins.fasta)
columns=$(echo $filename | grep -o '\.' | wc -l)
proteinsn=$(echo $filename | cut -f1-${columns} -d '.')

cp proteins.fasta ${proteinsn}.renamed.fasta
map_fasta_ids hsap_contig.map ${proteinsn}.renamed.fasta

filename=$(basename transcripts.fasta)
columns=$(echo $filename | grep -o '\.' | wc -l)
transcriptsn=$(echo $filename | cut -f1-${columns} -d '.')

cp transcripts.fasta ${transcriptsn}.renamed.fasta
map_fasta_ids hsap_contig.map ${transcriptsn}.renamed.fasta
map_data_ids hsap_contig.map output.iprscan
map_data_ids hsap_contig.map output.blastp

maker_functional_gff ${base} output.blastp ${gffn}.renamed.gff3 > ${gffn}.renamed.function.gff
maker_functional_fasta ${base} output.blastp ${proteinsn}.renamed.fasta > ${proteinsn}.renamed.function.fasta
maker_functional_fasta ${base} output.blastp ${transcriptsn}.renamed.fasta > ${transcriptsn}.renamed.function.fasta
ipr_update_gff ${gffn}.renamed.function.gff output.iprscan > ${gffn}_contig.renamed.putative_function.domain_added.gff
rm ${gffn}.renamed.gff3
rm ${proteinsn}.renamed.fasta
rm ${transcriptsn}.renamed.fasta
rm ${gffn}.renamed.function.gff
rm hsap_contig.map
rm output.blastp
rm output.iprscan

