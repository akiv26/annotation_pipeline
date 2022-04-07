#!/bin/bash
genome=        #указать название файла с геномом в формате fasta
est=           #указать название файла с последовательностью в формате fasta
maker -CTL
touch protein.maker.gff
sed -i 's/^model_org=all/model_org=/g' maker_opts.ctl 
sed -i 's/^est2genome=0/est2genome=1/g' maker_opts.ctl       
sed -i 's/^genome=/genome=dpp_contig.fasta/g' maker_opts.ctl  
sed -i 's/^est=/est=dpp_est.fasta/g' maker_opts.ctl
maker -base protein
gff3_merge -d ./protein.maker.output/protein_master_datastore_index.log
cat protein.all.gff | grep "       maker   "
cat protein.all.gff | sed -n 2,12p >> protein.maker.gff
rm -rf protein.maker.output
rm -rf protein.all.gff
rm *.ctl
