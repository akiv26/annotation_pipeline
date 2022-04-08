#!/bin/bash
genome=$1       
file=$2 
f=$3      
maker -CTL
prefix=$(basename ${genome} | sed "s/.fasta//g")
sed -i 's/^model_org=all/model_org=/g' maker_opts.ctl 
sed -i "s/^genome=/genome=${genome}/g" maker_opts.ctl 
if [[ ${f} == protein ]]
	then   
		sed -i 's/^protein2genome=0/protein2genome=1/g' maker_opts.ctl       
		sed -i "s/^protein=/protein=${file}/g" maker_opts.ctl
		maker -base protein_${prefix} -fix_nucleotides
		gff3_merge -d ./protein_${prefix}.maker.output/protein_${prefix}_master_datastore_index.log
		cat protein_${prefix}.all.gff | grep "	maker	" >  protein_${prefix}.maker.gff
		rm -rf protein.maker.output
		rm -rf protein.all.gff
		rm *.ctl
elif [[ ${f} == est ]]
	then 
		sed -i 's/^est2genome=0/est2genome=1/g' maker_opts.ctl       
		sed -i "s/^est=/est=${file}/g" maker_opts.ctl
		maker -base est_${prefix} -fix_nucleotides
		gff3_merge -d ./est_${prefix}.maker.output/est_${prefix}_master_datastore_index.log
		cat est_${prefix}.all.gff | grep "	maker	" >  est_${prefix}.maker.gff
		rm -rf est.maker.output
		rm -rf est.all.gff
		rm *.ctl
elif [[ ${f} == est_protein ]]
	then
		est=$(echo $2 | cut -f1 -d ',')
		protein=$(echo $2 | cut -f2 -d ',')   
		sed -i 's/^protein2genome=0/protein2genome=1/g' maker_opts.ctl       
		sed -i "s/^protein=/protein=${protein}/g" maker_opts.ctl
		maker -base protein_${prefix} -fix_nucleotides
		gff3_merge -d ./protein_${prefix}.maker.output/protein_${prefix}_master_datastore_index.log
		cat protein_${prefix}.all.gff | grep "	maker	" >  protein_${prefix}.maker.gff
		rm -rf protein.maker.output
		rm -rf protein.all.gff
		sed -i 's/^est2genome=0/est2genome=1/g' maker_opts.ctl       
		sed -i "s/^est=/est=${est}/g" maker_opts.ctl
		maker -base est_${prefix} -fix_nucleotides
		gff3_merge -d ./est_${prefix}.maker.output/est_${prefix}_master_datastore_index.log
		cat est_${prefix}.all.gff | grep "	maker	" >  est_${prefix}.maker.gff
		rm -rf est.maker.output
		rm -rf est.all.gff
		rm *.ctl
fi
