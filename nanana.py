#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 19:24:23 2022

@author: vika
"""
import pandas as pd
import os
from Bio import SeqIO
import argparse

liste=["seqid","source","type","start","end","score","strand","frame","attribute"]

def genename(g):
    global way
    global right_scaffolds
    way=os.path.dirname(g)
    gene_info=pd.read_csv(g,sep='\t',names=liste)
    df = gene_info[gene_info['type'] == 'mRNA']
    df.reset_index(drop=True, inplace=True)
    df = df[['attribute', 'seqid']]
    right_scaffolds = list()
    for i in range(len(df)):
            if pattern in df["attribute"][i]:
                right_scaffolds.append(df['seqid'][i])

    right_scaffolds = list(set(right_scaffolds))
    for a in right_scaffolds:
        os.mkdir(str(way)+'/scaf_'+ str(a))
        kj = str(way)+'/scaf_'+ str(a)+'/'
        tmg_frame = gene_info[gene_info['seqid'] == a]
        tmg_frame.to_csv(kj+'ewm_scaf_'+str(a)+'.gff', index=False)
                         

def gffki(n):
    info=pd.read_csv(n,sep='\t',names=liste)
    nam=os.path.splitext(os.path.basename(n))[0]        
    for a in right_scaffolds:
        tmp_frame = info[info['seqid'] == a]
        kj = str(way)+'/scaf_'+ str(a)+'/'
        tmp_frame.to_csv(kj+nam+'_scaf_'+str(a)+'.gff', index=False)

def fastas(f):
    fasta=os.path.abspath(f)
    record_dict = SeqIO.to_dict(SeqIO.parse(open(fasta), "fasta"))       
    for a in right_scaffolds:
        if str(a) in list(record_dict):
            print(a)
            fasta_dict = record_dict[str(a)]
            kj = '/home/vika/Загрузки/test/scaf_'+ str(a)+'/'
            SeqIO.write(fasta_dict, kj+'fasta_scaf'+str(a)+'.fasta', "fasta")


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--protein', type = str)
parser.add_argument('-t', '--transcript', type = str)
parser.add_argument('-g', '--gff', type = str)
parser.add_argument('-f', '--fasta', type = str)
parser.add_argument('-pr', '--predictions', type = str)
parser.add_argument('-pt', '--pattern', type = str)

parsed_arguments = parser.parse_args(); parsed_variables = vars(parsed_arguments)

protein = parsed_variables['protein']
transcript = parsed_variables['transcript']
gff = parsed_variables['gff']
fasta = parsed_variables['fasta']
pattern = parsed_variables['pattern']
predictions = parsed_variables['predictions']

genename(gff)
gffki(protein)
gffki(transcript)
gffki(predictions)
fastas(fasta)
    



