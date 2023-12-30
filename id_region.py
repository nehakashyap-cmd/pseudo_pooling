import click
import glob
import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import gzip
import itertools
import matplotlib.pyplot as plt
import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import json
import random

def seq_combos(curr_seq):
    ref_dict = {'Y': ['C', 'T'], 'M': ['C', 'A']}
    built_seq = ''
    seq_dict = {}
    seq_dict_key = 0
    for base in curr_seq:
        if base not in ref_dict.keys():
            built_seq = built_seq + base
        else:
            pot_bases = ref_dict[base]
            for pot_base in pot_bases:
                seq_dict_key = seq_dict_key + 1




    return seq_arr





fasta_dir = '/Users/jrabasc/Desktop/github/pseudo_pooling'
primer_data_dict = {}
for filename in os.listdir(fasta_dir):
    if (".fasta" in filename):
        f = os.path.join(fasta_dir, filename)
        for (record_one) in zip(SeqIO.parse(f, "fasta")):
            seq_id = record_one[0].id
            #temp_arr = str(primer_id).split('_')
            #primer_key = temp_arr[0]

            primer_530_f = 'GTGYCAGCMGCCGCGGTAA'
            primer_830_r = 'GGACTACNVGGGTWTCTAA'
            primer_one_dict = seq_combos('GTGYCAGCMGCCGCGGTAA')

            seq = record_one[0].seq
            if seq_id not in primer_data_dict.keys():
                primer_data_dict[seq_id] = [seq]
            else:
                primer_data_dict[seq_id].append(seq)


