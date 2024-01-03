import click
import glob
import os
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import scipy
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
import math

def seq_combos(curr_seq):
    symbol_dict = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'R': ['G', 'A'], 'Y': ['C', 'T'], 'K': ['G', 'T'],
                   'M': ['A', 'C'], 'S': ['G', 'C'], 'W': ['A', 'T'], 'B': ['G', 'T', 'C'], 'D': ['G', 'T', 'A'],
                   'H': ['A', 'T', 'C'],
                   'V': ['G', 'C', 'A'], 'N': ['G', 'C', 'A', 'T']}
    var_dict = {}
    dict_size = 1
    temp_storage = []
    indic_first = 0

    for bp in curr_seq:
        if indic_first == 0:
            temp_storage = symbol_dict[bp]
        dict_size = dict_size * len(symbol_dict[bp])
        indic_first = indic_first + 1
    for x in range(0, dict_size):
        if len(temp_storage) > x:
            var_dict[x] = temp_storage[x]
        else:
            var_dict[x] = ""
    num_entries = 0
    not_first = 0
    for bp in curr_seq:
        if not_first > 0:
            temp_num_entries = num_entries
            for keyness in range(0, len(var_dict.keys())):
                if keyness <= temp_num_entries:
                    if len(symbol_dict[bp]) > 1:
                        indic = 0
                        temp_seq = var_dict[keyness]
                        for var_bp in symbol_dict[bp]:
                            if indic != 0:
                                num_entries = num_entries + 1
                                var_dict[num_entries] = temp_seq + var_bp
                            else:
                                var_dict[keyness] = temp_seq + var_bp
                            indic = indic + 1
                    else:
                        temp_seq = var_dict[keyness]
                        temp_seq = temp_seq + symbol_dict[bp][0]
                        var_dict[keyness] = temp_seq
        not_first = not_first + 1

    #converts dict of possible primer seqs to kmers where key is primer seq and value is list of kmers
    kmer_dict = {}
    for keyer in var_dict.keys():
        primer_seq = var_dict[keyer]
        kmer_list = generate_kmers(primer_seq, 8)
        kmer_dict[primer_seq] = kmer_list
    return kmer_dict

def generate_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    return kmers

def find_subseq(seq, subseq):
    start = seq.find(subseq)
    end = start + len(subseq) - 1
    return start, end


def ham_dist(seq_one, seq_two):
    if len(seq_one) != len(seq_two):
        raise ValueError("Seqs dont have the same length")

    dist = sum(bp_one != bp_two for bp_one, bp_two in zip(seq_one, seq_two))
    return dist

def sec_check(low_bound, up_bound, primer_seq, seq):
    if((low_bound-(len(primer_seq)-1)) < 0):
        lower_bound = 0
    else:
        lower_bound = low_bound - (len(primer_seq)-1)
    if((up_bound+(len(primer_seq)-1)) > len(seq)-1):
        upper_bound = len(seq)-1
    else:
        upper_bound = up_bound+(len(primer_seq)-1)
    narrow_seq = seq[lower_bound:upper_bound]
    kmer_list = generate_kmers(narrow_seq, len(primer_seq))
    min_dist=len(primer_seq)
    offset_val = 0
    temp = 0
    for kmer in kmer_list:
        dist = ham_dist(kmer, primer_seq)
        if min_dist >= dist:
            min_dist = dist
            offset_val = temp
        temp = temp + 1
    return (lower_bound+offset_val), (lower_bound+offset_val+len(primer_seq))

def histogram_graph(data, seq, name):
    # Create a histogram
    plt.hist(data, bins=range(0, len(seq) + 1), edgecolor='black')

    # Add labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(str(name))

    # Show the plot
    plt.show()


def find_dist(record_one, kmer_dict, primer_seq, rev_or_naw):
    seq_id = record_one[0].id
    seq = record_one[0].seq
    length = len(seq)
    score_arr = [0] * length
    temp_arr=[]

    # iterate through kmers for primer and if found add them to the score_arr
    for keyer in kmer_dict.keys():
        kmer_list = kmer_dict[keyer]
        for kmer in kmer_list:
            start, end = find_subseq(seq, kmer)
            if start != -1:
                for x in range(start, end + 1):
                    temp = score_arr[x]
                    score_arr[x] = temp + 1
                    temp_arr.append(x)
    peaks, _ =scipy.signal.find_peaks(score_arr)
    print(peaks)
    prominences = max(scipy.signal.peak_prominences(score_arr, peaks)[0])
    print(prominences)
    results_full = scipy.signal.peak_widths(score_arr, peaks, rel_height=1)
    print(results_full[0])
    print(*results_full[1:])

    max_index = score_arr.index(max(score_arr))
    temp_index = max_index
    temp_value = 1


    # find distribution of kmers around max value (low_bound)
    while (temp_value != 0 and temp_index != -1):
        temp_value = score_arr[temp_index]
        temp_index = temp_index - 1
    low_bound = temp_index + 2
    # find distribution of kmers around max value (upper_bound)
    temp_index = max_index
    temp_value = 1
    while (temp_value != 0 and temp_index != len(score_arr)):
        temp_value = score_arr[temp_index]
        temp_index = temp_index + 1
    up_bound = temp_index - 2

    if((up_bound - low_bound+1) != len(primer_seq)):
        if(rev_or_naw == True):
            test_val = int(math.sqrt(((up_bound - low_bound+1) - len(primer_seq))**2))
            if test_val < int(len(primer_seq)/2):
                low_bound, up_bound = sec_check(low_bound, up_bound, primer_seq, seq)
                return low_bound, up_bound, (up_bound - low_bound + 1), temp_arr
        return -99, -99, (up_bound - low_bound + 1), temp_arr
    else:
        return low_bound, up_bound, (up_bound - low_bound + 1), temp_arr

def find_dist_test(record_one, kmer_dict, primer_seq, rev_or_naw):
    seq_id = record_one[0].id
    seq = record_one[0].seq
    length = len(seq)
    score_arr = [0] * length
    temp_arr=[]

    # iterate through kmers for primer and if found add them to the score_arr
    for keyer in kmer_dict.keys():
        kmer_list = kmer_dict[keyer]
        for kmer in kmer_list:
            start, end = find_subseq(seq, kmer)
            if start != -1:
                for x in range(start, end + 1):
                    temp = score_arr[x]
                    score_arr[x] = temp + 1
                    temp_arr.append(x)
    peaks, _ = scipy.signal.find_peaks(score_arr)
    print(peaks)
    if len(peaks) != 0:
        top_prom_arr = scipy.signal.peak_prominences(score_arr, peaks)
        top_prom = max(top_prom_arr[0])
        max_indices = [index for index, value in enumerate(top_prom_arr[0]) if value == top_prom]
        #print(max_indices)
        prom_test_val = top_prom/len(kmer_dict.keys())
        print(prom_test_val)
        if(prom_test_val > 1):
            results_full = scipy.signal.peak_widths(score_arr, peaks, rel_height=1)
            for top_index in max_indices:
                print(top_index)

                print(*results_full[1:])
                if(len(peaks)==1):
                    low_bound = results_full[1:][top_index+1]
                    up_bound = results_full[1:][top_index+2]
                else:
                    low_bound = results_full[1:][1][top_index]
                    up_bound = results_full[1:][2][top_index]
                if(up_bound - (low_bound+1) > (len(primer_seq)/2)):
                    print(primer_seq)
                    print(seq[int(low_bound + 1):int(up_bound)])
                    return int(low_bound+1), int(up_bound), (up_bound - (low_bound + 1)), temp_arr
    return -99, -99, None, temp_arr




fasta_dir = '/Users/jrabasc/Desktop/github/pseudo_pooling'
primer_loc_dict = {}

#generate kmer_lists for fwd and rev primers
primer_530_f = 'GTGYCAGCMGCCGCGGTAA'
primer_830_r = 'GGACTACNVGGGTWTCTAAT'
fwd_kmer_dict = seq_combos(primer_530_f)
rev_kmer_dict = seq_combos(primer_830_r)
fwd_rev_comp_kmer_dict = seq_combos(Seq(primer_530_f).reverse_complement())
rev_rev_comp_kmer_dict = seq_combos(Seq(primer_830_r).reverse_complement())




for filename in os.listdir(fasta_dir):
    if ("silva_top_20.fasta" in filename):
        f = os.path.join(fasta_dir, filename)
        file_out_arr = filename.split('.fasta')
        file_out = file_out_arr[0] + '_V4_region_only.fasta'
        f_out = os.path.join(fasta_dir, file_out)
        with open(f_out, 'w') as read_output:
            for (record_one) in zip(SeqIO.parse(f, "fasta")):
                print(record_one[0].id)
                prim_one_low_bound, prim_one_up_bound, prim_one_len, temp_guy = find_dist_test(record_one, fwd_kmer_dict, primer_530_f, False)
                if(prim_one_low_bound == -99):
                    temp_seq = Seq(primer_530_f).reverse_complement()
                    prim_one_low_bound, prim_one_up_bound, prim_one_len, temp_guy = find_dist_test(record_one, fwd_rev_comp_kmer_dict, temp_seq, True)
                print("----------------------")
                prim_two_low_bound, prim_two_up_bound, prim_two_len, two_temp_guy = find_dist_test(record_one, rev_kmer_dict, primer_830_r, False)
                if(prim_two_low_bound == -99):
                    temp_seq = Seq(primer_830_r).reverse_complement()
                    prim_two_low_bound, prim_two_up_bound, prim_two_len, two_temp_guy = find_dist_test(record_one, rev_rev_comp_kmer_dict, temp_seq, True)
                #histogram_graph(two_temp_guy, record_one[0].seq, record_one[0].id)
                print(prim_one_low_bound)
                print(prim_one_up_bound)
                print(prim_two_low_bound)
                print(prim_two_up_bound)
                seq_low_bound = prim_one_low_bound
                seq_up_bound = prim_two_up_bound
                if(prim_one_up_bound > prim_two_low_bound):
                    seq_low_bound = prim_two_low_bound
                    seq_up_bound = prim_one_up_bound
                temp_seq = record_one[0].seq
                #print(seq_up_bound)
                #print(seq_low_bound)
                record_one[0].seq = temp_seq[seq_low_bound:seq_up_bound]
                #SeqIO.write(record_one, read_output, 'fasta')

