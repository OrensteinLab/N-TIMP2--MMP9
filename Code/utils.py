import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def translate_nuc_to_aa(seq):
    '''
    The function receives a nucleotide sequence and translates it to an amino acid sequence
    :param seq: string of nucleotide sequence
    :return: string of amino acid sequence
    '''
    dictionary = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'AGA': 'R', 'AGG': 'R', 'AAC': 'N', 'AAT': 'N',
        'GAC': 'D', 'GAT': 'D', 'TGC': 'C', 'TGT': 'C',
        'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
        'ATT': 'I', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L',
        'CTT': 'L', 'TTA': 'L', 'TTG': 'L', 'AAA': 'K',
        'AAG': 'K', 'ATG': 'M', 'TTC': 'F', 'TTT': 'F',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S',
        'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
        'ACG': 'T', 'ACT': 'T', 'TGG': 'W', 'TAC': 'Y',
        'TAT': 'Y', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
        'GTT': 'V', 'TAA': '*', 'TAG': '*', 'TGA': '*'
    }
    aa_seq = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            if codon not in dictionary:
                aa_seq += '?'
            else:
                aa_seq += dictionary[codon]
    return aa_seq


def fastq_reader(filename):
    '''
    The function receives a fastq file name and reads it into a dictionary
    :param filename: string - fastq file name
    :return: dictionary of the sequences in the file
    '''
    n = 4   # The number of repetitive lines in the file
    lines_keys = ['Name', 'Sequence', 'Optional', 'Quality']
    seq_dict = {}

    with open(filename, 'r') as fastq_file:
        lines = []
        for line in fastq_file:
            lines.append(line.rstrip())
            if len(lines) == n:
                d = dict(zip(lines_keys, lines))
                seq_dict[d['Name'].split(' ')[0]] = d['Sequence'] #the index of the read
                lines = []
    return seq_dict


def sort_seqs_to_mut(seq_dict, info):
    '''
    The function receives dictionary of sequences and returns dictionary sorted by the valid sequences into mutations
    :param  seq_dict: dictionary containing all the sequences
            info: list containing all the desired information:
                start_motif: string- after the motif the desired sequence will begin
                start_seq_pos: string- after this length of nucleotide the sequence starts
                nuc_length: int- the length of the desired sequence
                WTaa: string- The amino acid sequence of the WT
    :return: dictionary of sequences- the key: The name of the seq
                                      the values: aa_seq - amino acid seq, count_diff- number of differences,
                                       diff- string , mut_pos- list of mutation positions
    '''
    start_motif = info[0]
    start_seq_pos = info[1]
    nuc_length = info[2]
    WTaa = info[3]
    mutation_dict = dict()
    for name, seq in seq_dict.items():
        start_motif_pos = seq.find(start_motif)
        position = start_motif_pos + start_seq_pos
        if start_motif_pos != -1 and position + nuc_length <= len(seq):
            relevant_seq = seq[position: position + nuc_length]
            aa_seq = translate_nuc_to_aa(relevant_seq)
            stop_codon = aa_seq.find('*')
            unknown_codon = aa_seq.find('?')
            if stop_codon == -1 and unknown_codon == -1:    # If the seq is valid
                count_diff, diff, mut_pos = compare_seq(seq1=WTaa, seq2=aa_seq)
                mutation_dict[name] = aa_seq, count_diff, diff, mut_pos
    return mutation_dict


def sort_mut_by_number(mut_dict, num_of_mut, col_name):
    '''
    The function receives a dictionary of sequences, a number of mutations and a desired name
    and returns a series of read-counts for the mutations
    :param mut_dict: dictionary of sequences
    :param num_of_mut: int - number of mutation
    :param col_name: string - the desired name of the column
    :return: a series of read-counts for the mutations with the number of mutations
    '''
    df = pd.DataFrame.from_dict(mut_dict, orient='index', columns=[col_name, 'Number of Mutations',
                                                                   'Type of Mutation', 'Position'])
    if num_of_mut <= 7:
        ind_relevant_mut = df['Number of Mutations'] == num_of_mut
    else:
        ind_relevant_mut = df['Number of Mutations'] >= num_of_mut
    relevant_mut = df.loc[ind_relevant_mut, col_name].value_counts()
    return relevant_mut


def compare_seq(seq1, seq2):
    '''
    The function receives two sequences and compares them
    :param seq1: string of sequence
    :param seq2: string of sequence
    :return: the function returns 3 param: count_diff - int - the number of differences found between the 2 sequences,
                                           diff - string - the differences,
                                           position- list of int - the positions of the differences in the sequences
    '''
    count_diff = 0
    diff = ''
    position = []
    if len(seq1) == len(seq2):
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                count_diff += 1
                diff = diff + seq1[i] + str(i + 1) + seq2[i]
                position.append(i + 1)
        if count_diff == 0:
            diff = 'WT'
    return count_diff, diff, position


def is_interface_mut(seq, WTseq, interface_mut):
    '''
    The function receives a sequence, compares it to the WT sequence and returns True if the differences in interface positions
    :param seq: string of the sequence
    :param WTseq: string of the WT sequence
    :param interface_mut: list of int - the interface positions
    :return: True- if the differences are in interface positions
             False - if the differences are not in interface positions
    '''
    _, _, diff = compare_seq(seq1=WTseq, seq2=seq)
    if all(position in interface_mut for position in diff):
        return True
    return False


def oneHot(string):
    aa_dict = {"A": "0", "R": "1", "N": "2", "D": "3", "C": "4", "Q": "5", "E": "6", "G": "7", "H": "8",
               "I": "9", "L": "10", "K": "11", "M": "12", "F": "13", "P": "14", "S": "15", "T": "16", "W": "17",
               "Y": "18", "V": "19"}
    int_encode = [aa_dict[aa] for aa in string]
    one_hot_encode = []
    for i in int_encode:
        l = [0 for _ in range(20)]
        l[int(i)] = 1
        one_hot_encode.append(l)
    return one_hot_encode
