import pandas as pd
import numpy as np
from Code.utils import fastq_reader, sort_seqs_to_mut, sort_mut_by_number, is_interface_mut


# TIMP library NGS- variables
WTaa = 'CSCSPVHPQQAFCNADVVIRAKAVSEKEVDSGNDIYGNPIKRIQYEIKQIKMFKGPEKDIEFIYTAPSSAVCGVSLDVGGKKEYLIAGKAEGDGKMHIT'
WTnuc = 'TGCAGCTGCTCCCCGGTGCACCCGCAACAGGCGTTTTGCAATGCAGATGTAGTGATCAGGGCCAAAGCGGTCAGTGAGAAGGAAGTGGACTCTGGAAAC' \
        'GACATCTATGGCAACCCTATCAAGAGGATCCAGTATGAGATCAAGCAGATAAAGATGTTCAAAGGGCCTGAGAAGGATATAGAGTTTATCTACACGGCC' \
        'CCCTCCTCGGCAGTGTGTGGGGTCTCGCTGGACGTTGGAGGAAAGAAGGAATATCTCATTGCAGGAAAGGCCGAGGGGGACGGCAAGATGCACATCACC'
start_motif = 'AGAGATG'
start_motif_ala = 'AGAGAGC'
nuc_length = len(WTnuc)
wanted_position = [4, 35, 38, 68, 71, 97, 99]
aaList = 'ACDEFGHIKLMNPQRSTVWY'


def raw_data_to_mutations_dict(data_path, all_gates, mut_dict, info):
    '''
    The function receives a dictionary of the desired gates and mutations number and returns the dictionary full of found sequences
    :param data_path: string of data path
    :param all_gates: dictionary containing all the gates
    :param mut_dict: dictionary containing all the desired number of mutations
    :param info: list containing all the desired information:
                start_motif: string- after the motif the desired sequence will begin
                start_seq_pos: string- after this length of nucleotide the sequence starts
                nuc_length: int- the length of the desired sequence
                WTaa: string- The amino acid sequence of the WT
    :return: mut_dict: dictionary containing all the desired sequences
    '''
    for name in all_gates.keys():
        seq_dict = fastq_reader(data_path + name + '.fastq')
        all_gates[name] = sort_seqs_to_mut(seq_dict, info)
        # Creating a data frame of the reads count according to the number of mutations
        for num_of_mut in mut_dict.keys():
            mut_dict[num_of_mut] = pd.concat([mut_dict[num_of_mut],
                                              sort_mut_by_number(all_gates[name], num_of_mut, name)], axis=1, sort=True)


def remove_not_valid_variants(mut_dict, numerator_gate, denominator_gate):
    '''
        The function receives a dictionary of sequences and returns the dictionary without not valid sequences
        :param mut_dict: dictionary containing all the sequences
        :param numerator_gate: str of the desired numerator gate
        :param denominator_gate: str of the desired denominator gate
        :return: mut_dict: dictionary containing all the desired sequences
        '''

    for num_of_mut in mut_dict.keys():
        if denominator_gate == 'Gate3':
            mut_dict[num_of_mut]['Gate3'] = mut_dict[num_of_mut]['Gate3'].fillna(0) + mut_dict[num_of_mut]['Gate4'].fillna(0)
        WTnumerator_gate = int(mut_dict[0][numerator_gate])
        WTdenominator_gate = int(mut_dict[0][denominator_gate])
        mut_dict[num_of_mut]['read_count'] = mut_dict[num_of_mut].loc[:, 'Gate1'].fillna(0) + \
                                             mut_dict[num_of_mut].loc[:, 'Gate2'].fillna(0) + \
                                             mut_dict[num_of_mut].loc[:, 'Gate3'].fillna(0) + \
                                             mut_dict[num_of_mut].loc[:,  'sort1'].fillna(0)
        mut_dict[num_of_mut]['log2_read_count'] = np.log2(mut_dict[num_of_mut]['read_count'])
        mut_dict[num_of_mut]['full_seq'] = mut_dict[num_of_mut].index
        mut_dict[num_of_mut]['seq'] = mut_dict[num_of_mut]['full_seq'].str[3:4] + mut_dict[num_of_mut]['full_seq'].str[34:35]\
                                      + mut_dict[num_of_mut]['full_seq'].str[37:38] + mut_dict[num_of_mut]['full_seq'].str[67:68] + \
                                      mut_dict[num_of_mut]['full_seq'].str[70:71] + mut_dict[num_of_mut]['full_seq'].str[96:97] + \
                                      mut_dict[num_of_mut]['full_seq'].str[98:99]
        valid_ind1 = mut_dict[num_of_mut]['full_seq'].apply(is_interface_mut, WTseq=WTaa, interface_mut=wanted_position)
        valid_ind2 = mut_dict[num_of_mut].loc[:, numerator_gate].fillna(0) != 0
        valid_ind3 = mut_dict[num_of_mut].loc[:, denominator_gate].fillna(0) != 0
        not_valid = []
        for i in range(len(valid_ind1)):
            if (not valid_ind1.values[i]) or (not valid_ind2.values[i]) or (not valid_ind3.values[i]):
                not_valid.append(mut_dict[num_of_mut].index[i])
        mut_dict[num_of_mut].drop(not_valid, inplace=True)
        mut_dict[num_of_mut]['ER'] = (mut_dict[num_of_mut].loc[:, numerator_gate].fillna(0) / WTnumerator_gate) / \
                                     (mut_dict[num_of_mut].loc[:, denominator_gate].fillna(0) / WTdenominator_gate)
        mut_dict[num_of_mut]['log2_ER'] = np.log2(mut_dict[num_of_mut]['ER'])


def mutations_dict_to_df(mut_dict):
    '''
        The function receives a dictionary with valid sequences and returns a data frame containing all variants
        :param mut_dict: dictionary containing all the desired sequences
        :param numerator_gate: str of the desired numerator gate
        :param denominatorr_gate: str of the desired denominator gate
        :return: mut_total_data: dictionary containing all variants information
    '''
    mut_total_data = pd.DataFrame()
    for num_of_mut in mut_dict.keys():
        mut_dict[num_of_mut]['mutations number'] = num_of_mut
        mut_total = pd.concat([mut_dict[num_of_mut].loc[:, 'seq'],
                               mut_dict[num_of_mut].loc[:, 'mutations number'],
                               mut_dict[num_of_mut].loc[:, 'Gate1'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'Gate2'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'Gate3'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'sort1'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'read_count'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'log2_read_count'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'ER'].fillna(0),
                               mut_dict[num_of_mut].loc[:, 'log2_ER'].fillna(0)],
                              axis=1, sort=True)
        mut_total_data = pd.concat([mut_total_data, mut_total], sort=False)
    return mut_total_data



if __name__ == '__main__':
    data_path = 'Data/'
    all_gates_dict = {'Gate1': {}, 'Gate2': {}, 'Gate3': {}, 'Gate4': {}, 'sort1': {}}
    all_gates_ala_dict = {'Gate1': {}, 'Gate2': {}, 'Gate3': {}, 'Gate4': {}, 'sort1': {}}

    mut_dictionary = {0: pd.DataFrame(), 1: pd.DataFrame(), 2: pd.DataFrame(), 3: pd.DataFrame(), 4: pd.DataFrame(),
                      5: pd.DataFrame(), 6: pd.DataFrame(), 7: pd.DataFrame(), 8: pd.DataFrame()}  # 8-> 8 or bigger
    mut_dictionary_ala = {0: pd.DataFrame(), 1: pd.DataFrame(), 2: pd.DataFrame(), 3: pd.DataFrame(), 4: pd.DataFrame(),
                          5: pd.DataFrame(), 6: pd.DataFrame(), 7: pd.DataFrame(), 8: pd.DataFrame()}  # 8-> 8 or bigger

    raw_data_to_mutations_dict(data_path, all_gates_dict, mut_dictionary, [start_motif, len(start_motif) - 2 , nuc_length, WTaa])
    raw_data_to_mutations_dict(data_path, all_gates_ala_dict, mut_dictionary_ala,
                               [start_motif_ala, len(start_motif_ala) + 1, nuc_length, WTaa])
    remove_not_valid_variants(mut_dictionary, numerator_gate='Gate1', denominator_gate='Gate2')
    remove_not_valid_variants(mut_dictionary_ala, numerator_gate='Gate2', denominator_gate='Gate3')
    raw_data = mutations_dict_to_df(mut_dictionary)
    raw_data_ala = mutations_dict_to_df(mut_dictionary_ala)
    raw_data.to_csv(f'{data_path}/All_variant_no_ala.csv')
    raw_data_ala.to_csv(f'{data_path}/All_variant_ala.csv')

