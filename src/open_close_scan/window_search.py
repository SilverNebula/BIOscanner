import numpy as np
import io
import pandas as pd

source_seq = 'ATGGAGACCCCGTCCCAGCGGCGCGCCACCCGCAGCGGGGCGCAGGCCAGCTCCACTCCGCTGTCGCCCACCCGCATCACCCGGCTGCAGGAGAAGGAGGACCTGCAGGAGCTCAATGATCGCTTGGCGGTCTACATCGACCGTGTGCGCTCGCTGGAAACGGAGAACGCAGGGCTGCGCCTTCGCATCACCGAGTCTGAAGAGGTGGTCAGCCGCGAGGTGTCCGGCATCAAGGCCGCCTACGAGGCCGAGCTCGGGGATGCCCGCAAGACCCTTGACTCAGTAGCCAAGGAGCGCGCCCGCCTGCAGCTGGAGCTGAGCAAAGTGCGTGAGGAGTTTAAGGAGCTGAAAGCGCGCAATACCAAGAAGGAGGGTGACCTGATAGCTGCTCAGGCTCGGCTGAAGGACCTGGAGGCTCTGCTGAACTCCAAGGAGGCCGCACTGAGCACTGCTCTCAGTGAGAAGCGCACGCTGGAGGGCGAGCTGCATGATCTGCGGGGCCAGGTGGCCAAGCTTGAGGCAGCCCTAGGTGAGGCCAAGAAGCAACTTCAGGATGAGATGCTGCGGCGGGTGGATGCTGAGAACAGGCTGCAGACCATGAAGGAGGAACTGGACTTCCAGAAGAACATCTACAGTGAGGAGCTGCGTGAGACCAAGCGCCGTCATGAGACCCGACTGGTGGAGATTGACAATGGGAAGCAGCGTGAGTTTGAGAGCCGGCTGGCGGATGCGCTGCAGGAACTGCGGGCCCAGCATGAGGACCAGGTGGAGCAGTATAAGAAGGAGCTGGAGAAGACTTATTCTGCCAAGCTGGACAATGCCAGGCAGTCTGCTGAGAGGAACAGCAACCTGGTGGGGGCTGCCCACGAGGAGCTGCAGCAGTCGCGCATCCGCATCGACAGCCTCTCTGCCCAGCTCAGCCAGCTCCAGAAGCAGCTGGCAGCCAAGGAGGCGAAGCTTCGAGACCTGGAGGACTCACTGGCCCGTGAGCGGGACACCAGCCGGCGGCTGCTGGCGGAAAAGGAGCGGGAGATGGCCGAGATGCGGGCAAGGATGCAGCAGCAGCTGGACGAGTACCAGGAGCTTCTGGACATCAAGCTGGCCCTGGACATGGAGATCCACGCCTACCGCAAGCTCTTGGAGGGCGAGGAGGAGAGGCTACGCCTGTCCCCCAGCCCTACCTCGCAGCGCAGCCGTGGCCGTGCTTCCTCTCACTCATCCCAGACACAGGGTGGGGGCAGCGTCACCAAAAAGCGCAAACTGGAGTCCACTGAGAGCCGCAGCAGCTTCTCACAGCACGCACGCACTAGCGGGCGCGTGGCCGTGGAGGAGGTGGATGAGGAGGGCAAGTTTGTCCGGCTGCGCAACAAGTCCAATGAGGACCAGTCCATGGGCAATTGGCAGATCAAGCGCCAGAATGGAGATGATCCCTTGCTGACTTACCGGTTCCCACCAAAGTTCACCCTGAAGGCTGGGCAGGTGGTGACGATCTGGGCTGCAGGAGCTGGGGCCACCCACAGCCCCCCTACCGACCTGGTGTGGAAGGCACAGAACACCTGGGGCTGCGGGAACAGCCTGCGTACGGCTCTCATCAACTCCACTGGGGAAGAAGTGGCCATGCGCAAGCTGGTGCGCTCAGTGACTGTGGTTGAGGACGACGAGGATGAGGATGGAGATGACCTGCTCCATCACCACCACGGCTCCCACTGCAGCAGCTCGGGGGACCCCGCTGAGTACAACCTGCGCTCGCGCACCGTGCTGTGCGGGACCTGCGGGCAGCCTGCCGACAAGGCATCTGCCAGCGGCTCAGGAGCCCAGGTGGGCGGACCCATCTCCTCTGGCTCTTCTGCCTCCAGTGTCACGGTCACTCGCAGCTACCGCAGTGTGGGGGGCAGTGGGGGTGGCAGCTTCGGGGACAATCTGGTCACCCGCTCCTACCTCCTGGGCAACTCCAGCCCCCGAACCCAGAGCCCCCAGAACTGCAGCATCATGTAA'
char_map = {'A': 'T',
            'G': 'C',
            'C': 'G',
            'T': 'A'}


def get_rev(seq):
    res = list(seq)[::-1]
    for i in range(len(res)):
        res[i] = char_map[res[i]]
    res = ''.join(res)
    # print('before: {}  after: {}'.format(seq, res))
    return res


if __name__ == '__main__':
    input_file = 'D:/CODE/BIO/dataset/ASO/LMNA_gampers_from_IONIS.csv'
    df = pd.read_csv(input_file)
    df.columns = ['id', 'seq', 'inhibition', 'na']
    df = df[['id', 'seq', 'inhibition']]
    df = df[df['inhibition'] > 70]
    res = ''
    for idex, row in df.iterrows():
        rev = get_rev(row['seq'])
        fd = source_seq.find(rev)
        if fd == -1:
            continue
        res += '{}-{},'.format(fd, fd+20)
    print(res)