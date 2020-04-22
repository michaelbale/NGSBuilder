# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:57:14 2020

@author: Michael
"""

import argparse
import csv
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

def _getCode(row, threshold):
    base_dict = { 
        '0': 'A',
        '1': 'C',
        '2': 'G',
        '3': 'T',
        '01': 'M',
        '10': 'M',
        '02': 'R',
        '20': 'R',
        '03': 'W',
        '30': 'W',
        '12': 'S',
        '21': 'S',
        '13': 'Y',
        '31': 'Y',
        '23': 'K',
        '32': 'K',
        '012': 'H',
        '021': 'H',
        '102': 'H',
        '120': 'H',
        '201': 'H',
        '210': 'H',
        '013': 'V',
        '031': 'V',
        '130': 'V',
        '103': 'V',
        '310': 'V',
        '301': 'V',
        '023': 'D',
        '032': 'D',
        '320': 'D',
        '302': 'D',
        '203': 'D',
        '230': 'D',
        '123': 'B',
        '132': 'B',
        '321': 'B',
        '312': 'B',
        '231': 'B',
        '213': 'B'
    }
    row_indSort = np.argsort(row[:-1])[::-1]
    
    for ind in range(0, 3):
        if row[row_indSort[0:(ind+1)]].sum() >= threshold:
            return base_dict[''.join(row_indSort[0:(ind+1)].astype(str))]
    return('N')
    


def ask_row(
    row,
    N_read = 1,
    min_read = 15,
    max_read = 30,
    threshold = 0.8
):
    if row[4] < N_read:
        return '-'
    elif N_read <= row[4] < min_read:
        return 'N'
    elif min_read <= row[4] < max_read:
        return _getCode(row, threshold).lower()
    elif max_read <= row[4]:
        return _getCode(row, threshold).upper()
    else:
        raise Exception('Unexpected value in BFF' )
    


def main(args):
    
    parser = argparse.ArgumentParser(
        prog='generate_consensus.py',
        description="""
        Generater program for reading base frequency file from
        SHIVER and outputting a sequence file.
        """
    )
    parser.add_argument(
        '-i', '--bff', '--input',
        required=True,
        help="""
        Input file as base frequency file. Will assume csv-file.
        Format of file must have A counts, C counts, T counts,
        and G counts as columns 3 through 6 and header row.
        """
    )
    parser.add_argument(
        '-n', '--name',
        required=True,
        help="""
        Name of sequence to output.
        """            
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="""
        Name of output file to store sequence.
        """
    )
    parsed_args = vars(parser.parse_args(args))
    
    with open(parsed_args['bff'], 'r') as fh:
        reader_obj = csv.reader(fh)
        base_file = list(reader_obj)
        
    bases= [sublist[2:6] for sublist in base_file[1:]]
    bases_array = np.array(bases, dtype=int)
    
    np.seterr(all='ignore')
    read_counts = bases_array.sum(axis=1).reshape(bases_array.shape[0], 1)
    bases_freq = (bases_array.T/bases_array.sum(axis=1)).T
    bases_freq = np.nan_to_num(bases_freq)
    
    bases_freq =  np.append(bases_freq, read_counts, 1)
    seq=''
    
    for row in bases_freq:
        seq+= ask_row(row, threshold = 0.8)


if __name__ == "__main__":
    main(sys.argv[1:])
#/end if