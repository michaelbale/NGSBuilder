# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


#Import statements
import argparse
import sys
import warnings
import yaml

#Custom Library Import Statement
from Bio import SeqIO


def main(args):
    
    
    """
    Function to parse fasta and yaml file for renaming sequences based on yaml
    
    Params:
        --fasta, type=string: Fasta file of sequences to be renamed;
            N.B. input file style can be dynamic with @param --style
        --YAML, type=string: YAML-style file describing python dict of index
            names and sequence names for renaming. 
    
    Optional Params:
        --style, type=string: Denoting the style of the input file -- assumed
            to be fasta file. If @param --fasta is a fasta file, this flag is 
            not needed
        --output, type=string: Custom name of output file where renamed 
            sequences will be saved. Default behavior is to do 
            ${FASTA}_renamed.fas
            
    Logic:
        YAML[INDEX] = RENAME where
            RENAME will be the new name of input sequence 
            index_INDEX_shiverBuild_consensus
        if sequences.id =~ INDEX
            sequences.id = YAML[INDEX] //Renamed
    
    Warnings:
        Loop keeps track of number of renamed sequences; if the number of
            renames is less than the number of supplied sequences, a warning
            message will be displayed
    """
    #Begin parser
    parser = argparse.ArgumentParser(
        prog='rename_seqs.py',
        description='''
        Renames sequences in supplied Fasta file based on YAML-style dictionary
        from supplied names file.
        ''',
        epilog='''
        Any questions can be directed to michaeljbale93@gmail.com or 
        michael.bale@nih.gov
        ''',
        add_help=True
    )
    parser.add_argument(
        '-i',
        '--fasta',
        '-f',
        required=True,
        help='Input fasta file'
    )
    parser.add_argument(
        '-n',
        '-y',
        '--YAML',
        '--yaml',
        '--names',
        required=True,
        help='Input yaml file containing dictionary of indexes and new names'
    )
    parser.add_argument(
        '-s',
        '--style',
        '--format',
        required=False,
        default="fasta",
        help="Optional style for input file; output style will always be fasta"
    )
    parser.add_argument(
        '-o',
        '--output',
        '--outfile',
        required=False,
        help='''
        Optional Arg for naming output file. Default behavior will be to have
        output file as ${INPUTFILE}_renamed.fas
        '''
    )
    #/end Parser
    
    #Parse and handle supplied args
    parsed_args = vars(parser.parse_args(args))
    
    parsed_args['output'] = \
        parsed_args['fasta'] + "_renamed.fas" if \
            parsed_args['output'] is None \
        else parsed_args['output']
    
    #IO Block
    try:
        with open(parsed_args['YAML'], 'r') as yml:
            names_dict = yaml.load(yml, Loader=yaml.FullLoader)
    except FileNotFoundError as err:
        print("Please supply valid YAML File; {0} does not exist".format(err.args))
    #/end try-except
    
    try:
        sequences = list(
            SeqIO.parse(parsed_args['fasta'], parsed_args['style'])
        )
        
    except FileNotFoundError as err:
        print("Please supply valid Sequence File; {0} does not exist".format(err.args))    
    #/end try-except
    
    count  = 0
    
    for key in names_dict.keys():
        for seq in sequences:
            if key in seq.id:
                seq.id = names_dict[key]
                seq.name=""
                seq.description=""
                count += 1
            #/end if
        #/next seq
    #/next key
                
    if count != len(sequences):
        warnings.warn('Not all Sequences were renamed')
        print(str(count))
    #/end if
        
    SeqIO.write(sequences, parsed_args['output'], 'fasta')
#/end def main()




if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
#/end rename_seqs.py