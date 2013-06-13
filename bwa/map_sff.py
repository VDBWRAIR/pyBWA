from Bio import SeqIO

from argparse import ArgumentParser
import os
import sys
import os.path
import glob

def main():
    pass

def parse_args( ):
    parser = ArgumentParser( description='Map sff files using bwa mem' )

    parser.add_argument( dest='input', help='Sff file or directory of sff files to concat and map' )
    parser.add_argument( dest='index', help='Reference index to map against' )

def sff_to_fastq( sfffiles ):
    '''
        Convert sff files to fastq
    '''

def index_file( filetoindex ):
    '''
        Run bwa index on filetoindex
    '''

def map_sffs( sffs, ref ):
    '''
        Given a list of sffs, concat them into a single fastq and
        then run bwa mem on them
        Index ref as well
    '''
