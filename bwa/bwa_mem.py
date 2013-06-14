from argparse import ArgumentParser
from subprocess import check_output

import bwa
import seqio

import logging
import os.path
import sys

logging.basicConfig( level=logging.DEBUG )
logger = logging.getLogger( os.path.splitext( __file__ )[0] )

def compile_reads( reads ):
    '''
        Compile all given reads from directory of reads or just return reads if it is fastq
        If reads is sff file then convert to fastq

        @param reads - Directory/file of .fastq or .sff
        @return fastq with all reads from reads
    '''
    return seqio.sffs_to_fastq( seqio.get_reads( reads ) )

def main():
    args = parse_args().__dict__

    ref_file = args['index']
    del args['index']
    read_path = compile_reads( args['reads'] )
    del args['reads']
    mates_path = args['mates']
    del args['mates']
    output_file = args['output']
    del args['output']
    args['bwa_path'] = bwa.which_bwa()

    ret = 1
    try:
        ret = bwa.BWAIndex( ref_file, bwa_path=bwa.which_bwa() ).run()
    except ValueError as e:
        logger.error( e )

    if ret != 0:
        logger.error( "Error running bwa index on {}".format( ref_file ) )
        sys.exit( ret )
    else:
        logger.info( "bwa index ran on {}".format(ref_file) )
        
    ret = 1
    try:
        if mates_path:
            ret = bwa.BWAMem( ref_file, read_path, mates_path, **args ).run( output_file )
        else:
            ret = bwa.BWAMem( ref_file, read_path, **args ).run( output_file )
    except ValueError as e:
        logger.error( str(e) )

    if ret != 0:
        logger.error( "Error running bwa mem" )
        sys.exit( ret )
    else:
        logger.info( "Finished running bwa" )
        logger.info( "Ref: {} Input: {} Mates: {}".format(ref_file,read_path,mates_path) )
        logger.info( "Options: {}".format(args) )

def parse_args( ):
    parser = ArgumentParser( epilog='Python wrapper around bwa mem' )

    parser.add_argument( '-t', type=int, default=1, help='number of threads' )
    parser.add_argument( '-k', type=int, default=19, help='minimum seed length' )
    parser.add_argument( '-w', type=int, default=100, help='band width for banded alignment' )
    parser.add_argument( '-d', type=int, default=100, help='off-diagonal X-dropoff' )
    parser.add_argument( '-r', type=float, default=1.5, help='look for internal seeds inside a seed longer than {-k} * FLOAT' )
    parser.add_argument( '-c', type=int, default=10000, help='skip seeds with more than INT occurrences' )
    parser.add_argument( '-S', type=bool, default=False, help='skip mate rescue' )
    parser.add_argument( '-P', type=bool, default=False, help='skip pairing; mate rescue performed unless -S also in use' )
    parser.add_argument( '-A', type=int, default=1, help='score for a sequence match' )
    parser.add_argument( '-B', type=int, default=4, help='penalty for a mismatch' )
    parser.add_argument( '-O', type=int, default=6, help='gap open penalty' )
    parser.add_argument( '-E', type=int, default=1, help='gap extension penalty; a gap of size k cost {-O} + {-E}*k' )
    parser.add_argument( '-L', type=int, default=5, help='penalty for clipping' )
    parser.add_argument( '-U', type=int, default=17, help='penalty for an unpaired read pair' )
    parser.add_argument( '-p', type=bool, default=False, help='first query file consists of interleaved paired-end sequences' )
    parser.add_argument( '-R', type=str, default=None, help='read group header line such as \'@RG\tID:foo\tSM:bar\'' )
    parser.add_argument( '-v', type=int, default=3, help='verbose level: 1=error, 2=warning, 3=message, 4+=debugging' )
    parser.add_argument( '-T', type=int, default=30, help='minimum score to output' )
    parser.add_argument( '-a', type=bool, default=False, help='output all alignments for SE or unpaired PE' )
    parser.add_argument( '-C', type=bool, default=False, help='append FASTA/FASTQ comment to SAM output' )
    parser.add_argument( '-H', type=bool, default=False, help='hard clipping' )
    parser.add_argument( '-M', type=bool, default=False, help='mark shorter split hits as secondary (for Picard/GATK compatibility)' )
    parser.add_argument( '--output', metavar='output_file', default='bwa.sai', help='Output file to put sam output in[Default:bwa.sai]' )

    parser.add_argument( dest='index', help='Reference location' )
    parser.add_argument( dest='reads', help='Read or directory of reads to be mapped(.fastq and .sff supported)' )
    parser.add_argument( dest='mates', nargs='?', help='mates file' )

    return parser.parse_args()

if __name__ == '__main__':
    main()
