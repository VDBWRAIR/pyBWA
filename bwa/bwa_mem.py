from argparse import ArgumentParser
from subprocess import check_output

import bwa
import seqio

import logging
import os.path
import sys
import fnmatch
import glob

logging.basicConfig( level=logging.DEBUG )
logger = logging.getLogger( os.path.basename( os.path.splitext( __file__ )[0] ) )

def compile_reads( reads ):
    '''
        Compile all given reads from directory of reads or just return reads if it is fastq
        If reads is sff file then convert to fastq

        @param reads - Directory/file of .fastq or .sff
        @return fastq with all reads from reads
    '''
    if os.path.isdir( reads ):
        reads = seqio.get_reads( reads )

    logger.info( "Concatting and Converting {} to fastq".format(reads) )
    return seqio.sffs_to_fastq( reads )

def compile_refs( refs ):
    '''
        Compile all given refs into a single file to be indexed

        @param refs - Directory/file of fasta formatted files
        @return path to concatted indexed reference file
    '''
    ref_files = []
    ref_extensions = ('.fa', '.fasta', '.fna', '.fas')

    if os.path.isdir( refs ):
        logger.info( "Compiling and concatting refs inside of {}".format(refs) )
        files = glob.glob( os.path.join( refs, '*' ) )
        logger.debug( "All files inside of {}: {}".format( files, refs ) )
        ref_files = [f for f in files if os.path.splitext(f)[1] in ref_extensions]
        logger.debug( "Filtering files down to only files with extensions in {}".format(ref_extensions) )
        logger.debug( "Filtered files to concat: {}".format( ref_files ) )
        try:
            seqio.concat_files( ref_files, 'reference.fa' )
        except (OSError,IOError,ValueError) as e:
            logger.error( "There was an error with the references in {}".format(refs) )
            logger.error( str( e ) )
            sys.exit(1)
        return 'reference.fa'
    else:
        return refs

def main():
    args = parse_args().__dict__

    ref_file = compile_refs( args['index'] )
    del args['index']

    read_path = compile_reads( args['reads'] )
    del args['reads']

    mates_path = args['mates']
    del args['mates']
    output_file = args['output']

    del args['output']
    args['bwa_path'] = bwa.which_bwa()

    ret = 1
    ret = bwa.index_ref( ref_file )

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
        logger.debug( "Ref: {} Input: {} Mates: {}".format(ref_file,read_path,mates_path) )
        logger.debug( "Options: {}".format(args) )

def parse_args( ):
    parser = ArgumentParser( epilog='Python wrapper around bwa mem' )

    parser.add_argument( '-t', help='number of threads' )
    parser.add_argument( '-k', help='minimum seed length' )
    parser.add_argument( '-w', help='band width for banded alignment' )
    parser.add_argument( '-d', help='off-diagonal X-dropoff' )
    parser.add_argument( '-r', help='look for internal seeds inside a seed longer than {-k} * FLOAT' )
    parser.add_argument( '-c', help='skip seeds with more than INT occurrences' )
    parser.add_argument( '-S', help='skip mate rescue' )
    parser.add_argument( '-P', help='skip pairing; mate rescue performed unless -S also in use' )
    parser.add_argument( '-A', help='score for a sequence match' )
    parser.add_argument( '-B', help='penalty for a mismatch' )
    parser.add_argument( '-O', help='gap open penalty' )
    parser.add_argument( '-E', help='gap extension penalty; a gap of size k cost {-O} + {-E}*k' )
    parser.add_argument( '-L', help='penalty for clipping' )
    parser.add_argument( '-U', help='penalty for an unpaired read pair' )
    parser.add_argument( '-p', help='first query file consists of interleaved paired-end sequences' )
    parser.add_argument( '-R', help='read group header line such as \'@RG\tID:foo\tSM:bar\'' )
    parser.add_argument( '-v', help='verbose level: 1=error, 2=warning, 3=message, 4+=debugging' )
    parser.add_argument( '-T', help='minimum score to output' )
    parser.add_argument( '-a', help='output all alignments for SE or unpaired PE' )
    parser.add_argument( '-C', help='append FASTA/FASTQ comment to SAM output' )
    parser.add_argument( '-H', help='hard clipping' )
    parser.add_argument( '-M', help='mark shorter split hits as secondary (for Picard/GATK compatibility)' )
    parser.add_argument( '--output', metavar='output_file', default='bwa.sai', help='Output file to put sam output in[Default:bwa.sai]' )

    parser.add_argument( dest='index', help='Reference location' )
    parser.add_argument( dest='reads', help='Read or directory of reads to be mapped(.fastq and .sff supported)' )
    parser.add_argument( dest='mates', nargs='?', help='mates file' )

    return parser.parse_args()

if __name__ == '__main__':
    main()
