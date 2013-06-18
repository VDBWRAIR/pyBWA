from Bio import SeqIO

import logging
import re
from subprocess import Popen, PIPE, check_output
import tempfile
import os
import os.path
import sys

logger = logging.getLogger( __name__ )

def which_bwa( ):
    ''' Return output of which bwa '''
    return check_output( ['which', 'bwa'] ).strip()

def bwa_usage():
    return check_output( ['bwa', 'mem'] ).strip()

class BWA( object ):
    # Options that are required
    REQUIRED_OPTIONS = ['bwa_path', 'command']
    # regex to detect usage output
    USAGE_REGEX = re.compile( 'Usage:\s*bwa' )

    def __init__( self, *args, **kwargs ):
        '''
            Base class to run mem or aln options

            BWA Options:
                kwargs represent any option that has a dash before where the option
                    with the dash is the key and the value is the value

                args represent the required options(handled in subclasses)
                the first of these should be one of the main commands(mem, aln...)
            
            Class Options:
                bwa_path as a kwarg that specifies the path to the bwa executable
        '''
        # Save args, kwargs for parsing
        self.kwargs = kwargs
        self.args = list( args )
        # Options list
        self.options = []
        # Required options values. If you zip REQUIRED_OPTIONS and required_options_values you will get a
        #  mapping of k,v pairs
        self.required_options_values = []
        # Parse and remove required_options
        self.required_options()
        # Setup options from the rest of the kwargs
        self.compile_bwa_options()
        # This needs to be implemented in subclass
        self.required_args()

    def required_args( self ):
        ''' Sets self.args '''
        raise NotImplementedError( "This class is intended to be subclassed and not instantiated directly" )
    
    def required_options( self ):
        ''' Parse out REQUIRED_OPTIONS from kwargs and set them in self.required_options_values '''
        try:
            # Build up the values in order they appear in REQUIRED_OPTIONS
            for op in self.REQUIRED_OPTIONS:
                self.required_options_values.append( self.kwargs[op] )
                # No longer need this in kwargs
                del self.kwargs[op]
        except KeyError as e:
            # Detects if a parameter is missing
            raise ValueError( "{} is a required parameter".format(op) )

    def compile_bwa_options( self ):
        '''
            Convert kwargs to options list
            Assumes REQUIRED_OPTIONS are not part of this list
            @returns list of options aka [k1, v1, k2, v2...]
        '''
        # Build up self.options from kwargs
        for op, val in self.kwargs.items():
            # Append dash to option
            self.options.append( '-'+op )
            # Options should all be strings(just being passed to command line anyways)
            val = str(val)
            # True false values only have option
            if val.lower() not in ('true','false'):
                self.options.append( val )

    def bwa_return_code( self, output ):
        '''
            Parse stderr output to find if it executed without errors
            Since it seems that bwa does not set return codes we have to parse
            stderr output instead

            If the following regex is found then the Usage statement was printed which indicates a failure of 
             one of the options:
                ^Usage:\s+bwa
            
            Subclasses need to implement this as well and call this but they need to parse the rest of the output
            if this returns success in order to tell if the algorithm ran correctly or not
        
            @returns 0 if no usage was found, 1 if usage was found
        '''
        # Search the output
        m = self.USAGE_REGEX.search( output )

        # If there is a match return 1
        if m:
            logger.warning( "BWA Returned Usage help instead of running. This could indicate an error." )
            return 2
        # Otherwise return 0
        return 0

    def run( self, output_file='bwa.sai' ):
        ''' Call run_bwa lazily '''
        return self.run_bwa( self.required_options_values, self.options, self.args, output_file )

    def run_bwa( self, required_options, options_list, args_list, output_file='bwa.sai' ):
        '''
            required_options - Should correspond to self.REQUIRED_OPTIONS
            options_list - Full options for bwa as a list (ex. ['mem', '-t', '2'])
            args_list - Required arguments that come after options
            output_file - Output location for stdout

            @returns 0 for success, 2 if incorrect options

            Subclass implementation should return 1 for any other failures
        '''
        if not os.path.exists( required_options[0] ):
            raise ValueError( "{} is not a valid bwa path".format( required_options[0] ) )

        # Run bwa
        with open( output_file, 'wb' ) as fh:
            cmd = required_options + options_list + args_list
            logger.info( "Running {}".format( " ".join( cmd ) ) )
            p = Popen( cmd, stdout=fh, stderr=PIPE )

            # Get the output
            stdout, stderr = p.communicate()
        logger.debug( "STDERR: {}".format(stderr) )

        # Parse the status
        return self.bwa_return_code( stderr )

    def validate_indexed_fasta( self, fastapath ):
        ''' Ensure fastapath is valid path and is bwa index'd '''
        if not os.path.exists( fastapath + '.bwt' ):
            raise ValueError( "{} does not have an index".format(fastapath) )
        if not os.path.exists( fastapath ):
            raise ValueError( "{} does not exist".format(fastapath) )

    def validate_input( self, inputpath ):
        if not os.path.exists( inputpath ):
            raise ValueError( "{} is not a valid input file".format(inputpath) )

    def reads_in_file( self, filename ):
        ftype = 'fasta'
        with open( filename ) as fh:
            firstline = fh.readline()
            if firstline.startswith( '>' ):
                ftype = 'fasta'
            elif firstline.startswith( '@' ):
                ftype = 'fastq'
        return sum( [1 for seq in SeqIO.parse( filename, ftype )] )

class BWAIndex( BWA ):
    def __init__( self, *args, **kwargs ):
        ''' Injects index command and runs super '''
        kwargs['command'] = 'index'
        super( BWAIndex, self ).__init__( *args, **kwargs )

    def required_args( self ):
        '''
            Index only requires an input fasta file to index
            Validate that it is an actual fasta file
        '''
        if len( self.args ) != 1:
            raise ValueError( "bwa index needs only 1 parameter" )
        self.validate_input( self.args[0] )
        if self.reads_in_file( self.args[0] ) == 0:
            raise ValueError( "{} is not a valid file to index".format(self.args[0]) )

    def bwa_return_code( self, stderr ):
        ''' 
            Missing file:
                [bwa_index] fail to open file 'bob'. Abort!

            bwa index runs successfully pretty much no matter what
        '''
        if '[bwa_index] fail to open file' in stderr:
            return 1
        return 0

    def run( self ):
        '''
            Call super and then remove output file
            hackish
        '''
        ret = super( BWAIndex, self ).run( 'removeme.sai' )
        os.unlink( 'removeme.sai' )
        return ret

class BWAMem( BWA ):
    def __init__( self, *args, **kwargs ):
        ''' Injects mem command and runs super '''
        kwargs['command'] = 'mem'
        super( BWAMem, self ).__init__( *args, **kwargs )

    def required_args( self ):
        '''
            Mem requires 2 args
                db.prefix - Indexed reference genome
                reads.fq - fastq file to map reads from

            It also has one optional argument
                mates.fq - mates fastq file

            Ensures args are valid
        '''
        # Validate args
        # First argument has to be a valid indexed fasta
        self.validate_args( self.args )

    def validate_args( self, args ):
        if len( args ) > 3:
            raise ValueError( "Too many arguments supplied to BWAMem" )
        elif len( args ) < 2:
            raise ValueError( "Too few arguments supplied to BWAMem" )
        else:
            self.validate_indexed_fasta( self.args[0] )
            self.validate_input( self.args[1] )
            if len( args ) == 3:
                self.validate_input( self.args[2] )

    def bwa_return_code( self, output ):
        '''
            Just make sure bwa output has the following regex and make sure the read \d counts
            up to how many sequence lines there are
            Should end with Version line

            Example Line:
                [M::main_mem] read 100 sequences (111350 bp)...
                [main] Version: 0.7.4-r385
        '''
        read_line_pat = '\[M::main_mem\] read (\d+) sequences \((\d+) bp\)...'
        cpat = re.compile( read_line_pat )

        total_reads = 0
        total_bp = 0
        counts = cpat.findall( output )
        for reads, bps in counts:
            total_reads += int( reads )
            total_bp += int( bps )

        # Count num of read sequences
        expected_reads = self.reads_in_file( self.args[1] )
        # If mates file was given count them too
        if len( self.args ) == 3:
            expected_reads += self.reads_in_file( self.args[2] )

        # No lines found in input file?
        if expected_reads == 0:
            return 1

        if total_reads != expected_reads:
            logger.warning( "Expecting BWA to process {} reads but processed {}".format(expected_reads, total_reads) )
            return 1

        return super( BWAMem, self ).bwa_return_code( output )
