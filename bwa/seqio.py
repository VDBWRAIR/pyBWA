from Bio import SeqIO

import os
import sys
import os.path
import glob

def sffs_to_fastq( sffs, output='sff.fastq' ):
    '''
        Given a list of sffs, concat them into a single fastq
        Nothing created if empty list given

        @raises ValueError if invalid sff file is encountered or output is invalid path
        @param sffs - List of sff file paths to convert and concat into fastq
        @param output - Output fastq file path[Default: sff.fastq]
        @return Path to fastq file created or None if empty list given
    '''
    # Has to be a list
    if not isinstance( sffs, list ):
        raise ValueError( "{} is not a list" )

    # Empty list gets ignored
    if not sffs:
        return

    try:
        # Concat all sequences to output file
        with open( output, 'w' ) as fh:
            for sff in sffs:
                    SeqIO.write( SeqIO.parse( sff, 'sff' ), fh, 'fastq' )
    except (OSError, IOError) as e:
        raise ValueError( "{} is not a valid output file".format(output) )
    except ValueError as e:
        raise ValueError( "{} is not a valid sff file".format(sff) )
    
    return output

def get_reads( dir_path ):
    '''
        Return a list of sff and fastq files in a given dir_path

        @raises ValueError if invalid path given
        @param dir_path - Path to directory of fastq and sff files
        @return list of fastq and sff files found in dir_path. Each has dir_path prefixed to them. Empty list if none found
    '''
    if not os.path.isdir( dir_path ):
        raise ValueError( "{} is not a valid directory".format(dir_path) )
    return glob.glob( os.path.join( dir_path, '*.sff' ) ) + glob.glob( os.path.join( dir_path, '*.fastq' ) ) 
