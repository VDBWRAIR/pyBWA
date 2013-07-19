# pyBWA

Python libraries to interact with BWA mapper

The goal is to make it very easy to run bwa from any python script


## Simple Example

```python
import bwa
import sys

reference_path = sys.argv[1]
read_path = sys.argv[2]

# Ensure reference is indexed
bwa.index_ref( reference_path )

# Setup and run bwa mem
mem = bwa.BWAMem( reference_path, read_path )
retstat = mem.run( 'myoutput.sai' )

# Check return status
if retstat != 0:
    sys.stderr.write( "Error running bwa" )
```

## Multiple reads and references

```python
import bwa
import sys

reference_path = sys.argv[1]
read_path = sys.argv[2]

# If read_path or reference_path are a directory
# they can contain multiple files that will be concatted
# together into a single file.
# Reads can also be .sff files that will be converted to fastq
reads = bwa.compile_reads( read_path )
refs = bwa.compile_refs( reference_path )

# Ensure reference is indexed
bwa.index_ref( refs )

# Setup and run bwa mem
mem = bwa.BWAMem( refs, reads )
retstat = mem.run( 'myoutput.sai' )

# Check return status
if retstat != 0:
    sys.stderr.write( "Error running bwa" )
```
