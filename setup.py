import os
import os.path
from setuptools import setup, find_packages

from fnmatch import fnmatch
import subprocess
import sys
import glob
from bwa.install import install_bwa

from bwa._version import __version__

# Install bwa into bin directory so it will be copied with all of the other
# scripts inside of bin
install_bwa( 'bin/' )

# Utility function to read the README file.
# Used for the long_description. It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def scripts( ):
    return [os.path.join( 'bin', f ) for f in os.listdir( 'bin' ) \
        if not fnmatch( f, '*.swp' ) and not fnmatch( f, '*.pyc' )]

def git_branch():
    ''' Return the current checked out branch name '''
    try:
        output = subprocess.check_output( ['git', 'branch'] ).splitlines()
    except:
        print "unable to get git branch"
        return ""

    # Get the line that the astriks is in
    branch = [x for x in output if '*' in x][0]
    branch = branch.replace( '*', '' ).strip()
    # Only return branches other than master
    if branch != 'master':
        return branch
    else:
        return ''

setup(
    name = "pyBWA",
    version = __version__,
    author = "Tyghe Vallard",
    author_email = "vallardt@gmail.com",
    description = ("Python wrapper for bwa mapper"),
    keywords = "bwa walter reed research python library",
    url = "https://github.com/VDBWRAIR/pyBWA",
    packages = find_packages(),
    scripts = scripts(),
    data_files = [
    ],
    setup_requires = [
    ],
    install_requires = [
        'biopython'
    ],
    tests_require = [
        'nose',
        'mock',
    ]
)

