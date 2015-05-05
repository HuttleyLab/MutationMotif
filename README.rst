##############
Mutation Motif
##############

This library controls log-linear analysis of neighbourhood base influences on mutation and provides a sequence logo like representation of influences.

*****
Usage
*****

The primary tool is installed as a command line executable, ``mutation_analysis``. It requires a counts table where the table contains counts for a specified flank size (maximum of 2 bases, presumed to be either side of the mutated base). It assumes the counts all reflect a specific mutation direction (e.g. A to G) and that counts from a control distribution are also included.

::

    $ mutation_analysis --help

will display options.