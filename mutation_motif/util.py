import os, gzip, bz2, sys, platform, json
from traceback import format_exc

from pandas import read_json

from cogent import DNA, LoadTable
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import Alignment, DenseAlignment

def spectra_table(table, group_label):
    """returns a table with columns without position information"""
    assert 'direction' in table.Header
    if 'mut' in table.Header:
        # remove redundant category (counts of M == U)
        table = table.filtered("mut=='M'")
    
    columns = ['count', 'direction', group_label]
    table = table.getColumns(columns)
    # so we have a table with counts per direction
    results = []
    group_categories = table.getDistinctValues(group_label)
    filter_template = "direction=='%(direction)s' and %(label)s=='%(category)s'"
    for direction in table.getDistinctValues('direction'):
        start = direction[0]
        for group_category in group_categories:
            condition = dict(direction=direction, label=group_label,
                            category=group_category)
            sub_table = table.filtered(filter_template % condition)
            total = sub_table.summed('count')
            results.append([total, start, direction, group_category])
    result = LoadTable(header=['count', 'start', 'direction', group_label],
                rows=results)
    return result

def get_subtables(table, group_label='direction'):
    """returns [(group, subtable),...] for distinct values of group_label"""
    groups = table.getDistinctValues(group_label)
    tables = []
    for group in groups:
        subtable = table.filtered(lambda x: x == group, columns=group_label)
        tables.append((group, subtable))
    return tables

def dump_loglin_stats(data, outfile_path):
    '''save data in json format to outfile_path'''
    # convert all pandas data frames to json
    saveable = {}
    for position_set in data:
        curr = {}
        new_key = str(position_set)
        for key, value in data[position_set].items():
            if key == 'stats':
                value = data[position_set][key].to_json()
            curr[key] = value
                
        saveable[new_key] = curr
    
    with open(outfile_path, mode="w") as outfile:
        json.dump(saveable, outfile)

def load_loglin_stats(infile_path):
    '''read in data in json format'''
    # convert all 'stats' to pandas data frames
    with open(infile_path) as infile:
        data = json.load(infile)
    
    new_data = {}
    for position_set in data:
        try:
            new_key = eval(position_set)
        except NameError:
            new_key = position_set
        
        new_data[new_key] = {}
        for key, value in data[position_set].items():
            if key == 'stats':
                value = read_json(value)
            new_data[new_key][key] = value
    return new_data


def is_valid(data):
    """returns True if all elements satisfy 0 <= e < 4"""
    return (data >= 0).all() and (data < 4).all()

def load_from_fasta(filename):
    infile = open_(filename)
    parser = MinimalFastaParser(infile)
    seqs = [(n, s) for n, s in parser]
    return DenseAlignment(data=seqs, MolType=DNA)

def array_to_str(data):
    """convert numpy array back to DNA sequence"""
    return [''.join(DNA.Alphabet.fromIndices(v)) for v in data]

def seqs_to_array(d_aln):
    """get input Dense alignment data, transfer them into a numpy array matrix,
    whith accordant numbers for DNA bases
    also
    filter sequences and save just_nuc sequences.
    """
    expect = set([str(s) for s in d_aln.values() if set(s) <= set(DNA)])
    just_bases = just_nucs(d_aln.ArraySeqs)
    return just_bases

def gaps_omitted_alignment(aln):
    """get rid of gaps from alignment sequences"""
    aln = Alignment(aln)
    aln = aln.filtered(lambda x: '-' not in x)
    return aln

def just_nucs(seqs):
    """eliminate sequences containing gaps/Ns,
    along the 1st axis, match each base in seqs to <= 3 element-wise,
    give the indices of those just_nucs seq idx.
    """
    indices, = (seqs <= 3).all(axis=1).nonzero()
    just_bases = seqs.take(indices, axis=0)
    return just_bases

def open_(filename, mode='r'):
    """handles different compression"""
    op = {'gz': gzip.open, 'bz2': bz2.BZ2File}.get(filename.split('.')[-1], open)
    return op(filename, mode)

def abspath(path):
    """returns an expanded, absolute path"""
    return os.path.abspath(os.path.expanduser(path))

def create_path(path):
    """creates dir path"""
    try:
        os.makedirs(path)
    except OSError, e:
        pass

def sibling_path(ref_path, new_dir):
    """returns a new path string with of dirname(ref_path)/basename(new_dir)"""
    ref_path = abspath(ref_path)
    new_dir = abspath(new_dir)
    dirname = os.path.dirname(ref_path)
    basename = os.path.basename(new_dir)
    return os.path.join(dirname, basename)
