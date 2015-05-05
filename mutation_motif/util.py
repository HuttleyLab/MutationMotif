import os, gzip, bz2, sys, platform
from traceback import format_exc
from socket import gethostname
import logging, hashlib

from cogent import DNA
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import Alignment, DenseAlignment

class CachingLogger(object):
    """stores log messages until a log filename is provided"""
    def __init__(self, log_file_path=None):
        super(CachingLogger, self).__init__()
        self.log_file_path = log_file_path
        self._started = False
        self._messages = []
    
    def write(self, msg):
        """writes a log message"""
        
        if not self._started:
            self._messages.append(msg)
        else:
            logging.info(msg)
        
        if not self._started and self.log_file_path:
            # start the logger and flush the message cache
            set_logger(self.log_file_path)
            for m in self._messages:
                logging.info(m)
            
            self._messages = []
            self._started = True
            
    

def set_logger(logfile_name, level=logging.DEBUG):
    """setup logging"""
    logging.basicConfig(filename=logfile_name, filemode='w', level=level,
                format='%(asctime)s\t%(levelname)s\t%(message)s',
                datefmt="%Y-%m-%d %H:%M:%S")
    logging.info('system_details: system=%s:python=%s' % (platform.version(),
                                                  platform.python_version()))

def get_file_hexdigest(filename):
    '''returns the md5 hexadecimal checksum of the file'''
    # from http://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python
    with open(filename) as infile:
        md5 = hashlib.md5()
        while True:
            data = infile.read(128)
            if not data:
                break
            
            md5.update(data)
    return md5.hexdigest()




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
