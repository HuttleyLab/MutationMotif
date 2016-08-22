from cogent import LoadTable, DNA
from cogent.core.moltype import IUPAC_DNA_ambiguities_complements as COMPLEMENTS

MUTATION_COMPLEMENTS = {'CtoG': 'GtoC',
                        'CtoA': 'GtoT',
                        'AtoT': 'TtoA',
                        'CtoT': 'GtoA',
                        'AtoC': 'TtoG',
                        'AtoG': 'TtoC'}

def _reverse_complement(table):
    '''returns a table with sequences reverse complemented'''
    pos_indices = [i for i, c in enumerate(table.Header) if c.startswith('pos')]
    
    rows = table.getRawData()
    for row in rows:
        # we use the cogent DnaSeq object to do reverse complementing
        seq = DNA.makeSequence(''.join(row[i] for i in pos_indices))
        seq = list(seq.rc())
        for i, index in enumerate(pos_indices):
            row[index] = seq[i]
    if rows:
        new = LoadTable(header=table.Header, rows=rows)
    else:
        new = None
    return new

def add_strand_column(rows, strand):
    for row in rows:
        row.append(strand)
    return rows

def make_strand_symmetric_table(table):
    '''takes a combined counts table and returns a table with reverse complemented seqs
    
    Uses MUTATION_COMPLEMENTS'''
    
    new_data = []
    direction_index = [i for i in range(len(table.Header))
                            if table.Header[i] == 'direction'][0]
    for plus, minus in list(MUTATION_COMPLEMENTS.items()):
        plus_table = table.filtered('direction=="%s"' % plus)
        plus_data = add_strand_column(plus_table.getRawData(), '+')
        new_data.extend(plus_data)
        
        minus_table = table.filtered('direction=="%s"' % minus)
        if minus_table.Shape[0] == 0:
            continue
        minus_table = _reverse_complement(minus_table)
        minus_data = minus_table.getRawData()
        for row in minus_data:
            row[direction_index] = plus
        minus_data = add_strand_column(minus_data, '-')
        new_data.extend(minus_data)
    
    return LoadTable(header=table.Header[:] + ['strand'], rows=new_data)
