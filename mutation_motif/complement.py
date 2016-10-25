from cogent3 import LoadTable, DNA

MUTATION_COMPLEMENTS = {'CtoG': 'GtoC',
                        'CtoA': 'GtoT',
                        'AtoT': 'TtoA',
                        'CtoT': 'GtoA',
                        'AtoC': 'TtoG',
                        'AtoG': 'TtoC'}


def _reverse_complement(table):
    '''returns a table with sequences reverse complemented'''
    pos_indices = [i for i, c in enumerate(
        table.header) if c.startswith('pos')]

    rows = table.tolist()
    for row in rows:
        # we use the cogent3 DnaSeq object to do reverse complementing
        seq = DNA.make_seq(''.join(row[i] for i in pos_indices))
        seq = list(seq.rc())
        for i, index in enumerate(pos_indices):
            row[index] = seq[i]
    if rows:
        new = LoadTable(header=table.header, rows=rows)
    else:
        new = None
    return new


def add_strand_column(rows, strand):
    for row in rows:
        row.append(strand)
    return rows


def make_strand_symmetric_table(table):
    '''takes a combined counts table and returns a table with reverse
    complemented seqs

    Uses MUTATION_COMPLEMENTS'''

    new_data = []
    direction_index = [i for i in range(len(table.header))
                       if table.header[i] == 'direction'][0]
    for plus, minus in list(MUTATION_COMPLEMENTS.items()):
        plus_table = table.filtered('direction=="%s"' % plus)
        plus_data = add_strand_column(plus_table.tolist(), '+')
        new_data.extend(plus_data)

        minus_table = table.filtered('direction=="%s"' % minus)
        if minus_table.shape[0] == 0:
            continue
        minus_table = _reverse_complement(minus_table)
        minus_data = minus_table.tolist()
        for row in minus_data:
            row[direction_index] = plus
        minus_data = add_strand_column(minus_data, '-')
        new_data.extend(minus_data)

    return LoadTable(header=table.header[:] + ['strand'], rows=new_data)
