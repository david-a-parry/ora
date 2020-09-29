import re
from Bio.SubsMat.MatrixInfo import pam250

cigar_re = re.compile('(\d+)?([\D])')


def cigar_to_align_string(seq, cigar):
    ctups = ((1, x[1]) if x[0] == '' else (int(x[0]), x[1]) for x in
             cigar_re.findall(cigar))
    i = 0
    align = ''
    for ct in ctups:
        if ct[1] == 'D':
            align += '-' * ct[0]
        elif ct[1] == 'M':
            align += seq[i:i+ct[0]]
            i += ct[0]
        else:
            raise ValueError("Unexpected cigar operator '{}'".format(ct[1]))
    return align


def get_align_pos(seq, p):
    ''' Return the position within alignment of given residue number'''
    x = 0
    for i in range(len(seq)):
        if seq[i] == '-':
            continue
        x += 1
        if x == p:
            return i
    return -1


def align_pos_to_amino_acid(seq, i):
    ''' Returns position within protein and the amino acid residue'''
    if seq[i] == '-':  # no residue at this position in ortholog
        return -1, '-'
    p = seq[:i + 1].replace('-', '')
    return len(p), p[-1]


def seq_and_pos_from_results(results):
    pairs = dict()
    for res in (x for x in results if 'query_seq' in x):
        align_pos = get_align_pos(res['query_seq'], res['query_pos'])
        k = (res['query_gene'], res['homolog_gene'])
        if k not in pairs:
            pairs[k] = dict(align_positions=list(), query_positions=list(),
                            query_aas=list())
            for x in ('query_species', 'species', 'query_protein', 'query_seq',
                      'homolog_protein', 'homolog_seq'):
                pairs[k][x] = res[x]
        pairs[k]['align_positions'].append(align_pos)
        pairs[k]['query_positions'].append(res['query_pos'])
        pairs[k]['query_aas'].append(res['query_aa'])
    return pairs


def get_conservation_symbol(aa1, aa2):
    if aa1 == aa2:
        return '*' if aa1 != '-' else ' '
    s = pam250.get((aa1, aa2), pam250.get((aa2, aa1), -99))
    if s > 0.5:
        return ':'
    if s >= 0.0:
        return '.'
    return ' '


def conservation_status(seq1, seq2):
    return ''.join(get_conservation_symbol(x, y) for x, y in zip(seq1, seq2))


def arrange_labels(label_tups, line_length, l_margin):
    '''labels is a list of tuples of labels and postitions on the line'''
    label_tups = sorted(label_tups, key=lambda x: x[1])
    labels, positions = zip(*label_tups)
    written = set()
    label_lines = []
    tag_pos = []  # position of | chars tagging label to seq
    for i in range(len(labels)):
        if labels[i] in written:
            continue
        new_tag_pos = []
        s = ' ' * (line_length + l_margin)
        span_i = positions[i] + len(labels[i]) - 1
        s = s[:positions[i]] + labels[i] + s[span_i:]
        new_tag_pos.append(positions[i])
        written.add(labels[i])
        for j in range(i + 1, len(labels)):
            if labels[j] in written:
                continue
            if positions[i] <= positions[j] <= span_i:
                continue
            span_j = positions[j] + len(labels[j]) - 1
            s = s[:positions[j]] + labels[j] + s[span_j:]
            written.add(labels[j])
            new_tag_pos.append(positions[j])
        for t in tag_pos:
            # we add | tags AFTER adding labels because each added label would
            # shift the position of the | tag if we added tags first
            s = s[:t] + '|' + s[t+1:]
        tag_pos.extend(new_tag_pos)
        label_lines.append(s.rstrip())
    s = ' ' * l_margin + '_' * line_length
    for t in tag_pos:
        s = s[:t] + '|' + s[t+1:]
    label_lines.append(s)
    return label_lines


def write_alignments(results, fh, gene2symbol, linelen=60):
    for query_hom, res in seq_and_pos_from_results(results).items():
        query_gene, hom_gene = query_hom
        qsymbol = gene2symbol.get(query_gene, '?')
        hsymbol = gene2symbol.get(hom_gene, '?')
        header = "|{} {} vs {} {} position ".format(
            qsymbol, res['query_species'], hsymbol, res['species']) + \
            ",".join(str(x) for x in sorted(set(res['query_positions']))) + "|"
        lmargin = max(len(res['query_protein']), len(res['homolog_protein']))
        fh.write('-' * len(header) + "\n")
        fh.write(header + "\n")
        fh.write('-' * len(header) + "\n")
        qlen = len(res['query_seq'])
        for i in range(0, qlen, linelen):
            # conservation status
            pos_overlaps = [(x, y, z) for x, y, z in
                            zip(res['align_positions'],
                                res['query_positions'],
                                res['query_aas']) if i <= x < i + linelen]
            if pos_overlaps:
                l_pad = lmargin + 2
                labels = set()
                for align_pos, qpos, qaa in pos_overlaps:
                    j = align_pos - i + l_pad
                    labels.add(("{}{}".format(qaa, qpos), j))
                label_text = arrange_labels(labels, linelen, l_pad)
                fh.write('\n'.join(label_text) + '\n')
            else:
                fh.write("\n")
            # query sequence
            fh.write("{:>{fill}}: {}\n".format(res['query_protein'],
                                               res['query_seq'][i:i+linelen],
                                               fill=lmargin))
            # homolog sequence
            fh.write("{:>{fill}}: {}\n".format(res['homolog_protein'],
                                               res['homolog_seq'][i:i+linelen],
                                               fill=lmargin))
            fh.write("{:>{fill}}  ".format(' ', fill=lmargin))
            fh.write(conservation_status(res['query_seq'][i:i+linelen],
                                         res['homolog_seq'][i:i+linelen])
                     + "\n")
        fh.write("\n")
