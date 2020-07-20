#!/usr/bin/env python3
import sys
from collections import namedtuple
from ensembl_rest_queries import EnsemblRestQueries

ParalogLookup = namedtuple("ParalogLookup", "gene protein position")


def get_align_pos(seq, p):
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


def main(ensg, pos, paralog_lookups=True):
    try:
        pos = int(pos)
    except ValueError:
        sys.exit("ERROR: Position argument must be an integer\n")
    ens_rest = EnsemblRestQueries()
    data = ens_rest.get_homologies(ensg)
    if data is None:
        sys.exit("ERROR: Could not find ensembl gene '{}'".format(ensg))
    print("\t".join('''Query_Gene Query_Protein Query_Pos Query_Res
                       Homolog_Gene Homolog_Protein Homolog_Type Species %ID
                       Homolog_Pos Homolog_Res'''.split()))
    paralogs = parse_homology_data(data, pos)
    if paralog_lookups:
        for p in paralogs:
            try:
                pdata = ens_rest.get_homologies(p.gene)
            except Exception as e:
                sys.stderr.write("Error looking up paralog {}: {}".format(
                    p.gene,
                    e))
            _ = parse_homology_data(pdata,
                                    p.position,
                                    protein_id=p.protein,
                                    skip_paralogs=True)


def parse_homology_data(data, pos, protein_id=None, skip_paralogs=False):
    paralogs = []
    for i in range(len(data['data'][0]['homologies'])):
        s_id = data['data'][0]['homologies'][i]['source']['id']
        t_id = data['data'][0]['homologies'][i]['target']['id']
        s_protein = data['data'][0]['homologies'][i]['source']['protein_id']
        t_protein = data['data'][0]['homologies'][i]['target']['protein_id']
        if protein_id is not None and s_protein != protein_id:
            sys.stderr.write("Skipping paralog lookup for protein {}\n".format(
                protein_id))
            return paralogs
        hom_type = data['data'][0]['homologies'][i]['type']
        if skip_paralogs and 'paralog' in hom_type:
            continue
        s_seq = data['data'][0]['homologies'][i]['source']['align_seq']
        t_seq = data['data'][0]['homologies'][i]['target']['align_seq']
        p = get_align_pos(s_seq, pos)
        o, aa = align_pos_to_amino_acid(t_seq, p) if p > 0 else (-1, '-')
        if 'paralog' in hom_type and o > 0:
            paralogs.append(ParalogLookup(t_id, t_protein, o))
        print("\t".join(str(x) for x in (
            s_id,
            s_protein,
            pos,
            s_seq[p],
            t_id,
            t_protein,
            hom_type,
            data['data'][0]['homologies'][i]['target']['species'],
            data['data'][0]['homologies'][i]['target']['perc_id'],
            o,
            aa)))
    return paralogs


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: {} ENSG0123456789 100".format(sys.argv[0]))
    main(*sys.argv[1:])
