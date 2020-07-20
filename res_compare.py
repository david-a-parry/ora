#!/usr/bin/env python3
import sys
import logging
import argparse
from collections import namedtuple
from ensembl_rest_queries import EnsemblRestQueries
from id_parser import parse_id
from uniprot_lookups import get_uniprot_features

ParalogLookup = namedtuple("ParalogLookup", "gene protein position")
feat_fields = ['UniprotID', 'Start', 'Stop', 'Feature', 'Description']
ens_rest = EnsemblRestQueries()
ensg2symbol = dict()

logger = logging.getLogger("ORA")
#logger.setLevel(logging.INFO)
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)


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


def get_gene_details(x, species='human'):
    id_info = parse_id(x)
    logger.info("Interpretting input as {}".format(id_info['identifier_type']))
    if id_info['is_ensembl_id']:
        if id_info['is_transcript']:
            data = ens_rest.gene_from_enst(x)
        elif id_info['is_protein']:
            data = ens_rest.gene_from_ensp(x)
        else:
            return x  # is ensembl gene id
        if data:
            return data['id']
    else:
        data = ens_rest.get_via_xref(x, species, 'gene')
        if data:
            if len(data) > 1:
                logger.warn("Multiple genes identified for input '{}'"
                            .format(x))
                logger.warn("Using ID from first of multiple lookups ({})"
                            .format(data[0]['id']))
            return data[0]['id']
    sys.exit("FAILED: Could not retrieve Ensembl gene ID for input '{}'"
             .format(x))


def main(gene, pos, paralog_lookups=False, timeout=10.0, max_retries=2,
         all_homologs=False):
    try:
        pos = int(pos)
    except ValueError:
        sys.exit("ERROR: Position argument must be an integer\n")
    ensg = get_gene_details(gene)
    if ensg != gene:
        logger.info("Got Ensembl Gene ID '{}'".format(ensg))
    data = ens_rest.get_homologies(ensg)
    if data is None:
        sys.exit("ERROR: Could not find ensembl gene '{}'".format(ensg))
    print("\t".join('''Query_Symbol Query_Gene Query_Protein Query_Pos 
                       Query_Res Homolog_Symbol Homolog_Gene Homolog_Protein
                       Homolog_Type Species %ID Homolog_Pos Homolog_Res'''
                    .split() + feat_fields))
    paralogs = parse_homology_data(data, pos,
                                   output_all_orthologs=all_homologs)
    if paralog_lookups:
        for i in range(len(paralogs)):
            logger.info("Parsing paralog {:,} of {:,}".format(i + 1,
                                                              len(paralogs)))
            try:
                pdata = ens_rest.get_homologies(paralogs[i].gene)
            except Exception as e:
                sys.stderr.write("Error looking up paralog {}: {}".format(
                    paralogs[i].gene,
                    e))
            _ = parse_homology_data(pdata,
                                    paralogs[i].position,
                                    protein_id=paralogs[i].protein,
                                    skip_paralogs=True)


def lookup_symbol(ensg):
    data = ens_rest.lookup_id(ensg)
    ensg2symbol[ensg] = data.get('display_name', 'N/A') if data else None
    return ensg2symbol[ensg]


def parse_homology_data(data, pos, protein_id=None, skip_paralogs=False,
                        output_all_orthologs=False):
    paralogs = []
    uniprot_lookups = set()
    logger.info("Got {:,} homologs".format(len(data['data'][0]['homologies'])))
    for i in range(len(data['data'][0]['homologies'])):
        s_id = data['data'][0]['homologies'][i]['source']['id']
        t_id = data['data'][0]['homologies'][i]['target']['id']
        s_protein = data['data'][0]['homologies'][i]['source']['protein_id']
        t_protein = data['data'][0]['homologies'][i]['target']['protein_id']
        s_seq = data['data'][0]['homologies'][i]['source']['align_seq']
        t_seq = data['data'][0]['homologies'][i]['target']['align_seq']
        p = get_align_pos(s_seq, pos)
        s_symbol = ensg2symbol.get(s_id, lookup_symbol(s_id))
        if (s_protein, pos) not in uniprot_lookups:
            ufeats = get_uniprot_features(s_protein, pos, pos)
            if ufeats:
                for f in ufeats:
                    print("\t".join(str(x) for x in [
                            s_symbol,
                            s_id,
                            s_protein,
                            pos,
                            s_seq[p],
                            s_symbol,
                            s_id,
                            s_protein,
                            "self",
                            data['data'][0]['homologies'][i]['source']['species'],
                            100,
                            pos,
                            s_seq[p],
                        ] + [f[x] for x in feat_fields]))
            uniprot_lookups.add((s_protein, pos))
        if protein_id is not None and s_protein != protein_id:
            logger.info("Skipping paralog lookup for protein {}\n".format(
                protein_id))
            return paralogs
        hom_type = data['data'][0]['homologies'][i]['type']
        if skip_paralogs and 'paralog' in hom_type:
            continue
        o, aa = align_pos_to_amino_acid(t_seq, p) if p > 0 else (-1, '-')
        if 'paralog' in hom_type and o > 0:
            paralogs.append(ParalogLookup(t_id, t_protein, o))
        if o > 0 and (t_protein, o) not in uniprot_lookups:
            ufeats = get_uniprot_features(t_protein, o, o)
            if ufeats:
                t_symbol = ensg2symbol.get(t_id, lookup_symbol(t_id))
                for f in ufeats:
                    print("\t".join(str(x) for x in [
                        s_symbol,
                        s_id,
                        s_protein,
                        pos,
                        s_seq[p],
                        t_symbol,
                        t_id,
                        t_protein,
                        hom_type,
                        data['data'][0]['homologies'][i]['target']['species'],
                        data['data'][0]['homologies'][i]['target']['perc_id'],
                        o,
                        aa,
                    ] + [f[x] for x in feat_fields]))
            elif output_all_orthologs:
                t_symbol = ensg2symbol.get(t_id, lookup_symbol(t_id))
                print("\t".join(str(x) for x in [
                    s_symbol,
                    s_id,
                    s_protein,
                    pos,
                    s_seq[p],
                    t_symbol,
                    t_id,
                    t_protein,
                    hom_type,
                    data['data'][0]['homologies'][i]['target']['species'],
                    data['data'][0]['homologies'][i]['target']['perc_id'],
                    o,
                    aa,
                ] + ['-' for x in feat_fields]))
            uniprot_lookups.add((t_protein, o))
    if not skip_paralogs:
        logger.info("Got {} paralog sequences".format(len(paralogs)))
    return paralogs


def get_options():
    parser = argparse.ArgumentParser(
        description='''Search for protein features in orthologs''')
    parser.add_argument("-g", "--gene", required=True,
                        help='''Gene/protein ID to search with. Accession types
                        will be inferred from input and the corresponding
                        Ensembl gene ID will be identified.''')
    parser.add_argument("-p", "--pos", required=True, type=int,
                        help='''Position in protein to search. This must be
                        relative to the canonical Ensembl transcript as used by
                        Ensembl's Compara database.''')
    parser.add_argument("--paralog_lookups", action='store_true',
                        help='''Also perform homology lookups on paralogs
                        identified.''')
    parser.add_argument("--all_homologs", action='store_true',
                        help='''Output alignment information for all homologs
                        even if no features are found.''')
    parser.add_argument("--timeout", type=int, default=10,
                        help='''Lookup timeout (in seconds) for Ensembl REST
                        queries. Default=10''')
    parser.add_argument("--max_retries", type=int, default=2,
                        help='''Maximum retry attempts for Ensembl REST 
                        queries. Default=2''')
    return parser


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
