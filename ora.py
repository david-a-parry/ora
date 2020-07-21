#!/usr/bin/env python3
import sys
import logging
import argparse
from collections import namedtuple, defaultdict
from Bio.SubsMat.MatrixInfo import pam250
from ora.ensembl_rest_queries import EnsemblRestQueries
from ora.id_parser import parse_id
from ora.uniprot_lookups import get_uniprot_features
from ora.uniprot_lookups import logger as unipro_logger

ParalogLookup = namedtuple("ParalogLookup", "gene protein position")
feat_fields = ['UniprotID', 'Start', 'Stop', 'Feature', 'Description']
ens_rest = EnsemblRestQueries()
ensg2symbol = dict()
uniprot_lookups = set()

logger = logging.getLogger("ORA")
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
                logger.warning("Multiple genes identified for input '{}'"
                            .format(x))
                logger.warning("Using ID from first of multiple lookups ({})"
                            .format(data[0]['id']))
            return data[0]['id']
    sys.exit("FAILED: Could not retrieve Ensembl gene ID for input '{}'"
             .format(x))


def main(gene, pos, paralog_lookups=False, timeout=10.0, max_retries=2,
         all_homologs=False, output_alignments=None, quiet=False, debug=False,
         silent=False):
    if silent:
        logger.setLevel(logging.ERROR)
    elif quiet:
        logger.setLevel(logging.WARN)
    elif debug:
        logger.setLevel(logging.DEBUG)
    ens_rest.logger.setLevel(logger.level)
    unipro_logger.setLevel(logger.level)
    ensg = get_gene_details(gene)
    if ensg != gene:
        logger.info("Got Ensembl Gene ID '{}'".format(ensg))
    data = ens_rest.get_homologies(ensg)
    if data is None:
        sys.exit("ERROR: Could not find ensembl gene '{}'".format(ensg))
    if output_alignments:
        alignment_fh = open(output_alignments, 'wt')
    header_fields = '''Query_Symbol Query_Gene Query_Protein Query_Pos 
                       Query_AA Homolog_Symbol Homolog_Gene Homolog_Protein
                       Orthology_Type Species Percent_ID Percent_Pos
                       Homolog_Pos Homolog_AA'''.split()
    print("\t".join(header_fields + feat_fields))
    results, paralogs = parse_homology_data(data,
                                            pos,
                                            output_all_orthologs=all_homologs)
    if paralog_lookups:
        n_paralogs = len(paralogs.keys())
        i = 0
        for gene_id, paralog_list in paralogs.items():
            i += 1
            logger.info("Parsing paralog {:,} of {:,}".format(i, n_paralogs))
            try:
                pdata = ens_rest.get_homologies(gene_id)
            except Exception as e:
                logger.warning("Error looking up paralog {}: {}".format(
                    gene_id, e))
                continue
            pro2pos = defaultdict(list)
            for para in paralog_list:
                pro2pos[para.protein].append(para.position)
            for protein, positions in pro2pos.items():
                p_results, _ = parse_homology_data(pdata,
                                                   positions,
                                                   protein_id=protein,
                                                   skip_paralogs=True)
                results.extend(p_results)
    if not results:
        logger.info("No results for {} position {}".format(gene, pos))
        sys.exit()
    symbol_lookups = set(x['homolog_gene'] for x in results)
    symbol_lookups.update(x['query_gene'] for x in results)
    ensg2symbol = lookup_symbols(symbol_lookups)
    if output_alignments:
        write_alignments(results, alignment_fh, ensg2symbol)
        alignment_fh.close()
    for res in results:
        res['query_symbol'] = ensg2symbol.get(res['query_gene'], '-')
        res['homolog_symbol'] = ensg2symbol.get(res['homolog_gene'], '-')
        line = [str(res[x.lower()]) for x in header_fields] + \
               [str(res['features'][x]) for x in feat_fields]
        print("\t".join(line))


def seq_and_pos_from_results(results):
    pairs = dict()
    for res in (x for x in results if 'query_seq' in x):
        align_pos = get_align_pos(res['query_seq'], res['query_pos'])
        k = (res['query_gene'], res['homolog_gene'])
        if k not in pairs:
            pairs[k] = dict(align_positions=set(), query_positions=set())
            for x in ('query_species', 'species', 'query_protein', 'query_seq',
                      'homolog_protein', 'homolog_seq'):
                pairs[k][x] = res[x]
        pairs[k]['align_positions'].add(align_pos)
        pairs[k]['query_positions'].add(res['query_pos'])
    for k in pairs:
        for p in ('align_positions', 'query_positions'):
            pairs[k][p] = sorted(pairs[k][p])
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


def write_alignments(results, fh, gene2symbol, linelen=60):
    for query_hom, res in seq_and_pos_from_results(results).items():
        query_gene, hom_gene = query_hom
        qsymbol = gene2symbol.get(query_gene, '?')
        hsymbol = gene2symbol.get(hom_gene, '?')
        header = "|{} {} vs {} {} position ".format(
            qsymbol, res['query_species'], hsymbol, res['species']) + \
            ",".join(str(x) for x in res['query_positions']) + "|"
        lmargin = max(len(res['query_protein']), len(res['homolog_protein']))
        fh.write('-' * len(header) + "\n")
        fh.write(header + "\n")
        fh.write('-' * len(header) + "\n")
        qlen = len(res['query_seq'])
        for i in range(0, qlen, linelen):
            # query sequence
            fh.write("{:>{fill}}: {}\n".format(res['query_protein'],
                                               res['query_seq'][i:i+linelen],
                                               fill=lmargin))
            # homolog sequence
            fh.write("{:>{fill}}: {}\n".format(res['homolog_protein'],
                                               res['homolog_seq'][i:i+linelen],
                                               fill=lmargin))
            # conservation status
            fh.write("{:>{fill}}  ".format(' ', fill=lmargin))
            fh.write(conservation_status(res['query_seq'][i:i+linelen],
                                         res['homolog_seq'][i:i+linelen])
                     + "\n")
            for align_pos in res['align_positions']:
                if i <= align_pos < i + linelen:  # pos on current line
                    l_pos = align_pos - i
                    l_pad = lmargin + 2
                    r_pad = min(linelen, qlen - i) - l_pos - 1
                    fh.write(' ' * l_pad + '_' * l_pos + '^' + '_' * r_pad +
                             '\n')
            fh.write("\n")
        fh.write("\n")


def lookup_symbols(ensgs):
    logger.info("Looking up gene symbols for {:,} IDs".format(len(ensgs)))
    datastring = '{"ids": [' + ", ".join('"' + x + '"' for x in ensgs) + ' ] }'
    data = ens_rest.get_endpoint("/lookup/id", data=datastring)
    if data:
        return dict((k, v.get('display_name', '-')) for k, v in data.items())
    logger.warning("No gene symbols found by POST lookup")
    return dict()


def parse_homology_data(data, positions, protein_id=None, skip_paralogs=False,
                        output_all_orthologs=False):
    paralogs = defaultdict(list)
    results = []
    homs = data['data'][0]['homologies']
    logger.info("Got {:,} homologs".format(len(homs)))
    for pos in positions:
        for i in range(len(homs)):
            s_id = homs[i]['source']['id']
            t_id = homs[i]['target']['id']
            s_protein = homs[i]['source']['protein_id']
            t_protein = homs[i]['target']['protein_id']
            s_seq = homs[i]['source']['align_seq']
            t_seq = homs[i]['target']['align_seq']
            p = get_align_pos(s_seq, pos)
            if (s_protein, pos) not in uniprot_lookups:
                ufeats = get_uniprot_features(s_protein, pos, pos)
                if ufeats:
                    for f in ufeats:
                        result = dict(query_gene=s_id,
                                      query_protein=s_protein,
                                      query_pos=pos,
                                      query_aa=s_seq[p],
                                      homolog_gene=s_id,
                                      homolog_protein=s_protein,
                                      orthology_type="self",
                                      species=homs[i]['source']['species'],
                                      percent_id=100,
                                      percent_pos=100,
                                      homolog_pos=pos,
                                      homolog_aa=s_seq[p],
                                      features=f)
                        results.append(result)
                uniprot_lookups.add((s_protein, pos))
            if protein_id is not None and s_protein != protein_id:
                logger.info("Skipping paralog lookup for protein {}\n".format(
                    protein_id))
                return results, paralogs
            hom_type = homs[i]['type']
            if skip_paralogs and 'paralog' in hom_type:
                continue
            o, aa = align_pos_to_amino_acid(t_seq, p) if p > 0 else (-1, '-')
            if 'paralog' in hom_type and o > 0:
                paralogs[t_id].append(ParalogLookup(t_id, t_protein, o))
            if o > 0 and (t_protein, o) not in uniprot_lookups:
                result_template = dict(
                                    query_gene=s_id,
                                    query_protein=s_protein,
                                    query_pos=pos,
                                    query_aa=s_seq[p],
                                    homolog_gene=t_id,
                                    homolog_protein=t_protein,
                                    orthology_type="self",
                                    species=homs[i]['target']['species'],
                                    percent_id=homs[i]['target']['perc_id'],
                                    percent_pos=homs[i]['target']['perc_pos'],
                                    homolog_pos=o,
                                    query_seq=s_seq,
                                    homolog_seq=t_seq,
                                    query_species=homs[i]['source']['species'],
                                    homolog_aa=aa)
                ufeats = get_uniprot_features(t_protein, o, o)
                if ufeats:
                    for f in ufeats:
                        result = result_template.copy()
                        result['features'] = f
                        results.append(result)
                elif output_all_orthologs:
                    result = result_template.copy()
                    result['features'] = dict((x, '-') for x in feat_fields)
                    results.append(result)
                uniprot_lookups.add((t_protein, o))
    if not skip_paralogs:
        logger.info("Got {} paralog sequences".format(len(paralogs)))
    return results, paralogs


def get_options():
    parser = argparse.ArgumentParser(
        description='''Search for protein features in orthologs''')
    parser.add_argument("-g", "--gene", required=True,
                        help='''Gene/protein ID to search with. Accession types
                        will be inferred from input and the corresponding
                        Ensembl gene ID will be identified.''')
    parser.add_argument("-p", "--pos", required=True, type=int, nargs='+',
                        help='''Position(s) in protein to search. This must be
                        relative to the canonical Ensembl transcript as used by
                        Ensembl's Compara database.''')
    parser.add_argument("--paralog_lookups", action='store_true',
                        help='''Also perform homology lookups on paralogs
                        identified.''')
    parser.add_argument("--all_homologs", action='store_true',
                        help='''Output alignment information for all homologs
                        even if no features are found.''')
    parser.add_argument("--output_alignments", help='''Output alignments of
                        with hits to this file.''')
    parser.add_argument("--timeout", type=int, default=10,
                        help='''Lookup timeout (in seconds) for Ensembl REST
                        queries. Default=10''')
    parser.add_argument("--max_retries", type=int, default=2,
                        help='''Maximum retry attempts for Ensembl REST
                        queries. Default=2''')
    parser.add_argument("--debug", action='store_true',
                        help='Add debugging messages to output.')
    parser.add_argument("--quiet", action='store_true',
                        help='Only output warning logger messages.')
    parser.add_argument("--silent", action='store_true',
                        help='Only output error logger messages.')
    return parser


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
