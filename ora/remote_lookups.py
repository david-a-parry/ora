import sys
import logging
from collections import namedtuple, defaultdict
from ora.ensembl_rest_queries import EnsemblRestQueries
from ora.id_parser import parse_id
from ora.uniprot_lookups import feat_fields
from ora.uniprot_lookups import logger as unipro_logger
from ora.alignments import write_alignments
from ora.homology_parser import parse_homology_data

ensg2symbol = dict()
ens_rest = EnsemblRestQueries()
logger = logging.getLogger("ORA REST")
logger.setLevel(logging.INFO)


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
                logger.debug("Ids are: " + ", ".join(x['id'] for x in data))
                logger.warning("Using ID from first of multiple lookups ({})"
                               .format(data[0]['id']))
            return data[0]['id']
    sys.exit("FAILED: Could not retrieve Ensembl gene ID for input '{}'"
             .format(x))


def check_paralog_lookups(original_ensg, output_results, results):
    '''
        If we have results with a paralog of our original query gene as
        the query gene we should output an alignment of our original
        query gene and the paralog if we aren't already doing so. This
        returns results for our original query gene and these paralogs
        if not already in our output_results.
    '''
    # All pairwise comparisons in our output list
    pair_comps = set((r['query_gene'], r['homolog_gene']) for r in
                     output_results)
    # Pairwise comparisons where the query gene is paralogous to our query
    para_comps = set((r['query_gene'], r['homolog_gene']) for r in
                     output_results if r['query_gene'] != original_ensg)
    extra = list()
    for pair in para_comps:
        if (original_ensg, pair[0]) not in pair_comps:
            extra.extend(r for r in results if r['query_gene'] == original_ensg
                         and r['homolog_gene'] == pair[0])
    return extra


def lookup_symbols(ensgs):
    logger.info("Looking up gene symbols for {:,} IDs".format(len(ensgs)))
    datastring = '{"ids": [' + ", ".join('"' + x + '"' for x in ensgs) + ' ] }'
    data = ens_rest.get_endpoint("/lookup/id", data=datastring)
    if data:
        return dict((k, v.get('display_name', '-')) for k, v in data.items())
    logger.warning("No gene symbols found by POST lookup")
    return dict()


def remote_lookups(gene, pos, paralog_lookups=False, line_length=60,
                   timeout=10.0, max_retries=2, all_homologs=False,
                   output_alignments=None, quiet=False, debug=False,
                   silent=False):
    if silent:
        logger.setLevel(logging.ERROR)
    elif quiet:
        logger.setLevel(logging.WARN)
    elif debug:
        logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logger.level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    unipro_logger.setLevel(logger.level)
    ens_rest.logger.setLevel(logger.level)
    for handler in ens_rest.logger.handlers + unipro_logger.handlers:
        handler.setLevel(logger.level)
    ens_rest.timeout = timeout
    ens_rest.max_retries = max_retries
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
    results, paralogs = parse_homology_data(data['data'][0]['homologies'],
                                            pos,
                                            logger,
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
                p_results, _ = parse_homology_data(
                    pdata['data'][0]['homologies'],
                    positions,
                    protein_id=protein,
                    skip_paralogs=True)
                results.extend(p_results)
    output_results = [r for r in results if r['should_output']]
    if not output_results:
        logger.info("No results for {} position {}".format(gene, pos))
        sys.exit()
    symbol_lookups = set(x['homolog_gene'] for x in output_results)
    symbol_lookups.update(x['query_gene'] for x in output_results)
    ensg2symbol = lookup_symbols(symbol_lookups)
    if output_alignments:
        extra_alignments = check_paralog_lookups(ensg, output_results, results)
        write_alignments(output_results + extra_alignments, alignment_fh,
                         ensg2symbol,
                         linelen=line_length)
        alignment_fh.close()
    for res in output_results:
        res['query_symbol'] = ensg2symbol.get(res['query_gene'], '-')
        res['homolog_symbol'] = ensg2symbol.get(res['homolog_gene'], '-')
        line = [str(res[x.lower()]) for x in header_fields] + \
               [str(res['features'][x]) for x in feat_fields]
        print("\t".join(line))
