import sys
import logging
import sqlite3
from collections import defaultdict, namedtuple
from ora.alignments import write_alignments, cigar_to_align_string
from ora.homology_parser import parse_homology_data, check_paralog_lookups
from ora.homology_parser import header_fields as result_header_fields
from ora.id_parser import parse_id
from ora.uniprot_lookups import feat_fields
from ora.uniprot_lookups import logger as unipro_logger

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)

uniprot_lookups = set()
orth_blacklist = ('gene_split', 'alt_allele')
gene_fields = ['gene_member_id', 'stable_id', 'version', 'taxon_id',
               'biotype_group', 'canonical_member_id', 'display_label',
               'taxon_name']
ProteinSeq = namedtuple("ProteinSeq", "protein sequence length")


def get_gene_details(x, curr, species='human'):
    id_info = parse_id(x)
    logger.info("Interpretting input as {}".format(id_info['identifier_type']))
    if id_info['is_ensembl_id']:
        if id_info['is_transcript']:
            raise ValueError("Transcript lookups are not supported for local" +
                             " database lookups.")
        elif id_info['is_protein']:
            return ensp_lookup(x, curr)
        else:
            return ensg_lookup(x, curr)
    else:
        curr.execute('''SELECT taxon_id FROM ncbi_taxa_name WHERE name = ?''',
                     (species,))
        taxon_id = curr.fetchone()
        if not taxon_id:
            raise ValueError('Could not interpret provided species name ' +
                             '"{}"'.format(species))
        return symbol_lookup(x, curr, taxon_id[0])


def ensp_and_seq_from_seq_member(seq_member_id, curr):
    curr.execute('''SELECT stable_id, version, sequence_id from seq_member
                    WHERE seq_member_id = ?''', (str(seq_member_id),))
    result = curr.fetchone()
    curr.execute('''SELECT sequence, length from sequence WHERE
                    sequence_id = ?''', (str(result[2]),))
    seq_res = curr.fetchone()
    return result[0], seq_res[0], seq_res[1]


def parse_gene_details(gene_details, curr):
    d = dict((k, v) for k, v in zip(gene_fields, gene_details))
    protein, seq, lgth = ensp_and_seq_from_seq_member(d['canonical_member_id'],
                                                      curr)
    d['protein'] = ProteinSeq(protein, seq, lgth)
    return d


def ensg_lookup(ensg, curr):
    curr.execute('select * from gene_member WHERE stable_id = ?', (ensg,))
    results = curr.fetchone()
    if not results:
        raise ValueError("No results found for Ensembl gene '{}'".format(ensg))
    return parse_gene_details(results, curr)


def ensp_lookup(ensp, curr):
    curr.execute('SELECT seq_member_id FROM seq_member WHERE stable_id = ?',
                 (ensp,))
    ensgs = curr.fetchone()
    if not ensgs:
        raise ValueError("No results found for Ensembl protein '{}'".format(
            ensp))
    curr.execute('select * from gene_member WHERE canonical_member_id = ?',
                 (ensgs[0],))
    results = curr.fetchone()
    if not results:
        raise ValueError("No results found for Ensembl gene '{}' found".format(
                         ensgs) + " via lookup of Ensembp protein '{}'".format(
                         ensp))
    return parse_gene_details(results, curr)


def symbol_lookup(symbol, curr, taxon_id=9606):
    curr.execute('''select * from gene_member WHERE display_label = ? AND
                    taxon_id = ?''', (symbol, str(taxon_id)))
    results = curr.fetchall()
    if not results:
        raise ValueError("No results found for symbol '{}' and taxon ID {}."
                         .format(symbol, taxon_id))
    if len(results) != 1:
        results = [x for x in results if x[1].startswith("ENS")]  # ignore LRGs
        if len(results) != 1:
            raise ValueError("Multiple results found for symbol " +
                            "'{}' and taxon ID {}. ".format(symbol, taxon_id) +
                            "Try a gene ID instead.")
    return parse_gene_details(results[0], curr)


def combine_query_and_target(query, target, query_homology, target_homology):
    homology = dict(source=dict(), target=dict())
    for k, gene in zip(['source', 'target'], [query, target]):
        homology[k]['species'] = gene['taxon_name']
        homology[k]['id'] = gene['stable_id']
        homology[k]['symbol'] = gene['display_label']
        homology[k]['protein_id'] = gene['protein'].protein
        homology[k]['sequence'] = gene['protein'].sequence
    for k, hom in zip(['source', 'target'], [query_homology, target_homology]):
        homology[k]['perc_id'] = hom['perc_id']
        homology[k]['perc_pos'] = hom['perc_pos']
        homology[k]['cigar_line'] = hom['cigar_line']
        homology[k]['align_seq'] = cigar_to_align_string(
            homology[k]['sequence'],
            hom['cigar_line'])
    homology['type'] = hom['description']
    return homology


def get_homologies(gene_details, curr):
    homolog_fields = ['homology_id', 'gene_member_id', 'seq_member_id',
                      'cigar_line', 'perc_cov', 'perc_id', 'perc_pos',
                      'description', 'is_high_confidence']
    homologies = dict()
    target_homologies = dict()
    target_genes = dict()
    curr.execute('''SELECT * from homology_member WHERE gene_member_id = ?''',
                 (gene_details['gene_member_id'],))
    for row in curr:
        rd = dict((k, v) for k, v in zip(homolog_fields, row))
        if rd['description'] not in orth_blacklist:
            homologies[rd['homology_id']] = rd
    for h_id, gm_id in ((x['homology_id'], x['gene_member_id']) for x in
                        homologies.values()):
        curr.execute('''SELECT * from homology_member WHERE homology_id = ? AND
                     gene_member_id != ?''', (str(h_id), str(gm_id)))
        rd = dict((k, v) for k, v in zip(homolog_fields, curr.fetchone()))
        target_homologies[h_id] = rd
        t_gm_id = rd['gene_member_id']
        curr.execute('SELECT * from gene_member WHERE gene_member_id = ?',
                     (str(t_gm_id),))
        target_genes[t_gm_id] = parse_gene_details(curr.fetchone(), curr)
    hom_data = []
    for k in homologies:
        try:
            target = target_homologies[k]
        except KeyError:
            raise KeyError("Missing target homology for homology_id=" +
                           "{} and gene_member={}. ".format(
                               k, gene_details['gene_member_id']) +
                           "Please check your database is complete.")
        merged = combine_query_and_target(
            query=gene_details,
            target=target_genes[target['gene_member_id']],
            query_homology=homologies[k],
            target_homology=target)
        hom_data.append(merged)
    return hom_data


def local_lookups(gene, pos, db, paralog_lookups=False, line_length=60,
                  all_homologs=False, output_alignments=None, quiet=False,
                  debug=False, silent=False):
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
    for handler in unipro_logger.handlers:
        handler.setLevel(logger.level)
    conn = sqlite3.connect(db)
    curr = conn.cursor()
    lookup_result = get_gene_details(gene, curr)
    hom_data = get_homologies(lookup_result, curr)
    results, paralogs = parse_homology_data(hom_data,
                                            pos,
                                            logger,
                                            output_all_orthologs=all_homologs)
    if paralog_lookups:
        n_paralogs = len(paralogs.keys())
        i = 0
        for gene_id, paralog_list in paralogs.items():
            i += 1
            logger.info("Parsing paralog {:,} of {:,}".format(i, n_paralogs))
            paralog_details = ensg_lookup(gene_id, curr)
            pdata = get_homologies(paralog_details, curr)
            pro2pos = defaultdict(list)
            for para in paralog_list:
                pro2pos[para.protein].append(para.position)
            for protein, positions in pro2pos.items():
                p_results, _ = parse_homology_data(pdata,
                                                   positions,
                                                   logger,
                                                   protein_id=protein,
                                                   skip_paralogs=True)
                results.extend(p_results)
    output_results = [r for r in results if r['should_output']]
    if not output_results:
        logger.info("No results for {} position {}".format(gene, pos))
        sys.exit()
    print("\t".join(result_header_fields + feat_fields))
    if output_alignments:
        alignment_fh = open(output_alignments, 'wt')
        extra_alignments = check_paralog_lookups(lookup_result['stable_id'],
                                                 output_results,
                                                 results)
        write_alignments(output_results + extra_alignments,
                         alignment_fh,
                         linelen=line_length)
        alignment_fh.close()
    for res in output_results:
        line = [str(res[x.lower()]) for x in result_header_fields] + \
               [str(res['features'][x]) for x in feat_fields]
        print("\t".join(line))
