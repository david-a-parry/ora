import logging
import sqlite3
from collections import defaultdict, namedtuple
from ora.alignments import write_alignments, cigar_to_align_string
from ora.alignments import get_align_pos, align_pos_to_amino_acid
from ora.id_parser import parse_id
from ora.uniprot_lookups import get_uniprot_features, feat_fields
from ora.uniprot_lookups import logger as unipro_logger

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)

uniprot_lookups = set()
orth_types = ('ortholog_one2one', 'other_paralog', 'within_species_paralog')
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
        raise ValueError("Multiple results found for symbol " +
                         "'{}' and taxon ID {}. ".format(symbol, taxon_id) +
                         "Try a gene ID instead.")
    return parse_gene_details(results[0], curr)


def merge_and_parse_homologies(query, target, query_align, target_align):
    raise NotImplementedError("local homology parsing is not implmented yet!")


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
    homolog_fields = ['homology_id', 'gene_member_id', 'seq_member_id',
                      'cigar_line', 'perc_cov', 'perc_id', 'perc_pos',
                      'description', 'is_high_confidence']
    homologies = dict()
    target_homologies = dict()
    smember2prot = dict()  # type: dict[int, ProteinSeq]
    curr.execute('''SELECT * from homology_member WHERE gene_member_id = ?''',
                 (lookup_result['gene_member_id'],))
    for row in curr:
        rd = dict((k, v) for k, v in zip(homolog_fields, row))
        if rd['description'] in orth_types:
            homologies[rd['homology_id']] = rd
    for h_id, gm_id in ((x['homology_id'], x['gene_member_id']) for x in
                        homologies.values()):
        curr.execute('''SELECT * from homology_member WHERE homology_id = ? AND
                     gene_member_id != ?''', (str(h_id), str(gm_id)))
        rd = dict((k, v) for k, v in zip(homolog_fields, curr.fetchone()))
        target_homologies[h_id] = rd
        protein, seq, lngth = ensp_and_seq_from_seq_member(rd['seq_member_id'],
                                                           curr)
        smember2prot[rd['seq_member_id']] = ProteinSeq(protein, seq, lngth)
    q_seq = lookup_result['protein'][1]
    for k in homologies:
        try:
            target = target_homologies[k]
        except KeyError:
            raise KeyError("Missing target homology for homology_id=" +
                           "{} and gene_member={}. ".format(
                               k, lookup_result['gene_member_id']) +
                           "Please check your database is complete.")
        t_seq = smember2prot[target['seq_member_id']][1]
        q_align = cigar_to_align_string(q_seq, homologies[k]['cigar_line'])
        t_align = cigar_to_align_string(t_seq, target['cigar_line'])
        merge_and_parse_homologies(query=homologies[k], target=target,
                                   query_align=q_align, target_align=t_align)
#    if output_alignments:
#        alignment_fh = open(output_alignments, 'wt')
