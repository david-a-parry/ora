import logging
import sqlite3
import pysam
from collections import namedtuple
from ora import uniprot_lookups
from ora.alignments import align_pos_to_amino_acid, get_align_pos
from ora.local import ensg_lookup, get_homologies
from vase.vcf_reader import VcfReader

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)
ReadBuffer = namedtuple('ReadBuffer',
                        'record genes symbols proteins positions amino_acids')


def get_orthologies(gene, cursor, paralog_lookups=False):
    gene_lookup = ensg_lookup(gene, cursor)
    homology_data = get_homologies(gene_lookup, cursor)
    paralogs = dict()
    if paralog_lookups:
        for paralogy in (x for x in homology_data if 'paralog' in x['type']):
            p_lookup = ensg_lookup(paralogy['id'], cursor)
            paralogs[paralogy['id']] = get_homologies(p_lookup, cursor)
    return dict(homologies=homology_data, paralogs=paralogs)


def get_csqs(record, csq_types=['missense_variant', 'inframe_deletion',
                                'inframe_insertion']):
    return [x for x in record.CSQ for y in x['Consequence'].split('&')
            if y in csq_types]


def annotate_record(record, results):
    pass  # TODO!


def process_buffer(record_buffer, gene_orthologies):
    results = []
    for rb in record_buffer:
        for i in range(len(rb.genes)):
            source_features = uniprot_lookups.get_uniprot_features(
                rb.proteins[i])
            ref_aa, var_aa = rb.amino_acids[i].split('/')
            for f in source_features:
                res = dict(query_gene=rb.genes[i],
                           query_symbol=rb.symbols[i],
                           query_protein=rb.proteins[i],
                           query_pos=rb.positions[i],
                           query_aa=ref_aa,
                           homolog_gene=rb.genes[i],
                           homolog_symbol=rb.symbols[i],
                           homolog_protein=rb.proteins[i],
                           orthology_type="self",
                           species='human',
                           percent_id=100,
                           percent_pos=100,
                           homolog_pos=rb.positions[i],
                           homolog_aa=ref_aa,
                           query_species='human',
                           features=f)
                results.append(res)
            for orthology in (x for x in
                              gene_orthologies[rb.genes[i]]['homologies'] if
                              x['source']['protein_id'] == rb.proteins[i]):
                s_id = orthology['source']['id']
                t_id = orthology['target']['id']
                s_protein = orthology['source']['protein_id']
                t_protein = orthology['target']['protein_id']
                s_seq = orthology['source']['align_seq']
                t_seq = orthology['target']['align_seq']
                p = get_align_pos(s_seq, rb.positions[i])
                o, aa = align_pos_to_amino_acid(t_seq, p) if p > 0 else (-1,
                                                                         '-')
                if o < 1:
                    continue
                ufeats = uniprot_lookups.get_uniprot_features(t_protein, o, o)
                for f in ufeats:
                    res = dict(query_gene=s_id,
                               query_protein=s_protein,
                               query_symbol=orthology['source'].get('symbol'),
                               query_pos=rb.positions[i],
                               query_aa=s_seq[p],
                               homolog_gene=t_id,
                               homolog_protein=t_protein,
                               homolog_symbol=orthology['target'].get(
                                   'symbol'),
                               orthology_type=orthology['type'],
                               species=orthology['target']['species'],
                               percent_id=orthology['target']['perc_id'],
                               percent_pos=orthology['target']['perc_pos'],
                               homolog_pos=o,
                               query_seq=s_seq,
                               homolog_seq=t_seq,
                               query_species=orthology['source']['species'],
                               homolog_aa=aa,
                               features=f)
                    results.append(res)
        annotate_record(rb.record, results)
    return results


def vcf_annotator(args):
    if args.silent:
        logger.setLevel(logging.ERROR)
    elif args.quiet:
        logger.setLevel(logging.WARN)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
    ch = logging.StreamHandler()
    ch.setLevel(logger.level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    uniprot_lookups.logger.setLevel(logger.level)
    for handler in uniprot_lookups.logger.handlers:
        handler.setLevel(logger.level)
    logger.info("Connecting to local database '{}'".format(args.db))
    conn = sqlite3.connect(args.db)
    curr = conn.cursor()
    logger.info("Reading UniProt feature data")
    uniprot_lookups.initialize()
    record_buffer = []
    current_genes = set()
    gene_orthologies = dict()
    out = '-' if args.output is None else args.output
    logger.info("Beginning VCF processing")
    with VcfReader(args.input) as vcf:
        # TODO - add header information
        out_vcf = pysam.VariantFile(out, 'w', header=vcf.header.header)
        n = 0
        progress_interval = 10_000
        for record in vcf:
            # if does not match our consequence classes and record_buffer empty
            # just write variant and move on
            csqs = get_csqs(record)
            n += 1
            if n % progress_interval == 0:
                logger.info("Read {:,} VCF records".format(n))
            if not csqs and not record_buffer:
                out_vcf.write(record.record)
                continue
            these_genes = [x['Gene'] for x in csqs]
            if current_genes.isdisjoint(set(these_genes)):
                process_buffer(record_buffer, gene_orthologies)
                record_buffer = []
                current_genes.clear()
                gene_orthologies.clear()
            record_buffer.append(ReadBuffer(record=record,
                                            genes=these_genes,
                                            symbols=[x['SYMBOL'] for x in
                                                     csqs],
                                            proteins=[x['ENSP'] for x in csqs],
                                            positions=[x['Protein_position']
                                                       for x in csqs],
                                            amino_acids=[x['Amino_acids'] for x
                                                         in csqs]))
            current_genes.update(these_genes)
            for gene in (x for x in these_genes if x not in gene_orthologies):
                try:
                    gene_orthologies[gene] = get_orthologies(
                        gene, curr, args.paralog_lookups)
                except LookupError:
                    logger.warn("Gene ID '{}' not present in database".format(
                        gene))
                    gene_orthologies[gene] = None
    process_buffer(record_buffer, gene_orthologies)


#    default_args = [input, db, paralog_lookups=False, line_length=60,
#                    all_homologs=False, output_alignments=None, quiet=False,
#                    debug=False, silent=False]
