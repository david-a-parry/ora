import logging
import sqlite3
import sys
from collections import namedtuple
from ora import uniprot_lookups
from ora.alignments import align_pos_to_amino_acid, align_range_to_amino_acid
from ora.alignments import get_align_pos
from ora.local import ensg_lookup, get_homologies
from ora.homology_parser import header_fields as result_header_fields
from vase.vcf_reader import VcfReader

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)
ReadBuffer = namedtuple('ReadBuffer',
                        'record genes symbols proteins positions amino_acids')
variant_fields = ['chromosome', 'position', 'id', 'ref', 'alt']
header_fields = (variant_fields + result_header_fields +
                 uniprot_lookups.feat_fields)


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


def process_buffer(record_buffer, gene_orthologies):
    results = []
    for rb in record_buffer:
        for i in range(len(rb.genes)):
            positions = [int(x) for x in rb.positions[i].split('-')]
            if len(positions) > 1:
                start, stop = positions
            else:
                start, stop = positions[0], positions[0]
            ref_aa, var_aa = rb.amino_acids[i].split('/')
            source_features = uniprot_lookups.get_uniprot_features(
                rb.proteins[i], start, stop)
            if source_features is not None:
                for f in source_features:
                    res = dict(chromosome=rb.record.chrom,
                               position=rb.record.pos,
                               id=rb.record.id,
                               ref=rb.record.ref,
                               alt=rb.record.alt,
                               query_gene=rb.genes[i],
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
            if gene_orthologies[rb.genes[i]] is None:
                continue
            for orthology in (x for x in
                              gene_orthologies[rb.genes[i]]['homologies'] if
                              x['source']['protein_id'] == rb.proteins[i]):
                s_id = orthology['source']['id']
                t_id = orthology['target']['id']
                s_protein = orthology['source']['protein_id']
                t_protein = orthology['target']['protein_id']
                s_seq = orthology['source']['align_seq']
                t_seq = orthology['target']['align_seq']
                s_start = get_align_pos(s_seq, start)
                if s_start < 1:  # start is out of range of sequence
                    continue
                if start == stop:
                    s_stop = s_start
                    o_start, aa = align_pos_to_amino_acid(t_seq, s_start)
                    o_stop = o_start
                else:
                    s_stop = get_align_pos(s_seq, stop)
                    if s_stop < 1:  # out of range - treat as SNV
                        s_stop = s_start
                    o_start, o_stop, aa = align_range_to_amino_acid(t_seq,
                                                                    s_start,
                                                                    s_stop)
                if o_start < 1 and o_stop < 1:
                    continue
                elif o_stop < 1:  # only start position found, treat as SNV
                    o_stop = o_start
                elif o_start < 1:  # only stop position found, treat as SNV
                    o_start = o_stop
                ufeats = uniprot_lookups.get_uniprot_features(t_protein,
                                                              o_start,
                                                              o_stop)
                if ufeats is None:
                    continue
                for f in ufeats:
                    res = dict(chromosome=rb.record.chrom,
                               position=rb.record.pos,
                               id=rb.record.id,
                               ref=rb.record.ref,
                               alt=rb.record.alt,
                               query_gene=s_id,
                               query_protein=s_protein,
                               query_symbol=orthology['source'].get('symbol'),
                               query_pos=rb.positions[i],
                               query_aa=s_seq[s_start:s_stop + 1],
                               homolog_gene=t_id,
                               homolog_protein=t_protein,
                               homolog_symbol=orthology['target'].get(
                                   'symbol'),
                               orthology_type=orthology['type'],
                               species=orthology['target']['species'],
                               percent_id=orthology['target']['perc_id'],
                               percent_pos=orthology['target']['perc_pos'],
                               homolog_pos=o_start if o_start == o_stop
                               else "{}-{}".format(o_start, o_stop),
                               query_seq=s_seq,
                               homolog_seq=t_seq,
                               query_species=orthology['source']['species'],
                               homolog_aa=aa,
                               features=f)
                    results.append(res)
    return ["\t".join([str(res[x.lower()]) for x in variant_fields +
                       result_header_fields] +
                      [str(res['features'][x]) for x in
                      uniprot_lookups.feat_fields]) for res in results]


def annotate_variants(args):
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
    record_buffer = []
    current_genes = set()
    gene_orthologies = dict()
    with VcfReader(args.input) as vcf:  # check we can open VCF first
        logger.info("Connecting to local database '{}'".format(args.db))
        conn = sqlite3.connect(args.db)
        curr = conn.cursor()
        logger.info("Reading UniProt feature data")
        uniprot_lookups.initialize()
        out_fh = sys.stdout if args.output is None else open(args.output, 'wt')
        out_fh.write("\t".join(header_fields) + "\n")
        n = 0
        progress_interval = args.progress
        logger.info("Beginning VCF processing")
        for record in vcf:
            if progress_interval and n % progress_interval == 0 and n != 0:
                logger.info("Processed {:,} VCF records. At {}:{}".format(
                    n, record.chrom, record.pos))
            n += 1
            csqs = get_csqs(record)
            if not csqs and not record_buffer:
                continue
            these_genes = [x['Gene'] for x in csqs]
            if current_genes.isdisjoint(set(these_genes)):
                results = process_buffer(record_buffer, gene_orthologies)
                if results:
                    out_fh.write("\n".join(results) + "\n")
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
    results = process_buffer(record_buffer, gene_orthologies)
    if results:
        out_fh.write("\n".join(results) + "\n")
    logger.info("Finished. Processed {:,} VCF records".format(n))
