import csv
import gzip
import logging
import os
import re
import sqlite3
import sys
from Bio.Data.IUPACData import protein_letters_3to1
from collections import namedtuple
from ora import uniprot_lookups
from ora.alignments import align_pos_to_amino_acid, align_range_to_amino_acid
from ora.alignments import get_align_pos, write_alignments, score_alignment
from ora.local import ensg_lookup, get_homologies
from ora.homology_parser import header_fields as result_header_fields
from vase.vcf_reader import VcfReader

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)
ReadBuffer = namedtuple('ReadBuffer', 'record genes csqs')
variant_fields = ['Chromosome', 'Position', 'ID', 'Ref', 'Alt', 'Ref_AA',
                  'Alt_AA']
result_header_fields = [x for x in result_header_fields if x != 'Query_AA']
result_header_fields.insert(0, 'Variant_Description')
result_header_fields.insert(-1, 'Alt_Score')
header_fields = (variant_fields + result_header_fields +
                 uniprot_lookups.feat_fields)

ens_version_re = re.compile(r"\.\d+$")
snv_re = re.compile(r'''p.([A-z]{3})  # first amino acid
                        (\d+)         # position
                        ([A-z]{3})$   # AA if SNV, "del" or "dup" if 1 AA indel
                    ''', re.X)
indel_re = re.compile(r'''p.([A-z]{3})      # first amino acid
                            (\d+)           # start position
                            _([A-z]{3})     # end amino acid
                            (\d+)           # end position
                            (ins|del|dup)   # indel type
                            (\w+)?$         # optional amino acids or 'delins'
                      ''', re.X)
delins_re = re.compile(r'''p.([A-z]{3})     # first amino acid
                            (\d+)           # start position
                            del(ins\w+)$    # inserted amino acids
                       ''', re.X)
protein_letters_3to1['Ter'] = '*'
csq_types = ['missense_variant', 'inframe_deletion', 'inframe_insertion',
             'conservative_inframe_deletion', 'conservative_inframe_insertion',
             'disruptive_inframe_deletion', 'disruptive_inframe_insertion']
enst2ensp = None


class HgvspError(ValueError):
    pass


def get_orthologies(gene, cursor, paralog_lookups=False):
    gene_lookup = ensg_lookup(gene, cursor)
    homology_data = get_homologies(gene_lookup, cursor)
    paralogs = dict()
    if paralog_lookups:
        for paralogy in (x for x in homology_data if 'paralog' in x['type']):
            p_lookup = ensg_lookup(paralogy['target']['id'], cursor)
            paralogs[paralogy['target']['id']] = get_homologies(p_lookup,
                                                                cursor)
    return dict(homologies=homology_data, paralogs=paralogs)


def read_enst2ensp_data():
    ens_file = os.path.join(os.path.dirname(__file__),
                            "data",
                            "enst2ensp.txt.gz")
    t2p = dict()
    with gzip.open(ens_file, 'rt') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for f in ('Transcript stable ID', 'Protein stable ID'):
            if f not in reader.fieldnames:
                raise ValueError("Invalid header for ensembl transcript to " +
                                 "protein data file '{}'".format(ens_file))
            for row in reader:
                t2p[row['Transcript stable ID']] = row['Protein stable ID']
    return t2p


def snpeff2vep(ann, curr):
    csq = {'SYMBOL': ann['Gene_Name'], 'Gene': ann['Gene_ID']}
    enst = ann['Feature_ID'].split('.')[0]
    ensp = enst2ensp.get(enst, enst)
    csq['ENSP'] = ensp
    csq['HGVSp'] = "{}:{}".format(ensp, ann['HGVS.p'])
    indel_match = indel_re.match(ann['HGVS.p'])
    delins_match = delins_re.match(ann['HGVS.p'])
    refaa = '?'
    varaa = '?'
    if indel_match or delins_match:
        if indel_match:
            ref, start, ref_end, stop, indel, insertion = indel_match.groups()
        else:
            ref, start, insertion = delins_match.groups()
            stop = start
            indel = 'del'
        pos = "{}-{}".format(start, stop)
        if stop == start:
            refaa = protein_letters_3to1[ref]
        elif int(stop) == int(start) + 1:
            refaa = protein_letters_3to1[ref] + protein_letters_3to1[ref_end]
        else:
            curr.execute('''SELECT sequence from sequence JOIN seq_member ON
                            seq_member.sequence_id == sequence.sequence_id
                            WHERE seq_member.stable_id == ?''', (ensp,))
            seqs = curr.fetchone()
            if seqs:
                refaa = seqs[0][int(start)-1:int(stop)]
        if indel == 'dup':
            varaa = refaa + refaa
        elif indel == 'del':
            if insertion and insertion.startswith('ins'):
                varaa = "".join(protein_letters_3to1[insertion[i:i+3]] for i in
                                range(3, len(insertion), 3))
            else:
                varaa = '-'
        elif indel == 'ins':
            varaa = refaa + "".join(protein_letters_3to1[insertion[i:i+3]] for
                                    i in range(0, len(insertion), 3))
        else:
            raise ValueError("Could not parse indel type '{}' in HGVSp '{}'"
                             .format(indel, ann['HGVS.p']))
    else:
        snv_match = snv_re.match(ann['HGVS.p'])
        if snv_match:
            ref, pos, var = snv_match.groups()
            refaa = protein_letters_3to1[ref]
            if var == 'del':
                varaa = '-'
            elif var == 'dup':
                var = refaa + refaa
            else:
                varaa = protein_letters_3to1[var]
        else:
            raise HgvspError("Could not parse HGVSp annotation '{}'".format(
                ann['HGVS.p']))
    csq['Amino_acids'] = '{}/{}'.format(refaa, varaa)
    csq['Protein_position'] = pos
    return csq


def get_csqs_snpeff(record, curr):
    csqs = []
    for ann in (x for x in record.ANN for y in x['Annotation'].split('&') if y
                in csq_types):
        try:
            csqs.append(snpeff2vep(ann, curr))
        except HgvspError as e:
            logger.error('{} at {}:{}-{}/{}'.format(e,
                                                    record.chrom,
                                                    record.pos,
                                                    record.ref,
                                                    record.alt))
    return csqs


def get_csqs_vep(record, curr):
    return [x for x in record.CSQ for y in x['Consequence'].split('&') if y in
            csq_types]


def trim_ref_alt(ref, alt):
    if len(ref) >= len(alt):
        adjust = 0
        while len(alt) > 0:
            if ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                adjust += 1
            else:
                break
        if alt == '':
            alt = '-'
    elif len(ref) < len(alt):
        adjust = -1
        while len(ref) > 0:
            if ref[0] == alt[0]:
                ref = ref[1:]
                alt = alt[1:]
                adjust += 1
            else:
                break
        if ref == '':
            ref = '-'
    return ref, alt, adjust


def parse_amino_acids(csq):
    '''
        Return reference amino acid(s), variant amino acid(s), protein start
        position, protein end position and variant description from standard
        VEP fields. Used as a basic fallback when HGVSp information is
        unavailable.

        SNV - unchanged
        'K', 'E', 100, 100 => 'K', 'E', 100, 100, 'K100E'

        DELETION
        'KD', 'K', 100, 101 => 'D', '-', 101, 101, '101delD'

        DELETION
        'KDL', 'K', 100, 102 => 'DL', '-', 101, 102, '101_102_delDL'

        INSERTION
        'K', 'KQ', 100, 100 => 'K', 'Q', 100, 101, '100_101_insQ'

        INSERTION
        'K', 'KQR', 100, 100 => '-', 'QR', 100, 101, '100_101insQR'

        INSERTION
        'KQQ', 'KQQQR', 100, 100 => '-', 'QR', 102, 103, '102_103_insQR'

    '''
    ref, var = csq['Amino_acids'].split('/')
    positions = csq['Protein_position'].split('-')
    start = int(positions[0])
    if len(positions) == 1:
        stop = start
    else:
        stop = int(positions[1])
    if len(ref) != len(var):
        original_ref = ref
        ref, var, start_adjust = trim_ref_alt(ref, var)
        if start_adjust > 0:
            start += start_adjust
        if var == '-' or len(ref) < len(var):
            if start == stop:
                description = "{}del{}".format(start, ref)
            else:
                description = "{}_{}del{}".format(start, stop, ref)
        else:
            stop = start + 1
            description = "{}{}_{}ins{}".format(original_ref, start, stop, var)
    else:
        description = "{}{}{}".format(ref, start, var)
    return ref, var, start, stop, description


def parse_hgvsp(csq):
    '''
        Return protein ID, reference amino acid(s), variant amino acid(s),
        protein start and protein end positions extracted from HGVSp
        annotation from VEP
    '''
    protein, hgvs = csq['HGVSp'].split(':')
    ensp = ens_version_re.sub('', protein)
    try:
        indel_match = indel_re.match(hgvs)
        delins_match = delins_re.match(hgvs)
        if indel_match or delins_match:
            if indel_match:
                h_ref, start, h_ref_end, stop, indel, insertion = \
                    indel_match.groups()
            else:
                h_ref, start, insertion = delins_match.groups()
                stop = start
            try:
                ref, alt, _ = trim_ref_alt(*csq['Amino_acids'].split('/'))
                return (ensp, ref, alt, int(start), int(stop))
            except ValueError:
                raise HgvspError("Could not parse amino acids ({})".format(
                    csq['Amino_acids']) + "for HGVSp annotation '{}'".format(
                        csq['HGVSp']))
        match = snv_re.match(hgvs)
        if match:  # either SNV or single AA deletion or single AA duplication
            ref, start, var = match.groups()
            start = int(start)
            ref = protein_letters_3to1[ref]
            if var == 'del':
                var = '-'
                stop = start
            elif var == 'dup':
                var = ref + ref
                stop = start + 1
            else:
                var = protein_letters_3to1[var]
                stop = start
            return (ensp, ref, var, start, stop)
    except KeyError:
        raise HgvspError("Failed to parse amino acids in {}".format(
            csq['HGVSp']))
    raise HgvspError("Could not parse HGVSp annotation '{}'".format(
        csq['HGVSp']))


def parse_csq(csq):
    result = dict()
    result['protein'] = csq['ENSP']
    result['symbol'] = csq['SYMBOL']
    if csq['HGVSp']:
        hkeys = ['protein', 'ref_aa', 'var_aa', 'start', 'stop']
        try:
            result.update(dict((k, v) for k, v in
                               zip(hkeys, parse_hgvsp(csq))))
            result['description'] = csq['HGVSp']
            return result
        except HgvspError:
            logger.warn("Could not parse HGVSp annotation" +
                        "'{}' - falling back to  basic VEP annotations".format(
                            csq['HGVSp']))
    akeys = ['ref_aa', 'var_aa', 'start', 'stop', 'description']
    try:
        result.update(dict((k, v) for k, v in zip(akeys,
                                                  parse_amino_acids(csq))))
    except ValueError:
        logger.warn("Could not parse Amino_acids '{}' and position '{}'"
                    .format(csq['Amino_acids'], csq['Protein_position']))
        return None
    return result


def features_from_homology(homology, record, start, stop, pos, csq,
                           paralogs=None, skip_paralogs=False,
                           score_flanks=10):
    '''
        Args:
                homology:
                     single result generated get_orthologies method

                record:
                     corresponding VcfRecord

                start:
                     variant start position in protein

                stop:
                     end position of variant in protein

                pos: arbitrary label for variant position

                csq: VEP CSQ annotation for variant (i.e. a single entry
                     from VcfRecord.CSQ attribute)

                score_flanks:
                     get alignment score for this many residues either
                     side of the start/stop positions

    '''
    if skip_paralogs and 'paralog' in homology['type']:
        return []
    results = []
    s_id = homology['source']['id']
    t_id = homology['target']['id']
    s_protein = homology['source']['protein_id']
    t_protein = homology['target']['protein_id']
    s_seq = homology['source']['align_seq']
    t_seq = homology['target']['align_seq']
    s_symbol = homology['source'].get('symbol')
    t_symbol = homology['target'].get('symbol')
    align_start = get_align_pos(s_seq, start)
    if align_start < 1:  # start is out of range of sequence
        return []
    if start == stop:
        align_stop = align_start
        o_start, aa = align_pos_to_amino_acid(t_seq, align_start)
        o_stop = o_start
    else:
        align_stop = get_align_pos(s_seq, stop)
        if align_stop < 1:  # out of range - treat as SNV
            align_stop = align_start
        o_start, o_stop, aa = align_range_to_amino_acid(t_seq,
                                                        align_start,
                                                        align_stop)
    if o_start < 1 and o_stop < 1:
        return []
    elif o_stop < 1:  # only start position found, treat as SNV
        o_stop = o_start
    elif o_start < 1:  # only stop position found, treat as SNV
        o_start = o_stop
    ufeats = uniprot_lookups.get_uniprot_features(t_protein,
                                                  o_start,
                                                  o_stop)
    if paralogs and 'paralog' in homology['type']:
        logger.debug("Checking paralog for {}/{} for {}/{}".format(
            homology['target']['id'],
            homology['target']['protein_id'],
            homology['source']['id'],
            homology['source']['protein_id']))
        p_pos = o_start if o_start == o_stop else "{}-{}".format(o_start,
                                                                 o_stop)
        for paralogy in (x for x in paralogs if x['source']['protein_id'] ==
                         t_protein):
            paralog_results = features_from_homology(homology=paralogy,
                                                     record=record,
                                                     start=o_start,
                                                     stop=o_stop,
                                                     pos=p_pos,
                                                     csq=csq,
                                                     skip_paralogs=True)
            for para in paralog_results:
                para['query_symbol'] += ' ({} paralog)'.format(s_symbol)
                para['query_gene'] += ' ({} paralog)'.format(s_id)
            results.extend(paralog_results)
    if results and not ufeats:  # results from paralog lookups only
        # create dummy result for this orthology in order to output alignment
        # of paralog with original protein
        para_res = dict(query_gene=s_id,
                        query_protein=s_protein,
                        query_symbol=s_symbol,
                        query_pos=pos,
                        query_aa=s_seq[align_start:align_stop + 1],
                        ref_aa=s_seq[align_start:align_stop + 1],
                        alt_aa=csq['var_aa'],
                        homolog_gene=t_id,
                        homolog_protein=t_protein,
                        homolog_symbol=t_symbol,
                        orthology_type=homology['type'],
                        species=homology['target']['species'],
                        percent_id=homology['target']['perc_id'],
                        percent_pos=homology['target']['perc_pos'],
                        homolog_pos=o_start if o_start == o_stop
                        else "{}-{}".format(o_start, o_stop),
                        query_seq=s_seq,
                        homolog_seq=t_seq,
                        query_species=homology['source']['species'],
                        homolog_aa=aa,
                        variant_description=csq['description'],
                        variant_start=start,
                        variant_stop=stop,
                        features=None,
                        should_output=False)
        results.append(para_res)
    if ufeats is None:
        return results
    query_aa = s_seq[align_start:align_stop + 1]
    homolog_aa = t_seq[align_start:align_stop + 1]
    residue_score = score_alignment(query_aa, homolog_aa)
    if len(query_aa) == len(csq['var_aa']):
        alt_score = score_alignment(query_aa, csq['var_aa'])
    else:
        alt_score = None
    f_start = align_start - score_flanks
    f_stop = align_stop + score_flanks + 1
    flank_score = score_alignment(s_seq[f_start:f_stop],
                                  t_seq[f_start:f_stop])
    for f in ufeats:
        res = dict(chromosome=record.chrom,
                   position=record.pos,
                   id=record.id,
                   ref=record.ref,
                   alt=record.alt,
                   query_gene=s_id,
                   query_protein=s_protein,
                   query_symbol=s_symbol,
                   query_pos=pos,
                   query_aa=query_aa,
                   ref_aa=query_aa,
                   alt_aa=csq['var_aa'],
                   homolog_gene=t_id,
                   homolog_protein=t_protein,
                   homolog_symbol=t_symbol,
                   orthology_type=homology['type'],
                   species=homology['target']['species'],
                   percent_id=homology['target']['perc_id'],
                   percent_pos=homology['target']['perc_pos'],
                   homolog_pos=o_start if o_start == o_stop
                   else "{}-{}".format(o_start, o_stop),
                   query_seq=s_seq,
                   homolog_seq=t_seq,
                   query_species=homology['source']['species'],
                   homolog_aa=homolog_aa,
                   variant_description=csq['description'],
                   variant_start=start,
                   variant_stop=stop,
                   features=f,
                   residue_score=residue_score,
                   alt_score=alt_score,
                   flank_score=flank_score,
                   should_output=True)
        results.append(res)
    return results


def process_buffer(record_buffer, gene_orthologies, curr):
    results = []
    for rb in record_buffer:
        for i in range(len(rb.genes)):
            csq = parse_csq(rb.csqs[i])
            if csq is None:
                logger.warn("Skipping consequence {}".format(
                            rb.csqs[i]['HGVSp']))
                continue
            if csq['start'] == csq['stop']:
                pos = csq['start']
            else:
                pos = "{}-{}".format(csq['start'], csq['stop'])
            source_features = uniprot_lookups.get_uniprot_features(
                csq['protein'],
                csq['start'],
                csq['stop'])
            if source_features is not None:
                residue_score = score_alignment(csq['ref_aa'], csq['ref_aa'])
                if len(csq['ref_aa']) == len(csq['var_aa']):
                    alt_score = score_alignment(csq['ref_aa'], csq['var_aa'])
                else:
                    alt_score = None
                for f in source_features:
                    res = dict(chromosome=rb.record.chrom,
                               position=rb.record.pos,
                               id=rb.record.id,
                               ref=rb.record.ref,
                               alt=rb.record.alt,
                               query_gene=rb.genes[i],
                               query_symbol=csq['symbol'],
                               query_protein=csq['protein'],
                               query_pos=pos,
                               query_aa=csq['ref_aa'],
                               ref_aa=csq['ref_aa'],
                               alt_aa=csq['var_aa'],
                               homolog_gene=rb.genes[i],
                               homolog_symbol=csq['symbol'],
                               homolog_protein=csq['protein'],
                               orthology_type="self",
                               species='human',
                               percent_id=100,
                               percent_pos=100,
                               homolog_pos=pos,
                               homolog_aa=csq['ref_aa'],
                               query_species='human',
                               features=f,
                               residue_score=residue_score,
                               alt_score=alt_score,
                               variant_description=csq['description'],
                               flank_score=None,
                               should_output=True)
                    results.append(res)
            if gene_orthologies[rb.genes[i]] is None:
                continue
            for orthology in (x for x in
                              gene_orthologies[rb.genes[i]]['homologies'] if
                              x['source']['protein_id'] == csq['protein']):
                paralogs = gene_orthologies[rb.genes[i]]['paralogs'].get(
                    orthology['target']['id'])
                res = features_from_homology(homology=orthology,
                                             record=rb.record,
                                             start=csq['start'],
                                             stop=csq['stop'],
                                             pos=pos,
                                             csq=csq,
                                             paralogs=paralogs)
                results.extend(res)
    return results


def write_results_table(results, fh):
    fh.write("\n".join(["\t".join([str(res[x.lower()]) for x in variant_fields
                                   + result_header_fields] +
                                  [str(res['features'][x]) for x in
                                   uniprot_lookups.feat_fields]) for res in
                        results if res['should_output']]) + "\n")


def check_csq_annotations(vcf, logger):
    try:
        _ = vcf.header.csq_label
        logger.info("Found VEP annotations")
        return False
    except KeyError:
        try:
            _ = vcf.header.ann_label
            logger.info("Found SnpEff annotations")
            return True
        except KeyError:
            raise ValueError("No VEP or SnpEff annotations detected in " +
                             "input VCF header.")


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
        if args.snpeff_mode:
            snpeff_mode = True
        else:
            snpeff_mode = check_csq_annotations(vcf, logger)
        if snpeff_mode:
            get_csqs = get_csqs_snpeff
            global enst2ensp
            enst2ensp = read_enst2ensp_data()
        else:
            get_csqs = get_csqs_vep
        uniprot_lookups.initialize()
        if args.output is None:
            out_fh = sys.stdout
        elif args.output.endswith('.gz'):
            out_fh = gzip.open(args.output, 'wt')
        else:
            out_fh = open(args.output, 'wt')
        aln_fh = None
        if args.alignments is not None:
            aln_fh = gzip.open(args.alignments, 'wt') \
                if args.alignments.endswith('.gz') \
                else open(args.alignments, 'wt')
        out_fh.write("\t".join(header_fields) + "\n")
        n = 0
        progress_interval = args.progress
        logger.info("Beginning VCF processing")
        for record in vcf:
            if progress_interval and n % progress_interval == 0 and n != 0:
                logger.info("Processed {:,} VCF records. At {}:{}".format(
                    n, record.chrom, record.pos))
            n += 1
            csqs = get_csqs(record, curr)
            if csqs:
                these_genes = [x['Gene'] for x in csqs]
                if current_genes.isdisjoint(set(these_genes)):
                    results = process_buffer(record_buffer,
                                             gene_orthologies,
                                             curr)
                    if results:
                        write_results_table(results, out_fh)
                        if aln_fh is not None:
                            write_alignments(results, aln_fh)
                    record_buffer = []
                    current_genes.clear()
                    gene_orthologies.clear()
                record_buffer.append(ReadBuffer(record=record,
                                                genes=these_genes,
                                                csqs=csqs))
                current_genes.update(these_genes)
                for gene in (x for x in these_genes if x not in
                             gene_orthologies):
                    try:
                        gene_orthologies[gene] = get_orthologies(
                            gene, curr, args.paralog_lookups)
                    except LookupError:
                        logger.warn("Gene ID '{}' not present in database"
                                    .format(gene))
                        gene_orthologies[gene] = None
    results = process_buffer(record_buffer, gene_orthologies, curr)
    if results:
        write_results_table(results, out_fh)
        if aln_fh is not None:
            write_alignments(results, aln_fh)
    out_fh.close()
    logger.info("Finished. Processed {:,} VCF records".format(n))
