#!/usr/bin/env python3
import sys
import gzip
import re
import argparse
import logging
from Bio import SwissProt

logger = logging.getLogger("UniprotParser")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)

iso_re = re.compile(r'IsoId=(\w+)(-\d+); Sequence=(.*?);')
var_seq_re = re.compile(r'[A-Z]+ -> ([A-Z]+) ')
feature_outputs = {'CROSSLNK', 'LIPID', 'MOD_RES', 'MUTAGEN', 'SITE',
                   'VARIANT'}


def parse_record(record, featfile, varfile, ensfile, seqfile=None):
    var_seqs = get_var_seqs(record)
    canonical = None
    display_isoform = None
    iso_vars = dict()
    iso_seqs = dict()  # isoform to sequence inferred from var_seqs
    for match in (z for y in filter(None, (iso_re.findall(x)
                                           for x in record.comments))
                  for z in y):
        if match[2] == 'Displayed':
            canonical = match[0]
            display_isoform = match[0] + match[1]
            if seqfile is not None:
                iso_seqs[display_isoform] = record.sequence
        elif match[2] == 'External' or match[2] == 'Not described':
            continue
        else:
            iso = match[0] + match[1]
            iso_vars[iso] = [x.strip() for x in match[2].split(', ')]
            if seqfile is not None:
                iso_seqs[iso] = construct_seq(iso,
                                              record.sequence,
                                              match[2],
                                              var_seqs)
    if not canonical:
        canonical = record.accessions[0]
        display_isoform = canonical
        if seqfile is not None:
            iso_seqs[canonical] = record.sequence
    for xref in (x for x in record.cross_references if x[0] == 'Ensembl'):
        enst = xref[1]
        ensp = xref[2]
        ensgplus = xref[3].split('. ', maxsplit=1)
        ensg = ensgplus[0]
        if len(ensgplus) > 1:
            isos = ensgplus[1].strip('[]').split(', ')
        else:
            isos = [display_isoform]
        for iso in isos:
            variants = ",".join(iso_vars.get(iso, ''))
            ensfile.write("\t".join((canonical,
                                     ensg,
                                     enst,
                                     ensp,
                                     iso,
                                     str(int(iso == display_isoform)),
                                     variants)) + "\n")
    for iso, seq in iso_seqs.items():
        seqfile.write("\t".join((canonical,
                                 iso,
                                 str(int(iso == display_isoform)),
                                 seq)) + "\n")
    for ft in (x for x in record.features if x.type in feature_outputs):
        featfile.write("\t".join((canonical,
                                  ft.type,
                                  str(ft.location.start + 1),
                                  str(ft.location.end),
                                  ft.qualifiers['note'])) + "\n")
    for var_id, var in var_seqs.items():
        varfile.write("\t".join((var_id,
                                 str(var['start']),
                                 str(var['stop']),
                                 str(len(var['seq'])),
                                 var['seq'])) + "\n")


def get_var_seqs(record):
    varseqs = dict()
    for ft in [x for x in record.features if x.type == 'VAR_SEQ']:
        if ft.qualifiers['note'].startswith('Missing'):
            seq = ''
        else:
            seq_match = var_seq_re.match(ft.qualifiers['note'])
            if seq_match:
                seq = seq_match.group(1)
            else:
                logger.warn("Could not parse VAR_SEQ: {}".format(ft.id))
                continue
        if ft.id in varseqs:
            logger.warn("Duplicate VAR_SEQ ID '{}'".format(ft.id))
        else:
            varseqs[ft.id] = dict(start=int(ft.location.start),  # 0-based
                                  stop=int(ft.location.end),  # end inclusive
                                  seq=seq)
    return varseqs


def construct_seq(isoform, seq, mods, variants):
    vars = []
    for mod in mods.split(', '):
        var = variants.get(mod.strip())
        if var is None:
            logger.warn("Got unexpected VAR_SEQ: '{}' ".format(mod) +
                        "for isoform '{}'".format(isoform))
            continue
        vars.append(var)
    # Apply variants starting from end of sequence, else mod coordinates will
    # not be accurate for downstream modifications anymore
    vars.sort(key=lambda x: (x['start'], x['stop']), reverse=True)
    for var in vars:
        seq = seq[:var['start']] + var['seq'] + seq[var['stop']:]
    return seq


def main(data_file, output_prefix, progress_interval=10000, no_gzip=False,
         output_seq=False):
    open_func = open
    if data_file.endswith(".gz"):
        open_func = gzip.open
    write_func = open if no_gzip else gzip.open
    suffix = 'tab' if no_gzip else 'tab.gz'
    feat_fh = write_func(output_prefix + '.features.' + suffix, 'wt')
    var_fh = write_func(output_prefix + '.seq_vars.' + suffix, 'wt')
    ens_fh = write_func(output_prefix + '.ens_ids.' + suffix, 'wt')
    if output_seq:
        seq_fh = write_func(output_prefix + '.isoform_seqs.' + suffix, 'wt')
        seq_fh.write("\t".join("UniprotID Isoform Displayed Seq".split()) + "\n")
    else:
        seq_fh = None
    ens_fh.write("\t".join('''UniprotID Gene Transcript Protein
                              Isoform Displayed Variants'''.split()) + "\n")
    feat_fh.write("\t".join('''UniprotID Feature Start Stop
                               Description'''.split()) + "\n")
    var_fh.write("\t".join('Variation Start Stop Length Seq'.split()) + "\n")
    with open_func(data_file, 'rt') as filehandle:
        n = 0
        for record in SwissProt.parse(filehandle):
            parse_record(record,
                         featfile=feat_fh,
                         varfile=var_fh,
                         seqfile=seq_fh,
                         ensfile=ens_fh)
            n += 1
            if n % progress_interval == 0:
                logger.info("Processed {:,} records".format(n))
    for fh in (seq_fh, feat_fh, ens_fh, var_fh):
        if fh is not None:
            fh.close()
    logger.info("Finished - processed {:,} records".format(n))


def get_options():
    parser = argparse.ArgumentParser(
        description='''Parse Uniprot/Swissprot flat file into tables required
        by ORA''')
    parser.add_argument('data_file', help='''Uniprot flat file (gzipped is fine)
    to process (e.g. text file from https://www.uniprot.org/downloads).''')
    parser.add_argument('output_prefix', help='''Prefix for output files. This
    program will create files named <prefix>.features.tab.gz,
    <prefix>.seq_vars.tab.gz, <prefix>.ens_ids.tab.gz and optionally
    <prefix>.isoform_seqs.tab.gz.''')
    parser.add_argument('--output_seq', action='store_true', help='''Output an
    additional data file with the peptide sequences of all protein isoforms.''')
    parser.add_argument('--no_gzip', action='store_true', help='''Do not
    compress output files with gzip.''')
    return parser


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
