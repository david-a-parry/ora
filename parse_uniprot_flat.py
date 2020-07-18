#!/usr/bin/env python3
import sys
import gzip
import re
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
feature_outputs = {'CROSSLNK', 'LIPID', 'MOD_RES', 'MUTAGEN' 'SITE', 'VARIANT'}


def parse_record(record, featfile, seqfile, ensfile):
    var_seqs = get_var_seqs(record)
    canonical = None
    display_isoform = None
    isoforms = dict()
    for match in (z for y in filter(None, (iso_re.findall(x)
                                           for x in record.comments))
                  for z in y):
        if match[2] == 'Displayed':
            canonical = match[0]
            display_isoform = match[0] + match[1]
            isoforms[display_isoform] = record.sequence
        elif match[2] == 'External' or match[2] == 'Not described':
            continue
        else:
            iso = match[0] + match[1]
            isoforms[iso] = construct_seq(iso,
                                          record.sequence,
                                          match[2],
                                          var_seqs)
    if not canonical:
        canonical = record.accessions[0]
        display_isoform = canonical
        isoforms[canonical] = record.sequence
    for xref in (x for x in record.cross_references if x[0] == 'Ensembl'):
        enst = xref[1]
        ensp = xref[2]
        ensgplus = xref[3].split('. ', maxsplit=1)
        ensg = ensgplus[0]
        if len(ensgplus) > 1:
            isos = ensgplus[1].strip('[]').split(', ')
        else:
            isos = [canonical]
        for iso in isos:
            ensfile.write("\t".join((canonical,
                                     ensg,
                                     enst,
                                     ensp,
                                     iso)) + "\n")
    for iso, seq in isoforms.items():
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


def main(f, prefix, progress_interval=10000, use_gzip=True):
    open_func = open
    if f.endswith(".gz"):
        open_func = gzip.open
    write_func = gzip.open if use_gzip else open
    suffix = 'tab.gz' if use_gzip else 'tab'
    feat_fh = write_func(prefix + '.features.' + suffix, 'wt')
    seq_fh = write_func(prefix + '.isoform_seqs.' + suffix, 'wt')
    ens_fh = write_func(prefix + '.ens_ids.' + suffix, 'wt')
    seq_fh.write("\t".join("UniprotID Isoform Displayed Seq".split()) + "\n")
    ens_fh.write("\t".join('''UniprotID Gene Transcript Protein
                              Isoform'''.split()) + "\n")
    feat_fh.write("\t".join('''UniprotID Feature Start Stop
                               Description'''.split()) + "\n")
    with open_func(f, 'rt') as filehandle:
        n = 0
        for record in SwissProt.parse(filehandle):
            parse_record(record, feat_fh, seq_fh, ens_fh)
            n += 1
            if n % progress_interval == 0:
                logger.info("Processed {:,} records".format(n))
    for fh in (seq_fh, feat_fh, ens_fh):
        fh.close()
    logger.info("Finished - processed {:,} records".format(n))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit("Usage: {} <uniprot_sprot.dat.gz> <output_prefix>".format(
            sys.argv[0]))
    main(*sys.argv[1:])
