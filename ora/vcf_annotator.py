import logging
import sqlite3
from collections import namedtuple
from ora.local import ensg_lookup, get_homologies
from ora.uniprot_lookups import logger as unipro_logger
from vase.vcf_reader import VcfReader

logger = logging.getLogger("ORA")
logger.setLevel(logging.INFO)
ReadBuffer = namedtuple('record protein positions')


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
    unipro_logger.setLevel(logger.level)
    for handler in unipro_logger.handlers:
        handler.setLevel(logger.level)
    conn = sqlite3.connect(args.db)
    curr = conn.cursor()
    buffer = []
    current_genes = set()
    out = '-' if args.output is None else args.output
    with VcfReader(args.input) as vcf:
        out_vcf = pysam.VariantFile(out, 'w', header=vcf.header.header)
        for record in vcf:
            raise NotImplementedError("VCF annotations not implemented yet!")
            # if does not match our consequence classes and buffer empty
            # just write variant and move on
            #out_vcf.write(record)
            # if in homology class get homology details for gene and store them
            # store record, ENSPs and positions to lookup in ReadBuffer
            # store gene ID in current_genes

            # if current_genes and CSQ genes are disjoint, process and empty
            # buffer, clear current_genes

#    default_args = [input, db, paralog_lookups=False, line_length=60,
#                    all_homologs=False, output_alignments=None, quiet=False,
#                    debug=False, silent=False]
