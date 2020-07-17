#!/usr/bin/env python3
import sys
import gzip
import re
import logging
from Bio import SeqIO
logger = logging.getLogger("UniprotParser")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)



def main(f):
    open_func = open
    if f.endswith(".gz"):
        open_func = gzip.open
    data = defaultdict(list)
    current_data = ''
    current_key = ''
    with open_func(f, 'rt') as infile:
        for record in SeqIO.parse(handle, 'uniprot-xml'):
            parse_record(record)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} uniprot_sprot.dat.gz".format(sys.argv[0]))
    main(sys.argv[1])
