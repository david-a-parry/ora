#!/usr/bin/env python3
import sys
import gzip
import re
import logging
from collections import defaultdict

logger = logging.getLogger("UniprotParser")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)

acc_split =  re.compile(r';\s?')
n_aa_re = re.compile(r'.*\s(\d+)\s+AA\.')
iso_re = re.compile(r'IsoId=(\w+)(-\d+);.*Sequence=(.*);')

def parse_uniprot_flat(data):
    pids = []
    match = n_aa_re.match(data['ID'][0])
    if not match:
        logger.warning("Could not parse data starting {}".format(data[0]))
        return
    else:
        p_length = match.group(1)
    for acc in data['AC']:
        accessions = [x for x in acc_split.split(match.group(1)) if x]
    if not accessions:
        logger.warning("Could not parse data starting {}".format(data['ID']))
        return
    seq = ''.join(x.replace(' ', '') for x in data['SQ'][0].split("\n")[1:])
    canonical = None
    isoforms = dict()
    if len(accessions) > 1:
        for match in filter(None, (iso_re.match(x) for x in data['CC'])):
            if match.group(3) == 'Displayed':
                canonical = match.group(1)
                isoforms[match.group(1) + match.group(2)] = seq
            else:
                isoforms[match.group(1) + match.group(2)] = construct_seq(
                    seq, match.group(3))

    else:
        canonical = accessions[0]
    var_seqs = get_var_seqs(data)

    #TODO

def get_var_seqs(data):
    for ft in data['FT']:
        pass #TODO

def construct_seq(seq, mods):
    mods = mods.split(', ')
    #TODO


def main(f):
    open_func = open
    if f.endswith(".gz"):
        open_func = gzip.open
    data = defaultdict(list)
    current_data = ''
    current_key = ''
    with open_func(f, 'rt') as infile:
        n = 0
        for line in infile:
            n += 1
            line = line.rstrip()
            if line == '//':
                if not data:
                    raise ValueError("Encountered record separator without " +
                                     "any preceding data at line {}".format(n))
                if current_data:
                    data[current_key].append(current_data)
                parse_uniprot_flat(data)
                data.clear()
            elif line.startswith(' '):
                current_data += '\n{}'.format(line)
            else:
                if current_data:
                    data[current_key].append(current_data)
                current_key, current_data = line.split(maxsplit=1)
    if data:
        parse_uniprot_flat(data)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} uniprot_sprot.dat.gz".format(sys.argv[0]))
    main(sys.argv[1])
