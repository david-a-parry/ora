#!/usr/bin/env python3
import sys
from biomart import BiomartServer


def get_ens_ids(output):
    server = BiomartServer("http://www.ensembl.org/biomart")
    ensgenes = server.datasets['hsapiens_gene_ensembl']
    response = ensgenes.search({'filters': {'biotype': 'protein_coding'},
                                'attributes': ['ensembl_gene_id',
                                               'ensembl_transcript_id',
                                               'ensembl_peptide_id']},
                               header=1)
    if not response.ok:
        response.raise_for_status()
    with open(output, 'wt') as fh:
        for line in response.iter_lines(decode_unicode=True):
            cols = line.split('\t')
            if len(cols) == 3 and cols[2]:
                fh.write(line + "\n")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} <output_filename>".format(sys.argv[0]))
    get_ens_ids(sys.argv[1])
