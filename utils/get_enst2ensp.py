#!/usr/bin/env python3
import sys
from biomart import BiomartServer


def get_ens_ids(output):
    server = BiomartServer( "http://www.ensembl.org/biomart")
    ensgenes = server.datasets['hsapiens_gene_ensembl']
    response = ensgenes.search({'filters': {'biotype': 'protein_coding'},
                                'attributes': ['ensembl_gene_id',
                                               'ensembl_transcript_id',
                                               'ensembl_peptide_id']},
                               header=1)
    if not response.ok:
       response.raise_for_status()
    with open(output, 'wt') as fh:
        fh.write(response.text)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} <output_filename>".format(sys.argv[0]))
    get_ens_ids(sys.argv[1])
