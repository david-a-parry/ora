#!/usr/bin/env python3
import argparse


def get_options():
    parser = argparse.ArgumentParser(
        description='''Search for protein features in orthologs''')
    parser.add_argument("-g", "--gene", required=True,
                        help='''Gene/protein ID to search with. Accession types
                        will be inferred from input and the corresponding
                        Ensembl gene ID will be identified.''')
    parser.add_argument("-p", "--pos", required=True, type=int, nargs='+',
                        help='''Position(s) in protein to search. This must be
                        relative to the canonical Ensembl transcript as used by
                        Ensembl's Compara database.''')
    parser.add_argument("-d", "--db",
                        help='''Local sqlite3 database to use for lookups. If
                        not provided lookups will be performed via Ensembl REST
                        lookups.''')
    parser.add_argument("--paralog_lookups", action='store_true',
                        help='''Also perform homology lookups on paralogs
                        identified.''')
    parser.add_argument("--all_homologs", action='store_true',
                        help='''Output alignment information for all homologs
                        even if no features are found.''')
    parser.add_argument("--output_alignments", help='''Output alignments of
                        with hits to this file.''')
    parser.add_argument("--line_length", type=int, default=60,
                        help='''Length of sequence lines in alignment output.
                        Default=60.''')
    parser.add_argument("--timeout", type=int, default=10,
                        help='''Lookup timeout (in seconds) for Ensembl REST
                        queries. Default=10''')
    parser.add_argument("--max_retries", type=int, default=2,
                        help='''Maximum retry attempts for Ensembl REST
                        queries. Default=2''')
    parser.add_argument("--debug", action='store_true',
                        help='Add debugging messages to output.')
    parser.add_argument("--quiet", action='store_true',
                        help='Only output warning logger messages.')
    parser.add_argument("--silent", action='store_true',
                        help='Only output error logger messages.')
    return parser


def main(gene, pos, db=None, paralog_lookups=False, line_length=60,
         timeout=10.0, max_retries=2, all_homologs=False,
         output_alignments=None, quiet=False, debug=False, silent=False):
    if db is None:
        from ora.remote_lookups import remote_lookups
        remote_lookups(gene=gene, pos=pos, paralog_lookups=paralog_lookups,
                       line_length=line_length, timeout=timeout,
                       max_retries=max_retries, all_homologs=all_homologs,
                       output_alignments=output_alignments, quiet=quiet,
                       debug=debug, silent=silent)
    else:
        from ora.local_lookups import local_lookups
        local_lookups(gene=gene, pos=pos, db=db,
                      paralog_lookups=paralog_lookups, line_length=line_length,
                      all_homologs=all_homologs,
                      output_alignments=output_alignments, quiet=quiet,
                      debug=debug, silent=silent)


if __name__ == '__main__':
    argparser = get_options()
    args = argparser.parse_args()
    main(**vars(args))
