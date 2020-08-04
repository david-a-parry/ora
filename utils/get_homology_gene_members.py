#!/usr/bin/env python3
import os
import sys
import gzip


def main(dldir='.'):
    '''
        Get genes from homology members. Must be run AFTER
        get_human_relevant_homologies.py using the same download
        directory.
    '''
    gene_members = set()
    n = 0
    with gzip.open(os.path.join(dldir, "ora_homology_member.txt.gz"),
                   'rt') as infile:
        for line in infile:
            gene_members.add(line.split()[1])
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} genes\n".format(n))
    sys.stderr.write("Got {:,} seq members\n".format(len(gene_members)))
    n = 0
    p = 0
    infile = os.path.join(dldir, "gene_member.txt.gz")
    outfile = os.path.join(dldir, "ora_gene_member.txt.gz")
    with gzip.open(infile, 'rt') as fh, gzip.open(outfile, 'wt') as out:
        for line in fh:
            cols = line.rstrip().split("\t")
            if cols[4] == '9606' or cols[0] in gene_members:
                # `gene_member_id` integer NOT NULL PRIMARY KEY AUTOINCREMENT
                # ,  `stable_id` varchar(128) NOT NULL
                # ,  `version` integer DEFAULT 0
                # ,  `source_name` text  NOT NULL
                # ,  `taxon_id` integer  NOT NULL
                # ,  `genome_db_id` integer  DEFAULT NULL
                # ,  `biotype_group` text  NOT NULL DEFAULT 'coding'
                # ,  `canonical_member_id` integer  DEFAULT NULL
                # ,  `description` text DEFAULT NULL
                # ,  `dnafrag_id` integer  DEFAULT NULL
                # ,  `dnafrag_start` integer DEFAULT NULL
                # ,  `dnafrag_end` integer DEFAULT NULL
                # ,  `dnafrag_strand` integer DEFAULT NULL
                # ,  `display_label` varchar(128) DEFAULT NULL
                # ,  UNIQUE (`stable_id`)
                # we keep cols `gene_member_id`, `stable_id`, `version`
                #              `taxon_id`, `biotype_group`,
                #              `canonical_member_id`, `display_label`
                out.write("\t".join(cols[:3] + [cols[4], cols[6], cols[7],
                                                cols[13]]) + "\n")
                p += 1
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} sequences. ".format(n) +
                                 "Got {:,} homology gene_members.\n"
                                 .format(p))
    sys.stderr.write("Finised. Processed {:,} genes. ".format(n) +
                     "Got {:,} homology genes.\n".format(p))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} <download_dir>\n".format(sys.argv[0]))
    main(sys.argv[1])
