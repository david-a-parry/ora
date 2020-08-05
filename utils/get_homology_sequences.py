#!/usr/bin/env python3
import os
import sys
import gzip


def main(dldir='.'):
    '''
        Get sequences from homology members. Must be run AFTER
        get_human_relevant_homologies.py using the same download
        directory.
    '''
    seq_members = set()
    n = 0
    with gzip.open(os.path.join(dldir, "ora_homology_member.txt.gz"),
                   'rt') as infile:
        for line in infile:
            seq_members.add(line.rstrip().split("\t")[2])
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} genes\n".format(n))
    sys.stderr.write("Got {:,} seq members\n".format(len(seq_members)))
    n = 0
    p = 0
    seq_ids = set()
    memberfile = os.path.join(dldir, "seq_member.txt.gz")
    outfile = os.path.join(dldir, "ora_seq_member.txt.gz")
    with gzip.open(memberfile, 'rt') as fh, gzip.open(outfile, 'wt') as out:
        for line in fh:
            cols = line.rstrip().split("\t")
            if cols[0] in seq_members:
                seq_ids.add(cols[6])
                # | seq_member_id | stable_id | version | sequence_id |
                out.write("\t".join((cols[0], cols[1], cols[2], cols[6],)) +
                          "\n")
                p += 1
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} seq_members. ".format(n) +
                                 "Got {:,} homology seq_members.\n"
                                 .format(p))
    sys.stderr.write("Finised. Processed {:,} seq_members. ".format(n) +
                     "Got {:,} homology sequences members.\n".format(p))
    seqfile = os.path.join(dldir, "sequence.txt.gz")
    outfile = os.path.join(dldir, "ora_sequence.txt.gz")
    n = 0
    p = 0
    with gzip.open(seqfile, 'rt') as fh, gzip.open(outfile, 'wt') as out:
        for line in fh:
            cols = line.rstrip().split("\t")
            if cols[0] in seq_ids:
                p += 1
                out.write("\t".join((cols[0], cols[1], cols[3])) + "\n")
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} sequences. ".format(n) +
                                 "Got {:,} homology sequences.\n"
                                 .format(p))

    sys.stderr.write("Processed {:,} sequences. ".format(n) +
                     "Got {:,} homology sequences.\n"
                     .format(p))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} <download_dir>\n".format(sys.argv[0]))
    main(sys.argv[1])
