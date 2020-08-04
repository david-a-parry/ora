#!/usr/bin/env python3
import sys
import os
import gzip


class HomologyRecord(object):

    def __init__(self, hom_line, member1, member2):
        self.m1_cols = member1.split()
        self.m2_cols = member2.split()
        self.homology_id = self.m1_cols[0]
        assert(self.homology_id == self.m2_cols[0])
        self.gene_members = (self.m1_cols[1], self.m2_cols[1])
        h_cols = hom_line.split()
        assert(self.homology_id == h_cols[0])
        self.hom_columns = [h_cols[2], h_cols[14]]

    def __str__(self):
        return "\n".join(self.output_lines())

    def output_lines(self):
        l1 = "\t".join(self.m1_cols + self.hom_columns)
        l2 = "\t".join(self.m2_cols + self.hom_columns)
        return (l1, l2)


class HomologyIter(object):
    ''' Iterate over homology member and homology tables together '''

    def __init__(self, hom_file, member_file):
        self.hom_fh = gzip.open(hom_file, 'rt')
        self.member_fh = gzip.open(member_file, 'rt')

    def __iter__(self):
        return self

    def __next__(self):
        member1, member2 = None, None
        try:
            member1 = next(self.member_fh)
            member2 = next(self.member_fh)
        except StopIteration:
            if member1 is not None:
                raise ValueError("Member file ended before reaching second " +
                                 "record.")
            try:
                _ = next(self.hom_fh)
                raise ValueError("Member file ended before homology file")
            except StopIteration:
                self.hom_fh.close()
                self.member_fh.close()
                raise StopIteration
        try:
            homology = next(self.hom_fh)
        except StopIteration:
            raise ValueError("Homology file ended before member file")
        return HomologyRecord(homology, member1, member2)


def main(dldir='.'):
    '''Print homology_members for all human-centric homlogies.'''
    gene_members = set()
    n = 0
    with gzip.open(os.path.join(dldir, "gene_member.txt.gz"), 'rt') as infile:
        for line in infile:
            cols = line.split()
            if cols[4] == '9606':
                gene_members.add(cols[0])
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} genes\n".format(n))
    sys.stderr.write("Got {:,} gene members\n".format(len(gene_members)))
    n = 0
    p = 0
    hom_iter = HomologyIter(os.path.join(dldir, "homology.txt.gz"),
                            os.path.join(dldir, "homology_member.txt.gz"))
    with gzip.open(os.path.join(dldir, "ora_homology_member.txt.gz"),
                   'wt') as out:
        for homology in hom_iter:
            if any(x in gene_members for x in homology.gene_members):
                out.write(str(homology) + "\n")
                p += 1
            n += 1
            if n % 100000 == 0:
                sys.stderr.write("Processed {:,} homologies. ".format(n) +
                                 "Got {:,} human-centric relationships.\n"
                                 .format(p))
    sys.stderr.write("Finised. Processed {:,} homologies. ".format(n) +
                     "Got {:,} human-centric relationships.\n".format(p))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit("Usage: {} <download_dir>\n".format(sys.argv[0]))
    main(sys.argv[1])
