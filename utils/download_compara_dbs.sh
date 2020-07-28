#!/bin/bash

PASS=mysqlpassword
MYUSER=mysqluser
ENSURL=ftp://ftp.ensembl.org/pub/release-100/mysql/ensembl_compara_100

echo $(date) Retrieving tables
for TABLE in  sequence homology homology_member gene_member
do
    echo $(date) Getting $TABLE
    wget ${ENSURL}/${TABLE}.txt.gz
done

for TABLE in  sequence homology homology_member gene_member
do
    echo $(date) decompressing $TABLE
    gzip -d ${TABLE}.txt.gz
    echo $(date) importing $TABLE
    mysqlimport -u $MYUSER -p${PASS} --fields-terminated-by='\t' --fields-escaped-by=\\ ensembl_compara_100  -L ${TABLE}.txt -C
    echo $(date) re-compressing $TABLE
    gzip ${TABLE}.txt.gz
done

echo $(date) Done
