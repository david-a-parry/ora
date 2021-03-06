#!/bin/bash

if [ $# -gt 2 ] || [ $# -lt 1 ]
then
    echo "Usage: $0 <output_db_file> [download_directory]"
    exit 1
fi
DB=$1
DLDIR=$2

if [ -e "$DB" ]
then
    echo "Database '$DB' already exists - please delete or choose another filename"
    exit 2
fi

if [ "$DLDIR" == "" ]
then
    DLDIR=./
else
    mkdir -p $DLDIR
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
set -euo pipefail

echo $(date) Finding latest compara database
CURRENT_URL=ftp://ftp.ensembl.org/pub/current_mysql/ 
COMPARA=$(curl $CURRENT_URL | grep compara | rev | cut -f 1 -d ' ' | rev)
ENSURL=${CURRENT_URL}/${COMPARA}

echo $(date) Using URL: $ENSURL

echo $(date) Retrieving tables
for TABLE in seq_member sequence homology homology_member gene_member ncbi_taxa_name
do
    echo $(date) Getting $TABLE
    wget --directory-prefix=$DLDIR -c ${ENSURL}/${TABLE}.txt.gz
done

echo $(date) Identifying human-relevant homologies
python3 $DIR/get_human_relevant_homologies.py "$DLDIR"

echo $(date) Identifying homology-relevant genes
python3 $DIR/get_homology_gene_members.py "$DLDIR"

echo $(date) Identifying homology-relevant sequences
python3 $DIR/get_homology_sequences.py "$DLDIR"

echo "$(date) Creating output database '$DB'"
cat $DIR/create_tables.sql | sqlite3 $DB

echo $(date) Importing tables
for TABLE in seq_member sequence homology_member gene_member
do
    # create temporary init script
    commandfile=$(mktemp)
    cat <<EOF > $commandfile
.mode tab
.import /dev/stdin $TABLE
EOF
    # import
    echo $(date) Importing $TABLE
    gzip -d -c ${DLDIR}/ora_${TABLE}.txt.gz | sqlite3 --init $commandfile $DB
    rm $commandfile
done

commandfile=$(mktemp)
cat <<EOF > $commandfile
.mode tab
.import /dev/stdin ncbi_taxa_name
EOF

echo $(date) Importing ncbi_taxa_name
gzip -d -c ${DLDIR}/ncbi_taxa_name.txt.gz | \
    grep -e "scientific name" -e "genbank common name" | \
    sqlite3 --init $commandfile $DB
rm $commandfile

echo $(date) Finished imports - converting NULL fields
cat $DIR/fix_nulls.sql | sqlite3 $DB

echo $(date) Done
