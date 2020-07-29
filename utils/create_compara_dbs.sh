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
fi

ENSURL=ftp://ftp.ensembl.org/pub/release-100/mysql/ensembl_compara_100
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
TEMPLATE=${DIR}/compara_sqlite3_template.db

echo "$(date) Creating output database '$DB'"
cp -nv $TEMPLATE $DB

echo $(date) Retrieving tables
for TABLE in sequence homology homology_member gene_member
do
    echo $(date) Getting $TABLE
    wget --directory-prefix=$DLDIR -c ${ENSURL}/${TABLE}.txt.gz
    commandfile=$(mktemp)
    # create temporary init script
    cat <<EOF > $commandfile
.mode tab gene_member
.import /dev/stdin $TABLE
EOF
    # import
    echo $(date) Importing $TABLE
    gzip -d -c ${DLDIR}/${TABLE}.txt.gz | sqlite3 --init $commandfile $DB
done
echo $(date) Finished imports
echo $?
