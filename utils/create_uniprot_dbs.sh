#!/bin/bash

if [ $# -gt 3 ] || [ $# -lt 1 ]
then
    echo "Usage: $0 <output_directory> [download_directory] [ftp_url]"
    exit 1
fi
OUTDIR=$1
DLDIR=$2
FTPURL=$3 # e.g. ftp://ftp.expasy.org/databases/uniprot or ftp://ftp.ebi.ac.uk/pub/databases/uniprot

if [ "$DLDIR" == "" ]
then
    DLDIR=./
else
    mkdir -p $DLDIR
fi
mkdir -p $OUTDIR

if [ "$FTPURL" == "" ]
then
    FTPURL=ftp://ftp.uniprot.org/pub/databases/uniprot
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
set -euo pipefail

if [ -e "$DLDIR/uniprot_sprot.dat.gz" ]
then
    echo $(date) \| Data file "$DLDIR/uniprot_sprot.dat.gz" already exists
    echo $(date) \| Using pre-existing datafile. Delete file and run again if you want to download a new version
else
    echo $(date) \| Retrieving UniProt flat file
    wget --directory-prefix=$DLDIR $FTPURL/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
fi

echo $(date) \| Extracting data to $OUTDIR
python3 $DIR/parse_uniprot_flat.py $DLDIR/uniprot_sprot.dat.gz $OUTDIR/sprot

echo $(date) \| Finished
