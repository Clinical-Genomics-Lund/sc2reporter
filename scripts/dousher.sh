#!/bin/bash

source activate usher

today=`date +%y%m%d`

cd /data/bnf/sarscov2/usher
echo Downloading new database
./getdb.sh
echo Preparing metadata and fasta files
./prepusher.sh
echo Running UShER
./doushersampled.sh
cp ${today}/single-subtree.nh /data/bnf/sarscov2/results/usher/${today}-single-subtree.nh
cp meta_${today}.tsv /data/bnf/sarscov2/results/usher/${today}-meta.tsv

echo "cat ../results/usher/${today}-single-subtree.nh | xclip -se c"
echo "cat ../results/usher/${today}-meta.tsv | xclip -se c"

cd -
