mkdir blout
for f in signal-alignments/*.faa;
do
    gname=$(basename ${f} .faa);
    if [ ! -f blout/${gname}.blout.tsv ]
    then
        echo "diamond blastp -q ${f} -o blout/${gname}.blout.tsv -d taxon-compare/GTDB-r202.dmnd -f 6 -p 2 --sensitive"
        diamond blastp -q ${f} -o blout/${gname}.blout.tsv -d taxon-compare/GTDB-r202.dmnd -f 6  -p 2 --sensitive
    fi
done;

# Fill sequences
mkdir blout-full
for f in full-seqs/*.faa;
do
    gname=$(basename ${f} .faa);
    if [ ! -f blout-full/${gname}.blout.tsv ]
    then
        echo "diamond blastp -q ${f} -o blout-full/${gname}.blout.tsv -d taxon-compare/GTDB-r202.dmnd -f 6  -p 2 --sensitive"
        diamond blastp -q ${f} -o blout-full/${gname}.blout.tsv -d taxon-compare/GTDB-r202.dmnd -f 6  -p 2 --sensitive
    fi
done;