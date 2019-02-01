#!/bin/bash

cd ../genomes/curated-diffs
gdtools INTERSECT Ara-3_33000gen_CZB152.gd Ara-3_33000gen_CZB154.gd -o initLCA.gd
gdtools UNION initLCA.gd clpA-serW-mut.gd -o LCA.gd
rm initLCA.gd
gdtools APPLY -o LCA.gff3 -f GFF3 -r ../REL606.7.gbk LCA.gd
gdtools APPLY -o LCA.fasta -f FASTA -r ../REL606.7.gbk LCA.gd

## horrible hacks to make a genbank file from the LCA ggf3 file
## for the dice analysis.
python ../../src/external/gff_to_genbank.py LCA.gff3 LCA.fasta
sed -i '' -e 's/Name=/gene=/g' ../genomes/curated-diffs/LCA.gbk
sed -i '' -e 's/ID=/locus_tag=/g' ../genomes/curated-diffs/LCA.gbk
