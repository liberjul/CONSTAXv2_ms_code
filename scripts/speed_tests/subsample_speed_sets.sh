SILVA_DATA=../../data/classification/silva
SPEED_DATA=../../data/speed_tests
for q in 1 2 4
do
  mkdir $SPEED_DATA/seq_set"$q"k
  python ../make_quer_ref_files.py -i $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact.fasta -q "$q"000 -r 10000 -k 5 \
    --rname $SPEED_DATA/seq_set"$q"k/ref_iter --qname $SPEED_DATA/seq_set"$q"k/query_iter
done
for n in 500 1000 2000 4000 8000 16000
do
  python ../make_quer_ref_files.py -i $SILVA_DATA/SILVA_138_SSURef_tax_silva_bact.fasta -q 1 -r $r -k 5 \
    --rname $SPEED_DATA/silva_partition_query_sub_"$n"_gen_"$i" --qname $SPEED_DATA/train_quer_iter
done
