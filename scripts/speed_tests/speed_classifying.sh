SPEED_DATA=../../data/speed_tests
for q in 1 2 4
do
  for NTHREADS in 1 2 4 8 16 32 64 96
  do
    for i in {0..4}
    do
      echo "Started at" $(date +"%s.%3N") "|$NTHREADS|blast|$q"
      sh constax_w_timing.sh -i $SPEED_DATA/seq_set"$q"k/query_iter"$i".fasta -d $SPEED_DATA/seq_set"$q"k/ref_iter"$i".fasta \
      --mem 16000 -n $NTHREADS -b --conservative -f $SPEED_DATA/training_files_"$i"_blast -x $SPEED_DATA/seq_set"$q"k/tax_"$NTHREADS"_core -o $SPEED_DATA/seq_set"$q"k/out_"$NTHREADS"_core
      echo "Finished at" $(date +"%s.%3N")
      echo "Started at" $(date +"%s.%3N") "|$NTHREADS|utax|$q"
      sh constax_w_timing.sh -i $SPEED_DATA/seq_set"$q"k/query_iter"$i".fasta -d $SPEED_DATA/seq_set"$q"k/ref_iter"$i".fasta \
      --mem 16000 -n $NTHREADS --conservative -f $SPEED_DATA/training_files_"$i"_utax -x $SPEED_DATA/seq_set"$q"k/tax_"$NTHREADS"_core -o $SPEED_DATA/seq_set"$q"k/out_"$NTHREADS"_core
      echo "Finished at" $(date +"%s.%3N")
    done
  done
done
