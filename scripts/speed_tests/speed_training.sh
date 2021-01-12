SPEED_DATA=../../data/speed_tests
for n in 500 1000 2000 4000 8000 16000
do
  for i in {0..4}
  do
    echo "Started at" $(date +"%s.%3N") "Iteration $i blast|$n"
    sh constax_w_timing.sh -t -i $SPEED_DATA/oneseq.fasta -d $SPEED_DATA/silva_partition_query_sub_"$n"_gen_"$i".fasta --mem 16000 \
      -n 1 -b --conservative -f  $SPEED_DATA/training_files_"$i"_blast
    echo "Finished at" $(date +"%s.%3N")
    echo "Started at" $(date +"%s.%3N") "Iteration $i utax|$n"
    sh constax_w_timing.sh -t -i $SPEED_DATA/oneseq.fasta -d $SPEED_DATA/silva_partition_query_sub_"$n"_gen_"$i".fasta --mem 16000 \
      -n 1 --conservative -f $SPEED_DATA/training_files_"$i"_utax
    echo "Finished at" $(date +"%s.%3N")
  done
done
