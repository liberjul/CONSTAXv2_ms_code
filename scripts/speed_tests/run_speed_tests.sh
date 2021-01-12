#!/bin/sh
SPEED_DATA=../../data/speed_tests
sh subsample_speed_sets.sh
sh speed_training.sh > $SPEED_DATA/training_times.out
sh speed_classifying.sh > $SPEED_DATA/classifying_times.out

python run_time_extract_train_n.py $SPEED_DATA/training_times.out $SPEED_DATA/speed_training_t1_nvary.csv
python run_time_extract_classify_nthreads.py $SPEED_DATA/classifying_times.out $SPEED_DATA/speed_classify_comb.csv
