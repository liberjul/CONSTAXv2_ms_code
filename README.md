# CONSTAXv2 Manuscript Code

## How to run:
```
cd ./scripts/classification
sh classification_performance.sh
cd ../speed_tests
sh run_speed_tests.sh
```
Note: The scripts will take a VERY long time to run (on the order of weeks). They are meant to show the code used, but not necessarily run in serial. Dependent on machine architecture, tests for threads and memory may not work.

## Making plots and tables

`./scripts/speed_tests/speed_plots.R`, `./scripts/classification/format_class_table.py`,  and `./scripts/classification/runtime_and_partion_cv_plots.R` can be used to recreate plots and tables.
