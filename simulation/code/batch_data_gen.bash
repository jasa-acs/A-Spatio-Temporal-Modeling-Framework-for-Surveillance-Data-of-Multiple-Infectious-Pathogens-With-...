#!/bin/bash

nepi=100

for((r=1; r<=$nepi; r++))
do
        Rscript simulation/code/data_gen.R ${r}
	Rscript simulation/code/data_gen_2infweek.R ${r}
	Rscript simulation/code/lab_agg.R ${r}
done
