#!/bin/bash
SCRIPT_FIT=print_sla_sub_fastTopics.R
FILE_OUT=print_sla_sub_fastTopics
for K in 50
do
  Rscript ${SCRIPT_FIT} ${K} > ${FILE_OUT}_${K}.out
done

SCRIPT_FIT=print_sla_sub_snn.R
FILE_OUT=print_sla_sub_snn
for K in 50
do
  Rscript ${SCRIPT_FIT} ${K} > ${FILE_OUT}_${K}.out
done